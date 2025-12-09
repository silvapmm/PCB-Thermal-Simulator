# Module for the analytic thermal simulation model
# This module provides a simple analytic approximation for temperature estimation in PCBs

import numpy as np
import logging  # Fallback logger if none provided
import math      # For pi in via area calculations (if extended)
import wx        # For dialog and UI components
import matplotlib.pyplot as plt  # For chart generation
from io import BytesIO  # Save chart to memory buffer
import os       # For file paths in Save
import datetime # For timestamps in Save
import shutil # Import for high-level file operations (e.g., copying files or directories)
import config   # Import config values (sizes, colors, filenames)

class AnalyticSolver:
    """
    @brief Class for analytic thermal simulation model.
    @details Provides quick estimation using closed-form formulas.
             Supports 'basic' (total PCB area) and 'advanced' (heater area via shoelace) modes.
    """
    @staticmethod
    def simulate(json_config, logger=None):
        """
        @brief Perform analytic thermal simulation.
        @param json_config Dict with pcb_params, solver_mode, and elements.
        @param logger Optional logging.Logger instance.
        @return dict with estimated metrics (value and unit).
        """
        # Use module logger if none provided
        if logger is None:
            logger = logging.getLogger(__name__)

        # Extract PCB parameters
        pcb = json_config['pcb_params']
        T_amb = pcb['ambient_temp_global']  # Ambient temperature (°C)
        h_global = pcb['convection_h_global']  # Global convection coefficient (W/m²K)
        width = pcb['width']  # PCB width (mm)
        height = pcb['height']  # PCB height (mm)
        k_copper = pcb['copper_conductivity']  # Copper thermal conductivity (W/mK)
        thickness = pcb['copper_thickness'] / 1000 * pcb['layers']  # Total copper thickness (m)
        model_mode = pcb.get('model_mode', 'basic').lower()  # 'basic' or 'advanced'

        # Extract pcb_corners for pixel-to-mm scale factor
        pcb_corners = pcb.get('pcb_corners', [[0, 0], [width, height]])  # Fallback if no corners
        x_min, y_min = pcb_corners[0]
        x_max, y_max = pcb_corners[1]
        pcb_width_px = abs(x_max - x_min)
        scale_factor = width / pcb_width_px if pcb_width_px > 0 else 1.0  # mm per pixel

        # Get all elements and filter heaters
        elements = json_config.get('elements', [])
        heaters = [e for e in elements if e['type'] == 'Heater']

        # Initialize power and heater area
        total_power = 0.0
        heater_area_m2 = 0.0

        # Process each heater: sum power and calculate area
        for heater in heaters:
            total_power += heater.get('power', 0.0)  # Add heater power (W)

            # Calculate polygon area using Shoelace formula if points exist
            if 'points' in heater and len(heater['points']) > 2:
                x_coords = [(p[0] - x_min) * scale_factor for p in heater['points']]  # Convert to mm
                y_coords = [(p[1] - y_min) * scale_factor for p in heater['points']]  # Convert to mm

                # Shoelace formula for polygon area in mm²
                area_mm2 = 0.5 * abs(sum(
                    x_coords[i] * y_coords[(i + 1) % len(x_coords)] -
                    y_coords[i] * x_coords[(i + 1) % len(x_coords)]
                    for i in range(len(x_coords))
                ))
                heater_area_m2 += area_mm2 * 1e-6  # Convert mm² to m²

        # Calculate effective area based on model mode
        pcb_area_m2 = (width / 1000) * (height / 1000)  # Total PCB area in m²
        effective_area = pcb_area_m2 if model_mode == 'basic' else (heater_area_m2 or pcb_area_m2)

        # Log effective area used
        logger.debug(f"Analytic: Using effective area = {effective_area:.6f} m² (mode: {model_mode})")

        # Thermal resistance: convection (R_th_conv = 1 / (h * A))
        R_th_conv = 1 / (h_global * effective_area)  # K/W

        # Thermal resistance: conduction (R_th_cond = L / (k * A * t))
        L_avg = min(width, height) / 2000  # Approximate conduction path (m)
        R_th_cond = L_avg / (k_copper * effective_area * thickness)  # K/W

        # Total thermal resistance (series)
        R_th_total = R_th_conv + R_th_cond

        # Max temperature: T_max = T_amb + P * R_th_total
        max_temp = T_amb + total_power * R_th_total

        # Average temperature: halfway between ambient and max
        avg_temp = T_amb + (total_power * R_th_total) / 2

        # Delta T: temperature rise
        delta_t = max_temp - T_amb

        # Debug print for console
        print(f"Estimated max temperature using analytic model: {max_temp:.2f} °C")

        # Log all results
        logger.info(f"Analytic Simulation Results: Max Temperature = {max_temp:.2f} °C")
        logger.info(f"Analytic Simulation Results: Average Temperature = {avg_temp:.2f} °C")
        logger.info(f"Analytic Simulation Results: Delta T = {delta_t:.2f} °C")
        logger.info(f"Analytic Simulation Results: Effective R_th = {R_th_total:.2f} K/W")
        logger.info(f"Analytic Simulation Results: Total Power = {total_power:.2f} W")
        logger.info(f"Analytic Simulation Results: Effective Area = {effective_area:.6f} m²")
        logger.info(f"Analytic Simulation Results: Note = Analytic approximation (mode: {model_mode})")

        # Return results in standard format
        results = {
            "max_temp": {"value": max_temp, "unit": "°C"},
            "avg_temp": {"value": avg_temp, "unit": "°C"},
            "delta_t": {"value": delta_t, "unit": "°C"},
            "effective_r_th": {"value": R_th_total, "unit": "K/W"},
            "total_power": {"value": total_power, "unit": "W"},
            "effective_area": {"value": effective_area, "unit": "m²"},
            "note": {"value": f"Analytic approximation - limited to simple PCBs (mode: {model_mode})", "unit": ""}
        }
        return results


class AnalyticResultsDialog(wx.Dialog):
    """
    @brief Dialog for displaying Analytic simulation results.
    @details Two-column layout: left chart, right text. Includes Save button.
             Uses config for styles. Buttons aligned right.
    """
    def __init__(self, parent, results, logger=None):
        """
        @brief Initialize the Analytic results dialog with hidden Save button.
        @details Creates the dialog with loading message, hides the Save button until results are available.
                 The Save button will be shown after simulation completes.
        @param parent Parent window.
        @param results Initial results (empty during loading).
        @param logger Optional logger.
        """
        super().__init__(parent, title="Analytic Simulation Results", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        # Set window size and center it on screen
        self.SetMinSize(config.DEFAULT_RESULTS_WINDOW_SIZE)
        self.SetSize(config.DEFAULT_RESULTS_WINDOW_SIZE)
        self.Center()

        self.results = results

        # Use passed logger or create new one
        self.logger = logger
        if logger is None:
            self.logger = logging.getLogger(__name__)

        # Bind close event to exit the application
        self.Bind(wx.EVT_CLOSE, self.on_close_window)

        # Main vertical sizer for the entire dialog
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Content sizer for chart and text (will be replaced after simulation)
        self.content_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.main_sizer.Add(self.content_sizer, 1, wx.ALL | wx.EXPAND, 5)

        # Button panel
        self.btn_panel = wx.Panel(self)
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.AddStretchSpacer()

        # Save button - hidden until results are ready
        self.save_btn = wx.Button(self.btn_panel, label="Save")
        self.save_btn.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_btn.Hide()  # Hide Save button during loading

        # Close button
        close_btn = wx.Button(self.btn_panel, label="Close")
        close_btn.Bind(wx.EVT_BUTTON, self.on_close)

        btn_sizer.Add(self.save_btn, 0, wx.ALL, 5)
        btn_sizer.Add(close_btn, 0, wx.ALL, 5)
        self.btn_panel.SetSizer(btn_sizer)
        self.main_sizer.Add(self.btn_panel, 0, wx.ALL | wx.ALIGN_RIGHT, 5)

        self.SetSizer(self.main_sizer)
        self.Layout()
    
    def show_loading(self):
        """
        @brief Show full-screen loading message (no tabs, no Save button).
        @details Clears the content sizer and displays a centered loading text.
        """
        # Clear any previous content in the content sizer
        self.content_sizer.Clear(True)

        # Create a new panel for the loading message
        loading_panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)

        # Create loading text label
        loading_text = wx.StaticText(
            loading_panel,
            label="Simulation in progress...\nPlease wait."
        )
        # Set font: size 12, Swiss family, normal style, normal weight (not bold)
        loading_text.SetFont(wx.Font(
            12,
            wx.FONTFAMILY_SWISS,
            wx.FONTSTYLE_NORMAL,
            wx.FONTWEIGHT_NORMAL
        ))

        # Add stretch spacers to center the text vertically
        sizer.AddStretchSpacer()
        sizer.Add(loading_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
        sizer.AddStretchSpacer()

        # Apply sizer to the loading panel
        loading_panel.SetSizer(sizer)

        # Add loading panel to the content sizer (full expand)
        self.content_sizer.Add(loading_panel, 1, wx.EXPAND)

        # Update layout to show the loading message immediately
        self.Layout()

    def _create_content_sizer(self):
        """
        @brief Build two-column content: chart (left) + text (right).
        @details Saves chart to config.DEFAULT_TEMP_DIR as 'analytic_summary_chart.png'.
                 Saves text to config.DEFAULT_TEMP_DIR as 'analytic_summary_report.txt'.
                 Files are copied on Save and deleted on close.
        @return wx.BoxSizer with full content.
        """
        content = wx.BoxSizer(wx.HORIZONTAL)

        # ------------------- LEFT: CHART -------------------
        left_panel = wx.Panel(self)
        chart_sizer = wx.BoxSizer(wx.VERTICAL)

        # Set chart font size
        plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})

        # Create figure
        fig, ax = plt.subplots(figsize=(4, 3))

        # Match dialog background
        bg = self.GetBackgroundColour()
        fig.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))
        ax.set_facecolor('white')

        # Add grids
        ax.grid(True, which='major', linestyle='--', alpha=0.7, color='gray', zorder=0)
        ax.grid(True, which='minor', linestyle=':', alpha=0.3, color='gray', zorder=0)
        ax.minorticks_on()

        # Prepare bar data
        vals = [
            self.results.get('max_temp', {}).get('value', 0.0),
            self.results.get('avg_temp', {}).get('value', 0.0),
            self.results.get('delta_t', {}).get('value', 0.0)
        ]
        labs = ['Max Temp', 'Avg Temp', 'Delta T']

        # Plot bars
        bars = ax.bar(labs, vals, color=config.DEFAULT_CHART_COLOR_SCHEME[0], zorder=3)

        # Rotate x-axis labels 90 degrees for better fit and coherence with other solvers
        ax.set_xticklabels(labs, rotation=90)

        # Y-axis margin from config
        max_val = max(vals)
        margin_percent = config.DEFAULT_CHART_Y_MARGIN_PERCENT / 100.0
        if margin_percent > 0:
            y_margin = max_val * margin_percent
            ax.set_ylim(0, max_val + y_margin)
            label_offset = y_margin * 0.3
        else:
            label_offset = max_val * 0.05

        # Title and label
        ax.set_title("Analytic Simulation Metrics")
        ax.set_ylabel("Temperature (°C)")

        # Add value labels on bars
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height() + label_offset,
                    f'{v:.2f}', ha='center', va='bottom')

        # Tight layout
        plt.tight_layout()

        # SAVE CHART TO TEMP DIRECTORY WITH FIXED NAME
        import os
        os.makedirs(config.DEFAULT_TEMP_DIR, exist_ok=True)
        chart_path = os.path.join(config.DEFAULT_TEMP_DIR, "analytic_summary_chart.png")
        fig.savefig(chart_path, format='png', dpi=100, bbox_inches='tight')
        plt.close(fig)

        # Load and display chart
        img = wx.Image(chart_path, wx.BITMAP_TYPE_PNG)
        bmp = wx.Bitmap(img) if img.IsOk() else wx.Bitmap(1, 1)
        chart_bitmap = wx.StaticBitmap(left_panel, bitmap=bmp)
        chart_bitmap.Refresh()
        chart_sizer.Add(chart_bitmap, 0, wx.ALL | wx.EXPAND, 5)
        left_panel.SetSizer(chart_sizer)
        content.Add(left_panel, 1, wx.ALL | wx.EXPAND, 5)

        # ------------------- RIGHT: TEXT -------------------
        right_panel = wx.Panel(self)
        txt_sizer = wx.BoxSizer(wx.VERTICAL)

        # Build text content
        text_content = (
            f"Max Temperature: {self.results.get('max_temp', {}).get('value', 0.0):.2f} {self.results.get('max_temp', {}).get('unit', '°C')}\n"
            f"Average Temperature: {self.results.get('avg_temp', {}).get('value', 0.0):.2f} {self.results.get('avg_temp', {}).get('unit', '°C')}\n"
            f"Delta T: {self.results.get('delta_t', {}).get('value', 0.0):.2f} {self.results.get('delta_t', {}).get('unit', '°C')}\n"
            f"Effective R_th: {self.results.get('effective_r_th', {}).get('value', 0.0):.2f} {self.results.get('effective_r_th', {}).get('unit', 'K/W')}\n"
            f"Total Power: {self.results.get('total_power', {}).get('value', 0.0):.2f} {self.results.get('total_power', {}).get('unit', 'W')}\n"
            f"Effective Area: {self.results.get('effective_area', {}).get('value', 0.0):.6f} {self.results.get('effective_area', {}).get('unit', 'm²')}\n\n"
            f"{config.EXPLANATION_TEXT['Analytic']}"
        )

        # SAVE TEXT TO TEMP DIRECTORY WITH FIXED NAME
        txt_path = os.path.join(config.DEFAULT_TEMP_DIR, "analytic_summary_report.txt")
        with open(txt_path, 'w', encoding='utf-8') as f:
            f.write(text_content)

        # Create text control
        ctrl = wx.TextCtrl(right_panel, value=text_content,
                          style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_BESTWRAP)
        ctrl.SetFont(wx.Font(config.DEFAULT_RESULTS_FONT_SIZE, wx.FONTFAMILY_SWISS,
                             wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        txt_sizer.Add(ctrl, 1, wx.ALL | wx.EXPAND, 5)
        right_panel.SetSizer(txt_sizer)
        content.Add(right_panel, 1, wx.ALL | wx.EXPAND, 5)

        return content

    def _delete_temp_dir(self):
        """
        @brief Delete the temp directory.
        @details Called after simulation to clean up.
                Removes all generated temporary files (charts, reports, GIFs, heatmaps).
                Uses shutil.rmtree for recursive deletion.
                Logs success or warning on failure.
        """
        try:
            # Check if the temporary directory exists
            if os.path.exists(config.DEFAULT_TEMP_DIR):
                # Recursively delete the entire directory and its contents
                shutil.rmtree(config.DEFAULT_TEMP_DIR)
                # Log successful cleanup
                self.logger.debug(f"Deleted temp directory: {config.DEFAULT_TEMP_DIR}")
        except Exception as e:
            # Log warning if deletion fails (e.g., file in use, permission denied)
            self.logger.warning(f"Failed to delete temp directory: {e}")

    def on_close(self, event):
        """Handle the Close button click."""
        # Delete temporary folder (commented out to persist temp dir until new simulation or app close)
        # self._delete_temp_dir()        
        # Destroy dialog first
        self.Destroy()
        # Exit app
        # wx.GetApp().ExitMainLoop()

    def on_close_window(self, event):
        """
        @brief Handle window close event.
        @details Deletes the entire temp directory (config.DEFAULT_TEMP_DIR).
                 Uses wx.CallAfter to avoid file-in-use errors.
        """
        # Delete temporary folder (commented out to persist temp dir until new simulation or app close)
        # self._delete_temp_dir() 
        # Destroy dialog first
        self.Destroy()
        # Exit app
        # wx.GetApp().ExitMainLoop()

    def on_save(self, event):
        """
        @brief Handle the Save button click.
        @details Copies all files from config.DEFAULT_TEMP_DIR to user-selected folder.
                Includes PNGs, GIF, TXT.
                Preserves file metadata (creation/modification time) using shutil.copy2.
                Shows message with number of files saved or error.
        """
        # Create directory selection dialog
        with wx.DirDialog(self, "Select save directory", style=wx.DD_DEFAULT_STYLE) as dlg:
            # Show dialog and wait for user selection
            if dlg.ShowModal() == wx.ID_CANCEL:
                # Exit function if user cancels
                return
            # Get selected destination directory
            dest_dir = dlg.GetPath()

        try:
            # Check if temporary directory exists
            if not os.path.exists(config.DEFAULT_TEMP_DIR):
                # Show warning if no files to save
                wx.MessageBox("No temporary files to save.", "Save", wx.OK | wx.ICON_WARNING)
                # Exit function
                return

            # Initialize counter for copied files
            copied = 0

            # Iterate through all files in temporary directory
            for item in os.listdir(config.DEFAULT_TEMP_DIR):
                # Build full source path
                src = os.path.join(config.DEFAULT_TEMP_DIR, item)
                # Build full destination path
                dst = os.path.join(dest_dir, item)
                # Copy file with metadata (creation/modification time)
                shutil.copy2(src, dst)
                # Increment counter
                copied += 1

            # Show success message with file count and destination
            wx.MessageBox(f"Saved {copied} file(s) to:\n{dest_dir}", "Save", wx.OK | wx.ICON_INFORMATION)

        except Exception as e:
            # Log error if copy operation fails
            self.logger.error(f"Save failed: {e}")
            # Show error dialog to user
            wx.MessageBox(f"Save failed: {e}", "Error", wx.OK | wx.ICON_ERROR)

