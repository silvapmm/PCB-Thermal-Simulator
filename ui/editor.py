# Module for the graphical editor using wxPython
# This module handles the UI for drawing PCB elements, editing parameters in tabs, and triggering simulations
import wx # Import wxPython for GUI components
import wx.lib.scrolledpanel as scrolled # Import for scrollable panel to handle canvas
import wx.grid # Import for grid table in Element Parameters tab (not used now, but kept for reference)
import math # Import for distance calculations in Delete Selected
import os # Import OS for path manipulation or system interaction
import json # Import for handling JSON data (e.g., loading/saving configuration or element data)
import shutil # Import for high-level file operations (e.g., copying files or directories)
import sys # Import for system-specific parameters and functions (e.g., accessing command line arguments or exiting the program)
import copy # Import for creating shallow and deep copies of objects
import logging # Import for logging events, status, and errors
from datetime import datetime # Import datetime functions
import config # Import config module from root directory (absolute import)
from core.solver import ThermalSolver # Import central solver module
from core.analytic_solver import AnalyticResultsDialog # For showing last analytic results
from core.lumped_solver import LumpedResultsDialog # For showing last analytic results
from core.fdm_solver import FDMResultsDialog # For showing last FDM results


class PCBEditor(wx.Frame):
    """
    @brief Graphical editor for PCB thermal simulation using wxPython.
    @details Allows drawing elements on a canvas with overlay image, editing parameters in tabs, and running simulation.
    """
    def __init__(self, parent):
        """
        @brief Initialize the editor with default values from config.py.
        @param parent Parent window (None for main window).
        """
       # Configure logging to write to a file in the parent directory and console
        log_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', config.LOG_FILENAME))
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        # Clear any existing handlers to avoid duplicates
        logger.handlers = []
        # File handler to write to log file, clearing it on each run
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
        # Console handler to output to console
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(console_handler)
        self.logger = logger

        # Initialize the wxPython frame with title from config.py and updated size
        wx.Frame.__init__(self, parent, title=f"{config.PRODUCT_NAME} V{config.APP_VERSION}", size=config.DEFAULT_MAIN_WINDOW_SIZE)

        # Set window icon using PyInstaller-compatible path
        icon_path = os.path.normpath(os.path.join(self.get_icon_base_path(), 'ico', 'pcb-board_128x128.ico'))
        if os.path.exists(icon_path):
            icon = wx.Icon(icon_path, wx.BITMAP_TYPE_ICO)
            if icon.IsOk():
                self.SetIcon(icon)
                self.logger.debug(f"Set window icon: {icon_path}")
            else:
                self.logger.warning(f"Failed to load icon: {icon_path}")
        else:
            self.logger.warning(f"Icon file not found: {icon_path}")

        # Initialize lists for drawing elements
        self.current_element = []  # List to store points of the current element being drawn
        self.elements = []  # List to store all completed elements (polygons and vias)
        # Current tool mode
        self.current_tool = None  # No default tool selected
        # Counters for element IDs
        self.element_counters = {"Zone": 0, "Heater": 0, "Via": 0, "Sink": 0}
        # Zoom and pan variables
        self.zoom_factor = 1.0  # Default zoom level (100%)
        self.pan_x = 0  # X offset for panning
        self.pan_y = 0  # Y offset for panning
        self.is_panning = False  # Flag to track panning state
        self.pan_start_x = 0  # Starting x position for panning
        self.pan_start_y = 0  # Starting y position for panning
        # Grid variables
        self.show_grid = False  # Grid visibility state
        self.grid_spacing = 20  # Default grid spacing in pixels
        self.show_solver_grid = False  # Solver grid visibility state
        # List to store PCB corner points for setting size
        self.pcb_corners = []  # Stores two points for PCB size
        # Scale factor for converting pixels to mm
        self.scale_factor = 1.0  # Default: 1 pixel = 1 mm (updated by Set PCB Size)
        # Element being duplicated
        self.duplicate_element = None  # Stores the element being duplicated
        # Current mouse position for preview
        self.preview_position = None  # Stores (x, y) for preview during duplication/moving
        # Element being moved
        self.moving_element = None  # Stores the element being moved
        # Starting position for moving
        self.move_start_pos = None  # Stores (x, y) of the initial click for moving
        # Element being edited
        self.editing_element = None  # Stores the element being edited
        # Index of the vertex being edited
        self.editing_vertex_index = None  # Stores the index of the vertex being moved
        # Original element for restoring on abort
        self.original_element = None  # Stores a copy of the element before editing/moving
        # Initialize PCB parameters with defaults from config.py
        self.pcb_params = {
            "width": config.DEFAULT_PCB_WIDTH,
            "height": config.DEFAULT_PCB_HEIGHT,
            "thickness": config.DEFAULT_PCB_THICKNESS,
            "layers": config.DEFAULT_PCB_LAYERS,
            "copper_thickness": config.DEFAULT_COPPER_THICKNESS_PER_LAYER,
            "copper_conductivity": config.DEFAULT_COPPER_CONDUCTIVITY,
            "fr4_conductivity": config.DEFAULT_FR4_CONDUCTIVITY,
            "ambient_temp_global": config.DEFAULT_AMBIENT_TEMP_GLOBAL,
            "convection_h_global": config.DEFAULT_CONVECTION_H_GLOBAL,
            "convection_h_top": config.DEFAULT_CONVECTION_H_TOP,
            "convection_h_bottom": config.DEFAULT_CONVECTION_H_BOTTOM,
            "grid_dx": config.DEFAULT_SOLVER_GRID_SIZE_MM,
            "grid_dy": config.DEFAULT_SOLVER_GRID_SIZE_MM,
            "tolerance": config.DEFAULT_FDM_TOLERANCE_CONVERGENCE,
            "max_iterations": config.DEFAULT_FDM_MAX_ITERATIONS,
            "over_relaxation_factor": config.DEFAULT_FDM_OVER_RELAXATION_FACTOR,
            "model_mode": config.DEFAULT_SOLVER_MODEL_MODE,
            "boundary_mode": config.DEFAULT_FDM_BOUNDARY_MODE,
            "hotspot_threshold": config.DEFAULT_FDM_HOTSPOT_THRESHOLD
        }
        # Dictionary to store global parameter entry widgets
        self.global_entries = {}  # Initialize global_entries
        # Set solver mode from config.py
        self.solver_mode = config.DEFAULT_SOLVER_METHOD
        # Initialize last simulation results
        self.last_results = None  # Store results of the last simulation
        self.last_solver_mode = None  # Store solver mode of the last simulation
        # Overlay variables
        self.overlay_image = None  # Stores the loaded overlay image (wx.Image)
        self.overlay_path = None  # Stores the original path of the loaded image
        self.overlay_visible = False  # Controls overlay visibility
        # Set up the user interface components
        self.setup_ui()

    def get_icon_base_path(self):
        """
        @brief Returns the absolute directory where the main.py script is located.
        @details For a Python script, it uses the directory of main.py. For a PyInstaller executable,
                it uses sys._MEIPASS to access the temporary resource directory.
        @return String representing the absolute path to the root directory containing main.py.
        """
        # Check if the code is running inside a PyInstaller bundle
        if getattr(sys, 'frozen', False):
            # When frozen, use sys._MEIPASS for bundled resources
            return getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
        else:
            # When running as a normal script, use the directory of main.py
            return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    def setup_ui(self):
        """
        @brief Setup the wxPython UI components: notebook with Drawing Area and Element Parameters tabs, surrounded by reserved areas.
        @details Creates panels, buttons, sliders, choices, and notebooks for tools, properties, and drawing.
                Includes updated element table with improved column names and new columns for coordinates/vertices.
                Adjusted column widths for better layout (Coordinates narrowed to 150px).
                Added bottom footer panel for copyright text from config.LEGAL_COPYRIGHT, centered horizontally.
                Added Copper Conductivity field (before FR4) with same decimals/increment as FR4.
                Corrected FR4 Thermal Conductivity unit to (W/mK).
                Updated Layer Filter dropdown to include dynamic options based on number of layers.
                Reordered Solver Parameters tab to Mode, Max Iterations, Solver Tolerance, Over Relaxation Factor, Grid X, Grid Y.
                Binds events for automatic parameter updates and removes Update Properties button.
                Disables Show/Hide Solver Grid button for non-FDM modes.
                Binds Simulation Results button to on_show_results.
                Removes Lock Overlay button and adds New button to reset to default state.
                Binds EVT_TEXT_ENTER to SpinCtrlDouble for manual value confirmation with Enter key.
                Ensures all FDM-specific fields (including Max Iterations) are enabled/disabled based on initial solver_mode.
        """
        # Create the main panel to hold all UI elements
        panel = wx.Panel(self)
        # Create a vertical box sizer for the main layout
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        # Create top reserved area with controls
        top_panel = wx.Panel(panel)
        top_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Add left-aligned controls
        btn_new = wx.Button(top_panel, label="New")
        btn_new.Bind(wx.EVT_BUTTON, self.new_project)
        top_sizer.Add(btn_new, 0, wx.ALL, 5)
        btn_load = wx.Button(top_panel, label="Load Overlay")
        btn_load.Bind(wx.EVT_BUTTON, self.load_overlay)
        top_sizer.Add(btn_load, 0, wx.ALL, 5)
        btn_toggle = wx.Button(top_panel, label="Show/Hide Overlay")
        btn_toggle.Bind(wx.EVT_BUTTON, self.toggle_overlay)
        top_sizer.Add(btn_toggle, 0, wx.ALL, 5)
        opacity_label = wx.StaticText(top_panel, label="Opacity:")
        top_sizer.Add(opacity_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.opacity_slider = wx.Slider(top_panel, value=100, minValue=0, maxValue=100, style=wx.SL_HORIZONTAL, size=(150, -1))
        self.opacity_slider.Bind(wx.EVT_SLIDER, self.update_opacity)
        top_sizer.Add(self.opacity_slider, 0, wx.ALL, 5)
        layer_label = wx.StaticText(top_panel, label="Layer Filter:")
        top_sizer.Add(layer_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        # Dynamic layer filter choices based on number of layers
        max_layers = int(self.pcb_params['layers'])
        layer_choices = ["All", f"1+{max_layers}"]
        if max_layers > 2:
            layer_choices.append("Inner")  # Renamed from Inner Layers for brevity
        layer_choices.extend([str(i) for i in range(1, max_layers + 1)])  # Add individual layers
        self.layer_choice = wx.Choice(top_panel, choices=layer_choices)
        self.layer_choice.SetSelection(0)  # Default to "All"
        self.layer_choice.Bind(wx.EVT_CHOICE, self.on_layer_choice)
        top_sizer.Add(self.layer_choice, 0, wx.ALL, 5)
        # Add grid controls
        btn_grid = wx.Button(top_panel, label="Show/Hide Grid")
        btn_grid.Bind(wx.EVT_BUTTON, self.on_toggle_grid)
        top_sizer.Add(btn_grid, 0, wx.ALL, 5)
        grid_spacing_label = wx.StaticText(top_panel, label="Grid Spacing:")
        top_sizer.Add(grid_spacing_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.grid_spacing_choice = wx.Choice(top_panel, choices=["10", "20", "50", "100"])
        self.grid_spacing_choice.SetSelection(1)  # Default to 20 pixels
        self.grid_spacing_choice.Bind(wx.EVT_CHOICE, self.on_grid_spacing_change)
        top_sizer.Add(self.grid_spacing_choice, 0, wx.ALL, 5)
        # Add stretch spacer to push remaining buttons to the right
        top_sizer.AddStretchSpacer()
        # Add right-aligned buttons
        btn_import = wx.Button(top_panel, label="Import JSON")
        btn_import.Bind(wx.EVT_BUTTON, self.import_json)
        top_sizer.Add(btn_import, 0, wx.ALL, 5)
        btn_export = wx.Button(top_panel, label="Export JSON")
        btn_export.Bind(wx.EVT_BUTTON, self.export_json)
        top_sizer.Add(btn_export, 0, wx.ALL, 5)
        top_panel.SetSizer(top_sizer)
        main_sizer.Add(top_panel, 0, wx.EXPAND | wx.ALL, 5)
        # Create horizontal sizer for left, notebook, and right areas
        horizontal_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Create left reserved area with tool buttons
        left_panel = wx.Panel(panel)
        left_sizer = wx.BoxSizer(wx.VERTICAL)
        # Add label for the tools section
        left_label = wx.StaticText(left_panel, label="Tools")
        left_sizer.Add(left_label, 0, wx.ALL | wx.ALIGN_LEFT, 5)
        # List of tool button labels
        tools = ["Set PCB Size", "Add Zone", "Add Heater", "Add Via", "Add Sink", "Move Selected", "Duplicate Selected", "Edit Selected", "Delete Selected", "Draw/Update"]
        # Add buttons to the tool panel with event bindings
        for tool in tools:
            btn = wx.Button(left_panel, label=tool)
            btn.Bind(wx.EVT_BUTTON, lambda e, t=tool: self.set_tool_and_show(t))
            left_sizer.Add(btn, 0, wx.ALL | wx.EXPAND, 5)
        left_panel.SetSizer(left_sizer)
        horizontal_sizer.Add(left_panel, 0, wx.EXPAND | wx.ALL, 5)
        # Create notebook for Drawing Area and Element Parameters
        notebook = wx.Notebook(panel)
        # Create Drawing Area tab
        drawing_panel = wx.Panel(notebook)
        drawing_sizer = wx.BoxSizer(wx.VERTICAL)
        # Create a horizontal sizer for zoom controls
        controls_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Add zoom label and buttons
        self.zoom_label = wx.StaticText(drawing_panel, label="Zoom: 100%")
        controls_sizer.Add(self.zoom_label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        zoom_in_button = wx.Button(drawing_panel, label="Zoom +")
        zoom_out_button = wx.Button(drawing_panel, label="Zoom -")
        reset_zoom_button = wx.Button(drawing_panel, label="100%")
        controls_sizer.Add(zoom_in_button, 0, wx.ALL, 5)
        controls_sizer.Add(zoom_out_button, 0, wx.ALL, 5)
        controls_sizer.Add(reset_zoom_button, 0, wx.ALL, 5)
        zoom_in_button.Bind(wx.EVT_BUTTON, self.on_zoom_in)
        zoom_out_button.Bind(wx.EVT_BUTTON, self.on_zoom_out)
        reset_zoom_button.Bind(wx.EVT_BUTTON, self.on_reset_zoom)
        drawing_sizer.Add(controls_sizer, 0, wx.EXPAND | wx.ALL, 5)
        # Create a scrolled panel for the drawing canvas
        self.canvas = scrolled.ScrolledPanel(drawing_panel, size=(500, 400))
        self.canvas.SetScrollbars(1, 1, 500, 400)
        # Bind events for painting, mouse clicks, panning, scrolling, and key presses
        self.canvas.Bind(wx.EVT_PAINT, self.on_paint)
        self.canvas.Bind(wx.EVT_LEFT_DOWN, self.on_canvas_click)
        self.canvas.Bind(wx.EVT_LEFT_DCLICK, self.on_canvas_double_click)
        self.canvas.Bind(wx.EVT_LEFT_UP, self.on_mouse_left_up)
        self.canvas.Bind(wx.EVT_MIDDLE_DOWN, self.on_middle_down)
        self.canvas.Bind(wx.EVT_MOTION, self.on_mouse_motion)
        self.canvas.Bind(wx.EVT_MIDDLE_UP, self.on_middle_up)
        self.canvas.Bind(wx.EVT_MOUSEWHEEL, self.on_mouse_wheel)
        self.canvas.Bind(wx.EVT_KEY_DOWN, self.on_key_down)
        self.Bind(wx.EVT_KEY_DOWN, self.on_key_down)
        drawing_sizer.Add(self.canvas, 1, wx.EXPAND | wx.ALL, 10)
        # Create a title label for instructions, placed below the canvas
        self.title_label = wx.StaticText(drawing_panel, label="Select a tool from the left panel to start drawing")
        drawing_sizer.Add(self.title_label, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        drawing_panel.SetSizer(drawing_sizer)
        notebook.AddPage(drawing_panel, "Drawing Area")
        # Create Element Parameters tab
        params_panel = wx.Panel(notebook)
        params_sizer = wx.BoxSizer(wx.VERTICAL)
        # Use wx.ListCtrl for a table-like structure
        self.element_table = wx.ListCtrl(params_panel, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        self.element_table.InsertColumn(0, "Element Type")
        self.element_table.InsertColumn(1, "Element ID")
        self.element_table.InsertColumn(2, "Coordinates")
        self.element_table.InsertColumn(3, "Num Vertices")
        self.element_table.InsertColumn(4, "Layer")
        self.element_table.InsertColumn(5, "Power (W)")
        self.element_table.InsertColumn(6, "Diameter (mm)")
        self.element_table.InsertColumn(7, "HTC Extra (W/m²K)")
        # Set column widths for better layout (Coordinates narrowed to 150px)
        self.element_table.SetColumnWidth(0, 100)   # Element Type
        self.element_table.SetColumnWidth(1, 80)    # Element ID
        self.element_table.SetColumnWidth(2, 150)   # Coordinates (narrowed for better space usage)
        self.element_table.SetColumnWidth(3, 100)   # Num Vertices
        self.element_table.SetColumnWidth(4, 60)    # Layer
        self.element_table.SetColumnWidth(5, 100)   # Power (W)
        self.element_table.SetColumnWidth(6, 100)   # Diameter (mm)
        self.element_table.SetColumnWidth(7, 120)   # HTC Extra (W/m²K)
        # Bind event for double-click editing
        self.element_table.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_list_item_activated)
        params_sizer.Add(self.element_table, 1, wx.EXPAND | wx.ALL, 10)
        params_panel.SetSizer(params_sizer)
        notebook.AddPage(params_panel, "Element Parameters")
        horizontal_sizer.Add(notebook, 1, wx.EXPAND | wx.ALL, 5)
        # Create right reserved area with property tabs
        prop_panel = wx.Panel(panel)
        prop_sizer = wx.BoxSizer(wx.VERTICAL)
        # Add label for the properties section
        prop_label = wx.StaticText(prop_panel, label="Properties")
        prop_sizer.Add(prop_label, 0, wx.ALL | wx.ALIGN_LEFT, 5)
        # Create notebook for property tabs
        prop_notebook = wx.Notebook(prop_panel)
        # Set minimum size for notebook to ensure visibility of long labels
        prop_notebook.SetMinSize((250, -1))
        # Define parameter labels for input fields
        param_labels = {
            "mode": "Solver Mode",
            "width": "PCB Width (mm)",
            "height": "PCB Height (mm)",
            "thickness": "PCB Thickness (mm)",
            "layers": "Number of Layers",
            "copper_thickness": "Copper Thickness per Layer (mm)",
            "copper_conductivity": "Copper Conductivity (W/mK)",
            "fr4_conductivity": "FR4 Thermal Conductivity (W/mK)",
            "ambient_temp_global": "Ambient Temperature (°C)",
            "convection_h_global": "Global Convection Coefficient (W/m²K)",
            "convection_h_top": "Top Convection Coefficient (W/m²K)",
            "convection_h_bottom": "Bottom Convection Coefficient (W/m²K)",
            "grid_dx": "Grid Resolution X (mm)",
            "grid_dy": "Grid Resolution Y (mm)",
            "tolerance": "Solver Tolerance (°C)",
            "max_iterations": "Max Iterations",
            "over_relaxation_factor": "Over Relaxation Factor",
            "model_mode": "Model Mode",  # New: Label for Model Mode dropdown
            "hotspot_threshold": "Hotspot Threshold (°C)"  # New: Label for Hotspot Threshold field
        }
        self.param_labels = param_labels  # Store param_labels for use in update_global_params
        # Create the PCB Parameters tab
        pcb_panel = wx.Panel(prop_notebook)
        pcb_sizer = wx.BoxSizer(wx.VERTICAL)
        pcb_entries = ['width', 'height', 'thickness', 'layers', 'copper_thickness', 'copper_conductivity', 'fr4_conductivity']
        for key in pcb_entries:
            label = wx.StaticText(pcb_panel, label=param_labels[key])
            pcb_sizer.Add(label, 0, wx.ALL, 5)
            if key == 'layers':
                entry = wx.Choice(pcb_panel, choices=[str(i) for i in range(1, config.DEFAULT_PCB_MAX_LAYERS + 1)])  # Options 1 to DEFAULT_PCB_MAX_LAYERS
                entry.SetSelection(self.pcb_params[key] - 1)  # Set default value
                entry.Bind(wx.EVT_CHOICE, self.on_param_change)  # Bind event for automatic updates
            else:
                # Set decimal places: 2 for width, height, thickness; 3 for copper_thickness, copper_conductivity, fr4_conductivity
                digits = 2 if key in ['width', 'height', 'thickness'] else 3
                # Set increments: 0.01 for width, height, thickness; 0.005 for copper_thickness, copper_conductivity, fr4_conductivity
                inc = 0.01 if key in ['width', 'height', 'thickness'] else 0.005
                entry = wx.SpinCtrlDouble(pcb_panel, min=0.0, max=1000.0, initial=self.pcb_params[key], inc=inc)
                entry.SetDigits(digits)  # Set decimal places
                entry.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_param_change)  # Bind event for automatic updates
                entry.Bind(wx.EVT_TEXT_ENTER, self.on_param_change)  # Bind Enter key for manual input confirmation
            pcb_sizer.Add(entry, 0, wx.ALL | wx.EXPAND, 5)
            self.global_entries[key] = entry
        pcb_panel.SetSizer(pcb_sizer)
        prop_notebook.AddPage(pcb_panel, "PCB Parameters")
        # Create the Simulation Parameters tab
        sim_panel = wx.Panel(prop_notebook)
        sim_sizer = wx.BoxSizer(wx.VERTICAL)
        sim_entries = ['ambient_temp_global', 'convection_h_global', 'convection_h_top', 'convection_h_bottom']
        for key in sim_entries:
            label = wx.StaticText(sim_panel, label=param_labels[key])
            sim_sizer.Add(label, 0, wx.ALL, 5)
            entry = wx.SpinCtrlDouble(sim_panel, min=0.0, max=1000.0, initial=self.pcb_params[key], inc=0.1)
            entry.SetDigits(1)  # Set to 1 decimal place
            entry.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_param_change)  # Bind event for automatic updates
            entry.Bind(wx.EVT_TEXT_ENTER, self.on_param_change)  # Bind Enter key for manual input confirmation
            sim_sizer.Add(entry, 0, wx.ALL | wx.EXPAND, 5)
            self.global_entries[key] = entry
        sim_panel.SetSizer(sim_sizer)
        prop_notebook.AddPage(sim_panel, "Simulation Parameters")
        # Create the Solver Parameters tab
        solver_panel = wx.Panel(prop_notebook)
        solver_sizer = wx.BoxSizer(wx.VERTICAL)
        solver_entries = ['max_iterations', 'tolerance', 'over_relaxation_factor', 'grid_dx', 'grid_dy']  # Reordered entries
        # Create Mode dropdown first to ensure solver_mode_choice exists
        mode_label = wx.StaticText(solver_panel, label=param_labels['mode'])
        solver_sizer.Add(mode_label, 0, wx.ALL, 5)
        self.solver_mode_choice = wx.Choice(solver_panel, choices=["Analytic", "Lumped", "FDM"])
        self.solver_mode_choice.SetSelection(self.solver_mode_choice.GetStrings().index(self.solver_mode))
        self.solver_mode_choice.Bind(wx.EVT_CHOICE, self.on_solver_mode_change)  # Bind event for mode change
        solver_sizer.Add(self.solver_mode_choice, 0, wx.ALL | wx.EXPAND, 5)
        # New: Add Model Mode dropdown after Solver Mode, enabled only for Analytic/Lumped
        model_mode_label = wx.StaticText(solver_panel, label=param_labels['model_mode'])
        solver_sizer.Add(model_mode_label, 0, wx.ALL, 5)
        self.model_mode_choice = wx.Choice(solver_panel, choices=["Basic", "Advanced"])
        self.model_mode_choice.SetSelection(self.model_mode_choice.GetStrings().index(self.pcb_params['model_mode']))
        self.model_mode_choice.Bind(wx.EVT_CHOICE, self.on_param_change)  # Bind event for automatic updates
        self.model_mode_choice.Enable(self.solver_mode in ["Analytic", "Lumped"])  # Enable only for Analytic/Lumped
        solver_sizer.Add(self.model_mode_choice, 0, wx.ALL | wx.EXPAND, 5)
        self.global_entries['model_mode'] = self.model_mode_choice
        
        # New: Add Boundary Mode dropdown right below Model Mode, enabled only for FDM
        boundary_label = wx.StaticText(solver_panel, label="Boundary Mode:")
        solver_sizer.Add(boundary_label, 0, wx.ALL, 5)
        self.boundary_choice = wx.Choice(solver_panel, choices=["Convective", "Insulated"])
        self.boundary_choice.SetSelection(0)  # Default to Convective
        self.boundary_choice.Bind(wx.EVT_CHOICE, self.on_param_change)  # Bind event for automatic updates
        self.boundary_choice.Enable(self.solver_mode == "FDM")  # Enable only for FDM
        solver_sizer.Add(self.boundary_choice, 0, wx.ALL | wx.EXPAND, 5)
        self.global_entries['boundary_mode'] = self.boundary_choice        
        
        # Add solver parameters in specified order (all FDM-specific, so enable/disable based on mode)
        for key in solver_entries:
            label = wx.StaticText(solver_panel, label=param_labels[key])
            solver_sizer.Add(label, 0, wx.ALL, 5)
            if key == 'max_iterations':
                entry = wx.SpinCtrl(solver_panel, min=1, max=100000, initial=int(self.pcb_params[key]))  # Integer input for iterations
                entry.Bind(wx.EVT_SPINCTRL, self.on_param_change)  # Bind event for automatic updates
                entry.Enable(self.solver_mode == "FDM")  # Enable only for FDM mode (new: added for max_iterations)
            elif key == 'over_relaxation_factor':
                entry = wx.SpinCtrlDouble(solver_panel, min=0.0, max=2.0, initial=self.pcb_params[key], inc=0.1)
                entry.SetDigits(1)  # Set to 1 decimal place
                entry.Enable(self.solver_mode == "FDM")  # Enable only for FDM mode
                entry.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_param_change)  # Bind event for automatic updates
                entry.Bind(wx.EVT_TEXT_ENTER, self.on_param_change)  # Bind Enter key for manual input confirmation
            else:
                entry = wx.SpinCtrlDouble(solver_panel, min=0.0, max=1000.0, initial=self.pcb_params[key], inc=0.01)
                entry.SetDigits(2)  # Allow up to 2 decimal places
                entry.Enable(self.solver_mode == "FDM")  # Enable only for FDM mode
                entry.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_param_change)  # Bind event for automatic updates
                entry.Bind(wx.EVT_TEXT_ENTER, self.on_param_change)  # Bind Enter key for manual input confirmation
            solver_sizer.Add(entry, 0, wx.ALL | wx.EXPAND, 5)
            self.global_entries[key] = entry
        # New: Add Hotspot Threshold after Grid Resolution Y, enabled only for FDM
        hotspot_label = wx.StaticText(solver_panel, label=param_labels['hotspot_threshold'])
        solver_sizer.Add(hotspot_label, 0, wx.ALL, 5)
        hotspot_entry = wx.SpinCtrlDouble(solver_panel, min=0.0, max=500.0, initial=self.pcb_params['hotspot_threshold'], inc=1.0)
        hotspot_entry.SetDigits(1)  # Set to 1 decimal place
        hotspot_entry.Enable(self.solver_mode == "FDM")  # Enable only for FDM mode
        hotspot_entry.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_param_change)  # Bind event for automatic updates
        hotspot_entry.Bind(wx.EVT_TEXT_ENTER, self.on_param_change)  # Bind Enter key for manual input confirmation
        solver_sizer.Add(hotspot_entry, 0, wx.ALL | wx.EXPAND, 5)
        self.global_entries['hotspot_threshold'] = hotspot_entry
        # Add spacer before Show/Hide Solver Grid button
        solver_sizer.AddSpacer(10)
        # Add Show/Hide Solver Grid button, disabled for non-FDM modes
        self.solver_grid_btn = wx.Button(solver_panel, label="Show/Hide Solver Grid")
        self.solver_grid_btn.SetMinSize((150, 30))  # Match size of other buttons
        self.solver_grid_btn.Enable(self.solver_mode == "FDM")  # Enable only for FDM mode
        self.solver_grid_btn.Bind(wx.EVT_BUTTON, self.on_toggle_solver_grid)
        solver_sizer.Add(self.solver_grid_btn, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        solver_panel.SetSizer(solver_sizer)
        prop_notebook.AddPage(solver_panel, "Solver Parameters")
        # Add notebook to sizer with expand to fill space
        prop_sizer.Add(prop_notebook, 1, wx.EXPAND | wx.ALL, 5)
        # Add the Run Simulation button
        run_btn = wx.Button(prop_panel, label="Run Simulation")
        run_btn.SetMinSize((150, 30))  # Set fixed size for consistent appearance
        run_btn.Bind(wx.EVT_BUTTON, self.run_simulation)
        prop_sizer.Add(run_btn, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        # Add the Simulation Results button (disabled by default)
        self.results_btn = wx.Button(prop_panel, label="Simulation Results")
        self.results_btn.SetMinSize((150, 30))  # Set same size as Run Simulation
        self.results_btn.Enable(False)
        self.results_btn.Bind(wx.EVT_BUTTON, self.on_show_results)  # Bind to show last results
        prop_sizer.Add(self.results_btn, 0, wx.ALL | wx.ALIGN_CENTER, 5)
        prop_panel.SetSizer(prop_sizer)
        prop_sizer.Fit(prop_panel)
        horizontal_sizer.Add(prop_panel, 0, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(horizontal_sizer, 1, wx.EXPAND | wx.ALL, 5)
        # Create bottom footer panel for copyright text
        bottom_panel = wx.Panel(panel)
        bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Add stretch spacer to center the text
        bottom_sizer.AddStretchSpacer()
        # Add copyright text from config (fallback to generic if not defined)
        copyright_text = getattr(config, 'LEGAL_COPYRIGHT', "© 2025 PCB Editor")
        copyright_label = wx.StaticText(bottom_panel, label=copyright_text)
        # Set smaller font size for a more compact footer
        copyright_label.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        bottom_sizer.Add(copyright_label, 0, wx.ALL | wx.ALIGN_CENTER, 2)
        # Add another stretch spacer for centering
        bottom_sizer.AddStretchSpacer()
        bottom_panel.SetSizer(bottom_sizer)
        main_sizer.Add(bottom_panel, 0, wx.EXPAND | wx.ALL, 5)
        # Set the main sizer for the panel
        panel.SetSizer(main_sizer)
        # Bind the close event to on_close_main
        self.Bind(wx.EVT_CLOSE, self.on_close_main)
        # Display the frame
        self.Show()

    def set_tool_and_show(self, tool):
        """
        @brief Set the active tool for drawing elements and update the instruction label.
        @details Resets current element and updates UI label based on tool type.
                For "Draw/Update", triggers renumbering of elements to remove ID gaps.
                For Add Zone/Heater/Sink/Via, shows parameter dialogs before allowing drawing.
        @param tool The selected tool (e.g., "Set PCB Size", "Add Zone", "Add Heater", "Add Via", "Add Sink", "Move Selected", "Duplicate Selected", "Edit Selected", "Delete Selected", "Draw/Update").
        """
        self.current_tool = tool
        self.current_element = []  # Reset current element when changing tools
        
        # Update the instruction label based on the selected tool
        if tool == "Set PCB Size":
            self.title_label.SetLabel("Set PCB Size - Click two opposite corners to define PCB size")
            self.pcb_corners = []  # Reset PCB corners only for Set PCB Size
        elif tool == "Add Via":
            # Show Via parameters dialog first
            if not self.show_via_dialog():
                self.current_tool = None  # Reset tool if dialog is cancelled
                self.title_label.SetLabel("Select a tool from the left panel to start drawing")
            else:
                self.title_label.SetLabel("Drawing Via - Click to add point, middle-click to pan")
        elif tool in ["Add Zone", "Add Heater", "Add Sink"]:
            # Show appropriate dialog for Zone/Heater/Sink before drawing
            success = False
            if tool == "Add Zone":
                success = self.show_zone_dialog()
            elif tool == "Add Heater":
                success = self.show_heater_dialog()
            elif tool == "Add Sink":
                success = self.show_sink_dialog()
            if not success:
                self.current_tool = None  # Reset tool if dialog is cancelled
                self.title_label.SetLabel("Select a tool from the left panel to start drawing")
            else:
                tool_type = tool.replace("Add ", "")
                self.title_label.SetLabel(f"Drawing {tool_type} - Click to add points, left double-click to close, middle-click to pan")
        elif tool == "Move Selected":
            self.title_label.SetLabel("Move Selected - Click to select element, click to place, esc to cancel")
            self.moving_element = None  # Reset moving element
            self.move_start_pos = None
            self.original_element = None
            self.canvas.SetFocus()  # Ensure canvas has focus for key events
        elif tool == "Duplicate Selected":
            self.title_label.SetLabel("Duplicate Selected - Click to select element, click to place copy")
            self.duplicate_element = None  # Reset duplicate element
        elif tool == "Edit Selected":
            self.title_label.SetLabel("Edit Selected - Click to select polygon, drag vertices, 'a' to add vertex, 'd' to delete vertex, double-click or esc to finish")
            self.editing_element = None  # Reset editing element
            self.editing_vertex_index = None
            self.original_element = None
            self.canvas.SetFocus()  # Ensure canvas has focus for key events
        elif tool == "Delete Selected":
            self.title_label.SetLabel("Delete Selected - Click to select element, middle-click to pan")
        elif tool == "Draw/Update":
            # Renumber elements to remove gaps and update UI
            self.renumber_elements()
            self.update_element_table()  # Refresh the element table with new IDs/labels
            self.canvas.Refresh()  # Redraw canvas to show updated labels
            self.current_tool = None  # Reset tool immediately after action (one-shot behavior)
            self.title_label.SetLabel("Elements renumbered. Select a tool from the left panel to start drawing")
        
        self.logger.debug(f"Selected tool: {tool}")  # Debug log
        self.canvas.Refresh()

    def show_zone_dialog(self):
        """
        @brief Show popup dialog for Zone parameters (Layer dropdown based on number of layers).
        @details Prompts user to select a layer for the Zone element.
        @return True if dialog was successful, False otherwise.
        """
        # Get maximum number of PCB layers from parameters
        max_layers = int(self.pcb_params['layers'])
        # Create list of layer options from 1 to max_layers
        layers = [str(i) for i in range(1, max_layers + 1)]
        # Create single-choice dialog with layer selection
        dialog = wx.SingleChoiceDialog(self, f"Select Layer for Zone (1=Top, {max_layers}=Bottom):", "Zone Parameters", layers)
        # Set default selection to first layer (index 0)
        dialog.SetSelection(0)
        # Show dialog and wait for user action
        if dialog.ShowModal() == wx.ID_OK:
            try:
                # Get selected layer as string and convert to integer
                layer = int(dialog.GetStringSelection())
                # Validate layer is within allowed range
                if 1 <= layer <= max_layers:
                    # Store selected layer for current element
                    self.current_layer = layer
                    # Destroy dialog to free resources
                    dialog.Destroy()
                    # Log successful layer selection
                    self.logger.info(f"Set Zone layer to {layer}")
                    # Return success
                    return True
                else:
                    # Log error for invalid layer number
                    self.logger.error(f"Layer must be between 1 and {max_layers}")
                    # Show error message to user
                    wx.MessageBox(f"Layer must be between 1 and {max_layers}.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            except ValueError:
                # Log error if conversion to int failed
                self.logger.error("Layer must be a number")
                # Show error message for non-numeric input
                wx.MessageBox("Layer must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
        # Destroy dialog if canceled or failed
        dialog.Destroy()
        # Return failure
        return False

    def show_heater_dialog(self):
        """
        @brief Show popup dialog for Heater parameters (Layer dropdown with extreme layers and Power input).
        @details Prompts user to select a layer (Top or Bottom) and enter power value.
        @return True if dialog was successful, False otherwise.
        """
        # Get maximum number of PCB layers
        max_layers = int(self.pcb_params['layers'])
        # Only allow Top (1) and Bottom (max_layers) for heaters
        layers = ['1', str(max_layers)]
        # Create layer selection dialog
        dialog = wx.SingleChoiceDialog(self, f"Select Layer for Heater (1=Top, {max_layers}=Bottom):", "Heater Parameters", layers)
        # Default to Top layer
        dialog.SetSelection(0)
        # Show layer dialog
        if dialog.ShowModal() == wx.ID_OK:
            # Parse selected layer
            layer = int(dialog.GetStringSelection())
            # Close layer dialog
            dialog.Destroy()
            # Create power input dialog with default value
            power_dialog = wx.TextEntryDialog(self, "Enter Power (W) for Heater:", "Heater Parameters", value="5.0")
            # Show power input dialog
            if power_dialog.ShowModal() == wx.ID_OK:
                try:
                    # Convert input to float
                    power = float(power_dialog.GetValue())
                    # Validate power is non-negative
                    if power >= 0:
                        # Store layer and power
                        self.current_layer = layer
                        self.current_power = power
                        # Close power dialog
                        power_dialog.Destroy()
                        # Log heater creation
                        self.logger.info(f"Added Heater with layer {layer} and power {power}")
                        # Return success
                        return True
                    else:
                        # Log negative power error
                        self.logger.error("Power must be non-negative")
                        # Show error message
                        wx.MessageBox("Power must be non-negative.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                except ValueError:
                    # Log invalid number error
                    self.logger.error("Power must be a number")
                    # Show error message
                    wx.MessageBox("Power must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                # Close power dialog on error
                power_dialog.Destroy()
            else:
                # User canceled power input
                power_dialog.Destroy()
        else:
            # User canceled layer selection
            dialog.Destroy()
        # Return failure
        return False

    def show_sink_dialog(self):
        """
        @brief Show popup dialog for Sink parameters (Layer dropdown with extreme layers and HTC Extra input).
        @details Prompts user to select a layer (Top or Bottom) and enter HTC extra value.
        @return True if dialog was successful, False otherwise.
        """
        # Get max layers from PCB parameters
        max_layers = int(self.pcb_params['layers'])
        # Allow only Top and Bottom layers for sinks
        layers = ['1', str(max_layers)]
        # Create layer selection dialog
        dialog = wx.SingleChoiceDialog(self, f"Select Layer for Sink (1=Top, {max_layers}=Bottom):", "Sink Parameters", layers)
        # Default to Top layer
        dialog.SetSelection(0)
        # Show layer dialog
        if dialog.ShowModal() == wx.ID_OK:
            # Parse selected layer
            layer = int(dialog.GetStringSelection())
            # Close layer dialog
            dialog.Destroy()
            # Create HTC input dialog
            htc_dialog = wx.TextEntryDialog(self, "Enter HTC Extra (W/m²K) for Sink:", "Sink Parameters", value="10.0")
            # Show HTC dialog
            if htc_dialog.ShowModal() == wx.ID_OK:
                try:
                    # Convert to float
                    htc_extra = float(htc_dialog.GetValue())
                    # Validate non-negative
                    if htc_extra >= 0:
                        # Store values
                        self.current_layer = layer
                        self.current_htc_extra = htc_extra
                        # Close dialog
                        htc_dialog.Destroy()
                        # Log sink creation
                        self.logger.info(f"Added Sink with layer {layer} and HTC Extra {htc_extra}")
                        # Return success
                        return True
                    else:
                        # Log negative HTC error
                        self.logger.error("HTC Extra must be non-negative")
                        # Show error
                        wx.MessageBox("HTC Extra must be non-negative.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                except ValueError:
                    # Log invalid number
                    self.logger.error("HTC Extra must be a number")
                    # Show error
                    wx.MessageBox("HTC Extra must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                # Close on error
                htc_dialog.Destroy()
            else:
                # User canceled
                htc_dialog.Destroy()
        else:
            # User canceled layer selection
            dialog.Destroy()
        # Return failure
        return False

    def show_via_dialog(self):
        """
        @brief Show popup dialog for Via parameters (Diameter input with default 0.5mm).
        @details Prompts user to enter diameter value for the Via.
        @return True if dialog was successful, False otherwise.
        """
        # Create diameter input dialog with default 0.5mm
        dialog = wx.TextEntryDialog(self, "Enter Diameter (mm) for Via:", "Via Parameters", value="0.5")
        # Show dialog
        if dialog.ShowModal() == wx.ID_OK:
            try:
                # Convert input to float
                diameter = float(dialog.GetValue())
                # Validate positive diameter
                if diameter > 0:
                    # Store diameter
                    self.current_diameter = diameter
                    # Close dialog
                    dialog.Destroy()
                    # Log via diameter
                    self.logger.info(f"Set Via diameter to {diameter} mm")
                    # Return success
                    return True
                else:
                    # Log non-positive error
                    self.logger.error("Diameter must be positive")
                    # Show error
                    wx.MessageBox("Diameter must be positive.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            except ValueError:
                # Log non-numeric error
                self.logger.error("Diameter must be a number")
                # Show error
                wx.MessageBox("Diameter must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
        # Close dialog if canceled
        dialog.Destroy()
        # Return failure
        return False

    def on_show_results(self, event):
        """
        @brief Handle the Simulation Results button click.
        @details Reopens the results dialog using stored last_results and last_solver_mode.
                 Reloads content from memory (last_results) and temp dir files (charts, GIFs, etc.).
                 If no results available, shows a message.
        @param event Button event.
        """
        # Check if there are stored results
        if not self.last_results or not self.last_solver_mode:
            wx.MessageBox("No simulation results available. Run a simulation first.", "No Results", wx.OK | wx.ICON_INFORMATION)
            return

        try:
            # Log attempt to show results for debugging
            self.logger.debug(f"Showing stored results for mode: {self.last_solver_mode}")

            # Get dialog class name based on last mode
            dialog_class_name = ThermalSolver.get_results_dialog_class_name(self.last_solver_mode)

            # Create dialog instance based on mode
            if dialog_class_name == 'AnalyticResultsDialog':
                # Analytic: Pass stored results directly
                dialog = AnalyticResultsDialog(None, self.last_results, logger=self.logger)
                # Show full content (not loading)
                content_sizer = dialog._create_content_sizer()
                dialog.content_sizer.Clear(True)
                dialog.content_sizer.Add(content_sizer, 1, wx.ALL | wx.EXPAND, 5)
                dialog.save_btn.Show()
                dialog.btn_panel.Layout()

            elif dialog_class_name == 'LumpedResultsDialog':
                # Lumped: Pass stored results directly
                dialog = LumpedResultsDialog(None, self.last_results, logger=self.logger)
                # Show full content (not loading)
                content_sizer = dialog._create_content_sizer()
                dialog.content_sizer.Clear(True)
                dialog.content_sizer.Add(content_sizer, 1, wx.ALL | wx.EXPAND, 5)
                dialog.save_btn.Show()
                dialog.btn_panel.Layout()

            elif dialog_class_name == 'FDMResultsDialog':
                # FDM: Pass stored results, T_history, and zone_polys (from last_results)
                T_history = self.last_results.get('T_history', [])
                zone_polys = self.last_results.get('zone_polys', {})  # Assuming stored in results
                dialog = FDMResultsDialog(None, self.last_results, T_history=T_history, zone_polys=zone_polys, logger=self.logger)
                
                # Clear loading and create notebook with tabs
                dialog.content_sizer.Clear(True)
                dialog.notebook = wx.Notebook(dialog)
                dialog.notebook.AddPage(dialog.create_summary_tab(), "Summary")
                dialog.notebook.AddPage(dialog.create_heatmap_tab(), "Heatmap")
                dialog.notebook.AddPage(dialog.create_evolution_tab(), "Evolution")
                dialog.notebook.AddPage(dialog.create_hotspots_tab(), "Hotspots")
                dialog.notebook.AddPage(dialog.create_gradients_tab(), "Gradients")
                dialog.notebook.AddPage(dialog.create_zones_tab(), "Zones")
                dialog.content_sizer.Add(dialog.notebook, 1, wx.EXPAND | wx.ALL, 5)
                
                # Show save button
                dialog.save_btn.Show()
                dialog.btn_panel.Layout()

            else:
                # Fallback to Analytic if unknown mode
                dialog = AnalyticResultsDialog(None, self.last_results, logger=self.logger)
                self.logger.warning(f"Unknown dialog class {dialog_class_name}, using Analytic")
                # Show full content (not loading)
                content_sizer = dialog._create_content_sizer()
                dialog.content_sizer.Clear(True)
                dialog.content_sizer.Add(content_sizer, 1, wx.ALL | wx.EXPAND, 5)
                dialog.save_btn.Show()
                dialog.btn_panel.Layout()

            # Final layout and show dialog
            dialog.Layout()
            dialog.Update()
            dialog.Show()
            self.logger.debug("Results dialog reopened with stored data")

        except Exception as e:
            # Log error and show message to user
            self.logger.error(f"Failed to show results: {str(e)}")
            wx.MessageBox(f"Failed to show results: {str(e)}", "Error", wx.OK | wx.ICON_ERROR)

        finally:
            # Allow event propagation
            event.Skip()
 
    def on_toggle_grid(self, event):
        """
        @brief Handle the Show/Hide Grid button click to toggle grid visibility.
        @details Toggles the show_grid flag and refreshes the canvas to redraw the grid.
        @param event Button event.
        """
        # Toggle grid visibility state
        self.show_grid = not self.show_grid
        # Log current grid state
        self.logger.debug(f"Grid visibility: {'On' if self.show_grid else 'Off'}")
        # Redraw canvas to apply change
        self.canvas.Refresh()

    def on_grid_spacing_change(self, event):
        """
        @brief Handle the grid spacing selection change.
        @details Updates the grid_spacing value from the choice selection and refreshes the canvas.
        @param event Choice event.
        """
        # Get selected spacing value as string and convert to integer
        self.grid_spacing = int(self.grid_spacing_choice.GetStringSelection())
        # Log new grid spacing
        self.logger.debug(f"Grid spacing set to: {self.grid_spacing} pixels")
        # Redraw canvas with new spacing
        self.canvas.Refresh()

    def on_toggle_solver_grid(self, event):
        """
        @brief Handle the Show/Hide Solver Grid button click to toggle solver grid visibility.
        @details Toggles the show_solver_grid flag and refreshes the canvas to redraw the solver grid.
        @param event Button event.
        """
        # Toggle solver grid visibility
        self.show_solver_grid = not self.show_solver_grid
        # Log current solver grid state
        self.logger.debug(f"Solver grid visibility: {'On' if self.show_solver_grid else 'Off'}")
        # Redraw canvas to apply change
        self.canvas.Refresh()

    def update_global_params(self, event=None):
        """
        @brief Update self.pcb_params from widget values
        @details Read values from controls, validate numeric fields only
        @param event Event (can be None)
        """
        # Update solver mode from dropdown selection
        self.solver_mode = self.solver_mode_choice.GetStringSelection()

        # List of numeric parameters that require validation
        numeric_keys = [
            'width', 'height', 'thickness', 'copper_thickness', 'copper_conductivity',
            'fr4_conductivity', 'ambient_temp_global', 'convection_h_global',
            'convection_h_top', 'convection_h_bottom', 'grid_dx', 'grid_dy',
            'tolerance', 'over_relaxation_factor', 'hotspot_threshold'
        ]

        # Loop through all global parameter entries
        for key, entry in self.global_entries.items():
            try:
                # Handle Choice widgets (dropdowns)
                if isinstance(entry, wx.Choice):
                    sel = entry.GetSelection()
                    if sel == wx.NOT_FOUND:
                        raise ValueError("Invalid selection")
                    # Convert selection index to layer number (0-based to 1-based)
                    value = sel + 1 if key == 'layers' else entry.GetString(sel)
                    # Special handling: model_mode and boundary_mode are strings, no numeric check
                    if key in ['model_mode', 'boundary_mode']:
                        self.pcb_params[key] = value
                        continue

                # Handle integer SpinCtrl (e.g., max_iterations)
                elif isinstance(entry, wx.SpinCtrl):
                    value = entry.GetValue()

                # Handle float SpinCtrlDouble
                elif isinstance(entry, wx.SpinCtrlDouble):
                    value = entry.GetValue()

                else:
                    # Fallback: try to convert to float
                    value = float(entry.GetValue())

                # Validate only numeric fields
                if key in numeric_keys:
                    if value <= 0:
                        raise ValueError(f"{key} must be positive")

                # Validate layer count range
                if key == 'layers' and not (1 <= value <= config.DEFAULT_PCB_MAX_LAYERS):
                    raise ValueError(f"Layers must be between 1 and {config.DEFAULT_PCB_MAX_LAYERS}")

                # Validate over-relaxation factor range
                if key == 'over_relaxation_factor' and not (0 <= value <= 2.0):
                    raise ValueError("Over Relaxation Factor must be between 0.0 and 2.0")

                # Update pcb_params with valid value
                self.pcb_params[key] = value
                # Log successful update
                self.logger.debug(f"Parameter {key} updated to {value}")

            except Exception as exc:
                # Log warning but do not revert user input
                self.logger.warning(f"Invalid input for {key}: {exc}")
                # Show error dialog with parameter label
                wx.MessageBox(
                    f"Invalid value for {self.param_labels.get(key, key)}:\n{exc}",
                    "Input Error",
                    wx.OK | wx.ICON_ERROR,
                )

        # Check if number of layers changed
        old_layers = getattr(self, '_last_layer_count', self.pcb_params.get('layers', config.DEFAULT_PCB_LAYERS))
        if self.pcb_params.get('layers') != old_layers:
            # Update stored layer count
            self._last_layer_count = self.pcb_params['layers']
            max_layers = int(self.pcb_params['layers'])
            # Rebuild layer filter choices
            choices = ["All", f"1+{max_layers}"]
            if max_layers > 2:
                choices.append("Inner")
            choices.extend(str(i) for i in range(1, max_layers + 1))
            # Preserve current selection if possible
            cur = self.layer_choice.GetSelection()
            self.layer_choice.SetItems(choices)
            self.layer_choice.SetSelection(min(cur, len(choices) - 1))
            # Log layer filter update
            self.logger.debug(f"Updated layer filter for {max_layers} layers")

        # Redraw canvas to reflect any visual changes
        self.canvas.Refresh()

        # Propagate event if provided
        if event:
            event.Skip()

    def on_param_change(self, event):
        """
        @brief Handle parameter control changes
        @details Always update pcb_params when any control changes
        @param event Event from SpinCtrl, SpinCtrlDouble, or Choice
        """
        # Sync all parameters immediately
        self.update_global_params(event)
        # Log parameter change
        self.logger.debug("Parameter changed - auto-updating pcb_params")

    def on_solver_mode_change(self, event):
        """
        @brief Handle solver mode dropdown changes
        @details Update solver_mode and enable/disable FDM controls
        @param event Choice event from solver_mode_choice
        """
        # Get selected solver mode
        selected_mode = self.solver_mode_choice.GetString(self.solver_mode_choice.GetSelection())
        self.solver_mode = selected_mode
        # Log mode change
        self.logger.debug(f"Solver mode changed to: {self.solver_mode}")

        # Determine if FDM is selected
        is_fdm = self.solver_mode == "FDM"
        # Enable/disable FDM-specific controls
        self.global_entries['max_iterations'].Enable(is_fdm)
        self.global_entries['tolerance'].Enable(is_fdm)
        self.global_entries['grid_dx'].Enable(is_fdm)
        self.global_entries['grid_dy'].Enable(is_fdm)
        self.global_entries['over_relaxation_factor'].Enable(is_fdm)
        self.solver_grid_btn.Enable(is_fdm)
        self.global_entries['hotspot_threshold'].Enable(is_fdm)
        self.boundary_choice.Enable(is_fdm)

        # Enable Model Mode only for Analytic and Lumped
        self.model_mode_choice.Enable(self.solver_mode in ["Analytic", "Lumped"])

        # Propagate event
        event.Skip()

    def on_zoom_in(self, event):
        """
        @brief Handle Zoom In button click.
        @details Increases zoom_factor by 20%, caps at 1000%, updates label, and refreshes canvas.
        @param event Button event.
        """
        # Increase zoom by 20%
        self.zoom_factor *= 1.2
        # Cap zoom at 1000%
        self.zoom_factor = min(self.zoom_factor, 10.0)
        # Update zoom label
        self.zoom_label.SetLabel(f"Zoom: {int(self.zoom_factor * 100)}%")
        # Log zoom level
        self.logger.debug(f"Zoomed in to: {self.zoom_factor:.2f}x")
        # Redraw canvas
        self.canvas.Refresh()

    def on_zoom_out(self, event):
        """
        @brief Handle Zoom Out button click.
        @details Decreases zoom_factor by 20%, caps at 10%, updates label, and refreshes canvas.
        @param event Button event.
        """
        # Decrease zoom by 20%
        self.zoom_factor /= 1.2
        # Floor zoom at 10%
        self.zoom_factor = max(0.1, self.zoom_factor)
        # Update zoom label
        self.zoom_label.SetLabel(f"Zoom: {int(self.zoom_factor * 100)}%")
        # Log zoom level
        self.logger.debug(f"Zoomed out to: {self.zoom_factor:.2f}x")
        # Redraw canvas
        self.canvas.Refresh()

    def on_reset_zoom(self, event):
        """
        @brief Handle Reset Zoom button click to set zoom to 100% and reset pan.
        @details Resets zoom_factor to 1.0, pan_x/y to 0, updates label, and refreshes canvas.
        @param event Button event.
        """
        # Reset zoom to 100%
        self.zoom_factor = 1.0
        # Reset horizontal pan
        self.pan_x = 0
        # Reset vertical pan
        self.pan_y = 0
        # Update zoom label
        self.zoom_label.SetLabel("Zoom: 100%")
        # Log reset
        self.logger.debug("Zoom and pan reset")
        # Redraw canvas
        self.canvas.Refresh()
    
    def on_list_item_activated(self, event):
        """
        @brief Handle double-click on ListCtrl item to edit allowed fields directly.
        @details Opens appropriate dialog based on element type (Zone, Heater, Sink, Via) and updates element and table.
                Adjusted for new column indices: Layer now col 4, Power col 5, Diameter col 6, HTC col 7.
        @param event List event.
        """
        # Get index of double-clicked item
        item = event.GetIndex()
        # Get element type from column 0
        element_type = self.element_table.GetItem(item, 0).GetText()
        # Get element ID from column 1
        element_id = self.element_table.GetItem(item, 1).GetText()
        # Get maximum number of PCB layers
        max_layers = int(self.pcb_params['layers'])

        # Edit Zone: change layer only
        if element_type == "Zone":
            # Create list of all layer numbers as strings
            layers = [str(i) for i in range(1, max_layers + 1)]
            # Create dialog for layer selection
            dialog = wx.SingleChoiceDialog(self, f"Select Layer for Zone (1=Top, {max_layers}=Bottom):", "Edit Zone Layer", layers)
            # Get current layer from element (default to '1' if missing)
            current_layer = str(self.elements[item].get('layer', '1'))
            # Set current selection if valid, else default to first
            dialog.SetSelection(layers.index(current_layer) if current_layer in layers else 0)
            # Show dialog and wait for OK
            if dialog.ShowModal() == wx.ID_OK:
                try:
                    # Convert selected string to integer
                    layer = int(dialog.GetStringSelection())
                    # Validate layer is within range
                    if 1 <= layer <= max_layers:
                        # Update element's layer
                        self.elements[item]['layer'] = layer
                        # Update table cell in column 4
                        self.element_table.SetItem(item, 4, str(layer))
                        # Log successful update
                        self.logger.info(f"Updated {element_type} {element_id} - Layer: {layer}")
                    else:
                        # Log error for out-of-range layer
                        self.logger.error(f"Layer must be between 1 and {max_layers}")
                        # Show error dialog
                        wx.MessageBox(f"Layer must be between 1 and {max_layers}.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                except ValueError:
                    # Log error for non-numeric input
                    self.logger.error("Layer must be a number")
                    # Show error dialog
                    wx.MessageBox("Layer must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            # Destroy dialog to free resources
            dialog.Destroy()

        # Edit Heater: change layer (top/bottom) and power
        elif element_type == "Heater":
            # Allow only top and bottom layers
            layers = ['1', str(max_layers)]
            # Create custom dialog for combined editing
            dialog = wx.Dialog(self, title="Edit Heater Parameters")
            sizer = wx.BoxSizer(wx.VERTICAL)
            # Layer selection label
            layer_label = wx.StaticText(dialog, label=f"Layer (1=Top, {max_layers}=Bottom):")
            # Layer dropdown
            layer_choice = wx.Choice(dialog, choices=layers)
            # Set current layer selection
            current_layer = str(self.elements[item].get('layer', '1'))
            layer_choice.SetSelection(layers.index(current_layer) if current_layer in layers else 0)
            sizer.Add(layer_label, 0, wx.ALL, 5)
            sizer.Add(layer_choice, 0, wx.ALL | wx.EXPAND, 5)
            # Power input label
            power_label = wx.StaticText(dialog, label="Power (W):")
            # Power input field with current value
            power_input = wx.TextCtrl(dialog, value=str(self.elements[item].get('power', '5.0')))
            sizer.Add(power_label, 0, wx.ALL, 5)
            sizer.Add(power_input, 0, wx.ALL | wx.EXPAND, 5)
            # Create OK/Cancel button sizer
            btn_sizer = wx.StdDialogButtonSizer()
            ok_button = wx.Button(dialog, wx.ID_OK)
            cancel_button = wx.Button(dialog, wx.ID_CANCEL)
            btn_sizer.AddButton(ok_button)
            btn_sizer.AddButton(cancel_button)
            btn_sizer.Realize()
            sizer.Add(btn_sizer, 0, wx.ALL | wx.ALIGN_CENTER, 5)
            # Apply sizer and fit dialog
            dialog.SetSizer(sizer)
            dialog.Fit()
            # Show dialog and check OK
            if dialog.ShowModal() == wx.ID_OK:
                try:
                    # Get selected layer
                    layer = int(layer_choice.GetStringSelection())
                    # Get power input
                    power = float(power_input.GetValue())
                    # Validate power and layer
                    if power >= 0 and layer in [1, max_layers]:
                        # Update element
                        self.elements[item]['layer'] = layer
                        self.elements[item]['power'] = power
                        # Update table cells
                        self.element_table.SetItem(item, 4, str(layer))
                        self.element_table.SetItem(item, 5, str(power))
                        # Log update
                        self.logger.info(f"Updated {element_type} {element_id} - Layer: {layer}, Power: {power}")
                    else:
                        # Log validation error
                        self.logger.error(f"Power must be non-negative and Layer must be 1 or {max_layers}")
                        # Show error
                        wx.MessageBox(f"Power must be non-negative and Layer must be 1 or {max_layers}.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                except ValueError:
                    # Log numeric error
                    self.logger.error("Power must be a number")
                    # Show error
                    wx.MessageBox("Power must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            # Destroy dialog
            dialog.Destroy()

        # Edit Sink: change layer (top/bottom) and HTC extra
        elif element_type == "Sink":
            # Allow only top and bottom layers
            layers = ['1', str(max_layers)]
            # Create custom dialog
            dialog = wx.Dialog(self, title="Edit Sink Parameters")
            sizer = wx.BoxSizer(wx.VERTICAL)
            # Layer label
            layer_label = wx.StaticText(dialog, label=f"Layer (1=Top, {max_layers}=Bottom):")
            # Layer dropdown
            layer_choice = wx.Choice(dialog, choices=layers)
            # Set current layer
            current_layer = str(self.elements[item].get('layer', '1'))
            layer_choice.SetSelection(layers.index(current_layer) if current_layer in layers else 0)
            sizer.Add(layer_label, 0, wx.ALL, 5)
            sizer.Add(layer_choice, 0, wx.ALL | wx.EXPAND, 5)
            # HTC input label
            htc_label = wx.StaticText(dialog, label="HTC Extra (W/m²K):")
            # HTC input field
            htc_input = wx.TextCtrl(dialog, value=str(self.elements[item].get('htc_extra', '10.0')))
            sizer.Add(htc_label, 0, wx.ALL, 5)
            sizer.Add(htc_input, 0, wx.ALL | wx.EXPAND, 5)
            # OK/Cancel buttons
            btn_sizer = wx.StdDialogButtonSizer()
            ok_button = wx.Button(dialog, wx.ID_OK)
            cancel_button = wx.Button(dialog, wx.ID_CANCEL)
            btn_sizer.AddButton(ok_button)
            btn_sizer.AddButton(cancel_button)
            btn_sizer.Realize()
            sizer.Add(btn_sizer, 0, wx.ALL | wx.ALIGN_CENTER, 5)
            # Apply sizer
            dialog.SetSizer(sizer)
            dialog.Fit()
            # Show dialog
            if dialog.ShowModal() == wx.ID_OK:
                try:
                    # Get layer
                    layer = int(layer_choice.GetStringSelection())
                    # Get HTC value
                    htc_extra = float(htc_input.GetValue())
                    # Validate
                    if htc_extra >= 0 and layer in [1, max_layers]:
                        # Update element
                        self.elements[item]['layer'] = layer
                        self.elements[item]['htc_extra'] = htc_extra
                        # Update table
                        self.element_table.SetItem(item, 4, str(layer))
                        self.element_table.SetItem(item, 7, str(htc_extra))
                        # Log update
                        self.logger.info(f"Updated {element_type} {element_id} - Layer: {layer}, HTC Extra: {htc_extra}")
                    else:
                        # Log validation error
                        self.logger.error(f"HTC Extra must be non-negative and Layer must be 1 or {max_layers}")
                        # Show error
                        wx.MessageBox(f"HTC Extra must be non-negative and Layer must be 1 or {max_layers}.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                except ValueError:
                    # Log numeric error
                    self.logger.error("HTC Extra must be a number")
                    # Show error
                    wx.MessageBox("HTC Extra must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            # Destroy dialog
            dialog.Destroy()

        # Edit Via: change diameter only
        elif element_type == "Via":
            # Create diameter input dialog
            dialog = wx.TextEntryDialog(self, "Enter Diameter (mm) for Via:", "Edit Via Diameter", value=str(self.elements[item].get('diameter', '0.5')))
            # Show dialog
            if dialog.ShowModal() == wx.ID_OK:
                try:
                    # Convert to float
                    diameter = float(dialog.GetValue())
                    # Validate positive
                    if diameter > 0:
                        # Update element
                        self.elements[item]['diameter'] = diameter
                        # Update table column 6
                        self.element_table.SetItem(item, 6, str(diameter))
                        # Log update
                        self.logger.info(f"Updated {element_type} {element_id} - Diameter: {diameter}")
                        # Refresh canvas to redraw via
                        self.canvas.Refresh()
                    else:
                        # Log non-positive error
                        self.logger.error("Diameter must be positive")
                        # Show error
                        wx.MessageBox("Diameter must be positive.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                except ValueError:
                    # Log non-numeric error
                    self.logger.error("Diameter must be a number")
                    # Show error
                    wx.MessageBox("Diameter must be a number.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            # Destroy dialog
            dialog.Destroy()

    def on_paint(self, event):
        """
        @brief Handle painting of the canvas to draw the overlay image, regular grid, solver grid, PCB outline, all elements, the current element, duplicate/move/edit element preview, and vertex markers.
        @details Uses BufferedPaintDC for smooth rendering, applies transformations, draws layers in order: overlay, background, grids, outline, elements, previews.
                Filters polygon elements (Zones, Heaters, Sinks) based on Layer Filter selection, but always draws Vias.
                Overlay is drawn at original size with zoom_factor and pan applied, preserving aspect ratio.
        @param event Paint event.
        """
        # Create buffered DC for flicker-free drawing
        dc = wx.BufferedPaintDC(self.canvas)
        # Prepare DC for scrolling
        self.canvas.DoPrepareDC(dc)
        # Apply zoom scaling
        dc.SetUserScale(self.zoom_factor, self.zoom_factor)
        # Apply pan offset
        dc.SetDeviceOrigin(self.pan_x, self.pan_y)
        # Set background brush
        dc.SetBackground(wx.Brush(self.canvas.GetBackgroundColour()))
        # Clear canvas
        dc.Clear()

        # Draw overlay if loaded and visible
        if hasattr(self, 'overlay_image') and self.overlay_image is not None and hasattr(self, 'overlay_visible') and self.overlay_visible:
            # Convert image to bitmap
            bitmap = wx.Bitmap(self.overlay_image)
            if bitmap.IsOk():
                # Get overlay position
                x, y = self.overlay_pos
                # Draw bitmap at position
                dc.DrawBitmap(bitmap, int(x), int(y))
                # Log overlay draw
                self.logger.debug(f"Drawing overlay at ({x}, {y}) with zoom {self.zoom_factor:.2f}, dimensions {self.overlay_image.GetWidth()}x{self.overlay_image.GetHeight()}px")
            else:
                # Log invalid bitmap
                self.logger.warning("Invalid bitmap for overlay image")

        # Create GraphicsContext for advanced drawing
        gc = wx.GraphicsContext.Create(dc)

        # Calculate opacity from slider (0-255)
        opacity_value = int((self.opacity_slider.GetValue() / 100.0) * 255)

        # Draw PCB outline if two corners defined
        if len(self.pcb_corners) >= 2:
            x1, y1 = self.pcb_corners[0]
            x2, y2 = self.pcb_corners[1]
            x = min(x1, x2)
            y = min(y1, y2)
            width = abs(x2 - x1)
            height = abs(y2 - y1)
            # Set dashed black pen
            gc.SetPen(gc.CreatePen(wx.Pen(wx.Colour(0, 0, 0), 3, wx.PENSTYLE_SHORT_DASH)))
            # Set semi-transparent white brush
            gc.SetBrush(gc.CreateBrush(wx.Brush(wx.Colour(255, 255, 255, opacity_value))))
            # Draw rectangle
            gc.DrawRectangle(int(x), int(y), int(width), int(height))
            # Log outline draw
            self.logger.debug(f"Drawing PCB outline at ({x}, {y}) with size {width}x{height}, opacity {opacity_value}")

        # Draw regular grid if enabled
        if self.show_grid:
            # Set light grey dotted pen
            dc.SetPen(wx.Pen(wx.Colour(200, 200, 200), 1, wx.PENSTYLE_DOT))
            # Get visible canvas size
            canvas_width, canvas_height = self.canvas.GetSize()
            canvas_width /= self.zoom_factor
            canvas_height /= self.zoom_factor
            spacing = self.grid_spacing
            # Draw vertical grid lines
            x = -self.pan_x / self.zoom_factor
            while x < canvas_width:
                dc.DrawLine(int(x), 0, int(x), int(canvas_height))
                x += spacing
            x = -self.pan_x / self.zoom_factor
            while x > -canvas_width:
                dc.DrawLine(int(x), 0, int(x), int(canvas_height))
                x -= spacing
            # Draw horizontal grid lines
            y = -self.pan_y / self.zoom_factor
            while y < canvas_height:
                dc.DrawLine(0, int(y), int(canvas_width), int(y))
                y += spacing
            y = -self.pan_y / self.zoom_factor
            while y > -canvas_height:
                dc.DrawLine(0, int(y), int(canvas_width), int(y))
                y -= spacing

        # Draw solver grid if enabled and PCB corners defined
        if self.show_solver_grid and len(self.pcb_corners) >= 2:
            # Set blue dotted pen
            dc.SetPen(wx.Pen(wx.Colour(0, 0, 255), 1, wx.PENSTYLE_DOT))
            x1, y1 = self.pcb_corners[0]
            x2, y2 = self.pcb_corners[1]
            x = min(x1, x2)
            y = min(y1, y2)
            width = abs(x2 - x1)
            height = abs(y2 - y1)
            # Scale grid spacing to pixels
            dx = max(self.pcb_params['grid_dx'] * self.scale_factor, 2.0)
            dy = max(self.pcb_params['grid_dy'] * self.scale_factor, 2.0)
            # Log solver grid parameters
            self.logger.debug(f"Drawing solver grid: dx={self.pcb_params['grid_dx']:.2f} mm, dy={self.pcb_params['grid_dy']:.2f} mm, scale_factor={self.scale_factor:.2f}, dx_pixels={dx:.2f}, dy_pixels={dy:.2f}")
            # Draw vertical lines
            px = x
            while px <= x + width + dx * 0.1:
                dc.DrawLine(int(px), int(y), int(px), int(y + height))
                px += dx
            # Draw horizontal lines
            py = y
            while py <= y + height + dy * 0.1:
                dc.DrawLine(int(x), int(py), int(x + width), int(py))
                py += dy

        # Get current layer filter
        layer_filter = self.layer_choice.GetStringSelection()
        max_layers = int(self.pcb_params['layers'])

        # Draw all completed elements
        for element in self.elements:
            # Skip elements being moved or edited
            if (self.current_tool == "Move Selected" and self.moving_element == element) or \
            (self.current_tool == "Edit Selected" and self.editing_element == element):
                continue

            # Determine if element should be drawn based on layer filter
            if element["type"] == "Via":
                should_draw = True
            else:
                element_layer = element.get("layer", 1)
                if layer_filter == "All":
                    should_draw = True
                elif layer_filter == f"1+{max_layers}":
                    should_draw = element_layer in [1, max_layers]
                elif layer_filter == "Inner":
                    should_draw = max_layers > 2 and element_layer in range(2, max_layers)
                else:
                    should_draw = str(element_layer) == layer_filter
            if not should_draw:
                continue

            # Set pen and brush based on element type
            if element["type"] == "Zone":
                pen = gc.CreatePen(wx.Pen(wx.Colour(210, 180, 140), 2))
                brush = gc.CreateBrush(wx.Brush(wx.Colour(210, 180, 140, opacity_value)))
            elif element["type"] == "Heater":
                pen = gc.CreatePen(wx.Pen(wx.Colour(255, 182, 193), 2))
                brush = gc.CreateBrush(wx.Brush(wx.Colour(255, 182, 193, opacity_value)))
            elif element["type"] == "Via":
                pen = gc.CreatePen(wx.Pen(wx.Colour(0, 255, 0), 1))
                brush = gc.CreateBrush(wx.Brush(wx.Colour(0, 255, 0)))
            elif element["type"] == "Sink":
                pen = gc.CreatePen(wx.Pen(wx.Colour(70, 130, 180), 2))
                brush = gc.CreateBrush(wx.Brush(wx.Colour(70, 130, 180, opacity_value)))

            gc.SetPen(pen)
            gc.SetBrush(brush)

            # Draw element
            if element["type"] == "Via":
                x, y = element["points"][0]
                # Draw fixed-size circle
                gc.DrawEllipse(int(x) - 5, int(y) - 5, 10, 10)
                # Draw label
                font = wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
                gc.SetFont(gc.CreateFont(font))
                gc.DrawText(element["label"], int(x) + 10, int(y) - 5)
            else:
                # Create path for polygon
                path = gc.CreatePath()
                for i, (px, py) in enumerate(element["points"]):
                    if i == 0:
                        path.MoveToPoint(int(px), int(py))
                    else:
                        path.AddLineToPoint(int(px), int(py))
                path.CloseSubpath()
                # Fill and stroke
                gc.FillPath(path)
                gc.StrokePath(path)
                # Draw label at centroid
                if len(element["points"]) >= 3:
                    points_list = element["points"][:-1]
                    x_sum = sum(p[0] for p in points_list) / len(points_list)
                    y_sum = sum(p[1] for p in points_list) / len(points_list)
                    font = wx.Font(10, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
                    gc.SetFont(gc.CreateFont(font))
                    gc.DrawText(element["label"], int(x_sum), int(y_sum))

        # Draw current element (dashed lines)
        if self.current_tool in ["Add Zone", "Add Heater", "Add Sink"] and len(self.current_element) >= 2:
            points = [(int(x), int(y)) for x, y in self.current_element]
            dc.SetPen(wx.Pen(wx.Colour(128, 128, 128), 1, wx.PENSTYLE_SHORT_DASH))
            dc.DrawLines(points)

        # Draw vertex markers for current element
        if self.current_tool in ["Add Zone", "Add Heater", "Add Sink"] and len(self.current_element) >= 1:
            dc.SetPen(wx.Pen(wx.Colour(255, 0, 0), 1))
            dc.SetBrush(wx.Brush(wx.Colour(255, 0, 0)))
            for x, y in self.current_element:
                dc.DrawCircle(int(x), int(y), 5)

        # Draw PCB corner markers
        if self.current_tool == "Set PCB Size" and self.pcb_corners:
            dc.SetPen(wx.Pen(wx.Colour(255, 0, 0), 1))
            dc.SetBrush(wx.Brush(wx.Colour(255, 0, 0)))
            for x, y in self.pcb_corners:
                dc.DrawCircle(int(x), int(y), 5)

        # Draw duplicate preview
        if self.current_tool == "Duplicate Selected" and self.duplicate_element is not None and self.preview_position is not None:
            mouse_x, mouse_y = self.preview_position
            if self.duplicate_element["type"] == "Via":
                dc.SetPen(wx.Pen(wx.Colour(128, 128, 128), 1, wx.PENSTYLE_SHORT_DASH))
                dc.SetBrush(wx.Brush(wx.Colour(0, 255, 0, 128)))
                dc.DrawEllipse(int(mouse_x) - 5, int(mouse_y) - 5, 10, 10)
            else:
                points = self.duplicate_element["points"][:-1]
                centroid_x = sum(p[0] for p in points) / len(points)
                centroid_y = sum(p[1] for p in points) / len(points)
                dx = mouse_x - centroid_x
                dy = mouse_y - centroid_y
                preview_points = [(int(p[0] + dx), int(p[1] + dy)) for p in points]
                dc.SetPen(wx.Pen(wx.Colour(128, 128, 128), 1, wx.PENSTYLE_SHORT_DASH))
                dc.SetBrush(wx.Brush(wx.Colour(255, 255, 255, 0)))
                dc.DrawLines(preview_points + [preview_points[0]])

        # Draw move preview
        if self.current_tool == "Move Selected" and self.moving_element is not None and self.preview_position is not None:
            mouse_x, mouse_y = self.preview_position
            if self.moving_element["type"] == "Via":
                dc.SetPen(wx.Pen(wx.Colour(128, 128, 128), 1, wx.PENSTYLE_SHORT_DASH))
                dc.SetBrush(wx.Brush(wx.Colour(0, 255, 0, 128)))
                dc.DrawEllipse(int(mouse_x) - 5, int(mouse_y) - 5, 10, 10)
            else:
                points = self.moving_element["points"][:-1]
                centroid_x = sum(p[0] for p in points) / len(points)
                centroid_y = sum(p[1] for p in points) / len(points)
                dx = mouse_x - centroid_x
                dy = mouse_y - centroid_y
                preview_points = [(int(p[0] + dx), int(p[1] + dy)) for p in points]
                dc.SetPen(wx.Pen(wx.Colour(128, 128, 128), 1, wx.PENSTYLE_SHORT_DASH))
                dc.SetBrush(wx.Brush(wx.Colour(255, 255, 255, 0)))
                dc.DrawLines(preview_points + [preview_points[0]])

        # Draw editing element with vertex markers
        if self.current_tool == "Edit Selected" and self.editing_element is not None:
            points = self.editing_element["points"][:-1]
            if self.editing_vertex_index is not None and self.preview_position is not None:
                points = points.copy()
                points[self.editing_vertex_index] = self.preview_position
            preview_points = [(int(x), int(y)) for x, y in points]
            dc.SetPen(wx.Pen(wx.Colour(128, 128, 128), 1, wx.PENSTYLE_SHORT_DASH))
            dc.SetBrush(wx.Brush(wx.Colour(255, 255, 255, 0)))
            dc.DrawLines(preview_points + [preview_points[0]])
            dc.SetPen(wx.Pen(wx.Colour(255, 0, 0), 1))
            dc.SetBrush(wx.Brush(wx.Colour(255, 0, 0)))
            for i, (x, y) in enumerate(points):
                dc.DrawCircle(int(x), int(y), 5)
                if self.editing_vertex_index == i:
                    dc.SetPen(wx.Pen(wx.Colour(255, 255, 0), 1))
                    dc.DrawCircle(int(x), int(y), 6)
                    dc.SetPen(wx.Pen(wx.Colour(255, 0, 0), 1))

        # Log render summary
        self.logger.debug(f"Rendering {len(self.elements)} elements, current tool: {self.current_tool}, zoom: {self.zoom_factor:.2f}, pan: ({self.pan_x:.2f}, {self.pan_y:.2f}), scale_factor: {self.scale_factor:.2f}")

    def on_canvas_click(self, event):
        """
        @brief Handle left mouse click on canvas for adding points to elements, deleting, duplicating, moving, or editing selected elements, or setting PCB size.
        @details Converts mouse position to adjusted coordinates, handles each tool's logic (e.g., find closest element, add points, validate), updates elements/table/canvas.
                Parameter dialogs (e.g., layer, diameter) are shown in set_tool_and_show, so only drawing logic here.
        @param event Mouse event with x, y coordinates.
        """
        # Ensure canvas has focus to receive key events (e.g., Esc, A, D)
        self.canvas.SetFocus()

        # Convert mouse position from scrolled to unscrolled coordinates
        x, y = self.canvas.CalcUnscrolledPosition(event.GetPosition())

        # Apply zoom and pan to get real drawing coordinates
        adjusted_x = (x - self.pan_x) / self.zoom_factor
        adjusted_y = (y - self.pan_y) / self.zoom_factor

        # Ignore click if no tool is selected
        if self.current_tool is None:
            # Log warning to user
            self.logger.warning("No tool selected. Please select a tool from the left panel.")
            # Exit function early
            return

        # Handle "Set PCB Size" tool
        if self.current_tool == "Set PCB Size":
            # Check if PCB dimensions are set
            if self.pcb_params['width'] <= 0 or self.pcb_params['height'] <= 0:
                # Log error
                self.logger.error("PCB width and height must be set in PCB Parameters before setting size")
                # Show error dialog
                wx.MessageBox("Please set PCB width and height in PCB Parameters first.", "Invalid Input", wx.OK | wx.ICON_ERROR)
                # Exit function
                return
            # Add clicked point as a corner
            self.pcb_corners.append((adjusted_x, adjusted_y))
            # Log new corner point
            self.logger.debug(f"Added PCB corner point ({adjusted_x:.2f}, {adjusted_y:.2f}) in pixels")
            # If two corners are set, calculate scale factor
            if len(self.pcb_corners) == 2:
                self.on_set_pcb_size()
            # Refresh canvas to show new marker
            self.canvas.Refresh()
            # Exit function
            return

        # Handle "Delete Selected" tool
        if self.current_tool == "Delete Selected":
            # Initialize variables to find closest element
            closest_element = None
            min_distance = float('inf')
            # Loop through all elements
            for element in self.elements:
                # Handle Via (single point)
                if element["type"] == "Via":
                    ex, ey = element["points"][0]
                    distance = math.sqrt((adjusted_x - ex) ** 2 + (adjusted_y - ey) ** 2)
                else:
                    # Handle polygon elements (Zone, Heater, Sink)
                    points = element["points"][:-1]  # Exclude closing point
                    if points:
                        # Calculate centroid
                        centroid_x = sum(p[0] for p in points) / len(points)
                        centroid_y = sum(p[1] for p in points) / len(points)
                        distance = math.sqrt((adjusted_x - centroid_x) ** 2 + (adjusted_y - centroid_y) ** 2)
                    else:
                        # Skip empty polygon
                        continue
                # Check if this element is closer and within threshold
                if distance < min_distance and distance < 20 / self.zoom_factor:
                    min_distance = distance
                    closest_element = element
            # If a close element was found
            if closest_element:
                # Remove element from list
                self.elements.remove(closest_element)
                # Log deletion
                self.logger.info(f"Deleted {closest_element['type']} {closest_element['id']}")
                # Update element table
                self.update_element_table()
                # Refresh canvas
                self.canvas.Refresh()
            else:
                # Log no element found
                self.logger.debug(f"No element found near ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels")

        # Handle "Duplicate Selected" tool
        elif self.current_tool == "Duplicate Selected":
            # First click: select element to duplicate
            if self.duplicate_element is None:
                closest_element = None
                min_distance = float('inf')
                for element in self.elements:
                    if element["type"] == "Via":
                        ex, ey = element["points"][0]
                        distance = math.sqrt((adjusted_x - ex) ** 2 + (adjusted_y - ey) ** 2)
                    else:
                        points = element["points"][:-1]
                        if points:
                            centroid_x = sum(p[0] for p in points) / len(points)
                            centroid_y = sum(p[1] for p in points) / len(points)
                            distance = math.sqrt((adjusted_x - centroid_x) ** 2 + (adjusted_y - centroid_y) ** 2)
                        else:
                            continue
                    if distance < min_distance and distance < 20 / self.zoom_factor:
                        min_distance = distance
                        closest_element = element
                if closest_element:
                    # Deep copy selected element
                    import copy
                    self.duplicate_element = copy.deepcopy(closest_element)
                    # Log selection
                    self.logger.info(f"Selected {closest_element['type']} {closest_element['id']} for duplication")
                    # Update instruction label
                    self.title_label.SetLabel("Duplicate Selected - Click to place the copied element")
                else:
                    # Log no element found
                    self.logger.debug(f"No element found near ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels for duplication")
            # Second click: place the copy
            else:
                if self.duplicate_element["type"] == "Via":
                    # Update point position
                    self.duplicate_element["points"] = [(adjusted_x, adjusted_y)]
                    # Increment counter
                    self.element_counters["Via"] += 1
                    # Set new ID and label
                    self.duplicate_element["id"] = self.element_counters["Via"]
                    self.duplicate_element["label"] = f"Via {self.duplicate_element['id']}"
                else:
                    # Calculate centroid of original
                    points = self.duplicate_element["points"][:-1]
                    centroid_x = sum(p[0] for p in points) / len(points)
                    centroid_y = sum(p[1] for p in points) / len(points)
                    # Calculate offset
                    dx = adjusted_x - centroid_x
                    dy = adjusted_y - centroid_y
                    # Translate all points
                    self.duplicate_element["points"] = [(p[0] + dx, p[1] + dy) for p in self.duplicate_element["points"]]
                    # Increment counter
                    self.element_counters[self.duplicate_element["type"]] += 1
                    # Set new ID and label
                    self.duplicate_element["id"] = self.element_counters[self.duplicate_element["type"]]
                    self.duplicate_element["label"] = f"{self.duplicate_element['type']} {self.duplicate_element['id']}"
                # Add to elements list
                self.elements.append(self.duplicate_element)
                # Log duplication
                self.logger.info(f"Duplicated {self.duplicate_element['type']} {self.duplicate_element['id']} at ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels")
                # Reset duplication state
                self.duplicate_element = None
                # Update instruction label
                self.title_label.SetLabel("Duplicate Selected - Click to select element, click to place copy")
                # Update table and canvas
                self.update_element_table()
                self.canvas.Refresh()

        # Handle "Move Selected" tool
        elif self.current_tool == "Move Selected":
            # First click: select element to move
            if self.moving_element is None:
                closest_element = None
                min_distance = float('inf')
                for element in self.elements:
                    if element["type"] == "Via":
                        ex, ey = element["points"][0]
                        distance = math.sqrt((adjusted_x - ex) ** 2 + (adjusted_y - ey) ** 2)
                    else:
                        points = element["points"][:-1]
                        if points:
                            centroid_x = sum(p[0] for p in points) / len(points)
                            centroid_y = sum(p[1] for p in points) / len(points)
                            distance = math.sqrt((adjusted_x - centroid_x) ** 2 + (adjusted_y - centroid_y) ** 2)
                        else:
                            continue
                    if distance < min_distance and distance < 20 / self.zoom_factor:
                        min_distance = distance
                        closest_element = element
                if closest_element:
                    import copy
                    # Set moving element
                    self.moving_element = closest_element
                    # Store original for Esc cancel
                    self.original_element = copy.deepcopy(closest_element)
                    # Set preview position
                    self.preview_position = (adjusted_x, adjusted_y)
                    # Log selection
                    self.logger.info(f"Selected {closest_element['type']} {closest_element['id']} for moving")
                    # Update instruction
                    self.title_label.SetLabel("Move Selected - Click to place element, esc to cancel")
                    # Refresh to show preview
                    self.canvas.Refresh()
                else:
                    # Log no element found
                    self.logger.debug(f"No element found near ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels for moving")
            # Second click: place the element
            else:
                if self.moving_element["type"] == "Via":
                    # Update single point
                    self.moving_element["points"] = [(adjusted_x, adjusted_y)]
                else:
                    # Calculate centroid
                    points = self.moving_element["points"][:-1]
                    centroid_x = sum(p[0] for p in points) / len(points)
                    centroid_y = sum(p[1] for p in points) / len(points)
                    # Calculate offset
                    dx = adjusted_x - centroid_x
                    dy = adjusted_y - centroid_y
                    # Translate all points
                    self.moving_element["points"] = [(p[0] + dx, p[1] + dy) for p in self.moving_element["points"]]
                # Log move completion
                self.logger.info(f"Moved {self.moving_element['type']} {self.moving_element['id']} to ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels")
                # Reset move state
                self.moving_element = None
                self.original_element = None
                self.preview_position = None
                # Update instruction
                self.title_label.SetLabel("Move Selected - Click to select element, click to place, esc to cancel")
                # Update UI
                self.update_element_table()
                self.canvas.Refresh()

        # Handle "Edit Selected" tool
        elif self.current_tool == "Edit Selected":
            # First click: select polygon to edit
            if self.editing_element is None:
                closest_element = None
                min_distance = float('inf')
                for element in self.elements:
                    # Skip Vias (not editable)
                    if element["type"] == "Via":
                        continue
                    points = element["points"][:-1]
                    if points:
                        centroid_x = sum(p[0] for p in points) / len(points)
                        centroid_y = sum(p[1] for p in points) / len(points)
                        distance = math.sqrt((adjusted_x - centroid_x) ** 2 + (adjusted_y - centroid_y) ** 2)
                        if distance < min_distance and distance < 20 / self.zoom_factor:
                            min_distance = distance
                            closest_element = element
                if closest_element:
                    import copy
                    # Set editing element
                    self.editing_element = closest_element
                    # Store original for Esc
                    self.original_element = copy.deepcopy(closest_element)
                    # Log selection
                    self.logger.info(f"Selected {closest_element['type']} {closest_element['id']} for editing")
                    # Update instruction
                    self.title_label.SetLabel("Edit Selected - Click to select polygon, drag vertices, 'a' to add vertex, 'd' to delete vertex, double-click or esc to finish")
                    # Refresh to show markers
                    self.canvas.Refresh()
                else:
                    # Log no polygon found
                    self.logger.debug(f"No polygon found near ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels for editing")
            # Subsequent clicks: select vertex to move
            else:
                points = self.editing_element["points"][:-1]
                for i, (px, py) in enumerate(points):
                    distance = math.sqrt((adjusted_x - px) ** 2 + (adjusted_y - py) ** 2)
                    if distance < 5 / self.zoom_factor:
                        # Set vertex to move
                        self.editing_vertex_index = i
                        self.preview_position = (adjusted_x, adjusted_y)
                        # Log vertex selection
                        self.logger.debug(f"Selected vertex {i} at ({px:.2f}, {py:.2f}) for moving")
                        # Update instruction
                        self.title_label.SetLabel("Edit Selected - Drag vertex to move, release to place")
                        break
                else:
                    # Log no vertex found
                    self.logger.debug(f"No vertex found near ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels for moving")
                # Refresh to show preview
                self.canvas.Refresh()

        # Handle polygon drawing tools (Zone, Heater, Sink)
        elif self.current_tool in ["Add Zone", "Add Heater", "Add Sink"]:
            # Add clicked point to current element
            self.current_element.append((adjusted_x, adjusted_y))
            # Log point addition
            self.logger.debug(f"Added point ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels to current {self.current_tool}")
            # Refresh canvas to show new line
            self.canvas.Refresh()

        # Handle "Add Via" tool (single point)
        elif self.current_tool == "Add Via":
            # Increment Via counter
            self.element_counters["Via"] += 1
            element_id = self.element_counters["Via"]
            # Create new Via element
            self.elements.append({
                "type": "Via",
                "id": element_id,
                "points": [(adjusted_x, adjusted_y)],
                "diameter": self.current_diameter,
                "label": f"Via {element_id}"
            })
            # Log Via creation
            self.logger.info(f"Added Via {element_id} at ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels with diameter {self.current_diameter}")
            # Clear temporary diameter
            del self.current_diameter
            # Update table and canvas
            self.update_element_table()
            self.canvas.Refresh()

    def on_canvas_double_click(self, event):
        """
        @brief Handle double-click on canvas to close polygons or finish editing.
        @details For polygon tools, creates element from current_element and adds to elements; for Edit, finishes editing mode.
        @param event Mouse event with x, y coordinates.
        """
        # Ensure canvas has focus for key events (Esc, A, D)
        self.canvas.SetFocus()

        # Handle polygon drawing tools (Zone, Heater, Sink)
        if self.current_tool in ["Add Zone", "Add Heater", "Add Sink"]:
            # Require at least 3 points to form a valid polygon
            if len(self.current_element) >= 3:
                # Extract element type from tool name
                element_type = self.current_tool.replace("Add ", "")
                # Increment counter for this element type
                self.element_counters[element_type] += 1
                # Get new unique ID
                element_id = self.element_counters[element_type]
                # Create new element dictionary
                new_element = {
                    "type": element_type,
                    "id": element_id,
                    # Close the polygon by repeating first point at the end
                    "points": self.current_element + [self.current_element[0]],
                    "label": f"{element_type} {element_id}"
                }
                # Assign layer from previously selected value
                new_element["layer"] = self.current_layer
                # Add specific parameters based on type
                if element_type == "Heater":
                    new_element["power"] = self.current_power
                elif element_type == "Sink":
                    new_element["htc_extra"] = self.current_htc_extra
                # Add completed element to main list
                self.elements.append(new_element)
                # Log polygon closure
                self.logger.info(f"Closed {element_type} {element_id} with {len(self.current_element)} points")
                # Clean up temporary variables
                del self.current_layer
                if element_type == "Heater":
                    del self.current_power
                elif element_type == "Sink":
                    del self.current_htc_extra
                # Reset current element list
                self.current_element = []
                # Exit drawing mode
                self.current_tool = None
                # Restore default instruction
                self.title_label.SetLabel("Select a tool from the left panel to start drawing")
                # Update element table and canvas
                self.update_element_table()
                self.canvas.Refresh()

        # Handle "Edit Selected" tool - finish editing
        elif self.current_tool == "Edit Selected" and self.editing_element is not None:
            # Reset editing state
            self.editing_element = None
            self.editing_vertex_index = None
            self.original_element = None
            self.preview_position = None
            # Restore instruction label
            self.title_label.SetLabel("Edit Selected - Click to select polygon, drag vertices, 'a' to add vertex, 'd' to delete vertex, double-click or esc to finish")
            # Log completion
            self.logger.info("Finished editing polygon")
            # Update table and canvas
            self.update_element_table()
            self.canvas.Refresh()

    def on_set_pcb_size(self):
        """
        @brief Calculate scale factor based on two clicked corners and fixed PCB dimensions from PCB Parameters.
        @details Computes pixel-to-mm scale from corner distances and PCB params, logs, resets tool, refreshes canvas.
        """
        # Validate exactly two corners were clicked
        if len(self.pcb_corners) != 2:
            # Log error
            self.logger.error("Two corners are required to set PCB size")
            # Show error dialog
            wx.MessageBox("Please click two opposite corners to set PCB size.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            # Exit function
            return

        # Validate PCB dimensions are positive
        if self.pcb_params['width'] <= 0 or self.pcb_params['height'] <= 0:
            # Log error
            self.logger.error("PCB width and height must be set in PCB Parameters before setting size")
            # Show error dialog
            wx.MessageBox("Please set PCB width and height in PCB Parameters first.", "Invalid Input", wx.OK | wx.ICON_ERROR)
            # Reset state on error
            self.pcb_corners = []
            self.current_tool = None
            self.title_label.SetLabel("Select a tool from the left panel to start drawing")
            self.canvas.Refresh()
            # Exit function
            return

        # Extract corner coordinates
        x1, y1 = self.pcb_corners[0]
        x2, y2 = self.pcb_corners[1]
        # Calculate pixel distances
        pixel_width = abs(x2 - x1)
        pixel_height = abs(y2 - y1)

        # Get PCB dimensions in mm
        width_mm = self.pcb_params['width']
        height_mm = self.pcb_params['height']

        # Calculate scale factors (pixels per mm)
        scale_factor_x = pixel_width / width_mm
        scale_factor_y = pixel_height / height_mm
        # Use average scale factor for consistency
        self.scale_factor = (scale_factor_x + scale_factor_y) / 2

        # Log final scale factor
        self.logger.info(f"Set PCB scale factor to {self.scale_factor:.2f} pixels/mm for width={width_mm:.2f} mm, height={height_mm:.2f} mm")

        # Reset tool and UI
        self.current_tool = None
        self.title_label.SetLabel("Select a tool from the left panel to start drawing")
        # Refresh canvas to remove corner markers
        self.canvas.Refresh()

    def on_middle_down(self, event):
        """
        @brief Handle middle mouse button down to start panning.
        @details Sets panning flag and starting position, changes cursor to hand.
        @param event Mouse event with x, y coordinates.
        """
        # Set flag to indicate panning is active
        self.is_panning = True
        # Get current mouse position in unscrolled canvas coordinates
        self.pan_start_x, self.pan_start_y = self.canvas.CalcUnscrolledPosition(event.GetPosition())
        # Change cursor to hand icon to indicate panning mode
        self.canvas.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        # Allow event to propagate to other handlers
        event.Skip()

    def on_mouse_motion(self, event):
        """
        @brief Handle mouse motion on canvas for panning, previewing duplicated/moving elements, or editing vertices.
        @details Updates pan offsets during drag, sets preview_position for tools, refreshes canvas for real-time preview.
        @param event Mouse event with x, y coordinates.
        """
        # Convert mouse position to unscrolled canvas coordinates
        x, y = self.canvas.CalcUnscrolledPosition(event.GetPosition())
        # Apply current pan and zoom to get adjusted drawing coordinates
        adjusted_x = (x - self.pan_x) / self.zoom_factor
        adjusted_y = (y - self.pan_y) / self.zoom_factor

        # Handle panning when middle mouse button is held and dragging
        if event.Dragging() and event.MiddleIsDown():
            # Calculate movement delta from last recorded position
            dx = x - self.pan_start_x
            dy = y - self.pan_start_y
            # Update global pan offsets
            self.pan_x += dx
            self.pan_y += dy
            # Update starting point for next motion event
            self.pan_start_x = x
            self.pan_start_y = y
            # Refresh canvas to show updated view
            self.canvas.Refresh()
            # Log panning movement
            self.logger.debug(f"Panning: dx={dx:.2f}, dy={dy:.2f}, new pan: ({self.pan_x:.2f}, {self.pan_y:.2f})")

        # Handle preview for "Duplicate Selected" tool
        if self.current_tool == "Duplicate Selected" and self.duplicate_element is not None:
            # Update preview position with current mouse coordinates
            self.preview_position = (adjusted_x, adjusted_y)
            # Refresh canvas to show live preview
            self.canvas.Refresh()
            # Log preview position
            self.logger.debug(f"Previewing duplicate at ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels")

        # Handle preview for "Move Selected" tool
        if self.current_tool == "Move Selected" and self.moving_element is not None:
            # Update preview position
            self.preview_position = (adjusted_x, adjusted_y)
            # Refresh canvas
            self.canvas.Refresh()
            # Log preview position
            self.logger.debug(f"Previewing move at ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels")

        # Handle vertex drag preview in "Edit Selected" tool
        if self.current_tool == "Edit Selected" and self.editing_element is not None and self.editing_vertex_index is not None and event.LeftIsDown():
            # Update preview position during drag
            self.preview_position = (adjusted_x, adjusted_y)
            # Refresh canvas to show vertex movement
            self.canvas.Refresh()
            # Log vertex preview
            self.logger.debug(f"Previewing vertex move at ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels")

    def on_middle_up(self, event):
        """
        @brief Handle middle mouse button release to stop panning.
        @details Resets panning flag and cursor to default.
        @param event Mouse event.
        """
        # Reset panning flag
        self.is_panning = False
        # Restore default cursor (arrow)
        self.canvas.SetCursor(wx.NullCursor)
        # Allow event propagation
        event.Skip()

    def on_mouse_wheel(self, event):
        """
        @brief Handle mouse wheel scroll to adjust zoom level.
        @details Zooms in/out by 10% based on rotation, clamps bounds, updates label, refreshes canvas.
        @param event Mouse wheel event.
        """
        # Get wheel rotation direction and magnitude
        rotation = event.GetWheelRotation()
        # Zoom in on upward scroll
        if rotation > 0:
            self.zoom_factor *= 1.1
        # Zoom out on downward scroll
        elif rotation < 0:
            self.zoom_factor /= 1.1
        # Clamp zoom factor between 10% and 1000%
        self.zoom_factor = max(0.1, min(self.zoom_factor, 10.0))
        # Update zoom percentage label
        self.zoom_label.SetLabel(f"Zoom: {int(self.zoom_factor * 100)}%")
        # Refresh canvas to apply new zoom
        self.canvas.Refresh()
        # Allow event to propagate
        event.Skip()

    def on_mouse_left_up(self, event):
        """
        @brief Handle mouse left button release for moving vertices in Edit Selected.
        @details Updates vertex position in editing_element, resets indices, updates table/canvas.
        @param event Mouse event with x, y coordinates.
        """
        # Check if we are in Edit Selected mode and dragging a vertex
        if self.current_tool == "Edit Selected" and self.editing_element is not None and self.editing_vertex_index is not None:
            # Convert mouse position to unscrolled canvas coordinates
            x, y = self.canvas.CalcUnscrolledPosition(event.GetPosition())
            # Apply current pan and zoom to get real drawing coordinates
            adjusted_x = (x - self.pan_x) / self.zoom_factor
            adjusted_y = (y - self.pan_y) / self.zoom_factor
            # Get the points list of the editing element
            points = self.editing_element["points"]
            # Update the selected vertex with new position
            points[self.editing_vertex_index] = (adjusted_x, adjusted_y)
            # Ensure the closing point matches the first point
            points[-1] = points[0]
            # Log the vertex move completion
            self.logger.info(f"Moved vertex {self.editing_vertex_index} of {self.editing_element['type']} {self.editing_element['id']} to ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels")
            # Reset vertex editing state
            self.editing_vertex_index = None
            self.preview_position = None
            # Restore instruction label to default edit mode
            self.title_label.SetLabel("Edit Selected - Drag vertices, 'A' to add vertex, 'D' to delete vertex, double-click or Esc to finish")
            # Update element table to reflect changes
            self.update_element_table()
            # Refresh canvas to show final position
            self.canvas.Refresh()

    def on_key_down(self, event):
        """
        @brief Handle key press events for adding/deleting vertices or aborting operations.
        @details Processes Esc to abort, 'a/A' to add vertex near edge, 'd/D' to delete nearest vertex (min 3 points).
        @param event Key event.
        """
        # Get the key code of the pressed key
        keycode = event.GetKeyCode()
        # Log the key pressed for debugging
        self.logger.debug(f"Key pressed: {keycode}")

        # Handle Esc key to abort current operation
        if keycode == wx.WXK_ESCAPE:
            # Abort "Move Selected" tool
            if self.current_tool == "Move Selected" and self.moving_element is not None:
                # Restore original element points
                self.moving_element["points"] = self.original_element["points"]
                # Reset moving state
                self.moving_element = None
                self.original_element = None
                self.preview_position = None
                # Update instruction label
                self.title_label.SetLabel("Move Selected - Click to select element, click to place, esc to cancel")
                # Log abort
                self.logger.info("Aborted move operation")
                # Refresh canvas
                self.canvas.Refresh()
            # Abort "Edit Selected" tool
            elif self.current_tool == "Edit Selected" and self.editing_element is not None:
                # Restore original element points
                self.editing_element["points"] = self.original_element["points"]
                # Reset editing state
                self.editing_element = None
                self.editing_vertex_index = None
                self.original_element = None
                self.preview_position = None
                # Update instruction label
                self.title_label.SetLabel("Edit Selected - Click to select polygon, drag vertices, 'a' to add vertex, 'd' to delete vertex, double-click or esc to finish")
                # Log abort
                self.logger.info("Aborted edit operation")
                # Refresh canvas
                self.canvas.Refresh()
            # Exit function after handling Esc
            return

        # Handle key bindings only in "Edit Selected" mode
        if self.current_tool == "Edit Selected" and self.editing_element is not None:
            # Get current mouse position in unscrolled coordinates
            x, y = self.canvas.CalcUnscrolledPosition(self.canvas.ScreenToClient(wx.GetMousePosition()))
            # Apply pan and zoom to get adjusted coordinates
            adjusted_x = (x - self.pan_x) / self.zoom_factor
            adjusted_y = (y - self.pan_y) / self.zoom_factor

            # Handle 'a' or 'A' to add a vertex near the closest edge
            if keycode in [ord('a'), ord('A')]:
                # Get polygon points excluding closing point
                points = self.editing_element["points"][:-1]
                # Initialize variables to find closest edge
                min_distance = float('inf')
                closest_edge = None
                closest_point = None
                # Loop through all edges of the polygon
                for i in range(len(points)):
                    p1 = points[i]
                    p2 = points[(i + 1) % len(points)]
                    px, py = adjusted_x, adjusted_y
                    x1, y1 = p1
                    x2, y2 = p2
                    # Vector from p1 to p2
                    dx = x2 - x1
                    dy = y2 - y1
                    length_sq = dx * dx + dy * dy
                    # Skip degenerate edge
                    if length_sq == 0:
                        continue
                    # Calculate projection parameter t (0 to 1)
                    t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / length_sq))
                    # Calculate projection point on edge
                    proj_x = x1 + t * dx
                    proj_y = y1 + t * dy
                    # Calculate distance from mouse to projection point
                    distance = math.sqrt((px - proj_x) ** 2 + (py - proj_y) ** 2)
                    # Check if this is the closest edge within threshold
                    if distance < min_distance and distance < 10 / self.zoom_factor:
                        min_distance = distance
                        closest_edge = i
                        closest_point = (proj_x, proj_y)
                # If a close edge was found
                if closest_edge is not None:
                    # Insert new vertex after the closest edge
                    self.editing_element["points"].insert(closest_edge + 1, closest_point)
                    # Update closing point to match first
                    self.editing_element["points"][-1] = self.editing_element["points"][0]
                    # Log vertex addition
                    self.logger.info(f"Added vertex to {self.editing_element['type']} {self.editing_element['id']} at ({closest_point[0]:.2f}, {closest_point[1]:.2f}) pixels")
                    # Update table and canvas
                    self.update_element_table()
                    self.canvas.Refresh()
                else:
                    # Log no edge found
                    self.logger.debug(f"No edge found near ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels for adding vertex")

            # Handle 'd' or 'D' to delete the nearest vertex
            elif keycode in [ord('d'), ord('D')]:
                # Get polygon points excluding closing point
                points = self.editing_element["points"][:-1]
                # Prevent deletion if polygon would have less than 3 vertices
                if len(points) <= 3:
                    # Log warning
                    self.logger.warning("Cannot delete vertex: polygon must have at least 3 vertices")
                    # Show error dialog
                    wx.MessageBox("Polygon must have at least 3 vertices.", "Invalid Operation", wx.OK | wx.ICON_ERROR)
                    # Exit function
                    return
                # Initialize variables to find closest vertex
                min_distance = float('inf')
                closest_vertex = None
                # Loop through all vertices
                for i, (px, py) in enumerate(points):
                    distance = math.sqrt((adjusted_x - px) ** 2 + (adjusted_y - py) ** 2)
                    # Check if this is the closest vertex within threshold
                    if distance < min_distance and distance < 5 / self.zoom_factor:
                        min_distance = distance
                        closest_vertex = i
                # If a close vertex was found
                if closest_vertex is not None:
                    # Remove the vertex
                    self.editing_element["points"].pop(closest_vertex)
                    # Update closing point
                    self.editing_element["points"][-1] = self.editing_element["points"][0]
                    # Log deletion
                    self.logger.info(f"Deleted vertex {closest_vertex} from {self.editing_element['type']} {self.editing_element['id']}")
                    # Update table and canvas
                    self.update_element_table()
                    self.canvas.Refresh()
                else:
                    # Log no vertex found
                    self.logger.debug(f"No vertex found near ({adjusted_x:.2f}, {adjusted_y:.2f}) pixels for deletion")

    def update_element_table(self):
        """
        @brief Update the ListCtrl table with the current elements.
        @details Clears table, iterates elements, sets columns with formatted data:
                - Coordinates: List of (x,y) points (truncated if >3 for polygons; single for Vias).
                Formatted as "(x, y)" with space after comma for readability.
                - Num Vertices: Number of points for polygons; "1" for Vias; "-" if not applicable.
                - Empty fields: Display "-" for better visual perception.
                Logs the update count.
        """
        # Clear all existing items in the table to start fresh
        self.element_table.DeleteAllItems()

        # Iterate through all elements in the current design
        for element in self.elements:
            # Insert new row and set element type in column 0
            index = self.element_table.InsertItem(self.element_table.GetItemCount(), element["type"])
            # Set element ID in column 1
            self.element_table.SetItem(index, 1, str(element["id"]))

            # Handle Via elements (single point)
            if element["type"] == "Via":
                # Extract the single point coordinates
                x, y = element["points"][0]
                # Format coordinates as "(x, y)" with 1 decimal place
                coords_str = f"({x:.1f}, {y:.1f})"
                # Vias have 1 vertex
                num_vertices_str = "1"
            else:
                # Handle polygon elements (Zone, Heater, Sink)
                # Exclude the closing point (last entry repeats first)
                points = element["points"][:-1]
                # Count number of vertices
                num_vertices = len(points)
                # Show all points if 3 or fewer
                if num_vertices <= 3:
                    # Format each point as "(x, y)" with space after comma
                    coords_str = " ".join([f"({p[0]:.1f}, {p[1]:.1f})" for p in points])
                else:
                    # Truncate long polygons: show first 3 points + "..."
                    truncated = [f"({p[0]:.1f}, {p[1]:.1f})" for p in points[:3]]
                    coords_str = " ".join(truncated) + " ..."
                # Set number of vertices as string
                num_vertices_str = str(num_vertices)

            # Set coordinates in column 2
            self.element_table.SetItem(index, 2, coords_str)
            # Set number of vertices in column 3
            self.element_table.SetItem(index, 3, num_vertices_str)

            # Set layer in column 4: show "-" if missing
            layer_str = str(element.get("layer", ""))
            self.element_table.SetItem(index, 4, "-" if not layer_str else layer_str)

            # Set power in column 5: show "-" if missing (only for Heater)
            power_str = str(element.get("power", ""))
            self.element_table.SetItem(index, 5, "-" if not power_str else power_str)

            # Set diameter in column 6: show "-" if missing (only for Via)
            diameter_str = str(element.get("diameter", ""))
            self.element_table.SetItem(index, 6, "-" if not diameter_str else diameter_str)

            # Set HTC extra in column 7: show "-" if missing (only for Sink)
            htc_str = str(element.get("htc_extra", ""))
            self.element_table.SetItem(index, 7, "-" if not htc_str else htc_str)

        # Log total number of elements after update
        self.logger.debug(f"Updated element table with {len(self.elements)} elements")

    def renumber_elements(self):
        """
        @brief Renumber all elements sequentially by type to remove ID gaps.
        @details Groups elements by type (Zone, Heater, Via, Sink), reassigns IDs starting from 1,
                updates labels, and logs the changes for debugging.
        """
        # Create dictionary to group elements by type
        type_groups = {"Zone": [], "Heater": [], "Via": [], "Sink": []}
        # Populate groups from current elements
        for element in self.elements:
            if element["type"] in type_groups:
                type_groups[element["type"]].append(element)

        # Renumber each group sequentially starting from 1
        for elem_type, group in type_groups.items():
            for i, element in enumerate(group, start=1):
                # Store old ID for logging
                old_id = element["id"]
                # Assign new sequential ID
                element["id"] = i
                # Update label to match new ID
                element["label"] = f"{elem_type} {i}"
                # Log renumbering change
                self.logger.debug(f"Renumbered {elem_type} from ID {old_id} to {i}")

        # Update element counters to reflect highest ID in each group
        for elem_type in type_groups:
            if type_groups[elem_type]:
                # Set counter to max ID in group (for future additions)
                self.element_counters[elem_type] = max(e["id"] for e in type_groups[elem_type])
            else:
                # Reset counter if no elements of this type
                self.element_counters[elem_type] = 0

        # Log successful completion
        self.logger.info("All elements renumbered successfully")

    def load_overlay(self, event):
        """
        @brief Handle the Load Overlay button click to load a JPG/PNG image without resizing.
        @details Opens a file dialog to select an image, validates size and dimensions against config limits,
                loads it as wx.Image to preserve original dimensions, and sets overlay position to (0, 0).
                The image is not saved to DEFAULT_OUTPUT_DIR until export.
        @param event Button event.
        """
        # Create file dialog for selecting JPG or PNG files
        with wx.FileDialog(self, "Choose an image file", wildcard="Image files (*.jpg;*.png)|*.jpg;*.png",
                        style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as file_dialog:
            # Show dialog and wait for user selection
            if file_dialog.ShowModal() == wx.ID_OK:
                # Get full path of selected file
                file_path = file_dialog.GetPath()
                # Calculate maximum allowed size in bytes from config
                max_size_bytes = config.DEFAULT_MAX_IMAGE_SIZE_MB * 1024 * 1024
                # Check if file exceeds size limit
                if os.path.getsize(file_path) > max_size_bytes:
                    # Log error for oversized file
                    self.logger.error(f"Image file exceeds {config.DEFAULT_MAX_IMAGE_SIZE_MB}MB")
                    # Show error dialog to user
                    wx.MessageBox(f"Image file exceeds {config.DEFAULT_MAX_IMAGE_SIZE_MB}MB.", "Invalid File", wx.OK | wx.ICON_ERROR)
                    # Exit function
                    return
                # Load image using wx.Image to preserve original pixel data
                image = wx.Image(file_path)
                # Validate image loaded successfully
                if not image.IsOk():
                    # Log failure to load image
                    self.logger.error(f"Failed to load image: {file_path}")
                    # Show error dialog
                    wx.MessageBox("Failed to load image. Please select a valid JPG or PNG file.", "Invalid File", wx.OK | wx.ICON_ERROR)
                    # Exit function
                    return
                # Get maximum allowed dimensions from config
                max_width, max_height = config.DEFAULT_MAX_IMAGE_DIMENSIONS
                # Check if image exceeds dimension limits
                if image.GetWidth() > max_width or image.GetHeight() > max_height:
                    # Log dimension error
                    self.logger.error(f"Image dimensions exceed {max_width}x{max_height}px")
                    # Show error dialog
                    wx.MessageBox(f"Image dimensions exceed {max_width}x{max_height}px.", "Invalid File", wx.OK | wx.ICON_ERROR)
                    # Exit function
                    return
                # Store loaded image (original size, no scaling)
                self.overlay_image = image
                # Store original file path for later export
                self.overlay_path = file_path
                # Enable overlay visibility
                self.overlay_visible = True
                # Set overlay position to top-left corner
                self.overlay_pos = (0, 0)
                # Log successful load with dimensions
                self.logger.info(f"Loaded overlay image: {file_path}, dimensions {image.GetWidth()}x{image.GetHeight()}px")
                # Refresh canvas to display the image
                self.canvas.Refresh()
            else:
                # Log when user cancels the dialog
                self.logger.debug("Load Overlay cancelled")

    def toggle_overlay(self, event):
        """
        @brief Handle the Show/Hide Overlay button click to toggle overlay visibility.
        @details Toggles the overlay_visible flag and refreshes the canvas to show/hide the image.
        @param event Button event.
        """
        # Check if an overlay image is loaded
        if self.overlay_image is None:
            # Log warning if no image is available
            self.logger.warning("No overlay image loaded to show/hide")
            # Show warning dialog
            wx.MessageBox("No overlay image loaded.", "No Image", wx.OK | wx.ICON_WARNING)
            # Exit function
            return
        # Toggle visibility state
        self.overlay_visible = not self.overlay_visible
        # Log current visibility state
        self.logger.debug(f"Overlay visibility: {'On' if self.overlay_visible else 'Off'}")
        # Refresh canvas to apply visibility change
        self.canvas.Refresh()

    def update_opacity(self, event):
        """
        @brief Handle the opacity slider change and refresh canvas.
        @details Retrieves the current opacity value from the slider and forces a full canvas repaint
                to apply the updated alpha transparency to rendered elements (e.g., zones, heaters, sinks).
        @param event Slider event.
        """
        # Log current slider value (0-100) for debugging
        self.logger.debug(f"Opacity set to {self.opacity_slider.GetValue()}")
        # Force full repaint of canvas to update transparency of elements
        # True parameter erases background before repaint
        self.canvas.Refresh(True)

    def on_layer_choice(self, event):
        """
        @brief Handle the layer filter selection change.
        @details Updates the layer filter from the choice selection, logs the change, and refreshes the canvas to apply the filter.
        @param event Choice event.
        """
        # Get the currently selected layer filter string
        selected_filter = self.layer_choice.GetStringSelection()
        # Log the new filter selection
        self.logger.debug(f"Layer filter set to {selected_filter}")
        # Refresh canvas to immediately apply the new layer filter
        self.canvas.Refresh()

    def new_project(self, event):
        """
        @brief Handle the New button click to reset the editor to its default state.
        @details Clears all drawn elements, resets overlay and PCB parameters to defaults from config.py,
                updates the interface (canvas, parameter fields, layer filter), and refreshes the display.
                Ensures canvas is refreshed to remove overlay visually.
        @param event Button event.
        """
        # Clear current polygon being drawn
        self.current_element = []
        # Clear all completed elements
        self.elements = []
        # Reset element counters for all types
        self.element_counters = {"Zone": 0, "Heater": 0, "Via": 0, "Sink": 0}
        # Log element cleanup
        self.logger.debug("Cleared all drawn elements and counters")

        # Reset overlay image and related properties
        self.overlay_image = None
        self.overlay_path = None
        self.overlay_visible = False
        # Log overlay reset
        self.logger.debug("Reset overlay to default state")
        # Force canvas repaint to remove any visible overlay
        self.canvas.Refresh()

        # Reset zoom and pan to default values
        self.zoom_factor = 1.0
        self.pan_x = 0
        self.pan_y = 0
        self.is_panning = False
        self.pan_start_x = 0
        self.pan_start_y = 0
        # Log view reset
        self.logger.debug("Reset zoom and pan to default values")

        # Reset PCB corner points and scale factor
        self.pcb_corners = []
        self.scale_factor = 1.0
        # Log PCB setup reset
        self.logger.debug("Reset PCB corners and scale factor")

        # Reset tool-specific temporary states
        self.duplicate_element = None
        self.preview_position = None
        self.moving_element = None
        self.move_start_pos = None
        self.editing_element = None
        self.editing_vertex_index = None
        self.original_element = None
        # Log tool state reset
        self.logger.debug("Reset duplicate, move, and edit states")

        # Reset all PCB parameters to defaults from config.py
        self.pcb_params = {
            "width": config.DEFAULT_PCB_WIDTH,
            "height": config.DEFAULT_PCB_HEIGHT,
            "thickness": config.DEFAULT_PCB_THICKNESS,
            "layers": config.DEFAULT_PCB_LAYERS,
            "copper_thickness": config.DEFAULT_COPPER_THICKNESS_PER_LAYER,
            "copper_conductivity": config.DEFAULT_COPPER_CONDUCTIVITY,
            "fr4_conductivity": config.DEFAULT_FR4_CONDUCTIVITY,
            "ambient_temp_global": config.DEFAULT_AMBIENT_TEMP_GLOBAL,
            "convection_h_global": config.DEFAULT_CONVECTION_H_GLOBAL,
            "convection_h_top": config.DEFAULT_CONVECTION_H_TOP,
            "convection_h_bottom": config.DEFAULT_CONVECTION_H_BOTTOM,
            "grid_dx": config.DEFAULT_SOLVER_GRID_SIZE_MM,
            "grid_dy": config.DEFAULT_SOLVER_GRID_SIZE_MM,
            "tolerance": config.DEFAULT_FDM_TOLERANCE_CONVERGENCE,
            "max_iterations": config.DEFAULT_FDM_MAX_ITERATIONS,
            "over_relaxation_factor": config.DEFAULT_FDM_OVER_RELAXATION_FACTOR
        }
        # Log parameter reset
        self.logger.debug("Reset PCB parameters to default values from config")

        # Reset solver mode to default from config
        self.solver_mode = config.DEFAULT_SOLVER_METHOD
        # Log solver mode reset
        self.logger.debug(f"Reset solver mode to: {self.solver_mode}")

        # Clear last simulation results
        self.last_results = None
        self.last_solver_mode = None
        # Disable results button if it exists
        if hasattr(self, 'results_btn'):
            self.results_btn.Enable(False)
        # Log results cleanup
        self.logger.debug("Reset last simulation results and disabled Simulation Results button")

        # Update all parameter input fields with default values
        for key, entry in self.global_entries.items():
            if key == 'layers':
                # Set layer choice selection (0-based index)
                entry.SetSelection(self.pcb_params[key] - 1)
            elif key == 'max_iterations':
                # Set integer spin control
                entry.SetValue(int(self.pcb_params[key]))
            else:
                # Set float spin control
                entry.SetValue(self.pcb_params[key])
        # Log UI field update
        self.logger.debug("Updated global_entries with default values")

        # Update layer filter dropdown based on new layer count
        max_layers = int(self.pcb_params['layers'])
        layer_choices = ["All", f"1+{max_layers}"]
        if max_layers > 2:
            layer_choices.append("Inner")
        layer_choices.extend([str(i) for i in range(1, max_layers + 1)])
        current_selection = self.layer_choice.GetSelection()
        self.layer_choice.SetItems(layer_choices)
        self.layer_choice.SetSelection(min(current_selection, len(layer_choices) - 1))
        # Log layer filter update
        self.logger.debug(f"Updated layer filter dropdown to reflect {max_layers} layers")

        # Update solver mode dropdown and trigger field enable/disable
        self.solver_mode_choice.SetSelection(self.solver_mode_choice.GetStrings().index(self.solver_mode))
        self.on_solver_mode_change(None)
        # Log solver UI update
        self.logger.debug("Updated solver mode choice and field states")

        # Refresh canvas to clear all drawings
        self.canvas.Refresh()
        # Restore default instruction label
        self.title_label.SetLabel("Select a tool from the left panel to start drawing")
        # Log final UI refresh
        self.logger.debug("Refreshed canvas and title label")

        # Show confirmation dialog
        wx.MessageBox("New project created with default settings.", "New Project", wx.OK | wx.ICON_INFORMATION)

        # Allow event to propagate
        event.Skip()

    def import_json(self, event):
        """
        @brief Handle the Import JSON button click to load PCB configuration from a JSON file.
        @details Loads pcb_parameters (now including nested pcb_corners), solver (solver-specific params including mode),
                elements, and overlay from a JSON file.
                Extracts pcb_corners from pcb_parameters to self.pcb_corners after update.
                Updates self.pcb_params with both pcb_parameters and solver fields, sets self.solver_mode from solver['solver_mode'].
                Updates GUI elements from loaded params (handle different widget types). Logs warnings if values don't match choices.
                Triggers on_solver_mode_change after setting solver_mode to update enables/disables.
                Shows message on success or warning on errors.
        @param event Button event.
        """
        # Create file dialog for selecting JSON files
        with wx.FileDialog(self, "Select PCB configuration JSON file",
                        wildcard="JSON files (*.json)|*.json",
                        style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as file_dialog:
            # Show dialog and wait for user selection
            if file_dialog.ShowModal() == wx.ID_OK:
                # Get full path of selected file
                json_path = file_dialog.GetPath()
                # Log import attempt
                self.logger.info(f"Loading configuration from {json_path}")
                try:
                    # Open and read JSON file
                    with open(json_path, 'r', encoding='utf-8') as f:
                        loaded_data = json.load(f)

                    # Extract metadata for logging
                    metadata = loaded_data.get("metadata", {})
                    self.logger.info(f"Loaded config: {metadata.get('product_name', 'Unknown')} v{metadata.get('app_version', 'Unknown')}")

                    # Get PCB parameters from JSON
                    loaded_pcb_params = loaded_data.get("pcb_parameters", {})
                    self.logger.debug(f"Loaded PCB params: {loaded_pcb_params}")

                    # Extract pcb_corners before updating pcb_params
                    self.pcb_corners = loaded_pcb_params.pop('pcb_corners', [])
                    self.logger.debug(f"Extracted pcb_corners: {self.pcb_corners}")

                    # Update main PCB parameters
                    self.pcb_params.update(loaded_pcb_params)

                    # Get solver-specific parameters
                    loaded_solver_params = loaded_data.get("solver", {})
                    # Merge solver params into pcb_params for consistency
                    self.pcb_params.update(loaded_solver_params)
                    self.logger.debug(f"Updated solver params: {loaded_solver_params}")

                    # Extract solver mode
                    loaded_solver_mode = loaded_solver_params.get("solver_mode", config.DEFAULT_SOLVER_METHOD)
                    self.solver_mode = loaded_solver_mode

                    # Update solver mode dropdown if it exists
                    if hasattr(self, 'solver_mode_choice'):
                        choice_count = self.solver_mode_choice.GetCount()
                        choices = [self.solver_mode_choice.GetString(i) for i in range(choice_count)]
                        try:
                            # Set selection to loaded mode
                            index = choices.index(self.solver_mode)
                            self.solver_mode_choice.SetSelection(index)
                            # Create event to trigger UI updates
                            mode_event = wx.CommandEvent(wx.wxEVT_CHOICE, self.solver_mode_choice.GetId())
                            self.on_solver_mode_change(mode_event)
                            self.logger.debug(f"Set solver_mode_choice to index {index} ({self.solver_mode}) and triggered mode change")
                        except ValueError:
                            # Fallback to default if mode invalid
                            self.logger.warning(f"Invalid solver mode in JSON: {self.solver_mode}. Using default.")
                            self.solver_mode = config.DEFAULT_SOLVER_METHOD
                            default_index = choices.index(config.DEFAULT_SOLVER_METHOD)
                            self.solver_mode_choice.SetSelection(default_index)
                            mode_event = wx.CommandEvent(wx.wxEVT_CHOICE, self.solver_mode_choice.GetId())
                            self.on_solver_mode_change(mode_event)

                    # Load elements list
                    self.elements = loaded_data.get("elements", [])
                    # Update element counters based on loaded IDs
                    for elem in self.elements:
                        if elem["type"] in self.element_counters:
                            self.element_counters[elem["type"]] = max(self.element_counters[elem["type"]], int(elem.get("id", 0)) + 1)
                    self.logger.debug(f"Updated elements: {len(self.elements)} total, counters: {self.element_counters}")

                    # Handle overlay image
                    overlay_info = loaded_data.get("overlay", {})
                    overlay_file = overlay_info.get("file")
                    if overlay_file:
                        # Build expected path in DEFAULT_OUTPUT_DIR
                        overlay_path = os.path.join(config.DEFAULT_OUTPUT_DIR, overlay_file)
                        if os.path.exists(overlay_path):
                            try:
                                # Load image
                                image = wx.Image(overlay_path)
                                if image.IsOk():
                                    # Validate dimensions
                                    if image.GetWidth() > config.DEFAULT_MAX_IMAGE_DIMENSIONS[0] or \
                                    image.GetHeight() > config.DEFAULT_MAX_IMAGE_DIMENSIONS[1]:
                                        raise ValueError(f"Image dimensions exceed {config.DEFAULT_MAX_IMAGE_DIMENSIONS[0]}x{config.DEFAULT_MAX_IMAGE_DIMENSIONS[1]}px")
                                    # Validate size
                                    if os.path.getsize(overlay_path) > config.DEFAULT_MAX_IMAGE_SIZE_MB * 1024 * 1024:
                                        raise ValueError(f"Image file exceeds {config.DEFAULT_MAX_IMAGE_SIZE_MB}MB")
                                    # Store image and properties
                                    self.overlay_image = image
                                    self.overlay_path = overlay_path
                                    self.overlay_visible = True
                                    self.overlay_pos = (0, 0)
                                    self.logger.info(f"Loaded overlay image from {overlay_path}, dimensions {image.GetWidth()}x{image.GetHeight()}px")
                                else:
                                    raise ValueError("Invalid image")
                            except Exception as e:
                                self.logger.warning(f"Failed to load overlay from {overlay_path}: {str(e)}")
                        else:
                            # Prompt user to locate missing overlay
                            with wx.FileDialog(self, f"Overlay image {overlay_file} not found. Choose a replacement image",
                                            wildcard="Image files (*.jpg;*.png)|*.jpg;*.png",
                                            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as file_dialog_overlay:
                                if file_dialog_overlay.ShowModal() == wx.ID_OK:
                                    new_path = file_dialog_overlay.GetPath()
                                    try:
                                        if os.path.getsize(new_path) > config.DEFAULT_MAX_IMAGE_SIZE_MB * 1024 * 1024:
                                            raise ValueError(f"Image file exceeds {config.DEFAULT_MAX_IMAGE_SIZE_MB}MB")
                                        image = wx.Image(new_path)
                                        if not image.IsOk():
                                            raise ValueError("Invalid image")
                                        if image.GetWidth() > config.DEFAULT_MAX_IMAGE_DIMENSIONS[0] or \
                                        image.GetHeight() > config.DEFAULT_MAX_IMAGE_DIMENSIONS[1]:
                                            raise ValueError(f"Image dimensions exceed {config.DEFAULT_MAX_IMAGE_DIMENSIONS[0]}x{config.DEFAULT_MAX_IMAGE_DIMENSIONS[1]}px")
                                        # Copy to DEFAULT_OUTPUT_DIR
                                        os.makedirs(config.DEFAULT_OUTPUT_DIR, exist_ok=True)
                                        shutil.copy2(new_path, overlay_path)
                                        # Store new image
                                        self.overlay_image = image
                                        self.overlay_path = overlay_path
                                        self.overlay_visible = True
                                        self.overlay_pos = (0, 0)
                                        self.logger.info(f"Copied and loaded overlay image to {overlay_path}, dimensions {image.GetWidth()}x{image.GetHeight()}px")
                                    except Exception as e:
                                        self.logger.warning(f"Failed to load or copy overlay: {str(e)}")
                                        wx.MessageBox(f"Failed to load overlay image: {str(e)}", "Import Warning", wx.OK | wx.ICON_WARNING)
                                else:
                                    self.logger.debug("Overlay image selection cancelled")
                    else:
                        # No overlay in JSON
                        self.overlay_image = None
                        self.overlay_path = None
                        self.overlay_visible = False

                    # Update all GUI parameter fields
                    for key, entry in self.global_entries.items():
                        if key in self.pcb_params:
                            value = self.pcb_params[key]
                            try:
                                if isinstance(entry, wx.Choice):
                                    choice_count = entry.GetCount()
                                    choices = [entry.GetString(i) for i in range(choice_count)]
                                    str_value = str(value)
                                    if str_value in choices:
                                        index = choices.index(str_value)
                                        entry.SetSelection(index)
                                        self.logger.debug(f"Set Choice {key} to index {index} ({str_value})")
                                    else:
                                        self.logger.warning(f"Invalid choice value for {key}: {str_value}. Available: {choices}. Skipping.")
                                elif isinstance(entry, wx.SpinCtrl):
                                    entry.SetValue(int(value))
                                    self.logger.debug(f"Set SpinCtrl {key} to {int(value)}")
                                elif isinstance(entry, wx.SpinCtrlDouble):
                                    entry.SetValue(float(value))
                                    self.logger.debug(f"Set SpinCtrlDouble {key} to {float(value)}")
                                else:
                                    entry.SetValue(float(value))
                                    self.logger.debug(f"Set other widget {key} to {float(value)}")
                            except (ValueError, AttributeError) as e:
                                self.logger.warning(f"Failed to set value for {key}: {str(e)}")

                    # Update scale factor if corners exist
                    if self.pcb_corners and len(self.pcb_corners) == 2:
                        pcb_width_px = abs(self.pcb_corners[1][0] - self.pcb_corners[0][0])
                        self.scale_factor = self.pcb_params['width'] / pcb_width_px if pcb_width_px > 0 else 1.0

                    # Refresh UI
                    self.update_element_table()
                    self.canvas.Refresh()
                    # Log success
                    self.logger.info(f"Imported configuration from {json_path}")
                    wx.MessageBox(f"Imported {os.path.basename(json_path)}", "Import Successful", wx.OK | wx.ICON_INFORMATION)

                except json.JSONDecodeError as e:
                    # Log JSON parsing error
                    self.logger.error(f"Invalid JSON in {json_path}: {str(e)}")
                    wx.MessageBox(f"Invalid JSON file: {str(e)}", "Import Error", wx.OK | wx.ICON_ERROR)
                except Exception as e:
                    # Log general import error
                    self.logger.error(f"Failed to import JSON: {str(e)}")
                    wx.MessageBox(f"Failed to import: {str(e)}", "Import Error", wx.OK | wx.ICON_ERROR)
            else:
                # Log cancellation
                self.logger.debug("JSON import cancelled by user")

        # Allow event propagation
        event.Skip()

    def export_json(self, event):
        """
        @brief Handle the Export JSON button click to save PCB configuration and overlay to a JSON file.
        @details Opens a file dialog to let the user choose the export directory and filename, suggesting
                pcb_YYYYMMDDHHMMSS.json. Exports pcb_parameters (PCB params including pcb_corners nested),
                a 'solver' object (solver-specific params), elements, and overlay info to the chosen JSON file.
                The 'solver' object includes: grid_dx, grid_dy, tolerance, max_iterations, over_relaxation_factor,
                model_mode, hotspot_threshold, and solver_mode.
                Copies the overlay image (if loaded) to the same directory with a matching filename.
                Shows a message box with success or error. Ensures a new file is created with current timestamp.
        @param event Button event.
        """
        # Generate timestamp for suggested filename
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        suggested_filename = f"pcb_{timestamp}.json"
        output_dir = config.DEFAULT_OUTPUT_DIR
        default_path = os.path.join(output_dir, suggested_filename)

        # Create file dialog for saving
        with wx.FileDialog(
            self, "Save PCB Configuration", defaultDir=output_dir, defaultFile=suggested_filename,
            wildcard="JSON files (*.json)|*.json", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
        ) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                self.logger.info("Export JSON cancelled by user")
                return
            # Get the selected path and filename
            output_path = dlg.GetPath()
            self.logger.info(f"User selected export path: {output_path}")

        # Extract directory and filename from the chosen path
        output_dir = os.path.dirname(output_path)
        filename = os.path.basename(output_path)
        base_name = os.path.splitext(filename)[0]  # Remove .json extension

        # Define keys to include in pcb_parameters
        pcb_keys = [
            'width', 'height', 'thickness', 'layers', 'copper_thickness', 'copper_conductivity', 'fr4_conductivity',
            'ambient_temp_global', 'convection_h_global', 'convection_h_top', 'convection_h_bottom',
            'pcb_corners'
        ]
        # Filter and include pcb_corners
        filtered_pcb_params = {k: self.pcb_params[k] for k in pcb_keys if k in self.pcb_params}
        filtered_pcb_params['pcb_corners'] = self.pcb_corners
        self.logger.debug(f"Filtered PCB params for export (with corners): {filtered_pcb_params}")

        # Define solver-specific keys
        solver_keys = [
            'grid_dx', 'grid_dy', 'tolerance', 'max_iterations', 'over_relaxation_factor',
            'model_mode', 'boundary_mode', 'hotspot_threshold'
        ]
        # Build solver object
        solver_params = {k: self.pcb_params.get(k) for k in solver_keys if k in self.pcb_params}
        solver_params['solver_mode'] = self.solver_mode
        self.logger.debug(f"Solver params for export: {solver_params}")

        # Build main export dictionary
        export_data = {
            "metadata": {
                "app_version": config.APP_VERSION,
                "product_name": config.PRODUCT_NAME,
                "timestamp": timestamp,
                "company_name": config.COMPANY_NAME,
                "author_name": config.AUTHOR_NAME,
                "format_version": "1.0"
            },
            "pcb_parameters": filtered_pcb_params,
            "solver": solver_params,
            "elements": self.elements,
            "overlay": {}
        }

        # Handle overlay image export
        if self.overlay_image is not None and self.overlay_path is not None:
            overlay_filename = f"{base_name}.png"  # Match the base name with .png extension
            overlay_output_path = os.path.join(output_dir, overlay_filename)
            try:
                # Ensure output directory exists
                os.makedirs(output_dir, exist_ok=True)
                # Copy overlay image
                shutil.copy2(self.overlay_path, overlay_output_path)
                self.logger.info(f"Copied overlay image to {overlay_output_path}")
                # Add to export data
                export_data["overlay"] = {"file": overlay_filename}
            except Exception as e:
                self.logger.error(f"Failed to copy overlay image: {str(e)}")
                wx.MessageBox(f"Failed to copy overlay image: {str(e)}", "Export Error", wx.OK | wx.ICON_ERROR)
                return

        # Write JSON file
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(export_data, f, indent=2)
            # Log new fields
            if 'model_mode' in solver_params:
                self.logger.debug(f"Exported model_mode: {solver_params['model_mode']}")
            if 'hotspot_threshold' in solver_params:
                self.logger.debug(f"Exported hotspot_threshold: {solver_params['hotspot_threshold']}")
            # Log success
            self.logger.info(f"Successfully exported configuration to {output_path}")
            wx.MessageBox(f"Exported to {filename}", "Export Successful", wx.OK | wx.ICON_INFORMATION)
        except Exception as e:
            self.logger.error(f"Failed to export JSON: {str(e)}")
            wx.MessageBox(f"Failed to export JSON: {str(e)}", "Export Error", wx.OK | wx.ICON_ERROR)
            return

        # Refresh canvas for consistency
        self.canvas.Refresh()
        # Allow event propagation
        event.Skip()
    
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
                # Log attempt for debugging
                self.logger.info(f"Attempting to delete temp directory: {config.DEFAULT_TEMP_DIR}")
                # Recursively delete the entire directory and its contents
                shutil.rmtree(config.DEFAULT_TEMP_DIR)
                # Log successful cleanup
                self.logger.info(f"Successfully deleted temp directory: {config.DEFAULT_TEMP_DIR}")
        except Exception as e:
            # Log warning if deletion fails (e.g., file in use, permission denied)
            self.logger.warning(f"Failed to delete temp directory: {e}")

    def on_close_main(self, event):
        """
        @brief Handle main window close.
        @details Deletes temp directory, flushes logs, and exits app.
                 Calls delete synchronously after cleanup to ensure execution before app shutdown.
                 Adds logging.shutdown() to force log flush before exit.
        """
        # Log that close handler was called for debugging
        self.logger.info("Main window close handler called - starting cleanup")
        # Close all matplotlib figures to release any file handles (important for FDM plots and GIFs)
        self.logger.info("Closing all Matplotlib figures...")
        import matplotlib.pyplot as plt
        plt.close('all')
        self.logger.info("Matplotlib figures closed.")
        # Delete temp dir synchronously (before Destroy, to avoid event queue issues)
        self._delete_temp_dir()
        # Destroy dialog first
        self.Destroy()
        # Force log flush and close handlers to ensure logs are written before app exit
        # import logging
        # logging.shutdown()
        # Check if this is the last window and destroy the app
        if not wx.GetApp().GetTopWindow():
            # Exit app (no WakeUpMainThread needed)
            wx.GetApp().Destroy()

    def run_simulation(self, event):
        """
        @brief Run simulation with immediate loading UI
        @details Sync all parameters before copying to ensure
                latest control values are used by solver
        @param event Button event
        """
        try:
            # Step 1 - Sync all parameters from controls before copying
            # This ensures any recent changes are captured
            self.update_global_params()
            self.logger.debug("Parameters synchronized before simulation")

            # Step 2 - Deep copy mutable data to avoid side effects
            pcb_params = copy.deepcopy(self.pcb_params)
            pcb_params['pcb_corners'] = self.pcb_corners
            json_config = {
                'pcb_params': pcb_params,
                'solver_mode': self.solver_mode,
                'elements': copy.deepcopy(self.elements)
            }
            self.logger.debug(f"Prepared json_config for mode: {self.solver_mode}")

            # Step 3 - Validate configuration
            ThermalSolver.validate_config(json_config)

            # Delete Temporary folder
            self._delete_temp_dir()

            # Step 4 - Get dialog class name
            dialog_class_name = ThermalSolver.get_results_dialog_class_name(self.solver_mode)
            self.logger.debug(f"Selected dialog class: {dialog_class_name}")

            # Step 5 - Create dialog instance (empty)
            if dialog_class_name == 'AnalyticResultsDialog':
                dialog = AnalyticResultsDialog(None, {}, logger=self.logger)
            elif dialog_class_name == 'LumpedResultsDialog':
                dialog = LumpedResultsDialog(None, {}, logger=self.logger)
            elif dialog_class_name == 'FDMResultsDialog':
                dialog = FDMResultsDialog(None, {}, T_history=[], zone_polys={}, logger=self.logger)
            else:
                # Fallback to Analytic
                dialog = AnalyticResultsDialog(None, {})
                self.logger.warning(f"Unknown dialog class {dialog_class_name}, using Analytic")

            # Step 6 - Show loading screen (full screen, no Save button)
            dialog.Show()
            dialog.show_loading()  # Show loading message
            wx.Yield()  # Force UI update
            dialog.Update()  # Force redraw
            dialog.Refresh()  # Ensure visible
            self.logger.debug("Loading screen displayed")

            # Step 7 - Run simulation in main thread
            results = ThermalSolver.run_simulation(json_config, logger=self.logger)
            self.logger.info(f"Simulation completed: {self.solver_mode}")

            # Step 8 - Set results in dialog
            dialog.results = results
            # self.logger.debug(f"Dialog results keys: {list(results.keys()) if isinstance(results, dict) else 'None'}")

            # Step 9 - Show Save button for all solvers
            if hasattr(dialog, 'save_btn'):
                dialog.save_btn.Show()
                dialog.btn_panel.Layout()
                self.logger.debug("Save button shown after simulation")

            # Step 10 - Update dialog content with results
            if dialog_class_name == 'FDMResultsDialog':
                # FDM: Create notebook with multiple tabs
                dialog.T_history = results.get('T_history', [])
                dialog.zone_polys = {}  # Initialize empty

                # Clear loading panel
                dialog.content_sizer.Clear(True)

                # Create notebook with tabs
                dialog.notebook = wx.Notebook(dialog)
                dialog.notebook.AddPage(dialog.create_summary_tab(), "Summary")
                dialog.notebook.AddPage(dialog.create_heatmap_tab(), "Heatmap")
                dialog.notebook.AddPage(dialog.create_evolution_tab(), "Evolution")
                dialog.notebook.AddPage(dialog.create_hotspots_tab(), "Hotspots")
                dialog.notebook.AddPage(dialog.create_gradients_tab(), "Gradients")
                dialog.notebook.AddPage(dialog.create_zones_tab(), "Zones")

                # Add notebook to content sizer
                dialog.content_sizer.Add(dialog.notebook, 1, wx.EXPAND | wx.ALL, 5)
            else:
                # Analytic/Lumped: Use single content sizer
                content_sizer = dialog._create_content_sizer()
                dialog.content_sizer.Clear(True)
                dialog.content_sizer.Add(content_sizer, 1, wx.ALL | wx.EXPAND, 5)

            # Final layout update
            dialog.Layout()
            dialog.Update()
            self.logger.debug("Dialog updated with results and content")

            # Step 11 - Store results for "Simulation Results" button
            self.last_results = results
            self.last_solver_mode = self.solver_mode
            self.results_btn.Enable(True)
            self.logger.debug("Results stored and Results button enabled")

        except Exception as e:
            # Error handling
            self.logger.error(f"Simulation failed: {str(e)}")
            if 'dialog' in locals():
                dialog.Destroy()
            wx.MessageBox(f"Simulation failed: {str(e)}", "Error", wx.OK | wx.ICON_ERROR)

        finally:
            event.Skip()


if __name__ == "__main__":
    # Create a wxPython application instance
    app = wx.App(False)
    # Create and show the main editor window
    frame = PCBEditor(None)
    frame.Show()
    # Start the application event loop
    app.MainLoop()