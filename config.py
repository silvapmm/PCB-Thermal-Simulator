# Configuration module for the PCB Thermal Simulator
# This module centralizes constants, defaults, and customizable parameters for the application
# It includes metadata, physical constants, simulation defaults, I/O settings, and UI/visualization styles

# APPLICATION METADATA
# Information about the tool, version, and author for reference and about dialog
APP_VERSION = "1.0.0"
COMPANY_NAME = "PMS Labs"
PRODUCT_NAME = "PCB Thermal Simulator Tool"
AUTHOR_NAME = "Pedro Silva"
FILE_DESCRIPTION = "2D PCB Thermal Simulation Tool."
INTERNAL_NAME = "pts_core_engine"
LEGAL_COPYRIGHT = "(C) 2025 Pedro Silva. All rights reserved."

# PHYSICAL CONSTANTS
# Default material properties (can be overridden in GUI/JSON for custom simulations)
DEFAULT_COPPER_CONDUCTIVITY = 385.0         # Thermal conductivity of Copper (W/mK)
DEFAULT_FR4_CONDUCTIVITY = 0.3              # Default thermal conductivity of FR4 (W/mK)
DEFAULT_VIA_DIAMETER = 0.6                  # Default diameter for vias (mm)
DEFAULT_HEATER_POWER = 5.0                  # Default power for heaters (W)
DEFAULT_SINK_HTC_EXTRA = 20.0               # Default extra heat transfer coefficient for sinks (W/m²K)

# PCB PHYSICAL DIMENSIONS
# Default PCB properties (can be overridden in GUI/JSON for custom simulations)
DEFAULT_PCB_WIDTH = 100.0                   # PCB Width (mm)
DEFAULT_PCB_HEIGHT = 100.0                  # PCB Height (mm)
DEFAULT_PCB_THICKNESS = 1.6                 # PCB Thickness (mm)
DEFAULT_PCB_LAYERS = 2                      # Number of Layers
DEFAULT_PCB_MAX_LAYERS = 8                  # Maximum number of PCB layers allowed
DEFAULT_COPPER_THICKNESS_PER_LAYER = 0.035  # Copper Thickness per Layer (mm)
DEFAULT_VIA_PLATING_THICKNESS = 0.025       # Thickness of via plating (mm)

# Ambient and Convection Conditions (can be overridden in GUI/JSON for custom simulations)
DEFAULT_AMBIENT_TEMP_GLOBAL = 25.0          # Ambient Temperature (°C)
DEFAULT_CONVECTION_H_GLOBAL = 10.0          # Global Convection Coefficient (W/m²K)
DEFAULT_CONVECTION_H_TOP = 10.0             # Top Convection Coefficient (W/m²K)
DEFAULT_CONVECTION_H_BOTTOM = 15.0          # Bottom Convection Coefficient (W/m²K)

# SIMULATION DEFAULT PARAMETERS (Solver)
# Numerical and convergence settings for the solvers (FDM, Lumped, Analytic)
DEFAULT_SOLVER_METHOD = "Analytic"                  # Default simulation mode (Analytic, Lumped, FDM)
DEFAULT_SOLVER_MODEL_MODE = "Basic"                 # Default model mode for Analytic and Lumped
DEFAULT_FDM_BOUNDARY_MODE = "Convective"            # Default boundary mode for FDM
DEFAULT_FDM_MAX_ITERATIONS = 50000                  # Safety limit for non-convergence in iterative solvers
DEFAULT_FDM_TOLERANCE_CONVERGENCE = 0.01            # Stopping criterion (max change in T, °C)
DEFAULT_FDM_OVER_RELAXATION_FACTOR = 1.5            # Omega for Successive Over-Relaxation (SOR) in FDM
DEFAULT_FDM_HOTSPOT_THRESHOLD = 100.0               # Default hotspot threshold for FDM solver (°C)
DEFAULT_SOLVER_GRID_SIZE_MM = 10.0                  # Default grid spacing for FDM (dx/dy in mm)


# I/O PARAMETERS (Input/Output)
# Paths, file names, and output formats for inputs and results
DEFAULT_OUTPUT_DIR = "output"                                   # Directory for output files (heatmaps, GIFs, reports, overlay)
DEFAULT_TEMP_DIR = "temp"                                       # Directory for temporary files (GIFs, PNGs, etc.)
LOG_FILENAME = "thermal_simulator.log"                          # Log file name for application
DEFAULT_PCB_CONFIG_FILENAME = "pcb_{timestamp}.json"            # PCB configuration file name format
DEFAULT_HEATMAP_FILENAME = "heatmap_{timestamp}.png"            # Heatmap file name format
DEFAULT_ITERATION_GIF_ENABLED = True                            # Option to generate GIF of iterations (for FDM)
DEFAULT_ITERATION_GIF_FILENAME = "iterations_{timestamp}.gif"   # GIF file name format
DEFAULT_PCB_OVERLAY_FILENAME = "overlay_{timestamp}.png"        # Default overlay file name format
DEFAULT_MAX_IMAGE_DIMENSIONS = (5000, 5000)                     # Max image dimensions (width, height) in pixels
DEFAULT_MAX_IMAGE_SIZE_MB = 10                                  # Max image file size in megabytes
DEFAULT_PCB_RESULTS_TEXT_FILENAME = "pcb_results_{solver}_{timestamp}.txt"  # Results text file name format
DEFAULT_PCB_RESULTS_CHART_FILENAME = "pcb_chart_{solver}_{timestamp}.png"   # Results chart file name format

# UI AND VISUALIZATION SETTINGS
# Default window sizes and font sizes for dialogs and UI elements
DEFAULT_MAIN_WINDOW_SIZE = (1300, 850)              # Main window dimensions (width, height) in pixels
DEFAULT_RESULTS_WINDOW_SIZE = (900, 600)            # Results window dimensions (width, height) in pixels
DEFAULT_RESULTS_FONT_SIZE = 9                       # Default font size for results text in ResultsDialog
DEFAULT_EXPLANATION_FONT_SIZE = 8                   # Default font size for explanation text in ResultsDialog
DEFAULT_CHART_FONT_SIZE = 8                         # Default font size for text in charts in ResultsDialog
DEFAULT_FONT_FAMILY = 'Arial'                       # Default font family for text and labels


# Visualization styles for charts and heatmaps
DEFAULT_BACKGROUND_COLOR = '#FFFFFF'                # Default background color for charts (white)
DEFAULT_CHART_COLOR_SCHEME = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # Default colors for bar charts
DEFAULT_CHART_Y_MARGIN_PERCENT = 15                 # 0 = no margin, 15 = 15% extra
DEFAULT_VISUALIZATION_TEMP_UNIT = "C"               # Temperature unit for display (C for Celsius, K for Kelvin)
# Heatmap Settings
DEFAULT_VISUALIZATION_HEAT_COLORMAP_STYLE = 'hot'   # Colormap style for matplotlib (e.g., 'hot', 'viridis', 'plasma')
DEFAULT_VISUALIZATION_TEMP_MIN_DEFAULT = 25.0       # Default minimum for color bar (Ambient temperature)
DEFAULT_VISUALIZATION_TEMP_MAX_DEFAULT = 175.0      # Default maximum for color bar (Safety limit)
# Hotspot Settings
DEFAULT_HOTSPOT_COLOR = "#E01EA6"                 # Default hotspot color
# Evolution Tab
DEFAULT_EVOLUTION_CAPTURE_STEP = 10                 # Capture frame every N iterations
DEFAULT_EVOLUTION_GIF_FRAME_DELAY_MS = 200          # ms per frame in GIF (5 FPS)
# Gradients Settings
DEFAULT_VISUALIZATION_GRADIENT_COLORMAP_STYLE = 'viridis'  # Default colormap for gradient plots

# Explanation text for simulation results by solver mode
EXPLANATION_TEXT = {
    "Analytic": (
        "Analytic Notes:\n"
        "- Max Temp: Maximum estimated temperature across the PCB (°C).\n"
        "- Avg Temp: Average estimated temperature across the PCB (°C).\n"
        "- Delta T: Temperature difference from ambient (max_temp - T_amb) (°C).\n"
        "- Effective R Th: Effective thermal resistance of the PCB (K/W).\n"
        "- Total Power: Total power dissipated by heaters (W).\n"
        "- Model Mode: Basic uses full PCB area for uniform dissipation; Advanced calculates real heater areas (Shoelace formula with pixel-to-mm scaling) for localized hotspots.\n"
        "- Analytic Note: Simplified model using closed-form equations, suitable for basic PCB layouts with limited accuracy.\n"
    ),
    "Lumped": (
        "Lumped Notes:\n"
        "- Max Temp: Maximum temperature across the PCB (°C).\n"
        "- Avg Temp: Average temperature per layer (°C).\n"
        "- T Top: Temperature at the top layer (°C).\n"
        "- T Bottom: Temperature at the bottom layer (°C).\n"
        "- T Inner: Temperature at inner layers (if applicable) (°C).\n"
        "- Effective R Th: Effective thermal resistance of the PCB (K/W).\n"
        "- R Th Copper Eff: Effective thermal resistance of copper layers (K/W).\n"
        "- R Th Vias Eff: Effective thermal resistance of vias (K/W).\n"
        "- Heat Flux: Heat dissipation per unit area (W/m²).\n"
        "- Delta T: Temperature difference from ambient (max_temp - T_amb) (°C).\n"
        "- Total Power: Total power dissipated by heaters (W).\n"
        "- Model Mode: Basic uses full PCB area; Advanced uses real areas of heaters and vias (Shoelace for polygons, approx for vias) for better multi-layer accuracy.\n"
        "- Lumped Note: Approximation using a network of thermal resistances, suitable for multi-layer PCBs.\n"
    ),
    "FDM": (
        "FDM Notes:\n"
        "- Max Temp: Highest temperature in the grid (°C).\n"
        "- Avg Temp: Global average temperature across all nodes (°C).\n"
        "- Delta T: Maximum temperature rise above ambient (°C).\n"
        "- T Top: Average temperature on the top layer (°C).\n"
        "- T Bottom: Average temperature on the bottom layer (°C).\n"
        "- T Inner: Average temperature on inner layers (if >2 layers) (°C).\n"
        "- Effective R Th: Equivalent thermal resistance (K/W).\n"
        "- Total Power: Sum of all heater power inputs (W).\n"
        "- Effective Area: Total PCB area used in calculations (m²).\n"
        "- Iterations: Number of SOR iterations to converge.\n"
        "- Converged: Whether the simulation reached tolerance.\n"
        "- Max Delta: Largest temperature change in final iteration (°C).\n"
        "- FDM Note: Finite Difference Method with SOR solver - high spatial resolution, suitable for complex multi-layer PCBs.\n"
    )
}