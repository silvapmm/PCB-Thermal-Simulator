# Main entry point for the PCB Thermal Simulator project
# This script initializes and runs the graphical editor using wxPython
# All other modules are imported from here
# Import for console clearing
import os
# Import for cross-platform console clearing
import platform
# Import for the graphical user interface framework
import wx
# Import the editor class from ui module
from ui.editor import PCBEditor
# Import Config File
from config import PRODUCT_NAME, APP_VERSION

if __name__ == "__main__":
    # Clear the console screen (Windows: 'cls', Unix-based: 'clear')
    os.system('cls' if platform.system() == 'Windows' else 'clear')
    # Initialize wxPython application
    app = wx.App(False)
    # Create the main window instance with no parent
    frame = PCBEditor(None)
    # Set the window title using config values
    frame.SetTitle(f"{PRODUCT_NAME} V{APP_VERSION}")
    # Show the main window
    frame.Show()
    # Start the wxPython main event loop to run the UI
    app.MainLoop()