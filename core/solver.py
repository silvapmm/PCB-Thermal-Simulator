# Central module for thermal simulation solvers
# Routes calls to static methods without instantiating solver classes

from core.analytic_solver import AnalyticSolver
from core.lumped_solver import LumpedSolver
from core.fdm_solver import FDMSolver


class ThermalSolver:
    """
    @brief Central dispatcher for thermal solvers.
    @details Routes json_config to the correct static solve method.
             No instantiation required.
    """
    # Map solver mode to static simulate method
    _SOLVER_MAP = {
        'analytic': AnalyticSolver.simulate,
        'lumped': LumpedSolver.simulate,
        'fdm': FDMSolver.simulate
    }

    # Map solver mode to dialog class name
    _DIALOG_MAP = {
        'analytic': 'AnalyticResultsDialog',
        'lumped': 'LumpedResultsDialog',
        'fdm': 'FDMResultsDialog'
    }

    @staticmethod
    def validate_config(json_config):
        """
        @brief Validate input configuration.
        @param json_config Input dict.
        @raise ValueError if invalid.
        """
        # Check required top-level keys
        if 'pcb_params' not in json_config:
            raise ValueError("json_config must contain 'pcb_params'")
        if 'elements' not in json_config:
            raise ValueError("json_config must contain 'elements'")

        # Get pcb parameters
        pcb = json_config['pcb_params']

        # List of required pcb_params keys
        required = [
            'width', 'height', 'thickness', 'layers', 'copper_thickness',
            'copper_conductivity', 'fr4_conductivity', 'ambient_temp_global',
            'convection_h_global', 'convection_h_top', 'convection_h_bottom'
        ]

        # Validate each required key
        for key in required:
            if key not in pcb:
                raise ValueError(f"Missing required key in pcb_params: {key}")

        # Get and validate solver mode
        mode = json_config.get('solver_mode', 'fdm').lower()
        if mode not in ThermalSolver._SOLVER_MAP:
            raise ValueError(f"Invalid solver_mode: {mode}")

        # FDM-specific validation
        if mode == 'fdm':
            if 'grid_dx' not in pcb or 'grid_dy' not in pcb:
                raise ValueError("FDM requires 'grid_dx' and 'grid_dy'")

    @staticmethod
    def get_results_dialog_class_name(mode):
        """
        @brief Return dialog class name for given mode.
        @param mode Solver mode string.
        @return String with dialog class name.
        """
        # Return dialog class name or fallback to Analytic
        return ThermalSolver._DIALOG_MAP.get(mode.lower(), 'AnalyticResultsDialog')

    @staticmethod
    def run_simulation(json_config, logger=None):
        """
        @brief Run simulation by calling the correct static method.
        @param json_config Full input configuration.
        @param logger Optional logger.
        @return Results dictionary.
        """
        # Validate configuration first
        ThermalSolver.validate_config(json_config)

        # Get solver mode (default: fdm)
        mode = json_config.get('solver_mode', 'fdm').lower()

        # Get the correct simulate function
        solve_func = ThermalSolver._SOLVER_MAP[mode]

        # Run simulation and return results
        return solve_func(json_config, logger=logger)