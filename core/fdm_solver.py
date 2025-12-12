import numpy as np  # Import NumPy for fast array operations, grids, and numerical computations
import logging      # Import for fallback logger when no external logger is provided
import math         # Import for mathematical constants (e.g., pi for via area) and functions
import wx           # Import wxPython for GUI components (windows, dialogs, canvas, controls)
import wx.adv       # Import wx.adv for advanced widgets (e.g., animation, taskbar icon, etc.)
import matplotlib.pyplot as plt  # Import Matplotlib for generating charts and plots
from matplotlib.ticker import FixedLocator  # Import FixedLocator to control tick positions on axes
from matplotlib.animation import FuncAnimation # Import FuncAnimation to create animation gif
import shutil # Import for high-level file operations (e.g., copying files or directories)
from io import BytesIO  # Import BytesIO for in-memory binary streams (e.g., saving images to memory)
import os           # Import for operating system interface (file paths, directories, existence checks)
import datetime     # Import for handling timestamps and date/time formatting
import config       # Import local config module containing constants (sizes, colors, paths, defaults)

# Numba import and sor_update function (placed here for global access)
try:
    from numba import njit, float64, int64
    USE_NUMBA = True
except ImportError:
    USE_NUMBA = False
    logging.getLogger(__name__).warning("Numba not installed. Falling back to pure Python (slower). Install with 'pip install numba' for speedup.")
    # No-op decorator for fallback
    def njit(func):
        return func

@njit  # Effective only if Numba available; with explicit signature for typing safety
def sor_update(T: float64[:,:,:], heater_sources: float64[:,:,:], h_arr: float64[:,:,:], 
               cond_lat_x: float64[:,:,:], cond_lat_y: float64[:,:,:], cond_vert: float64[:,:,:], 
               ambient: float64, omega: float64, tol: float64, max_iter: int64, 
               boundary_mode_int: int64, h_side: float64, dx: float64, dy: float64, 
               effective_thickness: float64, gif_snapshot_interval: int64, 
               T_history_buffer: float64[:,:,:], T_history_count: int64[:]):
    """
    * @brief Numba-optimized SOR update function for thermal simulation.
    *        Modifies T in-place. For snapshots, uses a pre-allocated buffer and count.
    * @note If num_layers == 1, cond_vert.shape[0] == 0 (handled by skipping vertical loops).
    * @param T: 3D np.array (num_layers, ny, nx) temperatures
    * @param heater_sources: 3D np.array heat sources (W/node)
    * @param h_arr: 3D np.array convection coeffs (W/mm² K)
    * @param cond_lat_x: 3D np.array (num_layers, ny, nx-1) horizontal conductances
    * @param cond_lat_y: 3D np.array (num_layers, ny-1, nx) vertical conductances
    * @param cond_vert: 3D np.array (num_layers-1, ny, nx) interlayer conductances
    * @param ambient: float ambient temperature (C)
    * @param omega: float over-relaxation factor
    * @param tol: float convergence tolerance
    * @param max_iter: int maximum iterations
    * @param boundary_mode_int: int 1 for convective, 0 for insulated
    * @param h_side: float side convection coeff (W/mm² K)
    * @param dx, dy: float grid spacings (mm)
    * @param effective_thickness: float effective thickness for sides (mm)
    * @param gif_snapshot_interval: int snapshot every N iterations
    * @param T_history_buffer: 3D np.array (max_snapshots, ny, nx) pre-allocated for snapshots
    * @param T_history_count: 1D np.array [int] to track number of snapshots (modified in-place)
    * @return tuple (iterations: int64, converged: bool, max_delta: float64)
    """
    num_layers, ny, nx = T.shape
    converged = False
    iterations = 0
    max_delta = 0.0
    boundary_convective = boundary_mode_int == 1  # Scalar bool (safe in Numba)
    count = T_history_count[0]  # Mutable count for snapshots (via array ref)
    
    # Main iteration loop for SOR solver
    for iteration in range(max_iter):
        local_max_delta = 0.0
        # Loop over all layers
        for l in range(num_layers):
            # Loop over all grid points in y and x directions
            for i in range(ny):
                for j in range(nx):
                    sum_cond_T = 0.0  # Sum of (conductance * neighbor_temp)
                    sum_cond = 0.0    # Sum of conductances
                    
                    # Lateral neighbors: left (if not on left boundary)
                    if j > 0:
                        cond_left = cond_lat_x[l, i, j - 1]
                        sum_cond_T += cond_left * T[l, i, j - 1]
                        sum_cond += cond_left
                    
                    # Lateral neighbors: right (if not on right boundary)
                    if j < nx - 1:
                        cond_right = cond_lat_x[l, i, j]
                        sum_cond_T += cond_right * T[l, i, j + 1]
                        sum_cond += cond_right
                    
                    # Lateral neighbors: down (if not on bottom boundary)
                    if i > 0:
                        cond_down = cond_lat_y[l, i - 1, j]
                        sum_cond_T += cond_down * T[l, i - 1, j]
                        sum_cond += cond_down
                    
                    # Lateral neighbors: up (if not on top boundary)
                    if i < ny - 1:
                        cond_up = cond_lat_y[l, i, j]
                        sum_cond_T += cond_up * T[l, i + 1, j]
                        sum_cond += cond_up
                    
                    # Vertical neighbors: above (only if multi-layer and not top layer)
                    if l > 0 and num_layers > 1:
                        cond_above = cond_vert[l - 1, i, j]
                        sum_cond_T += cond_above * T[l - 1, i, j]
                        sum_cond += cond_above
                    
                    # Vertical neighbors: below (only if multi-layer and not bottom layer)
                    if l < num_layers - 1 and num_layers > 1:
                        cond_below = cond_vert[l, i, j]
                        sum_cond_T += cond_below * T[l + 1, i, j]
                        sum_cond += cond_below
                    
                    # Convection term (top/bottom/inner surfaces)
                    h_node = h_arr[l, i, j]
                    area = dx * dy  # Node area in mm²
                    sum_cond += h_node * area
                    sum_cond_T += h_node * area * ambient
                    
                    # Boundary conditions for sides if convective mode enabled
                    if boundary_convective:
                        # Left/right borders: vertical face convection
                        if j == 0 or j == nx - 1:
                            cond_side_vertical = h_side * dy * effective_thickness
                            sum_cond += cond_side_vertical
                            sum_cond_T += cond_side_vertical * ambient
                        # Bottom/top borders: horizontal face convection
                        if i == 0 or i == ny - 1:
                            cond_side_horizontal = h_side * dx * effective_thickness
                            sum_cond += cond_side_horizontal
                            sum_cond_T += cond_side_horizontal * ambient
                    
                    # Heat source at this node
                    source = heater_sources[l, i, j]
                    
                    # Compute new temperature using energy balance
                    if sum_cond > 0:
                        T_new = (sum_cond_T + source) / sum_cond
                    else:
                        T_new = ambient  # Fallback to ambient if no conduction
                    
                    # Apply over-relaxation (SOR method)
                    delta = T_new - T[l, i, j]
                    T[l, i, j] += omega * delta  # Update temperature
                    if abs(delta) > local_max_delta:
                        local_max_delta = abs(delta)  # Track max change for convergence
        
        max_delta = local_max_delta
        iterations = iteration + 1
        
        # Snapshot for history (copy top layer to buffer if interval matches)
        if (iteration % gif_snapshot_interval == 0 or iteration == max_iter - 1) and count < T_history_buffer.shape[0]:
            T_history_buffer[count] = T[0].copy()  # Copy top layer
            count += 1
        
        # Check convergence based on max delta
        if local_max_delta < tol:
            converged = True
            break
    
    T_history_count[0] = count  # Update count
    return iterations, converged, max_delta

class FDMSolver:
    """
    * @class FDMSolver
    * @brief Main class for simulating thermal distribution on a PCB using FDM with SOR.
    * Handles loading processed JSON config, grid setup, element processing, simulation, and results logging.
    * No visualization exports; focuses on console output and result dictionary.
    * Expects input as dict with 'pcb_params', 'zones', 'heaters', 'vias', 'sinks' keys.
    * Usage: FDMSolver.run_simulation(json_config, logger=logger)
    """

    @classmethod
    def _point_in_polygon(cls, point, poly):
        """
        * @brief Checks if a point is inside a polygon using ray casting algorithm.
        * @param point Tuple (x, y) in mm.
        * @param poly List of tuples [(x1,y1), ..., (xn,yn)] defining the polygon in mm.
        * @return True if point is inside, False otherwise.
        """
        x, y = point
        n = len(poly)
        inside = False
        p1x, p1y = poly[0]
        # Ray casting: count intersections with polygon edges
        for i in range(n + 1):
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    @classmethod
    def _polygon_area(cls, poly):
        """
        * @brief Computes the area of a polygon using the shoelace formula.
        * @param poly List of tuples [(x1,y1), ..., (xn,yn)] in mm.
        * @return Area in mm².
        """
        n = len(poly)
        area = 0.0
        # Shoelace formula: sum of (xi*yi+1 - xi+1*yi) / 2
        for i in range(n):
            j = (i + 1) % n
            area += poly[i][0] * poly[j][1]
            area -= poly[j][0] * poly[i][1]
        return abs(area) / 2.0

    @classmethod
    def _pixel_to_mm(cls, px, py, x_min, y_min, x_max, y_max, width_mm, height_mm):
        """
        * @brief Converts pixel coordinates to real mm coordinates.
        * @param px Pixel x.
        * @param py Pixel y.
        * @param x_min, y_min, x_max, y_max: Corner pixels for scaling.
        * @param width_mm, height_mm: PCB dimensions in mm.
        * @return Tuple (mm_x, mm_y). Note: y is inverted so y=0 is at bottom.
        """
        # Linear scaling from pixel space to mm space
        mx = (px - x_min) / (x_max - x_min) * width_mm
        # Invert y-axis to match PCB coordinate system (y=0 at bottom)
        my = height_mm - (py - y_min) / (y_max - y_min) * height_mm
        return mx, my

    @classmethod
    def _process_elements(cls, elements, corners, pcb, nx, ny, dx, dy, num_layers):
        """
        * @brief Processes elements from JSON to fill masks, sources, extras, vias, and zone_polys.
        *        Uses scalar loops for polygon checks to avoid external dependencies.
        * @param elements: List of element dicts (zones, heaters, vias, sinks).
        * @param corners: PCB corners for pixel-to-mm scaling.
        * @param pcb: PCB parameters dict (for dimensions).
        * @param nx, ny: Grid dimensions.
        * @param dx, dy: Grid spacings.
        * @param num_layers: Number of layers.
        * @return Tuple (copper_masks, heater_sources, h_extras, vias, zone_polys).
        """
        # Extract corners for scaling
        x_min, y_min = corners[0]
        x_max, y_max = corners[1]
        width_mm = pcb['width']
        height_mm = pcb['height']

        # Initialize data structures for masks, sources, etc.
        copper_masks = np.zeros((num_layers, ny, nx), dtype=bool)  # Boolean mask for copper areas per layer
        heater_sources = np.zeros((num_layers, ny, nx))  # W per node
        h_extras = np.zeros((num_layers, ny, nx))        # extra h in W/m² K (converted later)
        vias = []  # list of (x_mm, y_mm, diam_mm)
        zone_polys = {}  # key: (id, layer), value: poly_mm

        # Loop over all elements
        for elem in elements:
            typ = elem['type']
            points_px = elem['points']
            # Convert all points from pixels to mm
            points_mm = [cls._pixel_to_mm(p[0], p[1], x_min, y_min, x_max, y_max, width_mm, height_mm) for p in points_px]
            
            if typ in ['Zone', 'Heater', 'Sink']:
                poly_mm = points_mm
                layer = elem.get('layer', 1) - 1  # 0-based, default to layer 1 (0)
                area_mm2 = cls._polygon_area(poly_mm)  # Compute polygon area for power density
                
                # Loop over all grid nodes to check containment
                for i in range(ny):
                    for j in range(nx):
                        x = j * dx
                        y = i * dy
                        inside = cls._point_in_polygon((x, y), poly_mm)
                        
                        if inside:
                            if typ == 'Zone':
                                copper_masks[layer, i, j] = True  # Mark as copper area
                                # Store polygon for zone analysis
                                zone_id = elem['id']
                                zone_polys[(zone_id, layer)] = poly_mm
                            
                            elif typ == 'Heater':
                                if area_mm2 > 0:
                                    pd = elem['power'] / area_mm2  # W/mm² power density
                                    heater_sources[layer, i, j] += pd * dx * dy  # Distribute power to node
                            
                            elif typ == 'Sink':
                                h_extra = elem['htc_extra']  # Extra heat transfer coefficient
                                h_extras[layer, i, j] += h_extra
            
            elif typ == 'Via':
                # Vias go through all layers; no layer specified
                x_mm, y_mm = points_mm[0]
                diam = elem['diameter']
                vias.append((x_mm, y_mm, diam))  # Store for vertical conductance boost

        return copper_masks, heater_sources, h_extras, vias, zone_polys

    @classmethod
    def _precompute_conductances(cls, copper_masks, vias, pcb, dx, dy, dz, num_layers, ny, nx):
        """
        * @brief Precomputes lateral and vertical conductances based on copper masks and vias.
        *        Uses vectorized operations for efficiency where possible.
        * @param copper_masks: 3D bool array for copper areas.
        * @param vias: List of via tuples (x_mm, y_mm, diam_mm).
        * @param pcb: PCB parameters dict (conductivities, thicknesses).
        * @param dx, dy, dz: Grid spacings (mm).
        * @param num_layers, ny, nx: Dimensions.
        * @return Tuple (cond_lat_x, cond_lat_y, cond_vert, h_arr) - base convection array.
        """
        k_cu = pcb['copper_conductivity'] / 1000.0  # W/mm K (convert from W/m K)
        k_fr4 = pcb['fr4_conductivity'] / 1000.0    # W/mm K (convert from W/m K)
        t_cu = pcb['copper_thickness']              # mm

        # Initialize conductance arrays
        cond_lat_x = np.zeros((num_layers, ny, nx - 1))  # cond to right (j+1), W/K
        cond_lat_y = np.zeros((num_layers, ny - 1, nx))  # cond to up (i+1), W/K
        cond_vert = np.zeros((num_layers - 1, ny, nx)) if num_layers > 1 else np.empty((0, ny, nx))
        conv_top = pcb['convection_h_top'] / 1e6        # W/mm² K (convert from W/m² K)
        conv_bottom = pcb['convection_h_bottom'] / 1e6  # W/mm² K (convert from W/m² K)
        h_arr = np.zeros((num_layers, ny, nx))  # Convection array (W/mm² K)

        # Lateral conductances: vectorized with np.roll for neighbor checks
        for l in range(num_layers):
            # Horizontal (x): between j and j+1 - conductance if both nodes are copper
            copper_shifted_x = np.roll(copper_masks[l], -1, axis=1)
            both_copper_x = copper_masks[l, :, :-1] & copper_shifted_x[:, :-1]
            cond_lat_x[l] = np.where(both_copper_x, k_cu * t_cu * dy / dx, 0.0)  # k * A / L formula
            
            # Vertical (y): between i and i+1 - similar for y-direction
            copper_shifted_y = np.roll(copper_masks[l], -1, axis=0)
            both_copper_y = copper_masks[l, :-1, :] & copper_shifted_y[:-1, :]
            cond_lat_y[l] = np.where(both_copper_y, k_cu * t_cu * dx / dy, 0.0)
        
        # Vertical conductances: base FR4 + vias (inter-layer)
        if num_layers > 1:
            # Base FR4 conductance for all nodes
            cond_vert.fill(k_fr4 * dx * dy / dz)
            # Boost with vias: add copper conductance for via area
            for x_mm, y_mm, d in vias:
                j = int(round(x_mm / dx))
                i = int(round(y_mm / dy))
                if 0 <= i < ny and 0 <= j < nx:
                    area_via = math.pi * (d / 2.0)**2  # Via cross-section area
                    cond_add = k_cu * area_via / dz    # Additional conductance
                    for ll in range(num_layers - 1):
                        cond_vert[ll, i, j] += cond_add
        
        # Base convection array: external layers only
        h_arr[0] = conv_top  # Top layer convection
        if num_layers == 1:
            h_arr[0] += conv_bottom  # Both sides for single layer
        else:
            h_arr[-1] = conv_bottom  # Bottom layer convection
        # Inner layers remain zero (h_extras added later)

        return cond_lat_x, cond_lat_y, cond_vert, h_arr

    @classmethod
    def _run_sor_simulation(cls, T, heater_sources, h_arr, cond_lat_x, cond_lat_y, cond_vert,
                            ambient, omega, tol, max_iter, boundary_mode_int, h_side,
                            dx, dy, effective_thickness, gif_snapshot_interval, use_numba):
        """
        @brief Runs the SOR simulation loop, using Numba JIT if available.
        @details Tracks convergence and collects temperature history snapshots.
                 Uses DEFAULT_EVOLUTION_CAPTURE_STEP from config.
        @param T: Initial temperature grid (modified in-place).
        @param heater_sources: Heat sources array.
        @param h_arr: Convection array.
        @param cond_lat_x, cond_lat_y, cond_vert: Conductance arrays.
        @param ambient: Ambient temperature.
        @param omega: Over-relaxation factor.
        @param tol: Tolerance.
        @param max_iter: Max iterations.
        @param boundary_mode_int: Boundary mode flag.
        @param h_side: Side convection coeff.
        @param dx, dy: Grid spacings.
        @param effective_thickness: Effective thickness.
        @param gif_snapshot_interval: Snapshot interval (from config).
        @param use_numba: Flag for Numba usage.
        @return Tuple (iterations, converged, max_delta, T_history).
        """
        # Extract grid dimensions from T (layer 0)
        ny, nx = T.shape[1], T.shape[2]

        # Estimate maximum number of snapshots needed
        max_snapshots = max_iter // gif_snapshot_interval + 2

        # Pre-allocate buffer for T_history (Numba-friendly)
        T_history_buffer = np.zeros((max_snapshots, ny, nx), dtype=np.float64)

        # Mutable counter for Numba compatibility
        T_history_count = np.array([0], dtype=np.int64)

        # Run SOR: Numba path if available, else pure Python fallback
        if use_numba:
            # Call Numba-compiled SOR function
            iterations, converged, max_delta = sor_update(
                T, heater_sources, h_arr, cond_lat_x, cond_lat_y, cond_vert,
                ambient, omega, tol, max_iter, boundary_mode_int, h_side,
                dx, dy, effective_thickness, gif_snapshot_interval,
                T_history_buffer, T_history_count
            )
        else:
            # Fallback to pure Python (if Numba not available)
            iterations, converged, max_delta = sor_update(
                T, heater_sources, h_arr, cond_lat_x, cond_lat_y, cond_vert,
                ambient, omega, tol, max_iter, boundary_mode_int, h_side,
                dx, dy, effective_thickness, gif_snapshot_interval,
                T_history_buffer, T_history_count
            )

        # Extract actual number of captured snapshots
        num_snaps = T_history_count[0]

        # Convert buffer to list of arrays (for UI)
        T_history = [T_history_buffer[k].copy() for k in range(num_snaps)]

        # Return simulation results and history
        return iterations, converged, max_delta, T_history
    
    @classmethod
    def _analyze_hotspots(cls, T, logger, hotspot_threshold=85.0, dx=1.0, dy=1.0):
        """
        * @brief Analyzes and logs hotspots (T > threshold) per layer.
        *        Samples up to 5 hotspots with positions if any exist.
        * @param T: 3D temperature grid.
        * @param logger: Logger instance.
        * @param hotspot_threshold: Temperature threshold in C.
        * @param dx, dy: Grid spacings for coordinates.
        * @return None (logs only).
        """
        num_layers = T.shape[0]
        # Loop over layers to find hotspots
        for l in range(num_layers):
            hotspots = np.where(T[l] > hotspot_threshold)  # Indices of hot nodes
            num_hotspots = len(hotspots[0])
            if num_hotspots == 0:
                logger.info(f"FDM Hotspots Analysis T > {hotspot_threshold:.2f} °C: Layer {l + 1} = {num_hotspots} hotspots")
            else:
                logger.info(f"FDM Hotspots Analysis T > {hotspot_threshold:.2f} °C: Layer {l + 1} = {num_hotspots} hotspots")
                samples = min(5, num_hotspots)  # Sample first 5
                for k in range(samples):
                    i, j = hotspots[0][k], hotspots[1][k]
                    x = j * dx  # Convert to mm
                    y = i * dy
                    t = T[l, i, j]
                    logger.info(f"  + Hotspot at ({x:.1f} mm, {y:.1f} mm): {t:.2f} °C")
                if num_hotspots > samples:
                    logger.info(f"  ... ({num_hotspots - samples} more hotspots)")

    @classmethod
    def _analyze_gradients(cls, T, logger, dz, dx, dy, num_layers):
        """
        * @brief Analyzes and logs max thermal gradients per layer (lateral and vertical).
        *        Uses np.gradient for lateral; direct diff for vertical.
        * @param T: 3D temperature grid.
        * @param logger: Logger instance.
        * @param dz: Vertical spacing.
        * @param dx, dy: Lateral spacings.
        * @param num_layers: Number of layers.
        * @return None (logs only).
        """
        # Loop over layers
        for l in range(num_layers):
            # Lateral gradients: compute magnitude using np.gradient
            grad_y, grad_x = np.gradient(T[l], dy, dx)
            grad_mag = np.sqrt(grad_x**2 + grad_y**2)  # Euclidean magnitude
            max_grad = np.max(grad_mag)
            max_idx = np.unravel_index(np.argmax(grad_mag), grad_mag.shape)  # Index of max
            max_x = max_idx[1] * dx  # Convert to mm
            max_y = max_idx[0] * dy
            logger.info(f"FDM Thermal Gradients Analysis: Layer {l + 1} Max gradient = {max_grad:.2f} °C/mm at ({max_x:.1f} mm, {max_y:.1f} mm)")
            
            # Vertical gradients (between layers)
            if num_layers > 1 and l < num_layers - 1:
                vert_grad = np.abs(T[l] - T[l+1]) / dz  # Absolute difference per dz
                max_vert_grad = np.max(vert_grad)
                max_vert_idx = np.unravel_index(np.argmax(vert_grad), vert_grad.shape)
                vert_x = max_vert_idx[1] * dx
                vert_y = max_vert_idx[0] * dy
                logger.info(f"FDM Thermal Gradients Analysis: Vertical to layer {l+2} Max gradient = {max_vert_grad:.2f} °C/mm at ({vert_x:.1f} mm, {vert_y:.1f} mm)")

    @classmethod
    def _analyze_heat_flow(cls, T, logger, h_arr, heater_sources, ambient, h_side, dx, dy, 
                           effective_thickness, boundary_mode, num_layers, ny, nx):
        """
        * @brief Analyzes total heat flow: input from heaters vs. dissipated (convection + sides).
        *        Sums convection over active nodes and side boundary losses.
        * @param T: 3D temperature grid.
        * @param logger: Logger instance.
        * @param h_arr: Convection array.
        * @param heater_sources: Heat sources.
        * @param ambient: Ambient temperature.
        * @param h_side: Side convection coeff.
        * @param dx, dy: Grid spacings.
        * @param effective_thickness: Effective thickness.
        * @param boundary_mode: Boundary mode string.
        * @param num_layers, ny, nx: Dimensions.
        * @return total_dissipated (float) - total heat loss.
        """
        total_power = np.sum(heater_sources)  # Total input power from all heaters
        
        # Convection losses: sum over nodes with h > 0
        total_convection_loss = 0.0
        for l in range(num_layers):
            h_active = h_arr[l] > 0  # Mask for active convection nodes
            if np.any(h_active):
                area = dx * dy  # Node area
                delta_t_active = T[l][h_active] - ambient  # Temperature rise
                total_convection_loss += np.sum(h_arr[l][h_active] * area * delta_t_active)  # Q = h * A * dT
        
        # Side losses: approximate boundary convection
        side_loss = 0.0
        if boundary_mode == 'convective':
            for l in range(num_layers):
                # Left/right borders: vertical faces
                left_right = T[l, :, 0] + T[l, :, -1]  # Temps on left + right
                side_loss += h_side * dy * effective_thickness * np.sum(left_right - 2 * ambient)
                # Bottom/top borders: horizontal faces
                bottom_top = T[l, 0, :] + T[l, -1, :]
                side_loss += h_side * dx * effective_thickness * np.sum(bottom_top - 2 * ambient)
        
        total_dissipated = total_convection_loss + side_loss  # Total output heat
        
        # Log results
        logger.info(f"Heat Flow Analysis: Total heater power input = {total_power:.2f} W")
        logger.info(f"Heat Flow Analysis: Total dissipated (convection + sides) = {total_dissipated:.2f} W")
        
        return total_dissipated

    @classmethod
    def _analyze_zone_temps(cls, T, logger, zone_polys, ny, nx, dx, dy):
        """
        * @brief Analyzes average temperatures in defined zones.
        *        Loops over grid nodes inside each polygon.
        * @param T: 3D temperature grid.
        * @param logger: Logger instance.
        * @param zone_polys: Dict of zone polygons { (id, layer): poly_mm }.
        * @param ny, nx: Grid dimensions.
        * @param dx, dy: Grid spacings.
        * @return None (logs only).
        """
        # Loop over all defined zones
        for (zone_id, layer), poly_mm in zone_polys.items():
            temps_in_zone = []
            # Check each grid node for containment
            for i in range(ny):
                y = i * dy
                for j in range(nx):
                    x = j * dx
                    if cls._point_in_polygon((x, y), poly_mm):
                        temps_in_zone.append(T[layer, i, j])  # Collect temps
            if temps_in_zone:
                avg_temp = np.mean(temps_in_zone)  # Average over nodes in zone
                logger.info(f"FDM Zone Temperatures Analysis: Zone {zone_id} (Layer {layer + 1}) Avg temp = {avg_temp:.2f} °C")
            else:
                logger.info(f"FDM Zone Temperatures Analysis: Zone {zone_id} (Layer {layer + 1}) = No nodes in zone")

    @classmethod
    def simulate(cls, json_config, logger=None):
        """
        @brief Main entry point for FDM simulation.
        @details Processes json_config with pcb_params (including pcb_corners), elements, and solver params.
                 Sets up grid, processes elements, runs SOR, analyzes results, and logs summaries.
                 Captures T_history every DEFAULT_EVOLUTION_CAPTURE_STEP iterations.
        @param json_config Dict with pcb_params, solver_mode, elements.
        @param logger Optional logging.Logger for output.
        @return Dict with simulation metrics and grids.
        """
        # Fallback to module logger if none provided
        if logger is None:
            logger = logging.getLogger(__name__)

        try:
            # Extract parameters from json_config
            data = json_config
            pcb = data['pcb_params']
            corners = pcb.get('pcb_corners', [[0, 0], [1000, 1000]])  # From provided example JSON
            x_min, y_min = corners[0]
            x_max, y_max = corners[1]
            width_mm = pcb['width']
            height_mm = pcb['height']
            hotspot_threshold= pcb.get('hotspot_threshold', 85.0)

            # Calculate effective area for compatibility with other solvers (total PCB area)
            effective_area = (width_mm / 1000) * (height_mm / 1000)  # m²

            # Reconstruct elements list from categorized lists (but use full 'elements' if available)
            elements = data.get('elements', [])  # Full list from json_config

            # Grid setup based on parameters
            dx = pcb['grid_dx']
            dy = pcb['grid_dy']
            nx = int(width_mm / dx) + 1
            ny = int(height_mm / dy) + 1
            num_layers = pcb['layers']

            # Unit conversions for consistency (mm, W, K)
            k_cu = pcb['copper_conductivity'] / 1000.0  # W/mm K
            k_fr4 = pcb['fr4_conductivity'] / 1000.0    # W/mm K
            t_cu = pcb['copper_thickness']              # mm
            dz = pcb['thickness'] / (num_layers - 1) if num_layers > 1 else 0.0  # mm between layers
            ambient = pcb['ambient_temp_global']         # C

            conv_global = pcb['convection_h_global'] / 1e6  # W/mm² K
            conv_top = pcb['convection_h_top'] / 1e6        # W/mm² K
            conv_bottom = pcb['convection_h_bottom'] / 1e6  # W/mm² K
            
            # New: Get boundary_mode from pcb_params (default to 'convective')
            boundary_mode = pcb.get('boundary_mode', 'Convective').lower()
            h_side = conv_global  # h for side boundaries in W/mm² K
            boundary_mode_int = 1 if boundary_mode == 'convective' else 0  # 1 for convective, 0 for insulated

            # Effective thickness for side convection approximation (per layer, focusing on copper)
            effective_thickness = t_cu

            # SOR parameters
            omega = pcb['over_relaxation_factor']
            tol = pcb['tolerance']
            max_iter = pcb['max_iterations']

            # Snapshot interval for temperature history (GIF animation)
            gif_snapshot_interval = config.DEFAULT_EVOLUTION_CAPTURE_STEP  # ← USE CONFIG

            # Check for Numba availability
            use_numba = USE_NUMBA

            # Step 1: Process elements
            copper_masks, heater_sources, h_extras, vias, zone_polys = cls._process_elements(
                elements, corners, pcb, nx, ny, dx, dy, num_layers
            )

            # Add h_extras to h_arr (post-conductance, as it's layer-specific)
            h_arr = np.zeros((num_layers, ny, nx))  # Re-init here for extras
            h_arr[0] += h_extras[0] / 1e6
            if num_layers > 1:
                h_arr[-1] += h_extras[-1] / 1e6

            # Step 2: Precompute conductances (base h_arr without extras already included above)
            cond_lat_x, cond_lat_y, cond_vert, base_h_arr = cls._precompute_conductances(
                copper_masks, vias, pcb, dx, dy, dz, num_layers, ny, nx
            )
            # Merge base_h_arr with extras (already done above, but ensure)
            h_arr += base_h_arr  # This adds conv_top/bottom if not already

            # Step 3: Initialize temperature grid
            T = np.full((num_layers, ny, nx), ambient, dtype=np.float64)  # Temperature grids (3D)

            # Step 4: Run SOR simulation
            iterations, converged, max_delta, T_history = cls._run_sor_simulation(
                T, heater_sources, h_arr, cond_lat_x, cond_lat_y, cond_vert,
                ambient, omega, tol, max_iter, boundary_mode_int, h_side,
                dx, dy, effective_thickness, gif_snapshot_interval, use_numba
            )

            # Log convergence with explicit prefix
            if converged:
                logger.info(f"FDM converged in {iterations} iterations with max delta {max_delta:.4f}")
            else:
                logger.warning(f"FDM did not converge after {max_iter} iterations, max delta {max_delta:.4f}")

            # Log min/max and average temps per layer with explicit format
            for l in range(num_layers):
                min_temp = np.min(T[l])
                max_temp = np.max(T[l])
                avg_temp = np.mean(T[l])
                logger.info(f"Layer {l + 1}: Min temp = {min_temp:.2f} °C, Max temp = {max_temp:.2f} °C, Avg temp = {avg_temp:.2f} °C")

            # Compute results
            global_min = np.min(T)
            global_max = np.max(T)
            global_avg = np.mean(T)
            avg_per_layer = [np.mean(T[l]) for l in range(num_layers)]
            t_top = avg_per_layer[0]
            t_bottom = avg_per_layer[-1] if num_layers > 1 else t_top
            # Compute inner average if more than 2 layers
            inner_avg = np.mean(avg_per_layer[1:-1]) if num_layers > 2 else None
            total_power = np.sum(heater_sources)
            delta_t = global_max - ambient
            effective_r_th = delta_t / total_power if total_power > 0 else 0.0
            area_m2 = (width_mm * height_mm) / 1e6  # mm² to m²
            heat_flux = total_power / area_m2 if area_m2 > 0 else 0.0

            # Step 5: Detailed analyses
            cls._analyze_hotspots(T, logger, hotspot_threshold, dx=dx, dy=dy)
            cls._analyze_gradients(T, logger, dz, dx, dy, num_layers)
            total_dissipated = cls._analyze_heat_flow(
                T, logger, h_arr, heater_sources, ambient, h_side, dx, dy,
                effective_thickness, boundary_mode, num_layers, ny, nx
            )
            cls._analyze_zone_temps(T, logger, zone_polys, ny, nx, dx, dy)

            # Log summary results with explicit prefix for each line (expanded for consistency)
            logger.info(f"FDM Simulation Results: Max Temperature = {global_max:.2f} °C")
            logger.info(f"FDM Simulation Results: Min Temperature = {global_min:.2f} °C")
            logger.info(f"FDM Simulation Results: Average Temperature = {global_avg:.2f} °C")
            logger.info(f"FDM Simulation Results: Delta T = {delta_t:.2f} °C")
            logger.info(f"FDM Simulation Results: Effective R_th = {effective_r_th:.2f} K/W")
            logger.info(f"FDM Simulation Results: Total Power = {total_power:.2f} W")
            logger.info(f"FDM Simulation Results: Effective Area = {effective_area:.6f} m²")
            logger.info(f"FDM Simulation Results: T Top = {t_top:.2f} °C")
            logger.info(f"FDM Simulation Results: T Bottom = {t_bottom:.2f} °C")
            if inner_avg is not None:
                logger.info(f"FDM Simulation Results: T Inner = {inner_avg:.2f} °C")
            logger.info(f"FDM Simulation Results: Avg Temp per Layer = {', '.join([f'{t:.2f}' for t in avg_per_layer])} °C")
            logger.info(f"FDM Simulation Results: Iterations = {iterations}")
            logger.info(f"FDM Simulation Results: Converged = {converged}")
            logger.info(f"FDM Simulation Results: Max Delta = {max_delta:.6f} °C")
            logger.info(f"FDM Simulation Results: Note = FDM with SOR - detailed spatial resolution on grid")

            # Prepare results dictionary
            results = {
                "max_temp": {"value": global_max, "unit": "°C"},
                "avg_temp": {"value": global_avg, "unit": "°C"},
                "delta_t": {"value": delta_t, "unit": "°C"},
                "effective_r_th": {"value": effective_r_th, "unit": "K/W"},
                "total_power": {"value": total_power, "unit": "W"},
                "effective_area": {"value": effective_area, "unit": "m²"},
                "t_top": {"value": t_top, "unit": "°C"},
                "t_bottom": {"value": t_bottom, "unit": "°C"},
                "t_inner": {"value": inner_avg, "unit": "°C"} if inner_avg is not None else None,
                "avg_temp_per_layer": avg_per_layer,
                "iterations": iterations,
                "converged": converged,
                "T_grid": T,
                "T_history": T_history,
                'dx': dx,
                'dy': dy,
                'width': width_mm,
                'height': height_mm,
                "hotspot_threshold": hotspot_threshold,
                "zone_polys": zone_polys,
                "note": {"value": "FDM with SOR - detailed spatial resolution on grid", "unit": ""}
            }
            return results

        except Exception as e:
            logger.error(f"Simulation failed: {str(e)}")
            raise

    

class FDMResultsDialog(wx.Dialog):
    """
    @brief Dialog for displaying FDM simulation results with multiple tabs.
    @details Uses wx.Notebook for tabbed interface. Includes Save/Close buttons aligned to the right.
             Tabs will be populated in subsequent steps.
    """
    def __init__(self, parent, results, T_history=None, zone_polys=None, logger=None):
        """
        @brief Initialize the FDM results dialog with correct size and loading state.
        @details Sets minimum and initial size, centers window, prepares content sizer.
        """
        super().__init__(parent, title="FDM Simulation Results", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        # Set size: minimum and initial
        self.SetMinSize(config.DEFAULT_RESULTS_WINDOW_SIZE)
        self.SetSize(config.DEFAULT_RESULTS_WINDOW_SIZE)  # ← Tamanho inicial
        self.Center()  # ← Centraliza na tela

        self.results = results
        self.T_history = T_history or []
        self.zone_polys = zone_polys or {}
        
        # Use module logger if none provided
        self.logger = logger
        if logger is None:
            self.logger = logging.getLogger(__name__)

        # Bind close event
        self.Bind(wx.EVT_CLOSE, self.on_close_window)

        # Main vertical sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Content sizer: loading panel or notebook
        self.content_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.main_sizer.Add(self.content_sizer, 1, wx.EXPAND | wx.ALL, 5)

        # Notebook: created after simulation
        self.notebook = None

        # Button panel
        self.btn_panel = wx.Panel(self)
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_sizer.AddStretchSpacer()

        # Save button (hidden until results)
        self.save_btn = wx.Button(self.btn_panel, label="Save")
        self.save_btn.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_btn.Hide()

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
        @brief Show full-screen loading message (no tabs).
        @details Clears current content, creates a centered loading panel with large text,
                and forces layout update to display immediately during simulation.
        """
        # Clear all existing widgets from content area (including any previous tabs)
        self.content_sizer.Clear(True)

        # Create a new panel to hold the loading message
        loading_panel = wx.Panel(self)

        # Create a vertical sizer to center the text vertically
        sizer = wx.BoxSizer(wx.VERTICAL)

        # Create static text with loading message
        loading_text = wx.StaticText(
            loading_panel,
            label="Simulation Running... Please wait."
        )

        # Set a larger, readable font for the loading message
        loading_text.SetFont(wx.Font(
            12,                              # Font size in points
            wx.FONTFAMILY_SWISS,             # Sans-serif font family
            wx.FONTSTYLE_NORMAL,             # Normal style (not italic)
            wx.FONTWEIGHT_NORMAL             # Normal weight (not bold)
        ))

        # Add flexible space above the text to push it down
        sizer.AddStretchSpacer()

        # Add the loading text: centered horizontally and vertically with 20px padding
        sizer.Add(loading_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)

        # Add flexible space below the text to push it up from bottom
        sizer.AddStretchSpacer()

        # Apply the sizer to the loading panel
        loading_panel.SetSizer(sizer)

        # Add the loading panel to the main content sizer with full expansion
        self.content_sizer.Add(loading_panel, 1, wx.EXPAND)

        # Force immediate layout update to show loading screen right away
        self.Layout()

    def create_summary_tab(self):
        """
        @brief Create the Summary tab for FDM simulation results.
        @details Displays bar chart and metrics.
                Saves chart to config.DEFAULT_TEMP_DIR as 'fdm_summary_chart.png'.
                Saves text content to config.DEFAULT_TEMP_DIR as 'fdm_summary_report.txt'.
                Files are copied on Save and deleted on close.
        @return wx.Panel with bar chart and text.
        """
        # Create main panel for the Summary tab
        panel = wx.Panel(self.notebook)
        # Create vertical sizer to organize content
        sizer = wx.BoxSizer(wx.VERTICAL)

        # Validate that results are available and in dictionary format
        if not isinstance(self.results, dict):
            # Create error message if results are missing or invalid
            error_text = wx.StaticText(panel, label="Error: No results.")
            # Add error text to sizer with padding and centering
            sizer.Add(error_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            # Apply sizer to panel
            panel.SetSizer(sizer)
            # Return early with error panel
            return panel

        # Helper function to safely extract nested values from results dictionary
        def safe_get(key, subkey=None, default=0.0):
            item = self.results.get(key, {})
            if subkey is None:
                return item if item is not None else default
            return item.get(subkey, default) if isinstance(item, dict) else default

        # Extract key temperature metrics using safe_get
        max_temp = safe_get('max_temp', 'value', 0.0)
        avg_temp = safe_get('avg_temp', 'value', 0.0)
        delta_t = safe_get('delta_t', 'value', 0.0)

        # Prepare data for bar chart
        vals = [max_temp, avg_temp, delta_t]
        labs = ['Max Temp', 'Avg Temp', 'Delta T']

        # Add T Top if available
        if 't_top' in self.results and isinstance(self.results['t_top'], dict):
            vals.append(self.results['t_top'].get('value', 0.0))
            labs.append('T Top')

        # Add T Bottom if available
        if 't_bottom' in self.results and isinstance(self.results['t_bottom'], dict):
            vals.append(self.results['t_bottom'].get('value', 0.0))
            labs.append('T Bottom')

        # Add T Inner if available and not None
        if 't_inner' in self.results and isinstance(self.results['t_inner'], dict) and self.results['t_inner'].get('value') is not None:
            vals.append(self.results['t_inner']['value'])
            labs.append('T Inner')

        # LEFT: Create panel for bar chart
        left_panel = wx.Panel(panel)
        chart_sizer = wx.BoxSizer(wx.VERTICAL)

        # Set global Matplotlib font size from config
        plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})

        # Create figure and axis with compact size
        fig, ax = plt.subplots(figsize=(4, 3))

        # Set figure background to match wxPython panel
        bg = self.GetBackgroundColour()
        fig.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))

        # Set plot area to white
        ax.set_facecolor('white')

        # Enable major grid with dashed lines
        ax.grid(True, which='major', linestyle='--', alpha=0.7, color='gray', zorder=0)
        # Enable minor grid with dotted lines
        ax.grid(True, which='minor', linestyle=':', alpha=0.3, color='gray', zorder=0)
        # Turn on minor ticks
        ax.minorticks_on()

        # Create bar chart using primary color from config
        bars = ax.bar(labs, vals, color=config.DEFAULT_CHART_COLOR_SCHEME[0], zorder=3)

        # Fix x-axis tick positions to avoid automatic rotation issues
        ax.xaxis.set_major_locator(FixedLocator(range(len(labs))))
        # Rotate x-axis labels 90 degrees for readability
        ax.set_xticklabels(labs, rotation=90)

        # Calculate Y-axis upper limit with margin
        max_val = max(vals) if vals else 0
        margin_percent = config.DEFAULT_CHART_Y_MARGIN_PERCENT / 100.0
        if margin_percent > 0:
            y_margin = max_val * margin_percent
            ax.set_ylim(0, max_val + y_margin)
            label_offset = y_margin * 0.3
        else:
            label_offset = max_val * 0.05

        # Set chart title and Y-axis label
        ax.set_title("FDM Simulation Temperatures")
        ax.set_ylabel("Temperature (°C)")

        # Add value labels on top of each bar
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2.,
                    bar.get_height() + label_offset,
                    f'{v:.2f}',
                    ha='center', va='bottom')

        # Optimize layout to prevent clipping
        plt.tight_layout()

        # Ensure temporary directory exists
        os.makedirs(config.DEFAULT_TEMP_DIR, exist_ok=True)
        # Define fixed filename for chart
        chart_path = os.path.join(config.DEFAULT_TEMP_DIR, "fdm_summary_chart.png")
        # Save chart as PNG with 100 DPI
        fig.savefig(chart_path, format='png', dpi=100, bbox_inches='tight')
        # Close figure to free memory
        plt.close(fig)

        # Load saved chart as wx.Image
        img = wx.Image(chart_path, wx.BITMAP_TYPE_PNG)
        # Convert to bitmap (fallback to 1x1 if invalid)
        bmp = wx.Bitmap(img) if img.IsOk() else wx.Bitmap(1, 1)
        # Create static bitmap control to display chart
        chart_bitmap = wx.StaticBitmap(left_panel, bitmap=bmp)
        # Add bitmap to sizer with padding
        chart_sizer.Add(chart_bitmap, 0, wx.ALL | wx.EXPAND, 5)
        # Apply sizer to left panel
        left_panel.SetSizer(chart_sizer)

        # RIGHT: Create panel for text summary
        right_panel = wx.Panel(panel)
        txt_sizer = wx.BoxSizer(wx.VERTICAL)

        # Build text content string with key metrics
        text_content = (
            f"Max Temperature: {max_temp:.2f} °C\n"
            f"Average Temperature: {avg_temp:.2f} °C\n"
            f"Delta T: {delta_t:.2f} °C\n"
            f"Effective R Th: {safe_get('effective_r_th', 'value', 0.0):.2f} K/W\n"
            f"Total Power: {safe_get('total_power', 'value', 0.0):.2f} W\n"
            f"Effective Area: {safe_get('effective_area', 'value', 0.0):.6f} m²\n"
            f"Iterations: {self.results.get('iterations', 'N/A')}\n"
            f"Converged: {self.results.get('converged', 'N/A')}\n"
            f"Max Delta: {self.results.get('max_delta', 0.0):.6f} °C\n\n"
        )

        # Add T Top if available
        if 't_top' in self.results and isinstance(self.results['t_top'], dict):
            text_content += f"T Top: {self.results['t_top'].get('value', 0.0):.2f} °C\n"

        # Add T Bottom if available
        if 't_bottom' in self.results and isinstance(self.results['t_bottom'], dict):
            text_content += f"T Bottom: {self.results['t_bottom'].get('value', 0.0):.2f} °C\n"

        # Add T Inner if available and not None
        if 't_inner' in self.results and isinstance(self.results['t_inner'], dict) and self.results['t_inner'].get('value') is not None:
            text_content += f"T Inner: {self.results['t_inner']['value']:.2f} °C\n"

        # Add explanation text from config
        text_content += f"\n{config.EXPLANATION_TEXT.get('FDM', 'FDM Notes not available.')}"

        # Define fixed filename for text report
        txt_path = os.path.join(config.DEFAULT_TEMP_DIR, "fdm_summary_report.txt")
        # Save text content to file
        with open(txt_path, 'w', encoding='utf-8') as f:
            f.write(text_content)

        # Create read-only multiline text control
        text_ctrl = wx.TextCtrl(right_panel, value=text_content, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_BESTWRAP)
        # Set font from config
        text_ctrl.SetFont(wx.Font(config.DEFAULT_RESULTS_FONT_SIZE, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        # Add text control to sizer
        txt_sizer.Add(text_ctrl, 1, wx.ALL | wx.EXPAND, 5)
        # Apply sizer to right panel
        right_panel.SetSizer(txt_sizer)

        # Create horizontal sizer to place chart and text side by side
        content_sizer = wx.BoxSizer(wx.HORIZONTAL)
        content_sizer.Add(left_panel, 1, wx.ALL | wx.EXPAND, 5)
        content_sizer.Add(right_panel, 1, wx.ALL | wx.EXPAND, 5)

        # Add horizontal layout to main vertical sizer
        sizer.Add(content_sizer, 1, wx.EXPAND)

        # Apply final sizer to panel
        panel.SetSizer(sizer)

        # Return completed panel
        return panel

    def create_heatmap_tab(self):
        """
        @brief Create the Heatmap tab for FDM simulation results.
        @details This tab shows:
                - One heatmap per layer using matplotlib's imshow
                - Temperature grid from self.results['T_grid']
                - Colorbar in °C with 'hot' colormap
                - FINE GRID LINES matching FDM dx/dy (full resolution)
                - TICK LABELS at 0%, 25%, 50%, 75%, 100% of the axis
                - CORRECT PHYSICAL SCALE: width = JSON 'width', height = JSON 'height'
                - COLORBAR TICKS ROUNDED TO 0 DECIMAL PLACES
                - DEBUG PRINTS: JSON width/height, grid size, expected nx/ny
                - Layer title ABOVE the plot (NO DUPLICATE BELOW)
                - Scrollable panel if multiple layers
                - SAVES EACH HEATMAP TO config.DEFAULT_TEMP_DIR as 'fdm_heatmap_layer_N.png'
                - Uses config for font size and colormap
        @return wx.Panel containing the Heatmap tab with scrollable heatmaps.
        """
        # Create main panel for the Heatmap tab
        panel = wx.Panel(self.notebook)
        # Create vertical sizer to organize content
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Ensure temporary directory exists for saving heatmaps
        os.makedirs(config.DEFAULT_TEMP_DIR, exist_ok=True)

        # Verify that temperature grid is present and is a NumPy array
        if 'T_grid' not in self.results or not isinstance(self.results['T_grid'], np.ndarray):
            # Log error if grid is missing or invalid
            self.logger.error("T_grid missing or invalid in results")
            # Create error message for user
            error_text = wx.StaticText(
                panel,
                label="Error: Temperature grid not available."
            )
            # Set bold font for emphasis
            error_text.SetFont(wx.Font(
                12, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD
            ))
            # Add error text to sizer with padding and centering
            main_sizer.Add(error_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            # Apply sizer and return error panel
            panel.SetSizer(main_sizer)
            return panel

        # Extract temperature grid (shape: num_layers, ny, nx)
        T_grid = self.results['T_grid']
        # Get grid dimensions
        num_layers, ny, nx = T_grid.shape

        # Extract physical PCB dimensions from results (defaults if missing)
        width_mm = self.results.get('width', 200.0)    # PCB width in mm
        height_mm = self.results.get('height', 100.0)  # PCB height in mm
        # Extract grid spacing
        dx = self.results.get('dx', 1.0)  # Grid spacing in X direction (mm)
        dy = self.results.get('dy', 1.0)  # Grid spacing in Y direction (mm)

        # Compute global temperature range across all layers
        global_min = np.min(T_grid)
        global_max = np.max(T_grid)
        
        # Create scrollable window with vertical scrollbar
        scroll = wx.ScrolledWindow(panel, style=wx.VSCROLL)
        scroll.SetScrollRate(20, 20)  # Set scroll speed
        scroll_sizer = wx.BoxSizer(wx.VERTICAL)

        
        # Loop through each layer and create heatmap
        for layer_idx in range(num_layers):
            # Create panel for current layer
            layer_panel = wx.Panel(scroll)
            layer_sizer = wx.BoxSizer(wx.VERTICAL)

            # Extract temperature data for this layer
            layer_data = T_grid[layer_idx]

            # Determine layer name based on position
            if num_layers == 1:
                layer_name = "Single Layer"
            elif layer_idx == 0:
                layer_name = "Top Layer"
            elif layer_idx == num_layers - 1:
                layer_name = "Bottom Layer"
            else:
                layer_name = f"Inner Layer {layer_idx}"

            # Create matplotlib figure and axis
            fig, ax = plt.subplots(figsize=(5, 4))
            # Set global font size from config
            plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})

            # Set figure background to match dialog
            bg = self.GetBackgroundColour()
            fig.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))
            # Set plot area to white
            ax.set_facecolor('white')

            # Display heatmap with correct physical extent and uniform color scale
            im = ax.imshow(
                layer_data,
                cmap=config.DEFAULT_VISUALIZATION_HEAT_COLORMAP_STYLE,  # Use 'hot' colormap from config
                origin='lower',  # Origin at bottom-left (matches FDM grid)
                extent=[0, width_mm, 0, height_mm],  # Physical scale in mm
                interpolation='nearest',  # No smoothing for grid accuracy
                vmin=global_min,  # Uniform min across all layers
                vmax=global_max   # Uniform max across all layers
            )

            # Add colorbar with rounded integer ticks
            cbar = fig.colorbar(im, ax=ax, label='Temperature (°C)')
            cbar.ax.tick_params(labelsize=config.DEFAULT_CHART_FONT_SIZE)
            # Create 5 evenly spaced ticks
            cbar_ticks = np.linspace(global_min, global_max, 5)
            # Round to nearest integer
            cbar_ticks_rounded = np.round(cbar_ticks).astype(int)
            # Set tick positions
            cbar.set_ticks(cbar_ticks)
            # Set tick labels as integers
            cbar.set_ticklabels([f"{int(t)}" for t in cbar_ticks_rounded])

            # Set title above the plot (with padding)
            ax.set_title(f"{layer_name} - Temperature Map", pad=15)

            # Set axis labels with units
            ax.set_xlabel("X (mm)")
            ax.set_ylabel("Y (mm)")

            # FINE GRID: Show lines at every dx/dy interval
            ax.set_xticks(np.arange(0, width_mm + dx, dx))
            ax.set_yticks(np.arange(0, height_mm + dy, dy))
            ax.grid(True, which='both', color='gray', linestyle='-', linewidth=0.3, alpha=0.4)

            # TICK LABELS: Show only at 0%, 25%, 50%, 75%, 100% — rounded to 0 decimals
            x_positions = [0, width_mm * 0.25, width_mm * 0.5, width_mm * 0.75, width_mm]
            y_positions = [0, height_mm * 0.25, height_mm * 0.5, height_mm * 0.75, height_mm]
            ax.set_xticks(x_positions)
            ax.set_yticks(y_positions)
            ax.set_xticklabels([f"{int(round(x))}" for x in x_positions])
            ax.set_yticklabels([f"{int(round(y))}" for y in y_positions])

            # Optimize layout to prevent clipping
            plt.tight_layout()

            # Define filename for this layer's heatmap
            heatmap_filename = f"fdm_heatmap_layer_{layer_idx + 1}.png"
            heatmap_path = os.path.join(config.DEFAULT_TEMP_DIR, heatmap_filename)
            # Save figure to temporary directory
            fig.savefig(heatmap_path, format='png', dpi=100, bbox_inches='tight')
            # Close figure to free memory
            plt.close(fig)

            # Load saved image as wx.Image
            img = wx.Image(heatmap_path, wx.BITMAP_TYPE_PNG)
            # Convert to bitmap (fallback to 1x1 if invalid)
            bmp = wx.Bitmap(img) if img.IsOk() else wx.Bitmap(1, 1)
            # Create static bitmap to display heatmap
            bitmap = wx.StaticBitmap(layer_panel, bitmap=bmp)
            # Add bitmap to layer sizer with padding
            layer_sizer.Add(bitmap, 0, wx.ALL | wx.ALIGN_CENTER, 10)
            # Apply sizer to layer panel
            layer_panel.SetSizer(layer_sizer)
            # Add layer panel to scroll sizer
            scroll_sizer.Add(layer_panel, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Apply sizer to scroll window
        scroll.SetSizer(scroll_sizer)
        scroll.FitInside()  # Adjust scroll area to content

        # Add scroll window to main sizer
        main_sizer.Add(scroll, 1, wx.EXPAND | wx.ALL, 5)
        # Apply final sizer to panel
        panel.SetSizer(main_sizer)

        # Return completed Heatmap tab panel
        return panel
    
    def create_evolution_tab(self):
        """
        @brief Create the Evolution tab with animated GIF.
        @details Uses T_history to generate GIF with EXACT same visual style as heatmap.
                 Saves to config.DEFAULT_TEMP_DIR as 'fdm_evolution_animation.gif'.
                 Displays with wx.adv.AnimationCtrl using SetAnimation().
                 File is deleted when dialog closes (via on_close_window).
        @return wx.Panel with animated GIF.
        """
        # Create main panel for the Evolution tab
        panel = wx.Panel(self.notebook)
        # Create vertical sizer to organize content
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Ensure temporary directory exists for saving the GIF
        import os
        os.makedirs(config.DEFAULT_TEMP_DIR, exist_ok=True)

        # Extract temperature history from results
        T_history = self.results.get('T_history')
        # Validate that history data exists and is not empty
        if not T_history or len(T_history) == 0:
            # Create placeholder message if no evolution data
            placeholder = wx.StaticText(panel, label="Evolution data not available.")
            # Add placeholder to sizer with padding and centering
            main_sizer.Add(placeholder, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            # Apply sizer and return panel
            panel.SetSizer(main_sizer)
            return panel

        # Extract grid spacing from results (defaults if missing)
        dx = self.results.get('dx', 1.0)  # Grid spacing in X direction (mm)
        dy = self.results.get('dy', 1.0)  # Grid spacing in Y direction (mm)
        # Get grid dimensions from first frame
        ny, nx = T_history[0].shape
        # Calculate physical PCB size from grid
        width_mm = (nx - 1) * dx
        height_mm = (ny - 1) * dy
        # Compute global temperature range across all frames
        global_min = min(np.min(f) for f in T_history)
        global_max = max(np.max(f) for f in T_history)
        # Create matplotlib figure and axis
        fig, ax = plt.subplots(figsize=(6, 4))
        # Initialize first frame of the animation
        im = ax.imshow(
            T_history[0],
            cmap=config.DEFAULT_VISUALIZATION_HEAT_COLORMAP_STYLE,  # Use 'hot' colormap from config
            origin='lower',  # Origin at bottom-left (matches FDM grid)
            extent=[0, width_mm, 0, height_mm],  # Physical scale in mm
            vmin=global_min,  # Uniform color scale min
            vmax=global_max,  # Uniform color scale max
            interpolation='nearest'  # No smoothing for grid accuracy
        )
        # Add colorbar with 5 ticks rounded to integers
        cbar = fig.colorbar(im, ax=ax, label='Temperature (°C)')
        cbar.set_ticks(np.linspace(global_min, global_max, 5))
        cbar.set_ticklabels([f"{int(round(t))}" for t in np.linspace(global_min, global_max, 5)])
        # Define update function for each animation frame
        def update(frame_idx):
            """
            @brief Update the plot for a given frame index.
            @param frame_idx Current iteration index.
            @return Updated image object.
            """
            # Update image data with current frame
            im.set_array(T_history[frame_idx])
            # Update title with iteration number
            ax.set_title(f"Thermal Evolution — Iteration {frame_idx}")
            # Return modified artist for blitting
            return im,
        # Create animation using FuncAnimation
        anim = FuncAnimation(
            fig, update, frames=len(T_history),
            interval=config.DEFAULT_EVOLUTION_GIF_FRAME_DELAY_MS,  # Delay between frames in ms
            blit=True,  # Enable blitting for performance
            repeat=True  # Loop animation
        )
        # Define fixed filename for the GIF
        gif_path = os.path.join(config.DEFAULT_TEMP_DIR, "fdm_evolution_animation.gif")
        # Save animation as GIF using Pillow writer
        anim.save(gif_path, writer='pillow', fps=1000 / config.DEFAULT_EVOLUTION_GIF_FRAME_DELAY_MS)
        # Close figure to free memory
        plt.close(fig)
        # Load and display the GIF
        try:
            # Create wx.Animation from saved GIF
            gif_anim = wx.adv.Animation(gif_path)
            # Create animation control
            anim_ctrl = wx.adv.AnimationCtrl(panel)
            # Assign animation to control
            anim_ctrl.SetAnimation(gif_anim)
            # Start playback
            anim_ctrl.Play()
            # Add control to sizer with expansion and padding — CENTRALIZADO
            main_sizer.AddStretchSpacer()
            main_sizer.Add(anim_ctrl, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            main_sizer.AddStretchSpacer()
        except Exception as e:
            # Log error if GIF fails to load
            self.logger.error(f"GIF failed: {e}")
            # Create fallback message
            fallback = wx.StaticText(panel, label="GIF failed to load.")
            # Add fallback to sizer
            main_sizer.Add(fallback, 0, wx.ALL | wx.ALIGN_CENTER, 20)
        # Apply sizer to panel
        panel.SetSizer(main_sizer)
        # Return completed Evolution tab panel
        return panel
    
    def create_hotspots_tab(self):
        """
        @brief Create the Hotspots tab for FDM simulation results.
        @details Displays two separate heatmaps for Top and Bottom layers only, stacked vertically.
                 Identifies hotspots above self.results['hotspot_threshold'] using T_grid.
                 Marks hotspots with config.DEFAULT_HOTSPOT_COLOR.
                 Adds legend at the bottom of each graph with a colored dot indicating hotspot color, positioned below the X-axis.
                 Saves each heatmap to config.DEFAULT_TEMP_DIR as 'fdm_hotspots_top.png' and 'fdm_hotspots_bottom.png'.
                 Uses scrollable panel with centered content, styled like the heatmap tab.
        @return wx.ScrolledWindow with hotspots visualization.
        """
        # Create main panel for the Hotspots tab
        panel = wx.Panel(self.notebook)
        # Create vertical sizer to organize content
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        # Ensure temporary directory exists for saving heatmaps
        os.makedirs(config.DEFAULT_TEMP_DIR, exist_ok=True)

        # Verify that temperature grid is present and is a NumPy array
        if 'T_grid' not in self.results or not isinstance(self.results['T_grid'], np.ndarray):
            # Log error if grid is missing or invalid
            self.logger.error("T_grid missing or invalid in results")
            # Create error message for user
            error_text = wx.StaticText(
                panel,
                label="Error: Temperature grid not available."
            )
            # Set bold font for emphasis
            error_text.SetFont(wx.Font(
                12, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD
            ))
            # Add error text to sizer with padding and centering
            main_sizer.Add(error_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            # Apply sizer and return error panel
            panel.SetSizer(main_sizer)
            return panel

        # Extract temperature grid (shape: num_layers, ny, nx)
        T_grid = self.results['T_grid']
        # Get grid dimensions
        num_layers = T_grid.shape[0] if T_grid.ndim == 3 else 1
        # Extract physical PCB dimensions from results (defaults if missing)
        width_mm = self.results.get('width', 200.0)    # PCB width in mm
        height_mm = self.results.get('height', 100.0)  # PCB height in mm
        # Extract grid spacing
        dx = self.results.get('dx', 1.0)  # Grid spacing in X direction (mm)
        dy = self.results.get('dy', 1.0)  # Grid spacing in Y direction (mm)

        # Use only Top (layer 0) and Bottom (last layer) for hotspots
        top_layer = T_grid[0]
        bottom_layer = T_grid[-1] if num_layers > 1 else top_layer
        ny, nx = top_layer.shape

        # Get hotspot threshold from results
        hotspot_threshold = self.results.get('hotspot_threshold', 85.0)

        # Create scrollable window with vertical scrollbar
        scroll = wx.ScrolledWindow(panel, style=wx.VSCROLL)
        scroll.SetScrollRate(20, 20)  # Set scroll speed
        scroll_sizer = wx.BoxSizer(wx.VERTICAL)

        # --- TOP LAYER HEATMAP ---
        # Create panel for Top layer
        layer_panel1 = wx.Panel(scroll)
        layer_sizer1 = wx.BoxSizer(wx.VERTICAL)
        # Create matplotlib figure and axis
        fig1, ax1 = plt.subplots(figsize=(5, 4))
        # Set global font size from config
        plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})
        # Set figure background to match dialog
        bg = self.GetBackgroundColour()
        fig1.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))
        # Set plot area to white
        ax1.set_facecolor('white')
        # Display heatmap with correct physical extent
        im1 = ax1.imshow(
            top_layer,
            cmap=config.DEFAULT_VISUALIZATION_HEAT_COLORMAP_STYLE,  # Use 'hot' colormap from config
            origin='lower',  # Origin at bottom-left (matches FDM grid)
            extent=[0, width_mm, 0, height_mm],  # Physical scale in mm
            interpolation='nearest'  # No smoothing for grid accuracy
        )
        # Find hotspots (above threshold from results)
        hotspots_top = np.where(top_layer > hotspot_threshold)
        hotspot_points1 = ax1.scatter([], [], c=config.DEFAULT_HOTSPOT_COLOR, marker='o', s=50, label='Hotspot')
        for y, x in zip(hotspots_top[0], hotspots_top[1]):
            ax1.plot(x * dx, y * dy, 'o', markersize=5, color=config.DEFAULT_HOTSPOT_COLOR)
        # Add colorbar with rounded integer ticks
        cbar1 = fig1.colorbar(im1, ax=ax1, label='Temperature (°C)')
        cbar1.ax.tick_params(labelsize=config.DEFAULT_CHART_FONT_SIZE)
        # Create 5 evenly spaced ticks
        cbar_ticks = np.linspace(np.min(top_layer), np.max(top_layer), 5)
        cbar_ticks_rounded = np.round(cbar_ticks).astype(int)
        cbar1.set_ticks(cbar_ticks)
        cbar1.set_ticklabels([f"{int(t)}" for t in cbar_ticks_rounded])
        # Set title above the plot (with padding)
        ax1.set_title("Top Layer - Hotspots", pad=15)
        # Set axis labels with units
        ax1.set_xlabel("X (mm)")
        ax1.set_ylabel("Y (mm)")
        # FINE GRID: Show lines at every dx/dy interval
        ax1.set_xticks(np.arange(0, width_mm + dx, dx))
        ax1.set_yticks(np.arange(0, height_mm + dy, dy))
        ax1.grid(True, which='both', color='gray', linestyle='-', linewidth=0.3, alpha=0.4)
        # TICK LABELS: Show only at 0%, 25%, 50%, 75%, 100% — rounded to 0 decimals
        x_positions = [0, width_mm * 0.25, width_mm * 0.5, width_mm * 0.75, width_mm]
        y_positions = [0, height_mm * 0.25, height_mm * 0.5, height_mm * 0.75, height_mm]
        ax1.set_xticks(x_positions)
        ax1.set_yticks(y_positions)
        ax1.set_xticklabels([f"{int(round(x))}" for x in x_positions])
        ax1.set_yticklabels([f"{int(round(y))}" for y in y_positions])
        # Optimize layout to prevent clipping
        plt.tight_layout()
        # Add legend with colored dot, positioned below X-axis
        ax1.legend([hotspot_points1], [f"Hotspots (> {hotspot_threshold}°C)"], loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1, fontsize=config.DEFAULT_CHART_FONT_SIZE)
        # Define filename for this layer's heatmap
        heatmap_filename1 = "fdm_hotspots_top.png"
        heatmap_path1 = os.path.join(config.DEFAULT_TEMP_DIR, heatmap_filename1)
        # Save figure to temporary directory with extra padding for legend
        fig1.savefig(heatmap_path1, format='png', dpi=100, bbox_inches='tight', pad_inches=0.2)
        # Close figure to free memory
        plt.close(fig1)
        # Load saved image as wx.Image
        img1 = wx.Image(heatmap_path1, wx.BITMAP_TYPE_PNG)
        # Convert to bitmap (fallback to 1x1 if invalid)
        bmp1 = wx.Bitmap(img1) if img1.IsOk() else wx.Bitmap(1, 1)
        # Create static bitmap to display heatmap
        bitmap1 = wx.StaticBitmap(layer_panel1, bitmap=bmp1)
        # Add bitmap to layer sizer with padding
        layer_sizer1.Add(bitmap1, 0, wx.ALL | wx.ALIGN_CENTER, 10)
        # Apply sizer to layer panel
        layer_panel1.SetSizer(layer_sizer1)
        # Add layer panel to scroll sizer
        scroll_sizer.Add(layer_panel1, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # --- BOTTOM LAYER HEATMAP ---
        # Create panel for Bottom layer
        layer_panel2 = wx.Panel(scroll)
        layer_sizer2 = wx.BoxSizer(wx.VERTICAL)
        # Create matplotlib figure and axis
        fig2, ax2 = plt.subplots(figsize=(5, 4))
        # Set global font size from config
        plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})
        # Set figure background to match dialog
        fig2.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))
        # Set plot area to white
        ax2.set_facecolor('white')
        # Display heatmap with correct physical extent
        im2 = ax2.imshow(
            bottom_layer,
            cmap=config.DEFAULT_VISUALIZATION_HEAT_COLORMAP_STYLE,  # Use 'hot' colormap from config
            origin='lower',  # Origin at bottom-left (matches FDM grid)
            extent=[0, width_mm, 0, height_mm],  # Physical scale in mm
            interpolation='nearest'  # No smoothing for grid accuracy
        )
        # Find hotspots (above threshold from results)
        hotspots_bottom = np.where(bottom_layer > hotspot_threshold)
        hotspot_points2 = ax2.scatter([], [], c=config.DEFAULT_HOTSPOT_COLOR, marker='o', s=50, label='Hotspot')
        for y, x in zip(hotspots_bottom[0], hotspots_bottom[1]):
            ax2.plot(x * dx, y * dy, 'o', markersize=5, color=config.DEFAULT_HOTSPOT_COLOR)
        # Add colorbar with rounded integer ticks
        cbar2 = fig2.colorbar(im2, ax=ax2, label='Temperature (°C)')
        cbar2.ax.tick_params(labelsize=config.DEFAULT_CHART_FONT_SIZE)
        # Create 5 evenly spaced ticks
        cbar_ticks = np.linspace(np.min(bottom_layer), np.max(bottom_layer), 5)
        cbar_ticks_rounded = np.round(cbar_ticks).astype(int)
        cbar2.set_ticks(cbar_ticks)
        cbar2.set_ticklabels([f"{int(t)}" for t in cbar_ticks_rounded])
        # Set title above the plot (with padding)
        ax2.set_title("Bottom Layer - Hotspots", pad=15)
        # Set axis labels with units
        ax2.set_xlabel("X (mm)")
        ax2.set_ylabel("Y (mm)")
        # FINE GRID: Show lines at every dx/dy interval
        ax2.set_xticks(np.arange(0, width_mm + dx, dx))
        ax2.set_yticks(np.arange(0, height_mm + dy, dy))
        ax2.grid(True, which='both', color='gray', linestyle='-', linewidth=0.3, alpha=0.4)
        # TICK LABELS: Show only at 0%, 25%, 50%, 75%, 100% — rounded to 0 decimals
        x_positions = [0, width_mm * 0.25, width_mm * 0.5, width_mm * 0.75, width_mm]
        y_positions = [0, height_mm * 0.25, height_mm * 0.5, height_mm * 0.75, height_mm]
        ax2.set_xticks(x_positions)
        ax2.set_yticks(y_positions)
        ax2.set_xticklabels([f"{int(round(x))}" for x in x_positions])
        ax2.set_yticklabels([f"{int(round(y))}" for y in y_positions])
        # Optimize layout to prevent clipping
        plt.tight_layout()
        # Add legend with colored dot, positioned below X-axis
        ax2.legend([hotspot_points2], [f"Hotspots (> {hotspot_threshold}°C)"], loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1, fontsize=config.DEFAULT_CHART_FONT_SIZE)
        # Define filename for this layer's heatmap
        heatmap_filename2 = "fdm_hotspots_bottom.png"
        heatmap_path2 = os.path.join(config.DEFAULT_TEMP_DIR, heatmap_filename2)
        # Save figure to temporary directory with extra padding for legend
        fig2.savefig(heatmap_path2, format='png', dpi=100, bbox_inches='tight', pad_inches=0.2)
        # Close figure to free memory
        plt.close(fig2)
        # Load saved image as wx.Image
        img2 = wx.Image(heatmap_path2, wx.BITMAP_TYPE_PNG)
        # Convert to bitmap (fallback to 1x1 if invalid)
        bmp2 = wx.Bitmap(img2) if img2.IsOk() else wx.Bitmap(1, 1)
        # Create static bitmap to display heatmap
        bitmap2 = wx.StaticBitmap(layer_panel2, bitmap=bmp2)
        # Add bitmap to layer sizer with padding
        layer_sizer2.Add(bitmap2, 0, wx.ALL | wx.ALIGN_CENTER, 10)
        # Apply sizer to layer panel
        layer_panel2.SetSizer(layer_sizer2)
        # Add layer panel to scroll sizer
        scroll_sizer.Add(layer_panel2, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Apply sizer to scroll window
        scroll.SetSizer(scroll_sizer)
        scroll.FitInside()  # Adjust scroll area to content
        # Add scroll window to main sizer
        main_sizer.Add(scroll, 1, wx.EXPAND | wx.ALL, 5)
        # Apply final sizer to panel
        panel.SetSizer(main_sizer)
        # Return completed Hotspots tab panel
        return panel
    
    def create_gradients_tab(self):
        """
        @brief Create the Gradients tab for FDM simulation results.
        @details Displays two separate gradient heatmaps for Top and Bottom layers only, stacked vertically.
                 Computes thermal gradients using np.gradient on T_grid.
                 Uses config.DEFAULT_VISUALIZATION_GRADIENT_COLORMAP_STYLE ('viridis') for visualization.
                 Saves each gradient map to config.DEFAULT_TEMP_DIR as 'fdm_gradients_top.png' and 'fdm_gradients_bottom.png'.
                 Adds colorbar with gradient magnitude in °C/mm, rounded to one decimal place.
                 Uses scrollable panel with centered content, styled like the heatmap tab.
        @return wx.ScrolledWindow with gradient visualization.
        """
        # Create main panel for the Gradients tab
        panel = wx.Panel(self.notebook)
        # Create vertical sizer to organize content
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        # Ensure temporary directory exists for saving heatmaps
        os.makedirs(config.DEFAULT_TEMP_DIR, exist_ok=True)

        # Verify that temperature grid is present and is a NumPy array
        if 'T_grid' not in self.results or not isinstance(self.results['T_grid'], np.ndarray):
            # Log error if grid is missing or invalid
            self.logger.error("T_grid missing or invalid in results")
            # Create error message for user
            error_text = wx.StaticText(
                panel,
                label="Error: Temperature grid not available."
            )
            # Set bold font for emphasis
            error_text.SetFont(wx.Font(
                12, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD
            ))
            # Add error text to sizer with padding and centering
            main_sizer.Add(error_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            # Apply sizer and return error panel
            panel.SetSizer(main_sizer)
            return panel

        # Extract temperature grid (shape: num_layers, ny, nx)
        T_grid = self.results['T_grid']
        # Get grid dimensions
        num_layers = T_grid.shape[0] if T_grid.ndim == 3 else 1
        # Extract physical PCB dimensions from results (defaults if missing)
        width_mm = self.results.get('width', 200.0)    # PCB width in mm
        height_mm = self.results.get('height', 100.0)  # PCB height in mm
        # Extract grid spacing
        dx = self.results.get('dx', 1.0)  # Grid spacing in X direction (mm)
        dy = self.results.get('dy', 1.0)  # Grid spacing in Y direction (mm)

        # Use only Top (layer 0) and Bottom (last layer) for gradients
        top_layer = T_grid[0]
        bottom_layer = T_grid[-1] if num_layers > 1 else top_layer
        ny, nx = top_layer.shape

        # Create scrollable window with vertical scrollbar
        scroll = wx.ScrolledWindow(panel, style=wx.VSCROLL)
        scroll.SetScrollRate(20, 20)  # Set scroll speed
        scroll_sizer = wx.BoxSizer(wx.VERTICAL)

        # --- TOP LAYER GRADIENT ---
        # Create panel for Top layer
        layer_panel1 = wx.Panel(scroll)
        layer_sizer1 = wx.BoxSizer(wx.VERTICAL)
        # Create matplotlib figure and axis
        fig1, ax1 = plt.subplots(figsize=(5, 4))
        # Set global font size from config
        plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})
        # Set figure background to match dialog
        bg = self.GetBackgroundColour()
        fig1.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))
        # Set plot area to white
        ax1.set_facecolor('white')
        # Compute gradients
        grad_y, grad_x = np.gradient(top_layer, dy, dx)
        gradient_magnitude = np.sqrt(grad_x**2 + grad_y**2)
        # Display gradient magnitude
        im1 = ax1.imshow(
            gradient_magnitude,
            cmap=config.DEFAULT_VISUALIZATION_GRADIENT_COLORMAP_STYLE,  # Use viridis colormap
            origin='lower',
            extent=[0, width_mm, 0, height_mm],
            interpolation='nearest'
        )
        # Add colorbar with gradient magnitude, rounded to one decimal place
        cbar1 = fig1.colorbar(im1, ax=ax1, label='Gradient Magnitude (°C/mm)')
        cbar1.ax.tick_params(labelsize=config.DEFAULT_CHART_FONT_SIZE)
        # Create 5 evenly spaced ticks
        cbar_ticks = np.linspace(np.min(gradient_magnitude), np.max(gradient_magnitude), 5)
        cbar_ticks_rounded = np.round(cbar_ticks, 1)  # Round to one decimal place
        cbar1.set_ticks(cbar_ticks)
        cbar1.set_ticklabels([f"{t:.1f}" for t in cbar_ticks_rounded])
        # Set title above the plot (with padding)
        ax1.set_title("Top Layer - Gradients", pad=15)
        # Set axis labels with units
        ax1.set_xlabel("X (mm)")
        ax1.set_ylabel("Y (mm)")
        # FINE GRID: Show lines at every dx/dy interval
        ax1.set_xticks(np.arange(0, width_mm + dx, dx))
        ax1.set_yticks(np.arange(0, height_mm + dy, dy))
        ax1.grid(True, which='both', color='gray', linestyle='-', linewidth=0.3, alpha=0.4)
        # TICK LABELS: Show only at 0%, 25%, 50%, 75%, 100% — rounded to 0 decimals
        x_positions = [0, width_mm * 0.25, width_mm * 0.5, width_mm * 0.75, width_mm]
        y_positions = [0, height_mm * 0.25, height_mm * 0.5, height_mm * 0.75, height_mm]
        ax1.set_xticks(x_positions)
        ax1.set_yticks(y_positions)
        ax1.set_xticklabels([f"{int(round(x))}" for x in x_positions])
        ax1.set_yticklabels([f"{int(round(y))}" for y in y_positions])
        # Optimize layout to prevent clipping
        plt.tight_layout()
        # Define filename for this layer's gradient
        gradient_filename1 = "fdm_gradients_top.png"
        gradient_path1 = os.path.join(config.DEFAULT_TEMP_DIR, gradient_filename1)
        # Save figure to temporary directory
        fig1.savefig(gradient_path1, format='png', dpi=100, bbox_inches='tight')
        # Close figure to free memory
        plt.close(fig1)
        # Load saved image as wx.Image
        img1 = wx.Image(gradient_path1, wx.BITMAP_TYPE_PNG)
        # Convert to bitmap (fallback to 1x1 if invalid)
        bmp1 = wx.Bitmap(img1) if img1.IsOk() else wx.Bitmap(1, 1)
        # Create static bitmap to display gradient
        bitmap1 = wx.StaticBitmap(layer_panel1, bitmap=bmp1)
        # Add bitmap to layer sizer with padding
        layer_sizer1.Add(bitmap1, 0, wx.ALL | wx.ALIGN_CENTER, 10)
        # Apply sizer to layer panel
        layer_panel1.SetSizer(layer_sizer1)
        # Add layer panel to scroll sizer
        scroll_sizer.Add(layer_panel1, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # --- BOTTOM LAYER GRADIENT ---
        # Create panel for Bottom layer
        layer_panel2 = wx.Panel(scroll)
        layer_sizer2 = wx.BoxSizer(wx.VERTICAL)
        # Create matplotlib figure and axis
        fig2, ax2 = plt.subplots(figsize=(5, 4))
        # Set global font size from config
        plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})
        # Set figure background to match dialog
        fig2.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))
        # Set plot area to white
        ax2.set_facecolor('white')
        # Compute gradients
        grad_y, grad_x = np.gradient(bottom_layer, dy, dx)
        gradient_magnitude = np.sqrt(grad_x**2 + grad_y**2)
        # Display gradient magnitude
        im2 = ax2.imshow(
            gradient_magnitude,
            cmap=config.DEFAULT_VISUALIZATION_GRADIENT_COLORMAP_STYLE,  # Use viridis colormap
            origin='lower',
            extent=[0, width_mm, 0, height_mm],
            interpolation='nearest'
        )
        # Add colorbar with gradient magnitude, rounded to one decimal place
        cbar2 = fig2.colorbar(im2, ax=ax2, label='Gradient Magnitude (°C/mm)')
        cbar2.ax.tick_params(labelsize=config.DEFAULT_CHART_FONT_SIZE)
        # Create 5 evenly spaced ticks
        cbar_ticks = np.linspace(np.min(gradient_magnitude), np.max(gradient_magnitude), 5)
        cbar_ticks_rounded = np.round(cbar_ticks, 1)  # Round to one decimal place
        cbar2.set_ticks(cbar_ticks)
        cbar2.set_ticklabels([f"{t:.1f}" for t in cbar_ticks_rounded])
        # Set title above the plot (with padding)
        ax2.set_title("Bottom Layer - Gradients", pad=15)
        # Set axis labels with units
        ax2.set_xlabel("X (mm)")
        ax2.set_ylabel("Y (mm)")
        # FINE GRID: Show lines at every dx/dy interval
        ax2.set_xticks(np.arange(0, width_mm + dx, dx))
        ax2.set_yticks(np.arange(0, height_mm + dy, dy))
        ax2.grid(True, which='both', color='gray', linestyle='-', linewidth=0.3, alpha=0.4)
        # TICK LABELS: Show only at 0%, 25%, 50%, 75%, 100% — rounded to 0 decimals
        x_positions = [0, width_mm * 0.25, width_mm * 0.5, width_mm * 0.75, width_mm]
        y_positions = [0, height_mm * 0.25, height_mm * 0.5, height_mm * 0.75, height_mm]
        ax2.set_xticks(x_positions)
        ax2.set_yticks(y_positions)
        ax2.set_xticklabels([f"{int(round(x))}" for x in x_positions])
        ax2.set_yticklabels([f"{int(round(y))}" for y in y_positions])
        # Optimize layout to prevent clipping
        plt.tight_layout()
        # Define filename for this layer's gradient
        gradient_filename2 = "fdm_gradients_bottom.png"
        gradient_path2 = os.path.join(config.DEFAULT_TEMP_DIR, gradient_filename2)
        # Save figure to temporary directory
        fig2.savefig(gradient_path2, format='png', dpi=100, bbox_inches='tight')
        # Close figure to free memory
        plt.close(fig2)
        # Load saved image as wx.Image
        img2 = wx.Image(gradient_path2, wx.BITMAP_TYPE_PNG)
        # Convert to bitmap (fallback to 1x1 if invalid)
        bmp2 = wx.Bitmap(img2) if img2.IsOk() else wx.Bitmap(1, 1)
        # Create static bitmap to display gradient
        bitmap2 = wx.StaticBitmap(layer_panel2, bitmap=bmp2)
        # Add bitmap to layer sizer with padding
        layer_sizer2.Add(bitmap2, 0, wx.ALL | wx.ALIGN_CENTER, 10)
        # Apply sizer to layer panel
        layer_panel2.SetSizer(layer_sizer2)
        # Add layer panel to scroll sizer
        scroll_sizer.Add(layer_panel2, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Apply sizer to scroll window
        scroll.SetSizer(scroll_sizer)
        scroll.FitInside()  # Adjust scroll area to content
        # Add scroll window to main sizer
        main_sizer.Add(scroll, 1, wx.EXPAND | wx.ALL, 5)
        # Apply final sizer to panel
        panel.SetSizer(main_sizer)
        # Return completed Gradients tab panel
        return panel
    
    def _point_in_polygon(self, point, poly):
        """
        @brief Checks if a point is inside a polygon using ray casting algorithm.
        @param point Tuple (x, y) in mm.
        @param poly List of tuples [(x1,y1), ..., (xn,yn)] defining the polygon in mm.
        @return True if point is inside, False otherwise.
        """
        x, y = point
        n = len(poly)
        inside = False
        p1x, p1y = poly[0]
        # Ray casting: count intersections with polygon edges
        for i in range(n + 1):
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    def create_zones_tab(self):
        """
        @brief Create the Zones tab for FDM simulation results.
        @details Displays two separate bar charts for Top and Bottom layers only, stacked vertically.
                 Computes average temperatures in zones defined by self.results['zone_polys'] using T_grid.
                 Uses style similar to create_summary_tab (bars with values, grid, Y-margin, named X-axis labels).
                 Saves each bar chart to config.DEFAULT_TEMP_DIR as 'fdm_zones_top.png' and 'fdm_zones_bottom.png'.
                 Adds colorbar with temperature in °C.
                 Uses scrollable panel with centered content, styled like the heatmap tab.
        @return wx.ScrolledWindow with zones visualization.
        """
        # Create main panel for the Zones tab
        panel = wx.Panel(self.notebook)
        # Create vertical sizer to organize content
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        # Ensure temporary directory exists for saving heatmaps
        os.makedirs(config.DEFAULT_TEMP_DIR, exist_ok=True)

        # Verify that temperature grid and zone_polys are present
        T_grid = self.results.get('T_grid', None)
        if T_grid is None or not isinstance(T_grid, np.ndarray):
            self.logger.error("Dialog updated with results and content: T_grid missing or invalid in results")
            error_text = wx.StaticText(
                panel,
                label="Error: Temperature grid not available."
            )
            error_text.SetFont(wx.Font(12, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
            main_sizer.Add(error_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            panel.SetSizer(main_sizer)
            return panel
        zone_polys = self.results.get('zone_polys', {})
        if not zone_polys:
            self.logger.error("Dialog updated with results and content: No zones defined in zone_polys")
            error_text = wx.StaticText(
                panel,
                label="Error: No zones defined."
            )
            error_text.SetFont(wx.Font(12, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
            main_sizer.Add(error_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
            panel.SetSizer(main_sizer)
            return panel

        # Extract temperature grid (shape: num_layers, ny, nx)
        num_layers = T_grid.shape[0] if T_grid.ndim == 3 else 1
        # Extract physical PCB dimensions from results (defaults if missing)
        width_mm = self.results.get('width', 200.0)    # PCB width in mm
        height_mm = self.results.get('height', 100.0)  # PCB height in mm
        # Extract grid spacing
        dx = self.results.get('dx', 1.0)  # Grid spacing in X direction (mm)
        dy = self.results.get('dy', 1.0)  # Grid spacing in Y direction (mm)

        # Use only Top (layer 0) and Bottom (last layer) for zones
        top_layer = T_grid[0]
        bottom_layer = T_grid[-1] if num_layers > 1 else top_layer
        ny, nx = top_layer.shape

        # Create scrollable window with vertical scrollbar
        scroll = wx.ScrolledWindow(panel, style=wx.VSCROLL)
        scroll.SetScrollRate(20, 20)  # Set scroll speed
        scroll_sizer = wx.BoxSizer(wx.VERTICAL)

        # Helper function to compute zone averages
        def compute_zone_averages(layer_data, layer_idx):
            zone_averages = {}
            for (zone_id, zone_layer), poly_mm in zone_polys.items():
                if zone_layer == layer_idx:
                    mask = np.zeros((ny, nx), dtype=bool)
                    for i in range(ny):
                        y = i * dy
                        for j in range(nx):
                            x = j * dx
                            if self._point_in_polygon((x, y), poly_mm):
                                mask[i, j] = True
                    temps_in_zone = layer_data[mask]
                    zone_averages[zone_id] = np.mean(temps_in_zone) if temps_in_zone.size > 0 else 0.0
            # Format log without np.float64, using float conversion
            formatted_averages = {k: round(float(v), 3) for k, v in zone_averages.items()}
            self.logger.debug("Zone averages for layer %d: %s", layer_idx, formatted_averages)
            return zone_averages

        # --- COMBINED ZONES ---
        # Create panel for combined zones
        layer_panel1 = wx.Panel(scroll)
        layer_sizer1 = wx.BoxSizer(wx.VERTICAL)
        # Create matplotlib figure and axis
        fig1, ax1 = plt.subplots(figsize=(4, 3))  # Same size as summary tab
        # Set global font size from config
        plt.rcParams.update({'font.size': config.DEFAULT_CHART_FONT_SIZE})
        # Match dialog background
        bg = self.GetBackgroundColour()
        fig1.patch.set_facecolor((bg.Red() / 255.0, bg.Green() / 255.0, bg.Blue() / 255.0))
        ax1.set_facecolor('white')
        # Compute zone averages for Top layer
        zone_averages_top = compute_zone_averages(top_layer, 0)
        zones_top = [f"Zone {zone_id}" for zone_id in zone_averages_top.keys()]  # Use zone IDs as labels
        averages_top = list(zone_averages_top.values())
        # Compute zone averages for Bottom layer
        zone_averages_bottom = compute_zone_averages(bottom_layer, num_layers - 1 if num_layers > 1 else 0)
        zones_bottom = [f"Zone {zone_id}" for zone_id in zone_averages_bottom.keys()]  # Use zone IDs as labels
        averages_bottom = list(zone_averages_bottom.values())
        # Combine top and bottom
        zones = zones_top + zones_bottom
        averages = averages_top + averages_bottom
        if not averages:
            self.logger.warning("Dialog updated with results and content: No valid zone averages")
            no_data_text = wx.StaticText(layer_panel1, label="No valid zone averages")
            no_data_text.SetFont(wx.Font(config.DEFAULT_RESULTS_FONT_SIZE, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
            layer_sizer1.Add(no_data_text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
        else:
            # Display bar chart
            bars1 = ax1.bar(zones, averages, color=config.DEFAULT_CHART_COLOR_SCHEME[0], zorder=3)
            # Add grids
            ax1.grid(True, which='major', linestyle='--', alpha=0.7, color='gray', zorder=0)
            ax1.grid(True, which='minor', linestyle=':', alpha=0.3, color='gray', zorder=0)
            ax1.minorticks_on()
            # Y-axis margin from config
            max_val = max(averages) if averages else 0
            margin_percent = config.DEFAULT_CHART_Y_MARGIN_PERCENT / 100.0
            if margin_percent > 0:
                y_margin = max_val * margin_percent
                ax1.set_ylim(0, max_val + y_margin)
                label_offset = y_margin * 0.3
            else:
                label_offset = max_val * 0.05
            # Rotate x-axis labels for readability
            ax1.set_xticklabels(zones, rotation=90)
            # Set title and label
            ax1.set_title("Zone Averages")
            ax1.set_ylabel("Temperature (°C)")
            # Add value labels on bars
            for bar, v in zip(bars1, averages):
                ax1.text(bar.get_x() + bar.get_width() / 2., bar.get_height() + label_offset,
                         f'{v:.2f}', ha='center', va='bottom')
            # Tight layout
            plt.tight_layout()
            # Define filename for this layer's zones
            zones_filename1 = "fdm_zones.png"
            zones_path1 = os.path.join(config.DEFAULT_TEMP_DIR, zones_filename1)
            # Save figure to temporary directory
            fig1.savefig(zones_path1, format='png', dpi=100, bbox_inches='tight')
            # Close figure to free memory
            plt.close(fig1)
            # Load saved image as wx.Image
            img1 = wx.Image(zones_path1, wx.BITMAP_TYPE_PNG)
            # Convert to bitmap (fallback to 1x1 if invalid)
            bmp1 = wx.Bitmap(img1) if img1.IsOk() else wx.Bitmap(1, 1)
            # Create static bitmap to display zones
            bitmap1 = wx.StaticBitmap(layer_panel1, bitmap=bmp1)
            # Add bitmap to layer sizer with padding
            layer_sizer1.Add(bitmap1, 0, wx.ALL | wx.ALIGN_CENTER, 10)
        # Apply sizer to layer panel
        layer_panel1.SetSizer(layer_sizer1)
        # Add layer panel to scroll sizer
        scroll_sizer.Add(layer_panel1, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Apply sizer to scroll window
        scroll.SetSizer(scroll_sizer)
        scroll.FitInside()  # Adjust scroll area to content
        # Add scroll window to main sizer
        main_sizer.Add(scroll, 1, wx.EXPAND | wx.ALL, 5)
        # Apply final sizer to panel
        panel.SetSizer(main_sizer)
        # Return completed Zones tab panel
        return panel
    
    def _create_content_sizer(self):
        """
        @brief Create loading content for the dialog (before results are ready).
        @return wx.BoxSizer with centered loading message.
        """
        content = wx.BoxSizer(wx.HORIZONTAL)
        loading = wx.StaticText(self, label="Simulation in progress. Please wait.")
        loading.SetFont(wx.Font(config.DEFAULT_RESULTS_FONT_SIZE,
                                wx.FONTFAMILY_SWISS,
                                wx.FONTSTYLE_NORMAL,
                                wx.FONTWEIGHT_NORMAL))
        content.Add(loading, 1, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL, 5)
        return content

    def _create_placeholder_tab(self, label):
        """
        @brief Create a placeholder panel for a tab.
        @param label Text to display.
        @return wx.Panel with centered text.
        """
        panel = wx.Panel(self.notebook)
        sizer = wx.BoxSizer(wx.VERTICAL)
        text = wx.StaticText(panel, label=f"[{label} tab content will be added here]")
        text.SetFont(wx.Font(10, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_NORMAL))
        sizer.Add(text, 0, wx.ALL | wx.ALIGN_CENTER, 20)
        panel.SetSizer(sizer)
        return panel
    
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
        # Stop any running animation and clean up Matplotlib figures
        if hasattr(self, 'notebook') and self.notebook:
            for page in self.notebook.GetChildren():
                if hasattr(page, 'anim_ctrl') and page.anim_ctrl:
                    page.anim_ctrl.Stop()
                    page.anim_ctrl.Destroy()
                # Fechar todas as figuras Matplotlib associadas a esta página
                import matplotlib.pyplot as plt
                plt.close('all')  # Fecha todas as figuras abertas
        # Destroy dialog first
        self.Destroy()
        
        # Delete temporary folder (commented out to persist temp dir until new simulation or app close)
        # wx.CallAfter(self._delete_temp_dir)

    def on_close_window(self, event):
        """
        @brief Handle window close event.
        @details Deletes the entire temp directory (config.DEFAULT_TEMP_DIR).
                 Uses wx.CallAfter to avoid file-in-use errors.
        """
        # Stop any running animation and clean up Matplotlib figures
        if hasattr(self, 'notebook') and self.notebook:
            for page in self.notebook.GetChildren():
                if hasattr(page, 'anim_ctrl') and page.anim_ctrl:
                    page.anim_ctrl.Stop()
                    page.anim_ctrl.Destroy()
                # Fechar todas as figuras Matplotlib associadas a esta página
                import matplotlib.pyplot as plt
                plt.close('all')  # Fecha todas as figuras abertas
        # Destroy dialog first
        self.Destroy()
        # Forçar o encerramento se for o último diálogo
        if not wx.GetApp().GetTopWindow():
            # Exit app
            wx.GetApp().Destroy()
        
        # Delete temporary folder (commented out to persist temp dir until new simulation or app close)
        # wx.CallAfter(self._delete_temp_dir)

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
    
    


