# 🔥 PCB Thermal Simulator Tool

## ⚙️ Overview

The **PCB Thermal Simulator Tool** is a graphical application developed for **simulating the thermal behavior** in printed circuit boards (**PCBs**).

It allows users to draw and configure PCB elements (such as **zones**, **heaters**, **vias**, and **sinks**), define detailed parameters, and run simulations using three distinct methodologies: **Analytic**, **Lumped Parameter Thermal Network (LPTN)**, and **Finite Difference Method (FDM)**.

The goal is to provide a robust and accessible tool for thermal analysis, complete with clear and detailed visual results.

***

### ✨ Technology Stack

| Technology | Purpose |
| :--- | :--- |
| **Python** | Primary development language. |
| **wxPython** | Used for the cross-platform Graphical User Interface (GUI). |
| **NumPy** | Utilized for all efficient numerical computations and operations. |
| **Matplotlib** | Used for data visualization and generating plots (heatmaps, etc.). |

### 📋 Project Information

| Detail | Value |
| :--- | :--- |
| **Version** | 1.0.0 |
| **Author** | Pedro Silva |
| **Company** | PMS Labs |
| **License** | (C) 2025 Pedro Silva. All rights reserved. (See `LICENSE` file for details) |

***

## 🚀 Key Features

The simulator is packed with features to cover the entire thermal analysis workflow, from modeling to result visualization.

### 🎨 Interactive Graphical Editor

* **Element Drawing:** Draw polygons for **zones**, **heaters**, and **sinks**, and place **vias** directly on the canvas.
* **Navigation:** Full support for **zooming**, **panning**, and grid overlays (including the solver grid).
* **Visual Reference:** Load and overlay PCB images for easier and more accurate modeling.
* **Complete Editing:** Functionalities to edit, duplicate, move, and delete elements quickly.
* **Geometric Setup:** Define exact PCB dimensions using corner points.

### 📝 Detailed Parameter Configuration

* **Global PCB Parameters:** Setup for width, height, thickness, number of layers, copper/FR4 conductivity, ambient temperature, and convection coefficients.
* **Element-Specific Parameters:** Power for heaters, diameter for vias, and extra Heat Transfer Coefficient (HTC) for sinks.
* **Solver Settings:** Choice of method (Analytic/Lumped/FDM), model mode (Basic/Advanced), grid size, tolerance, max iterations, and Successive Over-Relaxation (SOR) factor.

### 🔬 Simulation Solvers

| Solver | Description | Use Cases |
| :--- | :--- | :--- |
| **Analytic** | Quick, closed-form approximations for simple estimations. | Initial estimates and rapid validation. |
| **Lumped (LPTN)** | **Thermal Resistance Network** model. Supports multi-layer PCBs and basic/advanced area calculations. | Medium-complexity multi-layer PCB thermal analysis. |
| **FDM** | **Finite Difference Method** for high-resolution analysis. Uses **SOR** solver, including multi-layer support, boundary conditions, and convergence checks. | High-precision and high-resolution simulations. |

### 📊 Visualization and Results Export

* **Reports:** Generation of summary reports detailing temperatures, $\Delta T$, thermal resistances, and key metrics.
* **Visualization:** Display of **Heatmaps**, **Time Evolution GIFs** of temperature, **hotspot detection**, thermal gradients, and zone-specific charts.
* **Export:** Export results as **PNG** images, **TXT** reports, and **GIF** animations.

### 🛠️ Additional Tools

* **Configuration:** **JSON** import/export for PCB configurations.
* **Logging:** Robust logging system for debugging and error tracking.
* **File Management:** Temporary file management for outputs.

***

## 💻 Requirements and Installation

### 📚 Core Requirements

The following dependencies are essential for running the simulator:

* **Python 3.8+**
* `wxPython` (For the GUI)
* `NumPy` (For numerical computations)
* `Matplotlib` (For visualizations)

You can install the core dependencies using `pip`:

### ⚡ Optional Dependencies

These libraries are recommended for advanced features or performance enhancements:

* **Numba** (`pip install numba`): Highly recommended for **FDM performance optimization**.
* **SciPy** (Requires `pip install scipy`): Included for potential extension into advanced features in future versions.

