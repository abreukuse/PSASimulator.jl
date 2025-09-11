<img src="assets/logo.png" alt="PSA Process Diagram" width="300">

#

`PSASimulator.jl` is a Julia-based Pressure Swing Adsorption (PSA) simulation package for modeling gas separation processes, particularly CO₂ capture. This package is a Julia translation of the [original MATLAB-based simulator](https://github.com/PEESEgroup/PSA) at this moment, with additional features and improvements planned in the future.

If you find this package helpful, please cite: [Yin, Xiangyu, and Chrysanthos E. Gounaris. "Computational discovery of Metal–Organic Frameworks for sustainable energy systems: Open challenges." Computers & Chemical Engineering 167 (2022): 108022.](https://www.sciencedirect.com/science/article/pii/S0098135422003568)

## Installation

To install `PSASimulator.jl`, first ensure you have Julia 1.6 or higher installed. Then:

```julia
using Pkg

# Install directly from GitHub
Pkg.add(url="https://github.com/xyin-anl/PSASimulator.jl")

# Or for development (if you want to modify the code):
Pkg.develop(url="https://github.com/xyin-anl/PSASimulator.jl")

# Or if you have cloned the repository locally:
Pkg.develop(path="/path/to/local/PSASimulator.jl")
```

Dependencies:
- `LinearAlgebra`, `SparseArrays`, `Statistics` - Standard libraries
- `DifferentialEquations.jl` - Tools for differential equations
- `OrdinaryDiffEq.jl` - ODE solvers
- `DataFrames.jl` - Data manipulation
- `PrettyTables.jl` - Formatted output display

## Configuration

The simulation parameters are now managed through `config.yaml` in the project root. This file allows for easy adjustment of simulation settings without modifying the source code.

A typical `config.yaml` structure includes:

```yaml
simulation_settings:
  material_name: "Zeolite_13X"
  scenario_name: "Normal_Operation"
  N: 10
  max_iterations: 700
  optimization_scenario: "MaxPurity90"

process_variables:
  bed_length: 1.0
  intermediate_pressure: 10000.0

fault_injection:
  enabled: false
  type: None # e.g., Adsorbent_Degradation, Column_Leak, Stuck_Valve
  parameters: {}
```

## Basic Usage (Updated)

To run a simulation, configure your parameters in `config.yaml` and then execute `run.jl`:

```julia
julia run.jl
```

The script will automatically log the simulation progress and save all results in a structured output directory.

## Data Organization

Simulation results, including detailed cycle data and logs, are saved in the `simulation_output/` directory. The structure is as follows:

```
simulation_output/
└── <MaterialName>/
    └── <ScenarioName>/
        └── <Timestamp>/
            ├── simulation.log
            ├── cycle_1/
            │   ├── 1_Co-current_Pressurization.csv
            │   ├── 2_Adsorption.csv
            │   └── ...
            ├── cycle_2/
            │   ├── ...
            └── cycle_N/  (The final, converged cycle or max_iterations)
                ├── ...
```

Each `cycle_N` subdirectory contains CSV files for each step of that specific simulation cycle.

## Demo

A comprehensive demo is provided in the `demo/` directory that validates the simulator against literature data from [Yancy-Caballero et al. (2020)](https://pubs.rsc.org/en/content/articlelanding/2020/me/d0me00060d).

To run the demo, navigate to the `demo` directory:
```bash
cd demo
```

**Running All Scenarios**

To run all four validation scenarios, execute the script without arguments:
```bash
julia demo_psa_simulator.jl
```

**Running a Specific Scenario**

You can also run a single scenario by providing its number (1-4) as a command-line argument.

For example, to run only the third scenario (Economic - Max Productivity):
```bash
julia demo_psa_simulator.jl 3
```

Note: The demo script automatically activates the parent project to ensure all dependencies are available.

The demo tests four simulation scenarios:
1. Maximized purity with 90% CO₂ recovery constraint
2. Maximized purity with 95% CO₂ recovery constraint
3. Maximized CO₂ productivity [mol/kg/hr]
4. Minimized energy consumption [kWh/ton CO₂]

The demo includes 16 different adsorbent materials:
- Metal-Organic Frameworks (MOFs): Co-MOF-74, Cu-BTTri, Mg-MOF-74, MOF-177, etc.
- Zeolites: Zeolite 13X
- Other porous materials: ZIF-8, SIFSIX series

## Mathematical Model

`PSASimulator.jl` solves a system of Partial Differential Equations (PDEs) that describe the dynamics of gas separation in a packed bed column. The model uses dimensionless variables with spatial discretization via finite volumes and WENO schemes.

### Governing Equations

**Component Mass Balance (CO₂ mole fraction)**
```math
\begin{aligned}
\frac{\partial y}{\partial \tau} &= \frac{1}{\mathrm{Pe}}\left(\frac{\partial^2 y}{\partial z^2} + \frac{\partial y}{\partial z}\,\frac{\partial P}{\partial z}\,\frac{1}{P} - \frac{\partial y}{\partial z}\,\frac{\partial T}{\partial z}\,\frac{1}{T}\right)\\
&\quad -\frac{T}{P}\,\frac{\partial}{\partial z}\left(\frac{y P v}{T}\right)
+ \phi\,\frac{T}{P}\Bigl[(y-1)\frac{\partial x_1}{\partial \tau} + y\frac{\partial x_2}{\partial \tau}\Bigr]
\end{aligned}
```

**Total Mass Balance**
```math
\frac{\partial P}{\partial \tau} = -T\,\frac{\partial}{\partial z}\left(\frac{P v}{T}\right) - \phi T\left(\frac{\partial x_1}{\partial \tau} + \frac{\partial x_2}{\partial \tau}\right) + \frac{P}{T}\,\frac{\partial T}{\partial \tau}.
```

**Energy Balance**
```math
\begin{aligned}
\frac{\partial T}{\partial \tau} &= \frac{K_z}{v_0 L}\,\frac{1}{\zeta}\,\frac{\partial^2 T}{\partial z^2}\\
&\quad - \frac{\varepsilon C_{pg} P_0}{R T_0}\,\frac{1}{\zeta}\left[\frac{\partial (P v)}{\partial z} - T\,\frac{\partial}{\partial z}\left(\frac{P v}{T}\right)\right]\\\
&\quad + \frac{(1-\varepsilon) q_{s0}}{T_0}\,\frac{1}{\zeta}\Bigl[(-\Delta U_1 + R T_0 T)\frac{\partial x_1}{\partial \tau} + (-\Delta U_2 + R T_0 T)\frac{\partial x_2}{\partial \tau}\Bigr]
\end{aligned}
```
where \(\zeta = (1-\varepsilon)(\rho_s C_{ps} + q_{s0} C_{pa}) + \varepsilon\rho_g C_{pg}\) is the effective heat‐capacity term.

**Adsorption Kinetics (Linear Driving Force)**
```math
\frac{\partial x_i}{\partial \tau} = k_i\left(\frac{q_i^{\ast}}{q_{s0}} - x_i\right), \qquad i = 1,2.
```

**Momentum Balance (Ergun Equation)**
```math
-\frac{\partial P}{\partial z}\,\frac{P_0}{L} = \frac{150\,\mu(1-\varepsilon)^2}{4 r_p^{\,2} \varepsilon^{\,2}}\,v + \frac{1.75(1-\varepsilon)\,\rho_g\,\mathrm{MW}}{2 r_p \varepsilon}\,v^2.
```

### Equilibrium Model

The equilibrium loadings \(q_1^{\ast}\) and \(q_2^{\ast}\) follow a dual–site Langmuir isotherm:
```math
q_i = q_{s,b,i}\,\frac{B_i C_i}{1 + B_1 C_1 + B_2 C_2} + q_{s,d,i}\,\frac{D_i C_i}{1 + D_1 C_1 + D_2 C_2}.
```
with temperature–dependent adsorption constants
```math
B_i = b_i\,e^{-\Delta U_{b,i}/(R T)}, \qquad D_i = d_i\,e^{-\Delta U_{d,i}/(R T)}.
```

### Dimensionless Variables
- Pressure: $P = P_{\text{actual}}/P_0$
- Temperature: $T = T_{\text{actual}}/T_0$
- Velocity: $v = v_{\text{actual}}/v_0$
- Adsorbed amount: $x_i = q_i/q_{s0}$
- Spatial coordinate: $z = z_{\text{actual}}/L$
- Time: $\tau = t v_0 / L$

### Key Parameters
- $\mathrm{Pe} = v_0 L / D_l$ with $D_l = 0.7 D_m + v_0 r_p$ (axial dispersion)
- $\phi = R T_0 q_{s0} (1-\varepsilon) / (\varepsilon P_0)$ (capacity ratio)
- $k_i = k_{i,\text{LDF}} L / v_0$ (dimensionless LDF coefficient)

### Process Steps

The simulator models six process steps, each with specific boundary conditions:

- **Adsorption** - Feed gas at high pressure, CO₂ adsorption
- **Heavy Reflux** - Heavy product recycle for CO₂ recovery enhancement
- **Counter-current Depressurization** - Pressure reduction from feed end
- **Light Reflux** - Light product purge for purity improvement
- **Co-current Depressurization** - Pressure reduction from product end
- **Co-current Pressurization** - Re-pressurization with light product

### Numerical Solution

The PDEs are solved using:
- **Spatial discretization**: Finite volume method with WENO (Weighted Essentially Non-Oscillatory) schemes
- **Temporal integration**: Adaptive ODE solvers from `DifferentialEquations.jl`

The simulator iterates through cycles until the state variables at the beginning and end of a cycle converge, indicating cyclic steady state has been reached.

## Modifications and Enhancements

This repository, originally a Julia translation of the [MATLAB-based PSA simulator by PEESEgroup](https://github.com/PEESEgroup/PSA) (translated by [xyin-anl](https://github.com/xyin-anl/PSASimulator.jl)), has undergone several enhancements to improve its usability and extend its capabilities. These modifications were developed with the assistance of a large language model (Gemini).

Key enhancements include:

*   **Configuration Management:** Introduced `config.yaml` for externalizing simulation parameters, allowing for flexible scenario definition without code changes.
*   **Modular Data Handling:** Extracted all hardcoded simulation constants and optimization tables into a dedicated `config.jl` file, improving code cleanliness and maintainability.
*   **Comprehensive Logging:** Implemented detailed logging for each simulation run, capturing configuration details and progress in `simulation.log` files.
*   **Full Trajectory Data:** The simulator now saves data for every cycle, enabling analysis of transient behavior leading up to cyclic steady state.
*   **Fault Injection Framework (In Progress):** Prepared the codebase for configurable fault injection mechanisms, allowing for systematic study of system anomalies.
