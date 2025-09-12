# ===================================================================
# PSA SIMULATOR - STANDALONE RUN SCRIPT
# ===================================================================
#
# This script runs a single PSA simulation based on settings
# defined in the accompanying `config.yaml` file.
#
# ===================================================================

println("\n" * "="^60)
println(" PSA SIMULATOR - CONFIGURABLE RUN ")
println("="^60 * "\n")

# ===================================================================
# SETUP AND IMPORTS
# ===================================================================

using Pkg
# Activate the project environment and add YAML dependency if needed
Pkg.activate(@__DIR__)

using PSASimulator
using DataFrames
using CSV
using Dates
using YAML

# Include the configuration data module
include("config.jl")
using .PSAConfigData

println("✓ All modules loaded successfully")

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

function get_material_index(material_name)
    for material in MATERIALS_LIST
        if material.name == material_name
            return material.index
        end
    end
    error("Material '$(material_name)' not found in the dataset.")
end

function get_material_data(material_name)
    idx = get_material_index(material_name)
    properties = SIMULATION_PARAMETERS[idx, :]
    isotherm = ISOTHERM_PARAMETERS[idx, :]
    return properties, isotherm
end

function get_opt_vars(material_name, scenario_key)
    idx = get_material_index(material_name)
    if scenario_key == "MaxPurity90"
        return OPT_VARS_PURITY[idx, :]
    elseif scenario_key == "MaxPurity95"
        return OPT_VARS_RECOVERY[idx, :]
    elseif scenario_key == "MaxProductivity"
        return OPT_VARS_PRODUCTIVITY[idx, :]
    elseif scenario_key == "MinEnergy"
        return OPT_VARS_ENERGY[idx, :]
    else
        error("Optimization scenario key '$(scenario_key)' not recognized.")
    end
end

function save_simulation_data(traj, material_name, scenario_name, N, timestamp_dir)
    if traj === nothing
        println("⚠️ No trajectory data to save.")
        return
    end
    println("  Saving data to: $(timestamp_dir)")

    headers = ["Time"]
    append!(headers, ["P_Node$(i)" for i in 1:N+2])
    append!(headers, ["y_Node$(i)" for i in 1:N+2])
    append!(headers, ["x1_Node$(i)" for i in 1:N+2])
    append!(headers, ["x2_Node$(i)" for i in 1:N+2])
    append!(headers, ["T_Node$(i)" for i in 1:N+2])

    step_map = [
        ("a_storage", "t1_storage", "1_Co-current_Pressurization"),
        ("b_storage", "t2_storage", "2_Adsorption"),
        ("c_storage", "t3_storage", "3_Heavy_Reflux"),
        ("d_storage", "t4_storage", "4_Counter-current_Depressurization"),
        ("e_storage", "t5_storage", "5_Light_Reflux")
    ]

    if !haskey(traj, :a_storage) || isempty(traj[:a_storage])
        println("⚠️ Trajectory history not found in results. Nothing to save.")
        return
    end

    num_cycles = length(traj[:a_storage])
    println("  Found data for $(num_cycles) cycles.")

    for cycle_idx in 1:num_cycles
        cycle_dir = joinpath(timestamp_dir, "cycle_$(cycle_idx)")
        mkpath(cycle_dir)

        for (data_key, time_key, step_name) in step_map
            data_storage = traj[Symbol(data_key)]
            time_storage = traj[Symbol(time_key)]

            if length(data_storage) >= cycle_idx && length(time_storage) >= cycle_idx
                data = data_storage[cycle_idx]
                time = time_storage[cycle_idx]
                
                df = DataFrame(hcat(time, data), headers)
                filename = "$(step_name).csv"
                CSV.write(joinpath(cycle_dir, filename), df)
            end
        end
    end
end

# ===================================================================
# MAIN EXECUTION
# ===================================================================

function execute_simulation()
    # --- 1. LOAD CONFIGURATION ---
    config = YAML.load_file("config.yaml")
    sim_settings = config["simulation_settings"]
    proc_vars_config = config["process_variables"]
    fault_settings = config["fault_injection"]
    custom_proc_vars = get(config, "custom_process_variables", nothing)

    # Extract parameters from config
    N = sim_settings["N"]
    material_name = sim_settings["material_name"]
    scenario_name = sim_settings["scenario_name"]
    max_iterations = sim_settings["max_iterations"]
    run_type = "ProcessEvaluation" # This could also be in the config

    # --- LOAD FAULT SCENARIO ---
    fault_scenario_name = get(sim_settings, "fault_scenario", "None")
    overrides = Dict()
    if fault_scenario_name != "None" && fault_scenario_name != ""
        println("\n--- LOADING FAULT SCENARIO: $(fault_scenario_name) ---")
        fault_scenarios_file = "fault_scenarios.yaml"
        if isfile(fault_scenarios_file)
            all_faults = YAML.load_file(fault_scenarios_file)
            if haskey(all_faults, fault_scenario_name)
                fault_data = all_faults[fault_scenario_name]
                # The overrides are the fault data itself, minus the description.
                overrides = filter(p -> p.first != "description", fault_data)
                if !isempty(overrides)
                    println("  ✓ Applied parameter overrides for fault scenario.")
                else
                    println("  - No applicable parameter overrides found for this fault.")
                end
            else
                println("  ⚠️  WARNING: Fault scenario '$(fault_scenario_name)' not found in '$(fault_scenarios_file)'.")
            end
        else
            println("  ⚠️  WARNING: '$(fault_scenarios_file)' not found. Cannot apply fault.")
        end
    end


    # --- Create Timestamped Directory ---
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    base_dir = "simulation_output"
    material_dir = joinpath(base_dir, replace(material_name, r"[/: ]" => "_"))
    scenario_dir = joinpath(material_dir, replace(scenario_name, r"[/: ]" => "_"))
    timestamp_dir = joinpath(scenario_dir, timestamp)
    mkpath(timestamp_dir)

    log_file_path = joinpath(timestamp_dir, "simulation.log")
    open(log_file_path, "w") do log_file_stream
        original_stdout = stdout
        redirect_stdout(log_file_stream)
        try
            println("="^60)
            println(" PSA SIMULATOR RUN ")
            println("="^60)
            println("Timestamp: $(timestamp)")
            println("\n--- SIMULATION CONFIGURATION ---")
            for (key, value) in sim_settings
                println("  $(key): $(value)")
            end
            println("\n--- PROCESS VARIABLES (from config.yaml) ---")
            for (key, value) in proc_vars_config
                println("  $(key): $(value)")
            end
            if custom_proc_vars !== nothing
                println("\n--- CUSTOM PROCESS VARIABLES (from config.yaml) ---")
                for (key, value) in custom_proc_vars
                    println("  $(key): $(value)")
                end
            end

            # Log the applied fault overrides
            if !isempty(overrides)
                println("\n--- APPLIED FAULT PARAMETER OVERRIDES ---")
                for (key, value) in overrides
                    println("  $(key): $(value)")
                end
            end

            println("\n--- LEGACY FAULT INJECTION (from config.yaml) ---")
            for (key, value) in fault_settings
                println("  $(key): $(value)")
            end
            println("\n" * "="^60 * "\n")

            # --- 2. DETERMINE PROCESS VARIABLES TO USE ---
            local process_vars # Declare local to avoid scope issues
            local material_data # NEW: Declare material_data here to ensure it's always defined

            # Always get material_data, as it's needed by psacycle regardless of process_vars source
            material_properties, isotherm_params = get_material_data(material_name)
            material_data = (material_properties, isotherm_params)

            if custom_proc_vars !== nothing &&
               get(custom_proc_vars, "enabled", false) &&
               haskey(custom_proc_vars, "P_0") &&
               haskey(custom_proc_vars, "t_ads") &&
               haskey(custom_proc_vars, "alpha") &&
               haskey(custom_proc_vars, "beta") &&
               haskey(custom_proc_vars, "P_l")

                println("  Using custom process variables from config.yaml")
                # Construct process_vars from custom_proc_vars
                process_vars = [
                    proc_vars_config["bed_length"],            # L [m]
                    custom_proc_vars["P_0"],                   # P_0 [Pa]
                    custom_proc_vars["P_0"] * custom_proc_vars["beta"] / 8.314 / 313.15, # n_dot_0 [mol/s]
                    custom_proc_vars["t_ads"],                 # t_ads [s]
                    custom_proc_vars["alpha"],                 # alpha [-]
                    custom_proc_vars["beta"],                  # beta [-]
                    proc_vars_config["intermediate_pressure"], # P_I [Pa]
                    custom_proc_vars["P_l"]                    # P_l [Pa]
                ]
            else
                println("  Using process variables from optimization scenario lookup.")
                # Existing logic: Get optimization variables
                optimization_scenario = proc_vars_config["optimization_scenario"] # Get from new location
                opt_vars = get_opt_vars(material_name, optimization_scenario)

                process_vars = [
                    proc_vars_config["bed_length"],            # L [m]
                    opt_vars[1],                                # P_0 [Pa]
                    opt_vars[1] * opt_vars[4] / 8.314 / 313.15, # n_dot_0 [mol/s]
                    opt_vars[2],                                # t_ads [s]
                    opt_vars[3],                                # alpha [-]
                    opt_vars[5],                                # beta [-]
                    proc_vars_config["intermediate_pressure"],# P_I [Pa]
                    opt_vars[6]                                 # P_l [Pa]
                ]
            end

            # --- 3. RUN THE SIMULATION ---
            println("🚀 Running simulation... (This may take a moment)")
            result = PSASimulator.psacycle(process_vars, material_data; N=N, run_type=Symbol(run_type), it_disp=true, max_iters=max_iterations, overrides=overrides)

            # --- 4. LOG PARAMETERS ---
            if isdefined(result, :Params)
                println("\n" * "="^60)
                println(" LOGGING DIFFERENTIAL EQUATION PARAMETERS")
                println("="^60)

                Params = result.Params
                param_labels = [
                    "N", "ΔU₁", "ΔU₂", "ρ_s", "T₀", "ε", "r_p", "μ", "R", "v₀",
                    "q_s0", "C_pg", "C_pa", "C_ps", "D_m", "K_z", "P₀", "L", "MW_CO2", "MW_N2",
                    "k_CO2_LDF", "k_N2_LDF", "y₀", "τ", "P_l", "P_inlet", "y_LP", "T_LP", "ṅ_LP",
                    "α", "β", "P_I", "y_HR", "T_HR", "ṅ₀ * β", "y_LR_guess", "T_LR_guess", "ṅ₀", "feed_gas_type"
                ]

                println("\n--- DIFFERENTIAL EQUATION PARAMETERS (`Params` Vector) ---")
                for i in 1:length(Params)
                    println("  [$(i)] $(param_labels[i]): $(Params[i])")
                end

                # Calculate and log the final dimensionless parameters
                v0 = Params[10]
                L = Params[18]
                Dm = Params[15]
                rp = Params[7]
                R = Params[9]
                T0 = Params[5]
                q_s0 = Params[11]
                epsilon = Params[6]
                P0 = Params[17]
                k_CO2_LDF = Params[21]
                k_N2_LDF = Params[22]

                # Avoid division by zero if v0 is zero
                if v0 > 0
                    Pe = v0 * L / (0.7 * Dm + v0 * rp)
                    phi = R * T0 * q_s0 * (1 - epsilon) / (epsilon * P0)
                    k1 = k_CO2_LDF * L / v0
                    k2 = k_N2_LDF * L / v0

                    println("\n--- KEY DIMENSIONLESS PARAMETERS ---")
                    println("  Péclet Number (Pe): $(Pe)")
                    println("  Capacity Ratio (phi): $(phi)")
                    println("  Dim. LDF Coeff (k1 - CO2): $(k1)")
                    println("  Dim. LDF Coeff (k2 - N2): $(k2)")
                else
                    println("\n--- KEY DIMENSIONLESS PARAMETERS ---")
                    println("  v₀ is zero, dimensionless parameters cannot be calculated.")
                end
                println("\n" * "="^60)
            end


            # --- 5. SAVE THE RESULTS ---
            if result.traj !== nothing
                redirect_stdout(original_stdout) # Switch back to console for user messages
                println("\n💾 Simulation finished. Saving data...")
                save_simulation_data(result.traj, material_name, scenario_name, N, timestamp_dir)
                println("\n" * "="^60)
                println("✅ SIMULATION COMPLETE & DATA SAVED")
                println("="^60 * "\n")
            else
                redirect_stdout(original_stdout)
                println("❌ SIMULATION FAILED. No data to save.")
            end
        catch e
            redirect_stdout(original_stdout)
            println(stderr, "An error occurred during simulation: ", e)
            showerror(stderr, e, catch_backtrace())
        finally
            redirect_stdout(original_stdout) # Ensure it's always restored
        end
    end
    println("✓ Simulation log saved to: $(log_file_path)")
end

# --- RUN THE SCRIPT ---
execute_simulation()
