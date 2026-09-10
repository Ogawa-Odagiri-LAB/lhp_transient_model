module LHPDebugModel

using CSV, ComponentArrays, CoolProp, DataFrames, DifferentiationInterface, LinearAlgebra, NonlinearSolve, Parameters, Plots, Revise, StaticArrays, Interpolations, Statistics

export LHPGeometry, WickGeometry, GrvGeometry, TLGeometry, CondGeometry, CCGeometry, LHPTransientCondition, PropyleneLookups, AmmoniaLookups, FluidLookups
export init_propylene_lookups,init_ammonia_lookups, init_fluid_lookups
export initialize_transient_var, set_values_on_var!, initialize_T_wall, set_values_on_T_wall!, initialize_f, set_values_on_var_from_steady, set_values_on_T_wall_from_steady, set_values_on_lb!, set_values_on_ub!
export Euler_wall!, RK4_wall!
export fluid_wick, fluid_grv, fluid_vl!, fluid_cond!, fluid_ll!, fluid_cc
export save_T_wall_to_csv, save_var_to_csv, save_f_to_csv
export λ_plotter, show_sat
export restore_var, scale_var, f_transient_solve!

# ファイルを役割ごとに分割して読み込み
include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "props.jl"))
include(joinpath(@__DIR__, "wallfuncs.jl"))
include(joinpath(@__DIR__, "residualfuncs.jl"))
include(joinpath(@__DIR__, "fluidfuncs.jl"))
include(joinpath(@__DIR__, "savefuncs.jl"))
include(joinpath(@__DIR__, "solverfuncs.jl"))
include(joinpath(@__DIR__, "debugfuncs.jl"))

end