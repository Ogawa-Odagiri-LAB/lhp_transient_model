module LHPSteadyModel

using CSV, ComponentArrays, CoolProp, DataFrames, DifferentiationInterface, LinearAlgebra, NonlinearSolve, Parameters, Plots, Revise, StaticArrays, Interpolations

export LHPSteadyGeometry, WickGeometry, GrvGeometry, TLGeometry, CondGeometry, CCGeometry, LHPSteadyCondition
export init_fluid_lookups
export initialize_steady_vars, solve_steady_lhp!, f_cc_void_solve, save_steady_var_to_csv

# ファイルを役割ごとに分割して読み込み
include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "props.jl"))
include(joinpath(@__DIR__, "wallfuncs.jl"))
include(joinpath(@__DIR__, "residualfuncs.jl"))
include(joinpath(@__DIR__, "fluidfuncs.jl"))
include(joinpath(@__DIR__, "savefuncs.jl"))
include(joinpath(@__DIR__, "solverfuncs.jl"))
include(joinpath(@__DIR__, "steadyfuncs.jl"))
include(joinpath(@__DIR__, "debugfuncs.jl"))

end