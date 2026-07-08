module LHPTransientModel

using CSV, ComponentArrays, CoolProp, DataFrames, DifferentiationInterface, LinearAlgebra, NonlinearSolve, Parameters, Plots, Revise, StaticArrays, Interpolations, Statistics

export LHPGeometry, WickGeometry, GrvGeometry, TLGeometry, CondGeometry, CCGeometry, LHPTransientCondition, FluidLookups
export init_fluid_lookups
export solve_transient_lhp!

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