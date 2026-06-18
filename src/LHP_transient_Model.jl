module LHP_transient_Model

using CSV, ComponentArrays, CoolProp, DataFrames, DifferentiationInterface, LinearAlgebra, NonlinearSolve, Parameters, Plots, Revise, StaticArrays, Interpolations

export LHPGeometry, WickGeometry, GrvGeometry, TLGeometry, CondGeometry, CCGeometry, LHPTransientCondition, PropyleneLookups
export initialize_transient_vars, load_initial_vars_from_csv, prepare_f_res
export solve_transient_lhp!
export init_propylene_lookups

# ファイルを役割ごとに分割して読み込み
include("types.jl")       # 構造体の定義
include("physics.jl")     # 圧損や物性値の計算（旧metal_prop.jlの内容もここ）
include("transient_components.jl")  # 
include("utils.jl")
include("inits.jl")
include("prepare_props.jl")

end