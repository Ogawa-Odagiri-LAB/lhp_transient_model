module LHP_steady_Model

using CSV, ComponentArrays, CoolProp, DataFrames, DifferentiationInterface, LinearAlgebra, NonlinearSolve, Parameters, Plots, Revise, StaticArrays

export LHPGeometry, GrvGeometry, TLGeometry, CondGeometry, CCGeometry, LHPCondition
export initialize_steady_vars, solve_steady_lhp!, f_cc_void_solve, save_vars_to_csv

# ファイルを役割ごとに分割して読み込み
include("types.jl")       # 構造体の定義
include("physics.jl")     # 圧損や物性値の計算（旧metal_prop.jlの内容もここ）
include("steady_components.jl")  # calc_grv!, calc_vl! などのソルバー
include("utils.jl")
include("inits.jl")

end