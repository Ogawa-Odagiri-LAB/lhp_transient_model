using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート
using Revise                # 必ず最初のほうで Revise を読み込む
using ComponentArrays

# 自作モジュールをロードするための設定
# push!(LOAD_PATH, joinpath(@__DIR__, "../src"))

# 個別のファイル（types.jlなど）ではなく、モジュールのメインファイルを直接インクルード（Track）する
includet("../src/LHPDebugModel.jl") 
using .LHPDebugModel  # 現在の環境に読み込むため、ドット(.)をつけます

transient_condition = LHPTransientCondition(T_cond = 15.0 + 273.15, Q_apply_init = 50.0, Q_apply = 60.0, ΔT_startup = 0.0, dt = 0.1, t_init = 0.0, t_end = 0.01)

geo = LHPGeometry(
    WickGeometry(N=3),
    GrvGeometry(N=1),                                     # グルーブ
    TLGeometry(N=200, L=2, D_i=0.00435, D_o=0.00635),      # 蒸気管
    CondGeometry(N=100, L=1.0, D_i=0.00435, D_o=0.00635),  # 凝縮器 
    TLGeometry(N=200, L=2, D_i=0.00218, D_o=0.00318),      # 液管
    # TLGeometry(N=5, L=0.08, D_i=0.0018, D_o=0.0018),      # バイヨネット管
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037)         # CC
)
init_fluid_lookups(transient_condition.fluid, T_min=230.0, T_max=360.0, P_max=4.5e6)
var = initialize_transient_var(geo)
set_values_on_var!(var, transient_condition)
var_old = copy(var)
lb = copy(var)
ub = copy(var)
set_values_on_lb!(lb)
set_values_on_ub!(ub)

# T_wall = set_values_on_T_wall_from_steady(geo, transient_condition)
T_wall = initialize_T_wall(geo)
T_wall_old = copy(T_wall)

# scaling
scaled_var = scale_var(var, geo, transient_condition)
lb_scaled = scale_var(lb, geo, transient_condition)
ub_scaled = scale_var(ub, geo, transient_condition)

# var = set_values_on_var_from_steady(geo, 20)
# T_wall = set_values_on_T_wall_from_steady(geo, transient_condition)

save_T_wall_to_csv(T_wall, geo, "T_wall")
save_T_wall_to_csv(T_wall_old, geo, "T_wall_old")
save_var_to_csv(0.0, var, geo, "var")
save_var_to_csv(0.0, scaled_var, geo, "scaled_var")
save_var_to_csv(0.0, var_old, geo, "var_old")

# init_fluid_lookups(transient_condition.fluid, T_min=230.0, T_max=360.0, P_max=4.5e6)

##
set_values_on_var!(var, transient_condition)
set_values_on_T_wall!(T_wall, transient_condition)
# λ_plotter()

var_old = copy(var)
var_old.grv.T -= 0.5
T_wall_old = copy(T_wall)
dt = transient_condition.dt

T_wall_old .-= 0.5

@show T_wall.evap
m_evap = 0.00015
m_evap_old = 0.00015
save_T_wall_to_csv(T_wall_old, geo, "T_wall_old")
scaled_var = scale_var(var, geo, transient_condition)
params = (
        geo = geo,
        condition = transient_condition, # struct もそのまま格納できます
        m_evap = [0.00015],                    # 外側で初期値（0.0など）を入れておく
        T_wall = T_wall,
        var_old = var_old
    )
for t in 0.0:dt:0.01
    RK4_wall!(m_evap_old, var_old, T_wall, T_wall_old, geo, transient_condition)
    # f_transient_solve!(f, scaled_var, params)
    
    # m_evap, f.wick = fluid_wick(var, var_old, T_wall, geo, transient_condition)
    # f.grv.α, f.grv.m, f.grv.T, f.grv.P = fluid_grv(m_evap, var, var_old, T_wall, geo, transient_condition)
    # fluid_vl!(f.vl.α, f.vl.m, f.vl.T, f.vl.P, var, var_old, T_wall, geo, transient_condition)
    # fluid_cond!(f.cond.α, f.cond.m, f.cond.T, f.cond.P, var, var_old, T_wall, geo, transient_condition)
    # fluid_ll!(f.ll.α, f.ll.m, f.ll.T, f.ll.P, var, var_old, T_wall, geo, transient_condition)
    # f.cc.α, f.cc.m, f.cc.T, f.cc.P = fluid_cc(m_evap, var, var_old, T_wall, geo, transient_condition)
    # T_wall_old .= T_wall
end
save_T_wall_to_csv(T_wall, geo, "T_wall")
save_var_to_csv(m_evap, var, geo, "var")
save_var_to_csv(m_evap, scaled_var, geo, "scaled_var")
save_f_to_csv(f, geo, "residuals")
