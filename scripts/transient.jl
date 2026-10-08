using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# 自作モジュールをロードするための設定
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel

transient_condition = LHPTransientCondition(fluid = "ammonia", T_cond = 10.0 + 273.15, Q_apply_init = 20.0, Q_apply = 30.0, ΔT_startup = 0.0, dt = 0.1, t_init = 0.0, t_end = 1.0) # dt=0.1だと5.0e0, dt=0.01だと1.0e-3までいく

geo = LHPGeometry(
    WickGeometry(N=6),
    GrvGeometry(N=1),                                     # グルーブ
    TLGeometry(N=100, L=2, D_i=0.00435, D_o=0.00635),      # 蒸気管
    CondGeometry(N=300, L=1.0, D_i=0.00435, D_o=0.00635),  # 凝縮器 
    TLGeometry(N=100, L=2, D_i=0.00218, D_o=0.00318),      # 液管
    # TLGeometry(N=5, L=0.08, D_i=0.0018, D_o=0.0018),      # バイヨネット管
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037)         # CC
)

init_fluid_lookups(transient_condition.fluid, T_min=230.0, T_max=360.0, P_max=4.5e6)

solve_transient_lhp!(geo, transient_condition)
