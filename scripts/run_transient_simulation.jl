using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# 自作モジュールをロードするための設定
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHP_transient_Model

# 1. 幾何情報と計算条件を設定（ここでパラメータ変更が容易になる）
geo = LHPGeometry(
    WickGeometry(),                                     # ウィック
    GrvGeometry(),                                      # グルーブ
    TLGeometry(N=200, L=2, D_i=0.00435, D_o=0.00635),    # 蒸気管
    CondGeometry(N=200, L=1.0, D_i=0.00435, D_o=0.00635),  # 凝縮器 
    TLGeometry(N=200, L=2, D_i=0.00218, D_o=0.00318),    # 液管
    TLGeometry(N=100, L=0.08, D_i=0.0018, D_o=0.0018),    # バイヨネット管
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037)       # CC
)

# 条件設定
transient_condition = LHPTransientCondition(Q_apply_init = 50.0, Q_apply = 60.0, dt = 0.1, t_init = 0.0, t_end = 1.0)

init_propylene_lookups()

# 
solve_transient_lhp!(geo, transient_condition)
