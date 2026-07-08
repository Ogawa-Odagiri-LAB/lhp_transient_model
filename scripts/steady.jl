using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# 自作モジュールをロードするための設定
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPSteadyModel

# 1. 幾何情報と計算条件を設定（ここでパラメータ変更が容易になる）
geo = LHPSteadyGeometry(
    WickGeometry(),
    GrvGeometry(N=7),                                      # グルーブ
    TLGeometry(N=200, L=2, D_i=0.00435, D_o=0.00635),    # 蒸気管
    CondGeometry(N=5000, L=1.0, D_i=0.00435, D_o=0.00635),  # 凝縮器 
    TLGeometry(N=200, L=2, D_i=0.00218, D_o=0.00318),    # 液管
    TLGeometry(N=10, L=0.08, D_i=0.0018, D_o=0.0018),    # バイヨネット管
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037)       # CC
)
steady_condition = LHPSteadyCondition(fluid = "ammonia", Q_apply = 20.0, T_cond = 10.0 + 273.15) # 60Wで計算したい場合

init_fluid_lookups(steady_condition.fluid, T_min=230.0, T_max=360.0, P_max=11e6)
# α = calc_void_rate(302.5, 0.005, "propylene")
# @show α

# 2. 初期化と計算
vars = initialize_steady_vars(geo)
solve_steady_lhp!(vars, geo, steady_condition) # 解析を実行するメイン関数（※現状のコードに合わせてラップする）

# 3. CSVに保存
save_steady_var_to_csv(vars, geo, steady_condition)