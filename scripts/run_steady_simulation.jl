using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# 自作モジュールをロードするための設定
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHP_steady_Model

# 1. 幾何情報と計算条件を設定（ここでパラメータ変更が容易になる）
geo = LHPGeometry(
    WickGeometry(),
    GrvGeometry(),                                      # グルーブ
    TLGeometry(N=200, L=2, D_i=0.00435, D_o=0.00635),    # 蒸気管
    CondGeometry(N=200, L=1.0, D_i=0.00435, D_o=0.00635),  # 凝縮器 
    # TLGeometry(N=200, L=2, D_i=0.00218, D_o=0.00318),    # 液管
    TLGeometry(N=200, L=2, D_i=0.00435, D_o=0.00635),    # 液管（蒸気管と同じ太さ）
    TLGeometry(N=100, L=0.08, D_i=0.0018, D_o=0.0018),    # バイヨネット管
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037)       # CC
)
steady_condition = LHPSteadyCondition(Q_apply = 60.0) # 50Wで計算したい場合

init_propylene_lookups()

# 2. 初期化と計算
vars = initialize_steady_vars(geo)
solve_steady_lhp!(vars, geo, steady_condition) # 解析を実行するメイン関数（※現状のコードに合わせてラップする）

# 3. CSVに保存
save_vars_to_csv(vars, geo, steady_condition)