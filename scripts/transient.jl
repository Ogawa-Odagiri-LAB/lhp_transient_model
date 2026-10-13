using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# 自作モジュールをロードするための設定
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel

# --- 結果の保存先 ---
# 名前だけ指定すると ../res/<名前>/ に保存される。絶対パスを指定すると、そこに保存される。
result_name = "transient_ammonia"
result_dir = isabspath(result_name) ? result_name : joinpath(@__DIR__, "..", "res", result_name)

# --- 運転条件 ---
# Q_apply_init から Q_apply へ、t_step から τ_step かけて、滑らかにステップアップする
# (Q_apply_init = Q_apply ならステップなし)。時間刻みは、IDA が許容誤差で制御する。
condition = LHPTransientCondition(
    fluid = "ammonia",
    T_cond = 15.0 + 273.15,        # 凝縮器プレートの温度 [K]
    T_amb  = 25.0 + 273.15,        # 周囲温度 [K]
    T_init = 25.0 + 273.15,        # 初期の基準温度 [K]
    Q_apply_init = 5.0,           # ステップ前の印加熱量 [W]
    Q_apply = 10.0,                # ステップ後の印加熱量 [W]
    t_step = 50.0, τ_step = 1.0,   # t = 50 s から 1 s かけて 20 W → 40 W
    ΔT_startup = 0.0,              # 起動判定 (0 なら判定なし)
    t_init = 0.0, t_end = 200.0,
    saveat = 1.0,                  # 保存間隔 [s]。刻みへの影響はない (IDA の内部の刻みとは独立)
)

# --- 幾何 ---
geo = LHPGeometry(
    # ウィック: 固相の実効熱伝導率 k_eff は、SUS304 と液の熱伝導率、空隙率 φ から、
    # マクスウェルの複合則 (式53) で計算される (入力しない)
    WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0),
    GrvGeometry(N=1),                                     # グルーブ
    TLGeometry(N=20, L=2, D_i=0.00435, D_o=0.00635),      # 蒸気管
    CondGeometry(N=50, L=1.0, D_i=0.00435, D_o=0.00635),  # 凝縮器
    TLGeometry(N=20, L=2, D_i=0.00218, D_o=0.00318),      # 液管
    CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037),        # CC
)

# --- 作動流体の物性テーブル (初回のみ数分。以降はキャッシュから読み込む) ---
tables = init_fluid_lookups(condition.fluid, T_min=230.0, T_max=360.0, P_max=4.5e6)

# --- 計算 ---
# 既定: 線形ソルバー KLU (疎行列)、BDF の最大次数 2、reltol=1e-4、atol_frac_u=1e-7
# (solverfuncs.jl 冒頭に、これらの設定の根拠 (実測) を記載)
t_wall = @elapsed (sol, p, u0) = solve_transient_lhp!(geo, condition, tables)
println("計算時間: ", round(t_wall, digits = 1), " s (シミュレーション ", condition.t_end - condition.t_init, " s)")

# --- エネルギー収支 (蓄積の変化 = 印加熱 − 放熱。残差が入力の数 % を超えるなら要確認) ---
energy_balance(sol, p)

# --- 結果の保存 (CSV) ---
save_results(sol, p, tables, result_dir)

# --- プロット ---
# 保存した CSV を、plot_results.jl (スクリプト) で、計算をやり直さずに、いつでも眺められる。

# 例: 最終時刻の状態を、コンポーネントごとに取り出す
c = components_at(sol, p, length(sol.t))
@show c.T_wall.evap[1] c.cc.P[1] c.cc.T[1]