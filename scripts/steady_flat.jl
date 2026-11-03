using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel

# --- 結果の保存先 ---
result_name = "steady_flat_ammonia"
result_dir = isabspath(result_name) ? result_name : joinpath(@__DIR__, "..", "res", result_name)
# 定常状態の状態ベクトル (過渡計算の初期値に使う)
steady_file = joinpath(result_dir, "steady_state.jls")

# --- 条件 (実験条件に合わせる) ---
# 定常計算は、Q_apply_init の一定の熱負荷で行う (Q_apply, t_step などは使われない)。
# 熱負荷 0 W の定常状態を作る例: Q_apply_init = 0.0
condition = LHPTransientCondition(
    fluid = "ammonia",
    T_cond = 15.0 + 273.15,        # 凝縮器プレートの温度 [K]  ← 実験と同じ値
    T_amb  = 25.0 + 273.15,        # 周囲温度 [K]              ← 実験と同じ値
    T_init = 21.0 + 273.15,        # 擬似時間発展の出発点の基準温度 [K]。T_amb と T_cond の間 (実験の加熱前の壁温度) にすると速い。総質量にも効く
    Q_apply_init = 0.0,            # 定常計算の熱負荷 [W]
    Q_apply = 0.0,
    ΔT_startup = 0.0,
    t_init = 0.0, t_end = 1.0,     # 使われない (steady では t_chunk と t_max が決める)
)

# --- 幾何 (transient_flat.jl と同じにすること) ---
geo = FlatLHPGeometry(
    FlatWickGeometry(N=3, φ=0.65, ρ_wick=7929.0, Cp_wick=470.0),         # primary wick
    SecondaryWickGeometry(φ_s=0.50, ρ_swick=7929.0, Cp_swick=470.0),      # secondary wick
    FlatGrvGeometry(N=1),
    TLGeometry(N=20, L=2, D_i=0.00435, D_o=0.00635),
    CondGeometry(N=50, L=1.0, D_i=0.00435, D_o=0.00635),
    TLGeometry(N=20, L=2, D_i=0.00218, D_o=0.00318),
    FlatCCGeometry(N=1),
)

tables = init_fluid_lookups(condition.fluid, T_min=230.0, T_max=360.0, P_max=4.5e6)

# --- 定常計算 (擬似時間発展) ---
# t_chunk ごとに変化率を調べる。tol_dT: 壁温度の変化率 [K/s]、tol_dρ: コンポーネントごとの質量変化率 [1/s]
t_wall = @elapsed (sol, p, u_ss, info) = solve_steady_lhp!(geo, condition, tables;
    t_chunk = 50.0, t_max = 20000.0, tol_dT = 1.0e-3, tol_dρ = 1.0e-5)
println("計算時間: ", round(t_wall, digits = 1), " s (擬似時間 ", info.t, " s)")

# --- 保存 ---
save_steady_state(steady_file, u_ss, p, condition)
mkpath(result_dir)
using DataFrames, CSV
CSV.write(joinpath(result_dir, "steady_history.csv"), DataFrame(info.history))

# 最終状態の要約
using ComponentArrays
c = ComponentArray(copy(u_ss), p.axes)
@show c.T_wall.evap[1] c.T_wall.cc[1] c.swick.T[1] c.cc.P[1] c.cc.T[1] c.cc.x[1]