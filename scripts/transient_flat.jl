using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# 平板型蒸発器 (金属3Dプリンター製、secondary wick 付き) の過渡計算。
# transient.jl (円筒型) と同じ流れで、ジオメトリだけが違う。
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel

result_name = "transient_flat_ammonia"

condition = LHPTransientCondition(
    fluid = "ammonia",
    T_cond = 15.0 + 273.15,
    T_amb  = 25.0 + 273.15,
    T_init = 25.0 + 273.15,
    Q_apply_init = 5.0,
    Q_apply = 10.0,
    t_step = 50.0, τ_step = 1.0,
    ΔT_startup = 0.0,
    t_init = 0.0, t_end = 200.0,
    saveat = 1.0,
)

# --- 幾何 (断面図の寸法) ---
# primary wick と secondary wick は、空隙率・材料物性が別。
# φ, ρ_wick, Cp_wick, φ_s, ρ_swick, Cp_swick の値は仮置き (ご自身の値に置き換えてください)。
geo = FlatLHPGeometry(
    FlatWickGeometry(N=3, φ=0.30, ρ_wick=8900.0, Cp_wick=444.0),         # primary wick (高さ方向 N 層)
    SecondaryWickGeometry(φ_s=0.70, ρ_swick=8900.0, Cp_swick=444.0),      # secondary wick (1ノード、常に液で満たされる)
    FlatGrvGeometry(N=1),                                                 # グルーブ (平板)
    TLGeometry(N=20, L=1.17, D_i=0.00176, D_o=0.00318),                      # 蒸気管
    CondGeometry(N=50, L=0.56, D_i=0.00176, D_o=0.00318),                  # 凝縮器
    TLGeometry(N=20, L=1.15, D_i=0.00099, D_o=0.00159),                      # 液管
    FlatCCGeometry(N=1),                                                  # CC (半円筒 D_i=27 mm, L=25 mm)
)

tables = init_fluid_lookups(condition.fluid, T_min=230.0, T_max=360.0, P_max=4.5e6)

t_wall = @elapsed (sol, p, u0) = solve_transient_lhp!(geo, condition, tables)
println("計算時間: ", round(t_wall, digits = 1), " s (シミュレーション ", condition.t_end - condition.t_init, " s)")

energy_balance(sol, p)
save_results(sol, p, tables, joinpath(@__DIR__, "..", "res", result_name))

c = components_at(sol, p, length(sol.t))
@show c.T_wall.evap[1] c.T_wall.cc[1] c.swick.T[1] c.cc.P[1] c.cc.T[1]