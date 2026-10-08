# =============================================================================
# test_residuals_wick.jl
#
# residuals_wick.jl の動作確認。CoolPropのテーブルを実際に構築して、
#   1. 全ノード等温・ṁ=0 で全残差がゼロ
#   2. 対数温度分布 (円筒殻の定常伝導の厳密解) で、各セルの熱収支がゼロ
#   3. 蒸発器ケース・ウィック・グルーブ・CCの合計エネルギーが保存する
# を確認する。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
include(joinpath(SRC, "physics.jl"))      # k_sus304 (ウィックの k_eff の計算に使う)
include(joinpath(SRC, "smoothing.jl"))
include(joinpath(SRC, "types.jl"))
include(joinpath(SRC, "proptables.jl"))
include(joinpath(SRC, "propaccess.jl"))
include(joinpath(SRC, "residuals_wick.jl"))

fluid = "ammonia"
println("テーブル構築中 (初回は数分かかる場合があります。キャッシュされます)...")
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)

geo = WickGeometry(N=5, L=0.035, D_i=0.003, D_o=0.006,
                   φ=0.65, ρ_wick=8900.0, Cp_wick=444.0)
N = geo.N

P_grv, P_w, T_cc = 1.0e6, 1.05e6, 297.0
# 固相の実効熱伝導率 k_eff: SUS304 と液の熱伝導率から、マクスウェルの複合則 (式53) で計算
wp = wick_params(geo, k_sus304(T_cc), tables.k_l_sat(T_cc))
println("k_s = ", k_sus304(T_cc), ", k_l = ", tables.k_l_sat(T_cc), " W/mK → k_eff = ", wp.k_eff,
        " W/mK (φ=", geo.φ, ")")
T_sat_grv = Tsat(P_grv, tables)
ρ_l_ref = tables.ρ_l_sat(T_cc)

println("T_sat(P_grv) = ", T_sat_grv, " K")

function make_wick(T)
    E = [wp.V[i] * ((1 - wp.por) * wp.ρ_w * wp.C_w * T[i] +
                    wp.por * ρ_l_ref * liquid_props(T[i], P_w, tables).e) for i in 1:N]
    return ComponentArray(E = E, T = collect(T))
end

println("="^64)
println("テスト1: 全ノード等温 (=T_sat)・T_ec=T_sat (ṁ=0) → 全残差 ≈ 0")
println("="^64)
w = make_wick(fill(T_sat_grv, N))
du = ComponentArray(E = zeros(N), T = zeros(N))
resid = ComponentArray(E = zeros(N), T = zeros(N))
iface = wick_interface(w.T[N], P_grv, T_sat_grv, 1.0e4, 1.0e-3, wp.G[N+1], P_w, tables)
Q_leak = residual_wick!(resid, du, w, wp, iface, T_sat_grv, P_w, tables; ρ_l_ref = ρ_l_ref)
println("  ṁ = ", iface.ṁ, " (0 のはず)")
println("  max|resid.E| = ", maximum(abs, resid.E), ", max|resid.T| = ", maximum(abs, resid.T))
println("  Q_leak = ", Q_leak, " W (0 のはず)")

println()
println("="^64)
println("テスト2: 対数温度分布 (定常伝導の厳密解) → 各セルの熱収支 ≈ 0 (ṁ=0 とする)")
println("="^64)
Tin, Tout = T_cc, T_sat_grv + 3.0
Dn = wp.Dn
Tn = [Tin + (Tout - Tin) * log(Dn[k] / Dn[1]) / log(Dn[end] / Dn[1]) for k in 1:N+2]
w = make_wick(Tn[2:N+1])
iface0 = (ṁ = 0.0, T_out = Tn[end], h_out = liquid_props(Tn[end], P_w, tables).h,
          Q_ec = 0.0, Q_N = 0.0, h_N = 0.0, λ = 1.0, Q_extra = 0.0, g_s = 1.0)
Q_leak = residual_wick!(resid, du, w, wp, iface0, Tin, P_w, tables; ρ_l_ref = ρ_l_ref)
println("  max|resid.E| = ", maximum(abs, resid.E), " W  (≈0 のはず。-resid.E が各セルの正味受熱)")
println("  内表面からの伝導熱 Q_leak = ", Q_leak, " W")

println()
println("="^64)
println("テスト3: 4者 (蒸発器ケース・ウィック・グルーブ・CC) の合計エネルギー保存")
println("="^64)
for (title, Tcells, T_ec) in (
    ("セルNが界面に近い", range(T_cc + 0.5, T_sat_grv - 0.05, length = N), T_sat_grv + 8.0),
    ("過渡 (セルNが界面より10K低い)", range(T_cc - 1.0, T_sat_grv - 10.0, length = N), T_sat_grv + 8.0),
)
    for sens in (false, true)
        w = make_wick(collect(Tcells))
        iface = wick_interface(w.T[N], P_grv, T_ec, 1.0e4, 1.0e-3, wp.G[N+1], P_w, tables;
                               λ_includes_sensible = sens)
        Q_leak = residual_wick!(resid, du, w, wp, iface, T_cc, P_w, tables; ρ_l_ref = ρ_l_ref)
        dU_wick = -sum(resid.E)                           # du=0 での resid.E = -RHS なので、蓄熱変化率 = -sum(resid.E)
        dU_evap = -iface.Q_ec
        S_m, S_E = wick_groove_source(iface)
        ṁ = iface.ṁ
        wm = smooth_step(ṁ, 1e-8)
        h0 = liquid_props(T_cc, P_w, tables).h
        h1 = liquid_props(w.T[1], P_w, tables).h
        dU_cc = -ṁ * (wm * h0 + (1 - wm) * h1) + Q_leak   # CC: 面1の風上エンタルピーで流出
        total = dU_wick + dU_evap + S_E + dU_cc
        println("  ", title, ", λ_includes_sensible=", sens, ": 合計 = ", total, " W (≈0 のはず)")
    end
end

println()
println("="^64)
println("テスト4: 起動判定 (ΔT_startup>0) のゲートを入れても、エネルギーが保存する")
println("="^64)
Tcells = collect(range(T_cc - 1.0, T_sat_grv - 3.0, length = N))
for ΔTs in (0.0, 2.0, 10.0)
    w = make_wick(Tcells)
    iface = wick_interface(w.T[N], P_grv, T_sat_grv + 8.0, 1.0e4, 1.0e-3, wp.G[N+1], P_w, tables;
                           T_cc = T_cc, ΔT_startup = ΔTs)
    Q_leak = residual_wick!(resid, du, w, wp, iface, T_cc, P_w, tables; ρ_l_ref = ρ_l_ref)
    dU_wick = -sum(resid.E)
    dU_evap = -iface.Q_ec
    S_m, S_E = wick_groove_source(iface)
    ṁ = iface.ṁ
    wm = smooth_step(ṁ, 1e-8)
    h0 = liquid_props(T_cc, P_w, tables).h
    h1 = liquid_props(w.T[1], P_w, tables).h
    dU_cc = -ṁ * (wm * h0 + (1 - wm) * h1) + Q_leak
    total = dU_wick + dU_evap + S_E + dU_cc
    println("  ΔT_startup=", ΔTs, ": ゲート g_s=", round(iface.g_s, digits = 4), ", ṁ=", ṁ,
            ", Q_extra=", iface.Q_extra, " W, 合計 = ", total, " W (≈0 のはず)")
end
println("  (ΔT_startup=0 なら g_s=1, Q_extra=0 で従来と同じ。g_s<1 では ṁ が減り、",
        "蒸発に使われなかった熱が最外セルに渡される)")