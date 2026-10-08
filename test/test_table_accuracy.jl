# =============================================================================
# test_table_accuracy.jl
#
# (ρ,e) テーブルの精度を、領域 (圧縮液・二相・過熱蒸気) ごとに評価する。
#
# 評価するもの:
#   - P, T の絶対誤差 (テーブルの値 − CoolProp の厳密値)
#   - 圧力の傾き ∂P/∂ρ|e の比 (テーブルの傾き / 厳密な傾き)
#     液体は圧縮性が極端に小さく、わずかな密度変化 (相対 1e-7) で圧力が
#     数百Pa動く。液管の摩擦圧損 (数Pa/m) のような小さな圧力差が流れを
#     駆動するので、絶対値より傾きの精度が重要になる。
#
# 傾きの比が 1 から大きく外れる領域は、そのテーブルでは液体の圧力応答
# (音速、圧縮性) を正しく表せない。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
include(joinpath(SRC, "smoothing.jl"))
include(joinpath(SRC, "types.jl"))
include(joinpath(SRC, "proptables.jl"))
include(joinpath(SRC, "propaccess.jl"))

fluid = "ammonia"
println("テーブル構築中 (初回は数分かかる場合があります。キャッシュされます)...")
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)

println("格子: e ∈ [", first(tables.e_range), ", ", last(tables.e_range), "] (Δe = ",
        step(tables.e_range), " J/kg),  δ=ρ−ρ_sl(e) ∈ [", first(tables.δ_range), ", ",
        last(tables.δ_range), "] (Δδ = ", step(tables.δ_range), " kg/m3)")

function state_TP(T, P)
    ρ = CoolProp.PropsSI("Dmass", "T", T, "P", P, fluid)
    e = CoolProp.PropsSI("Umass", "T", T, "P", P, fluid)
    return ρ, e
end

function state_TQ(T, x)
    ρ = CoolProp.PropsSI("Dmass", "T", T, "Q", x, fluid)
    e = CoolProp.PropsSI("Umass", "T", T, "Q", x, fluid)
    return ρ, e
end

# 傾き ∂P/∂ρ|e を、CoolProp (厳密) と テーブルで比べる
function slope_ratio(ρ, e)
    dρ = 1.0e-4 * ρ / 600.0   # 液なら約1e-4 kg/m3 (相対 ~2e-7)。圧力が 100 Pa 動く程度
    dP_exact = CoolProp.PropsSI("P", "Dmass", ρ + dρ, "Umass", e, fluid) -
               CoolProp.PropsSI("P", "Dmass", ρ - dρ, "Umass", e, fluid)
    dP_table = tables.P_of_ρe(ρ + dρ, e) - tables.P_of_ρe(ρ - dρ, e)
    return dP_table / dP_exact
end

function report(label, ρ, e, P_exact, T_exact)
    P_t = tables.P_of_ρe(ρ, e)
    T_t = tables.T_of_ρe(ρ, e)
    sr = try
        slope_ratio(ρ, e)
    catch
        NaN
    end
    println(rpad(label, 34), "  P誤差 = ", lpad(round(P_t - P_exact, digits = 1), 12), " Pa (相対 ",
            round(abs(P_t - P_exact) / P_exact, sigdigits = 2), "),  T誤差 = ",
            round(T_t - T_exact, digits = 3), " K,  ∂P/∂ρ の比 = ", round(sr, digits = 3))
end

println()
println("="^100)
println("圧縮液 (∂P/∂ρ の比が 1 に近ければ良好。大きく外れると液体の圧縮性が狂う)")
println("="^100)
for (T, P) in ((250.0, 1.5e6), (270.0, 2.0e6), (290.0, 1.5e6), (300.0, 2.0e6), (320.0, 3.0e6))
    ρ, e = state_TP(T, P)
    P_sat_liq = CoolProp.PropsSI("P", "T", T, "Q", 0.0, fluid)   # 注意: Psat は propaccess.jl の関数名
    P > P_sat_liq || continue
    report("液 T=$(T)K P=$(P/1e6)MPa", ρ, e, P, T)
end

println()
println("="^100)
println("二相 (P = P_sat(T) なので、P,T の誤差は飽和線に対する誤差)")
println("="^100)
for (T, x) in ((270.0, 0.05), (297.0, 0.02), (297.0, 0.3), (320.0, 0.5))
    ρ, e = state_TQ(T, x)
    P = CoolProp.PropsSI("P", "T", T, "Q", x, fluid)
    report("二相 T=$(T)K x=$(x)", ρ, e, P, T)
end

println()
println("="^100)
println("過熱蒸気")
println("="^100)
for (T, P) in ((300.0, 0.8e6), (330.0, 1.5e6), (350.0, 1.0e6))
    P_sat_vap = CoolProp.PropsSI("P", "T", T, "Q", 1.0, fluid)
    P < P_sat_vap || continue
    ρ, e = state_TP(T, P)
    report("蒸気 T=$(T)K P=$(P/1e6)MPa", ρ, e, P, T)
end