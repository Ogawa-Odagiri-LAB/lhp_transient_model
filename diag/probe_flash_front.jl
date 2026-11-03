using Pkg; Pkg.activate(".")

# VL の前線セル (セル 12: 飽和液 ρ≈607 kg/m³, T≈22.4 °C, x≈0) まわりで、物性テーブルの P(ρ,e) を調べる。
# 診断のみ。モデルもテーブルも変更しない。
# 実行: julia --project diag/probe_flash_front.jl
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel
using Printf

tables = init_fluid_lookups("ammonia", T_min=230.0, T_max=360.0, P_max=4.5e6, Ne=1600, Δδ=0.5)

T0 = 22.4 + 273.15
ρ_l = tables.ρ_l_sat(T0)
e_sl = tables.h_l_sat(T0) - tables.Psat_T(T0) / ρ_l     # 飽和液の内部エネルギー
println("飽和液 (T = 22.4 °C): ρ_l = ", round(ρ_l, digits=3), " kg/m³,  e_sl = ", round(e_sl, digits=1), " J/kg,  Psat = ", round(tables.Psat_T(T0)/1e3, digits=2), " kPa")
println("(参考) 液の dP/dρ ≈ 1.4e6 Pa/(kg/m³) のはず")

f(ρ, e) = LHPTransientModel.flash_ρe(ρ, e, tables)

# ---- 1) ρ を飽和液の値に固定して e を動かす (飽和線を e 方向に横切る) ----
println("\n[1] ρ = ρ_l(22.4°C) 固定、e = e_sl + de。  P [kPa], T [°C], x, x_ps, dP/de [Pa/(J/kg)] (差分), dP/dρ [Pa/(kg/m³)] (差分)")
println("      de[J/kg]      P[kPa]     T[°C]        x      x_ps     dP/de      dP/dρ")
for de in -1500:100:1500
    e = e_sl + de
    r = f(ρ_l, e)
    dρ = 0.01; ee = 5.0
    dPde = (f(ρ_l, e + ee).P - f(ρ_l, e - ee).P) / (2ee)
    dPdρ = (f(ρ_l + dρ, e).P - f(ρ_l - dρ, e).P) / (2dρ)
    xps = LHPTransientModel.pseudo_quality_P(e, ρ_l, r.P, tables)
    @printf("  %10.0f  %10.3f  %8.3f  %8.5f  %8.5f  %10.3g  %10.3g\n", de, r.P/1e3, r.T-273.15, r.x, xps, dPde, dPdρ)
end

# ---- 2) e を飽和液の値に固定して ρ を動かす (密度が飽和液より小さくなる = 二相に入る) ----
println("\n[2] e = e_sl 固定、ρ = ρ_l + dρ。")
println("      dρ[kg/m³]     P[kPa]     T[°C]        x      x_ps     dP/dρ")
for dρ in (-5.0, -2.0, -1.0, -0.5, -0.2, -0.1, -0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
    ρ = ρ_l + dρ
    r = f(ρ, e_sl)
    h = 0.005
    dPdρ = (f(ρ + h, e_sl).P - f(ρ - h, e_sl).P) / (2h)
    xps = LHPTransientModel.pseudo_quality_P(e_sl, ρ, r.P, tables)
    @printf("  %10.3f  %10.3f  %8.3f  %8.5f  %8.5f  %10.3g\n", dρ, r.P/1e3, r.T-273.15, r.x, xps, dPdρ)
end

# ---- 3) 圧力の飽和圧力からのずれ: ΔP = P(ρ,e) − Psat(T(ρ,e)) ----
println("\n[3] 実際の前線セルの状態付近: (ρ, e) = 診断値 (ρ=606.85, T=22.397 °C, 前線の液) 周辺で、P − Psat(T) [Pa]")
for dρ in (-0.3, -0.1, 0.0, 0.1, 0.3), de in (-200.0, -50.0, 0.0, 50.0, 200.0, 500.0)
    ρ = 606.85 + dρ
    e = tables.h_l_sat(22.397 + 273.15) - tables.Psat_T(22.397 + 273.15) / tables.ρ_l_sat(22.397 + 273.15) + de
    r = f(ρ, e)
    @printf("  ρ=%8.3f  de=%6.0f  P=%9.3f kPa  T=%7.3f °C  x=%8.5f  P−Psat(T)=%9.1f Pa\n", ρ, de, r.P/1e3, r.T-273.15, r.x, r.P - tables.Psat_T(r.T))
end