using Pkg; Pkg.activate(".")
using CSV, DataFrames

# =============================================================================
# wick_dp_diag.jl … 保存済みの計算結果から、起動時にウィックの両側にかかる圧力差を、毛管圧の上限と比べる (計算のやり直しは不要)。
#   使い方: julia wick_dp_diag.jl      (この .jl は diag/ に置く。設定は下の 3 行をこのファイルで書き換える)
#   読むもの: res/<result_name>/summary.csv  (列: t, P_grv, P_cc, dP_ne, mdot, T_ec)
#   目的: ΔP_ne のスパイクと「ウィック界面の乾き・液の戻り」が関係するかを見る診断。物理は変えない。
#     ・ウィックの蒸気側 (グルーブ) の圧力 P_grv と、液側 (CC) の圧力 P_cc の差 ΔP_wick = P_grv − P_cc が、
#       メニスカスの保持できる最大の毛管圧 2σ/r を超えていないか。超えるなら、実機ではウィックに蒸気が入る (乾き・逆流) 可能性がある。
#       モデルのウィックは常に液で満たされ、毛管圧の上限の判定を持たないので、超えても何も起きない。
#     ・dP_ne は CC の蒸気側の過剰圧力 (液の流入による圧縮。モデルが再現しようとしているスパイク)。
# =============================================================================
result_name = "transient_flat_ammonia_v5"   # 見る結果のフォルダ名 (res/ の下)
t_heat = 10.0                                # モデルの加熱開始時刻 [s] (stage = :long)
σ_ref = 0.0209                               # アンモニアの表面張力 [N/m] (約 22.6 °C。概算)
r_list = (5.1e-6, 3.0e-6, 1.5e-6)            # 毛管圧の上限を計算する細孔半径 [m] (5.1 µm = 実測)

result_dir = joinpath(@__DIR__, "..", "res", result_name)
summ = CSV.read(joinpath(result_dir, "summary.csv"), DataFrame)
srow(t) = summ[argmin(abs.(summ.t .- t)), :]

println("毛管圧の上限 2σ/r [kPa]:  ", join(["r = $(round(r * 1e6, digits = 1)) µm → $(round(2 * σ_ref / r / 1e3, digits = 1))" for r in r_list], ",  "))
println("\n結果: ", result_name, "  (加熱開始からの時間 τ [s])")
println(rpad("  τ", 6), lpad("mdot[mg/s]", 12), lpad("P_grv[kPa]", 12), lpad("P_cc[kPa]", 12), lpad("P_grv-P_cc", 12), lpad("dP_ne[kPa]", 12))
for τ in (0, 10, 12, 14, 16, 18, 20, 21, 22, 23, 24, 26, 28, 30, 34, 38, 42, 44, 46, 48, 50, 55, 60, 70, 90)
    row = srow(t_heat + τ)
    println(lpad("+" * string(τ), 5), " ",
            lpad(round(1e6 * row.mdot, digits = 2), 11), " ",
            lpad(round(row.P_grv / 1e3, digits = 2), 11), " ",
            lpad(round(row.P_cc / 1e3, digits = 2), 11), " ",
            lpad(round((row.P_grv - row.P_cc) / 1e3, digits = 2), 11), " ",
            lpad(round(row.dP_ne / 1e3, digits = 2), 11))
end
dP = (summ.P_grv .- summ.P_cc) ./ 1e3
k = argmax(dP)
println("\nP_grv − P_cc の最大: ", round(dP[k], digits = 2), " kPa  (τ = ", round(summ.t[k] - t_heat, digits = 1), " s)")
println("毛管圧の上限を超える時間 (r = 5.1 µm: ", round(2 * σ_ref / 5.1e-6 / 1e3, digits = 1), " kPa): ",
        round(sum((dP .> 2 * σ_ref / 5.1e-6 / 1e3) .* [i == 1 ? 0.0 : summ.t[i] - summ.t[i-1] for i in eachindex(summ.t)]), digits = 1), " s")