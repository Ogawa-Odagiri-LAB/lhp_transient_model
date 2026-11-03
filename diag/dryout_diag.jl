using Pkg; Pkg.activate(".")
using CSV, DataFrames

# =============================================================================
# dryout_diag.jl … 保存済みの計算結果から、ウィック界面の「干上がり」(過熱した層の細孔が蒸気で満たされ、直後に液が戻る)が、
# 温度・圧力のスパイクにどの程度効きうるかを、上限で見積もる (計算のやり直しは不要。物理は変えない。逆流は含めない)。
#   使い方: julia dryout_diag.jl      (この .jl は diag/ に置く。設定は下の数行をこのファイルで書き換える)
#   読むもの: res/<result_name>/summary.csv (T_out, T_ec, mdot, dP_ne, Q_hb_ec), wick.csv (T, cell: 1..N = primary wick、N が界面側)
#
# 見積りの考え方:
#   (1) 過熱した層の顕熱 ρc_eff·θ で、その層の細孔 (体積分率 φ) のうち、飽和蒸気で満たせる割合
#         f = ρc_eff·θ / (φ ρ_v h_fg)    (θ = T_cell − T_sat(P_grv)。1 を超えたら 1)
#       液と蒸気の密度比 (約 83) が大きいので、θ ≈ 0.7 K で f = 1 (細孔が全部蒸気) になる。
#   (2) 干上がった細孔体積 V_dry = Σ f_i V_pore,i を、毛管力でウィックに液が戻る時間 (ダルシー則) で埋める:
#         t_refill = V_dry μ_l H / (K A ΔP_drive)
#       透過率 K と駆動圧 ΔP_drive は【概算】。K の実測値があれば差し替える (t_refill は 1/K に比例)。
#   (3) 干上がりの影響の上限:
#         温度: 干上がりの間は蒸発で熱が除かれないので、受熱面 P が Q_hb_ec·t_refill/C_P だけ余分に上がる
#               (センサの側面 S には、P–S の時定数 C_S/G_ps で遅れて、1 − exp(−t_refill/τ_PS) 倍だけ伝わる)。
#         圧力: 蒸気の発生が最大 t_refill だけ遅れるとして、ΔP_ne の立ち上がりの最大の傾き × t_refill。
# =============================================================================
result_name = "transient_flat_ammonia_v5"   # 見る結果のフォルダ名 (res/ の下)
t_heat = 10.0                                # モデルの加熱開始時刻 [s] (stage = :long)

# 物性・幾何 (概算。約 22.6 °C のアンモニア。ウィックの値はモデルと同じ)
φ, ρ_w, c_w = 0.30, 8900.0, 444.0
ρ_v, h_fg, ρ_l, c_l, μ_l = 7.3, 1.17e6, 607.0, 4.7e3, 1.4e-4
A_w, H_w, N_w = 0.030 * 0.030, 0.0032, 3     # primary wick の面積 [m²]・高さ [m]・層数
K_list = (1e-14, 3e-14, 1e-13, 3e-13, 1e-12) # 透過率 K [m²]。細孔半径 5.1 µm・空隙率 0.30 からの概算は 2e-13〜5e-13 (毛細管束 φ r²/16 ≈ 5e-13、コゼニー・カーマン ≈ 2e-13)。不確かさは 10 倍程度。下限は実験の定常の VL1−RES (≈0.03 K ≈ 0.9 kPa) から K ≳ 4e-15。実測があれば差し替える
ΔP_drive = 4.0e3                             # 再濡れの駆動圧 [Pa] (毛管圧 2σ/r ≈ 8 kPa の半分程度)
C_P, C_S, G_ps = 6.28, 3.48, 0.93            # EC の 2 節点 (transient_flat.jl と同じ)
τ_list = (16, 18, 19, 20, 21, 22, 23, 24, 26)  # 加熱開始からの時間 [s] (蒸発の開始の前後)

result_dir = joinpath(@__DIR__, "..", "res", result_name)
summ = CSV.read(joinpath(result_dir, "summary.csv"), DataFrame)
wick = CSV.read(joinpath(result_dir, "wick.csv"), DataFrame)
srow(t) = summ[argmin(abs.(summ.t .- t)), :]
tw = sort(unique(wick.t))
function wick_T(t)
    tt = tw[argmin(abs.(tw .- t))]
    d = sort(wick[(wick.t .== tt) .& (wick.cell .> 0), :], :cell)
    return Vector(d.T)
end

ρc_eff = (1 - φ) * ρ_w * c_w + φ * ρ_l * c_l
V_pore = φ * A_w * H_w / N_w                 # 1 層の細孔体積 [m³]
τ_PS = C_S / G_ps

println("結果: ", result_name, "   ρc_eff = ", round(ρc_eff / 1e6, digits = 2), "e6 J/(m³K),  細孔体積 (全層) = ", round(N_w * V_pore * 1e9), " mm³ (液 ",
        round(N_w * V_pore * ρ_l * 1e6), " mg)")
println("\n加熱開始からの時間 τ [s] ごとの、過熱した層の細孔を満たせる割合と、干上がり体積")
println(rpad("  τ", 6), lpad("P-Tsat[K]", 11), lpad("θ(界面側→内側)[K]", 24), lpad("f(界面側→内側)", 22), lpad("V_dry[mm³]", 12), lpad("蒸気[mg]", 10), lpad("エネルギー[J]", 14))
V_dry_max = 0.0
for τ in τ_list
    row = srow(t_heat + τ)
    Ts = row.T_out
    θ = reverse(wick_T(t_heat + τ)) .- Ts        # 界面側 (最外層) を先頭に
    f = clamp.(ρc_eff .* max.(θ, 0.0) ./ (φ * ρ_v * h_fg), 0.0, 1.0)
    V_dry = sum(f) * V_pore
    global V_dry_max = max(V_dry_max, V_dry)
    println(lpad("+" * string(τ), 5), lpad(round(row.T_ec - Ts, digits = 2), 11),
            lpad(join(round.(θ, digits = 2), " "), 24), lpad(join(round.(f, digits = 2), " "), 22),
            lpad(round(V_dry * 1e9), 12), lpad(round(V_dry * ρ_v * 1e6, digits = 1), 10), lpad(round(V_dry * ρ_v * h_fg, digits = 1), 14))
end

# 実際の蒸発開始の直後 (n が立ち上がるころ) の ∫ṁ dt (モデルの噴出量) と比べる
let ts = summ.t, m = summ.mdot
    sel = (ts .>= t_heat + 20) .& (ts .<= t_heat + 26)
    dt = diff(ts)
    mass = sum(m[2:end][sel[2:end]] .* dt[sel[2:end]]) * 1e6
    println("\nモデルが +20〜26 s に蒸発させた質量 ∫ṁ dt = ", round(mass, digits = 1), " mg  (干上がりの蒸気は最大 ", round(V_dry_max * ρ_v * 1e6, digits = 1), " mg)")
end

# 上限の見積り
dP = summ.dP_ne ./ 1e3
slope = maximum(diff(dP) ./ diff(summ.t))
amp = maximum(dP)
Q_hb = maximum(summ.Q_hb_ec[(summ.t .>= t_heat + 15) .& (summ.t .<= t_heat + 26)])
println("\nΔP_ne の最大 = ", round(amp, digits = 1), " kPa,  立ち上がりの最大の傾き = ", round(slope, digits = 2), " kPa/s,  Q_hb_ec (+15〜26 s の最大) = ", round(Q_hb, digits = 2), " W")
println("干上がりの影響の上限 (V_dry = ", round(V_dry_max * 1e9), " mm³、駆動圧 ", ΔP_drive / 1e3, " kPa):")
println(rpad("  K [m²]", 12), lpad("t_refill[s]", 13), lpad("ΔT_P 上限[K]", 14), lpad("ΔT_S 上限[K]", 14), lpad("ΔP ずれ[kPa]", 14), lpad("(ΔP_ne 最大に対し)", 20))
println("  ※ ΔP ずれ = 立ち上がりが t_refill だけ遅れたときの、最も急な所での ΔP のずれ。時間のずれであって、ピークの高さ・減衰は変わらない")
for K in K_list
    t_ref = V_dry_max * μ_l * H_w / (K * A_w * ΔP_drive)
    ΔT_P = Q_hb * t_ref / C_P
    ΔT_S = ΔT_P * (1 - exp(-t_ref / τ_PS))
    ΔP = slope * t_ref
    println(rpad("  " * string(K), 12), lpad(round(t_ref, digits = 2), 13), lpad(round(ΔT_P, digits = 3), 14), lpad(round(ΔT_S, digits = 3), 14),
            lpad(round(ΔP, digits = 2), 14), lpad(string(round(100 * ΔP / amp, digits = 1), " %"), 20))
end