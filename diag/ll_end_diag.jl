using Pkg; Pkg.activate(".")
using CSV, DataFrames

# =============================================================================
# ll_end_diag.jl … 保存済みの計算結果から、起動後 0〜70 s の液管の末端 3 セル (18〜20) と CC の圧力・温度を時間ごとに出す。計算のやり直しは不要。
#   実験の LL3 (CC から 7 cm) の +0.25 K の山 (τ≈19〜30 s、ピーク 23〜25 s) を、モデルが出すかを見る診断。
#   モデルの LL3 = 液管の最後のセル (セル 20) の壁温で代用。x > 0 は気相が残っている (二相)、T − Tsat は二相ならほぼ 0。
#   使い方: julia ll_end_diag.jl      (diag/ に置く。result_name を見たい結果に合わせる)
# =============================================================================
result_name = "transient_flat_ammonia_2h"
result_dir = joinpath(@__DIR__, "..", "res", result_name)
t_heat, t_exp_onset = 10.0, 303.8
K = 273.15

summ = CSV.read(joinpath(result_dir, "summary.csv"), DataFrame)
prof = CSV.read(joinpath(result_dir, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :rho, :u, :P, :T, :x, :T_wall])
ll = prof[(prof.comp .== "ll") .& (prof.cell .>= 18), :]
exp_file = joinpath(@__DIR__, "..", "data", "exp_5W_12min.csv")
ex = CSV.read(exp_file, DataFrame)
lin1(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
srow(t) = summ[argmin(abs.(summ.t .- t)), :]
cell(t, i) = (d = ll[(ll.t .== t) .& (ll.cell .== i), :]; nrow(d) == 0 ? nothing : d[1, :])
tl = sort(unique(ll.t))
nearest(t) = tl[argmin(abs.(tl .- t))]
f(x, n = 2) = lpad(string(round(x, digits = n)), 8)

P0 = srow(t_heat - 0.5).P_cc
println("液管の末端のセルと CC (加熱開始からの時間 τ)。圧力は τ=0 からの変化 [kPa]、温度は [°C]、x は乾き度 (質量基準)")
println(rpad("τ", 5), "  ΔP_cc   T_cc   |  セル20: ΔP     T   T−Tsat?     x     T_wall |  セル19: x    T_wall  |  セル18: x   T_wall | 実験LL3   モデルLL3(壁20)−初期  実験−初期")
T0w = cell(nearest(t_heat - 0.5), 20).T_wall - K
E0 = lin1(ex.t, ex.LL3, t_exp_onset - 0.5)
for τ in 0:2:70
    t = nearest(t_heat + τ)
    r = srow(t); c20 = cell(t, 20); c19 = cell(t, 19); c18 = cell(t, 18)
    (c20 === nothing || c19 === nothing || c18 === nothing) && continue
    e = lin1(ex.t, ex.LL3, t_exp_onset + τ)
    println(rpad(τ, 5), f((r.P_cc - P0) / 1e3, 1), f(r.T_cc - K), "  |", f((c20.P - P0) / 1e3, 1), f(c20.T - K), "  ", rpad(round(c20.x, sigdigits = 3), 9),
            f(c20.T_wall - K), " |", rpad(round(c19.x, sigdigits = 3), 9), f(c19.T_wall - K), " |", rpad(round(c18.x, sigdigits = 3), 9), f(c18.T_wall - K),
            " |", f(e), f(c20.T_wall - K - T0w, 3), f(e - E0, 3))
end
println("\n注: 実験の LL3 は τ≈18 s に上がり始め、23〜25 s に +0.24 K でピーク、その後 50 s までに初期値より 0.15 K 低い値へ戻る。実験の圧力 (P_EVA) は τ≈19〜21 s に +13 kPa 跳ね、P_RES は 21〜29 s に +7 → +22 kPa (圧力計の絶対値は不正確)。")