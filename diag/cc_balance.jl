using Pkg; Pkg.activate(".")
using CSV, DataFrames

# =============================================================================
# cc_balance.jl … 保存済みの計算結果 (res/<result_name>/summary.csv, profiles.csv) から、
#   CC の熱収支と「還流の過冷度で吸収できる熱 (G_sub = ṁ·c_p)」、ループ全体の蓄熱を時間ごとに出す。計算のやり直しは不要。
#   使い方: julia cc_balance.jl      (diag/ に置く)
# 出すもの (加熱開始からの時間 τ ごと):
#   T_cc, dT_cc/dt (±30 s の中心差分), ṁ,
#   Q_sub = ṁ (h_LL最後 − h_CC) … 還流液が CC から奪う熱 (負)。G_sub = −Q_sub/(T_cc − T_LL最後) [W/K] (= ṁ c_p の目安)
#   Q_ec_cc … 蒸発器ケース → CC ケース (summary の列)、Q_cc … CC 流体の自然対流の熱 (summary の列 nc.Q_fluid)
#   Q_amb_LL … 液管 (全セル) が周囲から受ける熱の目安 (G_AMB_CELL × Σ(T_amb − T_wall))。CC・蒸発器・HB は周囲と熱交換しない
#   Q_br … 蒸発器ケース → 液管 CC 側の熱橋 (summary の列)
#   蓄熱 = C_hb dT_hb/dt + C_ec dT_ec/dt + C_cc dT_cc/dt の目安 (印加熱 5 W に対する割合)
#   感度: 1/G_sub [K/W] と、CC に 0.01 W 余分に入ったときの T_cc の上昇の目安
# =============================================================================
result_name = "transient_flat_ammonia_v8e"   # 800 s の結果を見るときは "transient_flat_ammonia" + run_tag。2 h は "transient_flat_ammonia_2h" + run_tag
result_dir = joinpath(@__DIR__, "..", "res", result_name)
t_heat = 10.0
T_amb = 21.5
G_AMB_CELL = 3.3e-3        # 液管 1 セルの周囲との熱コンダクタンス [W/K] (目安。断熱材 ID10/OD30 の値)
C_hb, C_ec, C_cc = 47.2, 100.0, 80.0    # [J/K]。C_hb は transient_flat.jl の値、C_ec と C_cc は目安 (要確認)
K = 273.15

summ = CSV.read(joinpath(result_dir, "summary.csv"), DataFrame)
prof = CSV.read(joinpath(result_dir, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :rho, :Em, :P, :T, :T_wall])
nearest_s(t) = argmin(abs.(summ.t .- t))
tprof = sort(unique(prof.t))
nearest_p(t) = tprof[argmin(abs.(tprof .- t))]
function cell(comp, t, which)
    d = prof[(prof.comp .== comp) .& (prof.t .== t), :]
    d = sort(d, :cell)
    return which === :first ? d[1, :] : d[end, :]
end
hspec(r) = r.Em + r.P / r.rho
ddt(col, t) = (summ[nearest_s(t + 30), col] - summ[nearest_s(t - 30), col]) / (summ.t[nearest_s(t + 30)] - summ.t[nearest_s(t - 30)])
ndec(x, n = 3) = lpad(string(round(x, digits = n)), 9)

τs = [τ for τ in (60, 90, 120, 180, 240, 300, 400, 500, 600, 700, 780, 1200, 1800, 2400, 3600, 4800, 6000, 6800) if t_heat + τ + 30 <= summ.t[end]]
println("CC の熱収支 (モデル、加熱開始から τ [s])。熱は CC が受ける向きが正。")
println(rpad("τ", 6), "  T_cc[C]  dT_cc/dt[mK/s]  mdot[mg/s]  Q_sub[W]  G_sub[W/K]  Q_ec_cc[W]  Q_cc[W]  Q_amb_LL[W]  Q_br[W]  T_ec[C]  T_LL末[C]  蓄熱[W]  蓄熱/Q  1/G_sub[K/W]  ΔT_cc(0.01W)[K]")
for τ in τs
    t = t_heat + τ
    i = nearest_s(t)
    tp = nearest_p(t)
    rl = cell("ll", tp, :last)
    rc = cell("cc", tp, :first)
    ṁ = summ.mdot[i]
    if any(prof.comp .== "ccl")
        # 新しい CC (ccl): comp "cc" の行は気液混合 (x > 0 なら潜熱が混ざる) なので使わない。
        # 液の比熱を液管の全セルの hspec と T の傾きから求め、Q_sub = −ṁ c_p (T_l1 − T_LL最後) とする。
        dl = sort(prof[(prof.comp .== "ll") .& (prof.t .== tp), :], :cell)
        cp_l = (hspec(dl[end, :]) - hspec(dl[1, :])) / (dl.T[end] - dl.T[1])
        Tl1 = only(prof[(prof.comp .== "ccl") .& (prof.t .== tp) .& (prof.cell .== 1), :T])
        Q_sub = -ṁ * cp_l * (Tl1 - rl.T)
    else
        Q_sub = ṁ * (hspec(rl) - hspec(rc))
    end
    dT = summ.T_cc[i] - rl.T
    G_sub = abs(dT) > 1e-6 ? -Q_sub / dT : NaN
    dll = prof[(prof.comp .== "ll") .& (prof.t .== tp), :]
    Q_amb = G_AMB_CELL * sum((T_amb + K) .- dll.T_wall)
    stor = C_hb * ddt(:T_hb, t) + C_ec * ddt(:T_ec, t) + C_cc * ddt(:T_cc, t)
    println(rpad(τ, 6), ndec(summ.T_cc[i] - K, 2), ndec(1e3 * ddt(:T_cc, t), 2), "  ", ndec(1e6 * ṁ, 2), ndec(Q_sub, 4),
            ndec(G_sub, 4), ndec(summ.Q_ec_cc[i], 4), ndec(summ.Q_cc[i], 4), ndec(Q_amb, 4), ndec(summ.Q_br[i], 4), ndec(summ.T_ec[i] - K, 2), ndec(rl.T - K, 2),
            ndec(stor, 3), ndec(stor / max(summ.Q_apply[i], 1e-9), 3), ndec(1 / G_sub, 1), ndec(0.01 / G_sub, 3))
end
println("\n注: Q_sub = ṁ (h_LL最後 − h_CC) は還流液が CC を冷やす熱 (負)。G_sub ≈ ṁ c_p が小さいほど、CC に入る小さな熱でも T_cc が大きく動く。")
println("    蓄熱の C_ec, C_cc は目安 (transient_flat.jl の値に合わせて上の設定を直す)。")