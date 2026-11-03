using Pkg; Pkg.activate(".")
using CSV, DataFrames

# =============================================================================
# wall_heat_diag.jl … 壁面凝縮の熱 (Q_wall) が secondary wick (swick) に入った後、どこへ行くかを見る診断 (計算のやり直しは不要)。
#   使い方: julia wall_heat_diag.jl      (この .jl は diag/ に置く。設定は下の 3 行をこのファイルで書き換える)
#   読むもの: res/<result_name>/summary.csv  … 列 Q_wall, Q_cc, Q_cs, Q_leak, T_sw, T_cc (cc_wall_cond = true で再計算した結果。
#             savefuncs.jl の新しい版で保存した結果にだけ入っている)
#
#   swick の熱収支 (residual_swick!):  dE/dt = ṁ(h_in − h_out) + Q_cs − Q_nc + Q_leak ,   Q_nc = Q_cc − Q_wall
#     Q_wall … 壁面凝縮の熱 (CC の流体から swick へ)                       [W]
#     Q_cc   … swick → CC の流体の自然対流 (nc.Q_fluid。swick の温度が上がると増える) [W]
#     Q_cs   … CC ケース → swick の熱伝導 (swick の温度が上がると減る)           [W]
#     Q_leak … primary wick → swick (swick の温度が上がると減る)                [W]
#   Q_wall が入った後、基準時刻 τ_ref (凝縮の前) からの各項の変化を積分して、Q_wall の何割がどこへ出たかを出す。
#     流体へ戻る割合   = ∫(Q_cc − Q_cc,ref) dt / ∫Q_wall dt
#     CC ケースへ      = −∫(Q_cs − Q_cs,ref) dt / ∫Q_wall dt
#     primary wick へ  = −∫(Q_leak − Q_leak,ref) dt / ∫Q_wall dt
#     残り (swick に蓄積 + 流れ ṁ の持ち出し) = 1 − 上の 3 つの和
#   ※ 基準からの変化で見る近似 (他の原因による変化も入る)。
# =============================================================================
result_name = "transient_flat_ammonia_v7a"   # 見る結果のフォルダ名 (res/ の下)。再計算した結果
t_heat = 10.0                                # モデルの加熱開始時刻 [s] (stage = :long)
τ_ref = 15.0                                 # 基準の時刻 (加熱開始から。凝縮が始まる前)

result_dir = joinpath(@__DIR__, "..", "res", result_name)
summ = CSV.read(joinpath(result_dir, "summary.csv"), DataFrame)
for col in (:Q_wall, :Q_cc, :Q_cs, :Q_leak, :T_sw)
    String(col) in names(summ) || error("summary.csv に列 $(col) がありません。新しい savefuncs.jl で計算し直してください")
end
t = summ.t
iref = argmin(abs.(t .- (t_heat + τ_ref)))
ref = (Q_cc = summ.Q_cc[iref], Q_cs = summ.Q_cs[iref], Q_leak = summ.Q_leak[iref], T_sw = summ.T_sw[iref])
println("結果: ", result_name, "   基準 τ = ", τ_ref, " s:  Q_cc = ", round(ref.Q_cc, digits = 3), " W,  Q_cs = ", round(ref.Q_cs, digits = 3),
        " W,  Q_leak = ", round(ref.Q_leak, digits = 3), " W,  T_sw − T_cc = ", round(ref.T_sw - summ.T_cc[iref], digits = 3), " K")

# 台形則の積分 (基準時刻から τ まで)
function integ(y, i0, i1)
    s = 0.0
    for k in (i0 + 1):i1
        s += 0.5 * (y[k] + y[k-1]) * (t[k] - t[k-1])
    end
    return s
end

println("\n加熱開始からの時間 τ [s] ごと:  瞬時値 [W] と、基準 τ_ref からの積算 [J]")
println(rpad("  τ", 6), lpad("Q_wall", 9), lpad("ΔQ_cc", 9), lpad("ΔQ_cs", 9), lpad("ΔQ_leak", 9), lpad("T_sw-T_cc", 11),
        lpad("∫Q_wall", 10), lpad("→流体", 9), lpad("→CCケース", 11), lpad("→primary", 10), lpad("残り", 8), lpad("流体へ%", 9), lpad("ケースへ%", 10), lpad("wickへ%", 9))
for τ in (20, 22, 24, 26, 28, 30, 34, 38, 40, 42, 46, 50, 60)
    i1 = argmin(abs.(t .- (t_heat + τ)))
    i1 > iref || continue
    ∫wall = integ(summ.Q_wall, iref, i1)
    ∫nc   = integ(summ.Q_cc .- ref.Q_cc, iref, i1)
    ∫cs   = -integ(summ.Q_cs .- ref.Q_cs, iref, i1)
    ∫lk   = -integ(summ.Q_leak .- ref.Q_leak, iref, i1)
    rest  = ∫wall - ∫nc - ∫cs - ∫lk
    pc(x) = abs(∫wall) > 1e-6 ? string(round(100 * x / ∫wall, digits = 0)) : "NaN"
    println(lpad("+" * string(τ), 5), " ",
            lpad(round(summ.Q_wall[i1], digits = 3), 8), " ", lpad(round(summ.Q_cc[i1] - ref.Q_cc, digits = 3), 8), " ",
            lpad(round(summ.Q_cs[i1] - ref.Q_cs, digits = 3), 8), " ", lpad(round(summ.Q_leak[i1] - ref.Q_leak, digits = 3), 8),
            lpad(round(summ.T_sw[i1] - summ.T_cc[i1], digits = 3), 11),
            lpad(round(∫wall, digits = 2), 10), lpad(round(∫nc, digits = 2), 9), lpad(round(∫cs, digits = 2), 11), lpad(round(∫lk, digits = 2), 10),
            lpad(round(rest, digits = 2), 8), lpad(pc(∫nc), 9), lpad(pc(∫cs), 10), lpad(pc(∫lk), 9))
end