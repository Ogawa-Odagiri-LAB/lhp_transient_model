using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# tau_nuc_sweep.jl … τ_nuc (核生成ゲート n の時定数) を 0.5, 1.5, 3.0 s で続けて計算し、
#   加熱開始からの時間 τ = 14〜32 s の範囲で結果を比べる専用のスクリプト。
#   使い方: julia tau_nuc_sweep.jl      (scripts/ に置く。transient_flat.jl と同じフォルダ。実験の csv は data/ に置く)
#   しくみ: transient_flat.jl を 1 回ずつ include する。そのとき Main.SWEEP を定義しておくと、transient_flat.jl は
#           stage = :short (t_end = T_END) にして τ_nuc だけ置き換える (ほかの設定は transient_flat.jl のまま。
#           SWEEP が無い単独実行では何も変わらない)。結果は res/transient_flat_ammonia_tau<値>/ に保存される。
#   出すもの: 表 (n, Q_evap, ΔP_grv, EC の比較)、指標、res/tau_nuc_sweep/ の csv と図。
#   注: 計算時間はこのファイルを書いた時点では未測定 (起動の 45 s 分だけなので 1 本あたり数分〜十数分の見込み)。
#       RUN = false にすると、計算をやり直さず、保存済みの結果だけを比べる。
# =============================================================================

# ---------------------------- 設定 ----------------------------
sw_tau_list = (0.5, 1.5, 3.0)        # 試す τ_nuc [s]。この順に実行する (文字列にしたときのフォルダ名: _tau0.5, _tau1.5, _tau3.0)
sw_t_end    = 45.0                   # 計算の終わり [s] (加熱開始は 10 s なので、τ = 35 s まで)
sw_t_heat   = 10.0                   # transient_flat.jl の t_heat と同じ
sw_t_exp_onset = 303.8               # 実験の加熱開始時刻 [s]
sw_tau_tab  = collect(14.0:2.0:32.0) # 表に出す τ [s]
sw_RUN      = true                   # false: 計算せず、保存済みの結果だけ読んで比べる
sw_FIGURE   = true                   # 図を保存するか
# ----------------------------------------------------------------

using CSV, DataFrames, Statistics
import Plots

if sw_RUN
    for τn in sw_tau_list
        println("\n", "="^70, "\n  τ_nuc = ", τn, " s の計算を開始\n", "="^70)
        global SWEEP = (τ_nuc = Float64(τn), t_end = sw_t_end)
        try
            include(joinpath(@__DIR__, "transient_flat.jl"))
        catch err
            @warn "τ_nuc = $(τn) の実行が途中で失敗 (保存済みの結果があれば比較に使う)" err
        end
    end
end

# ---------------------------- 結果の読み込み ----------------------------
sw_dir(τn) = joinpath(@__DIR__, "..", "res", "transient_flat_ammonia_tau" * string(Float64(τn)))
sw_lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))

sw_res = Dict{Float64,Any}()
for τn in sw_tau_list
    d = sw_dir(τn)
    if !isfile(joinpath(d, "summary.csv")) || !isfile(joinpath(d, "profiles.csv"))
        @warn "結果がありません: $(d)"
        continue
    end
    summ = CSV.read(joinpath(d, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])
    sw_res[Float64(τn)] = (summ = summ, prof = prof, tp = sort(unique(prof.t)))
end
isempty(sw_res) && error("読める結果がありません。sw_RUN = true で計算してください")
sw_taus = [τn for τn in sw_tau_list if haskey(sw_res, Float64(τn))]

exp_file = joinpath(@__DIR__, "..", "data", "exp_5W_12min.csv")
isfile(exp_file) || error("exp_5W_12min.csv が見つかりません: ", exp_file)
sw_ex = CSV.read(exp_file, DataFrame)
sw_exv(cols, τ) = sum(sw_lin(sw_ex.t, Vector(sw_ex[!, c]), sw_t_exp_onset + τ) for c in cols) / length(cols)

K = 273.15
# モデルの値 (τ における): summary の列はそのまま (線形補間)、壁温度は profiles.csv の最も近い時刻の値
sw_col(r, col, τ) = sw_lin(r.summ.t, Vector(r.summ[!, col]), sw_t_heat + τ)
function sw_wall(r, comp, τ)
    t = r.tp[argmin(abs.(r.tp .- (sw_t_heat + τ)))]
    d = r.prof[(r.prof.comp .== comp) .& (r.prof.t .== t), :]
    nrow(d) == 0 && return NaN
    return mean(d.T_wall) - K
end
sw_P0(r) = sw_lin(r.summ.t, Vector(r.summ.P_grv), sw_t_heat - 0.5)
sw_dP(r, τ) = (sw_col(r, :P_grv, τ) - sw_P0(r)) / 1e3
function sw_EC(r, τ)               # EC のセンサ = 蒸発器ケースの側面の節点 S (comp = "eca")。なければ受熱面 P (summary の T_ec)
    v = sw_wall(r, "eca", τ)
    return isnan(v) ? sw_col(r, :T_ec, τ) - K : v
end
sw_tmax(r) = r.summ.t[end] - sw_t_heat

# ---------------------------- 表 ----------------------------
function sw_table(title, f; digits = 3)
    println("\n", title)
    println(rpad("τ [s]", 7), join([lpad("τ_nuc=" * string(τn), 13) for τn in sw_taus]))
    for τ in sw_tau_tab
        row = String[]
        for τn in sw_taus
            r = sw_res[Float64(τn)]
            push!(row, τ <= sw_tmax(r) ? lpad(string(round(f(r, τ), digits = digits)), 13) : lpad("-", 13))
        end
        println(rpad(τ, 7), join(row))
    end
end
sw_table("ゲート n (界面蒸発の活性度、0〜1)", (r, τ) -> sw_col(r, :n_nuc, τ))
sw_table("蒸発の熱 Q_evap [W]", (r, τ) -> sw_col(r, :Q_evap, τ), digits = 2)
sw_table("EC→ウィックの熱 Q_ec_wick [W]", (r, τ) -> sw_col(r, :Q_ec_wick, τ), digits = 2)
sw_table("グルーブ圧力の上昇 ΔP_grv [kPa]", (r, τ) -> sw_dP(r, τ), digits = 2)
sw_table("EC センサ (側面 S) [°C]", (r, τ) -> sw_EC(r, τ), digits = 3)
println("\n実験 EC (EC1, EC2 の平均) [°C] / 実験 HB / 実験 RES:")
println(rpad("τ [s]", 7), lpad("EC", 9), lpad("HB", 9), lpad("RES", 9))
for τ in sw_tau_tab
    println(rpad(τ, 7), lpad(string(round(sw_exv(["EC1", "EC2"], τ), digits = 3)), 9), lpad(string(round(sw_exv(["HB"], τ), digits = 3)), 9),
            lpad(string(round(sw_exv(["RES1", "RES2"], τ), digits = 3)), 9))
end
sw_table("EC の差 (モデル − 実験) [K]", (r, τ) -> sw_EC(r, τ) - sw_exv(["EC1", "EC2"], τ))
sw_table("HB の差 (モデル − 実験) [K]", (r, τ) -> sw_wall(r, "hb", τ) - sw_exv(["HB"], τ))
sw_table("RES の差 (モデル − 実験、CC の壁) [K]", (r, τ) -> sw_wall(r, "cc", τ) - sw_exv(["RES1", "RES2"], τ))

# ---------------------------- 指標 ----------------------------
println("\n指標 (τ = 14〜32 s の範囲):")
println(rpad("τ_nuc", 8), lpad("n=0.5 の時刻", 14), lpad("Q_evap 最大", 13), lpad("その時刻", 10), lpad("ΔP(26 s)", 11), lpad("EC の最低", 11), lpad("その時刻", 10))
sw_metrics = DataFrame(tau_nuc = Float64[], t_n50 = Float64[], Q_evap_max = Float64[], t_Q_evap_max = Float64[], dP_26 = Float64[], EC_min = Float64[], t_EC_min = Float64[])
τfine = collect(14.0:0.25:32.0)
for τn in sw_taus
    r = sw_res[Float64(τn)]
    τs = [τ for τ in τfine if τ <= sw_tmax(r)]
    nn = [sw_col(r, :n_nuc, τ) for τ in τs]
    i50 = findfirst(>=(0.5), nn)
    t50 = i50 === nothing ? NaN : (i50 == 1 ? τs[1] : τs[i50-1] + (0.5 - nn[i50-1]) / (nn[i50] - nn[i50-1]) * (τs[i50] - τs[i50-1]))
    qq = [sw_col(r, :Q_evap, τ) for τ in τs]
    iq = argmax(qq)
    # EC の最低: ピーク後 (Q_evap が最大になった後) の最小値
    ec = [sw_EC(r, τ) for τ in τs]
    ie = iq - 1 + argmin(ec[iq:end])
    dP26 = 26.0 <= sw_tmax(r) ? sw_dP(r, 26.0) : NaN
    push!(sw_metrics, (Float64(τn), t50, qq[iq], τs[iq], dP26, ec[ie], τs[ie]))
    println(rpad(τn, 8), lpad(round(t50, digits = 2), 14), lpad(round(qq[iq], digits = 2), 13), lpad(τs[iq], 10), lpad(round(dP26, digits = 2), 11),
            lpad(round(ec[ie], digits = 3), 11), lpad(τs[ie], 10))
end
println("\n見方: n=0.5 の時刻の差 = τ_nuc による蒸発の立ち上がりの遅れ。Q_evap 最大と EC の最低は、立ち上がり直後の蒸発の勢いと、EC の過熱度の崩れの大きさ。")

# ---------------------------- 保存 ----------------------------
sw_out = joinpath(@__DIR__, "..", "res", "tau_nuc_sweep")
mkpath(sw_out)
rows = DataFrame(tau_nuc = Float64[], tau = Float64[], n = Float64[], Q_evap = Float64[], Q_ec_wick = Float64[], dP_grv = Float64[], EC = Float64[], EC_exp = Float64[],
                 HB = Float64[], HB_exp = Float64[], RES = Float64[], RES_exp = Float64[])
for τn in sw_taus
    r = sw_res[Float64(τn)]
    for τ in τfine
        τ <= sw_tmax(r) || continue
        push!(rows, (Float64(τn), τ, sw_col(r, :n_nuc, τ), sw_col(r, :Q_evap, τ), sw_col(r, :Q_ec_wick, τ), sw_dP(r, τ), sw_EC(r, τ), sw_exv(["EC1", "EC2"], τ),
                     sw_wall(r, "hb", τ), sw_exv(["HB"], τ), sw_wall(r, "cc", τ), sw_exv(["RES1", "RES2"], τ)))
    end
end
CSV.write(joinpath(sw_out, "tau_nuc_sweep.csv"), rows)
CSV.write(joinpath(sw_out, "tau_nuc_metrics.csv"), sw_metrics)
println("\n保存しました: ", abspath(sw_out), " (tau_nuc_sweep.csv, tau_nuc_metrics.csv)")

if sw_FIGURE
    try
        cols = [:royalblue, :black, :orangered, :seagreen, :purple]
        panels = []
        function sw_panel(ycol, ylab; with_exp = nothing)
            p = Plots.plot(xlabel = "tau [s]", ylabel = ylab, legend = :best, xlims = (14, 32), framestyle = :box)
            for (j, τn) in enumerate(sw_taus)
                d = rows[rows.tau_nuc .== Float64(τn), :]
                Plots.plot!(p, d.tau, d[!, ycol], label = "tau_nuc = $(τn) s", color = cols[j], lw = 2)
            end
            if with_exp !== nothing
                d = rows[rows.tau_nuc .== Float64(sw_taus[1]), :]
                Plots.plot!(p, d.tau, d[!, with_exp], label = "experiment", color = :gray, ls = :dash, lw = 2)
            end
            return p
        end
        p1 = sw_panel(:n, "gate n [-]")
        p2 = sw_panel(:Q_evap, "Q_evap [W]")
        p3 = sw_panel(:dP_grv, "dP_grv [kPa]")
        p4 = sw_panel(:EC, "EC sensor [°C]"; with_exp = :EC_exp)
        p5 = sw_panel(:HB, "HB [°C]"; with_exp = :HB_exp)
        p6 = sw_panel(:RES, "RES (CC wall) [°C]"; with_exp = :RES_exp)
        fig = Plots.plot(p1, p2, p3, p4, p5, p6, layout = (3, 2), size = (1100, 1000))
        Plots.savefig(fig, joinpath(sw_out, "tau_nuc_sweep.png"))
        println("図を保存しました: ", abspath(joinpath(sw_out, "tau_nuc_sweep.png")))
    catch err
        @warn "図の保存に失敗 (表と csv は出力済み)" err
    end
end