using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# tamb_sweep.jl … 周囲温度 T_amb を 現状 (21.5 °C) / 22.3 °C にして、:long (816 s) を続けて 2 本計算し、実験と比べる
#   使い方: julia tamb_sweep.jl      (scripts/ に置く。transient_flat.jl・heat_budget_calc.jl と同じフォルダ。実験の csv は data/ に置く)
#   しくみ: f_shah_sweep.jl と同じ。Main.SWEEP を定義して transient_flat.jl を 1 回ずつ include する (T_amb は transient_flat.jl の sweep_get(:T_amb, …) で置き換わる)。
#   順序: base (T_amb = 21.5 °C のまま) → T_amb = 22.3 °C。結果は res/transient_flat_ammonia<ps_tag_prefix><名前>/ に保存 (budget.csv も出る)。
#   出すもの: 実験との差の表、基準との差の表、熱収支 (budget.csv) の比較、res/tamb_sweep/ の csv と図。
#   注: 1 本あたり :long で約 1400〜1600 s (前の記録)。2 本で 1 時間弱。ps_SKIP_EXISTING = true で保存済みのケースを飛ばす。ps_RUN = false なら比較だけ。
# =============================================================================

# ---------------------------- 設定 ----------------------------
ps_stage       = :long      # :long (816 s) / :short (t = ps_t_end まで) / :long2h
ps_t_end       = 816.0      # 表・図の範囲に使う終わり [s] (:long は 816 s)
ps_t_heat      = 10.0       # transient_flat.jl の t_heat と同じ
ps_t_exp_onset = 303.8      # 実験の加熱開始時刻 [s]
ps_tag_prefix = "_psTamb_"     # 結果フォルダ名の接尾辞の頭 (前の感度解析の _ps_ を上書きしないため)
ps_out_name   = "tamb_sweep"   # 表・csv・図の保存先 (res/<名前>/)
ps_base   = (T_amb = 21.5 + 273.15,)   # 基準値 (transient_flat.jl の現在の値 [K])。基準が 21.5 °C を兼ねるので、試す値は 22.3 °C の 1 つだけ
ps_values = (
    T_amb = [22.3 + 273.15],                          # [K]
)
# 比較にだけ使う保存済みの結果 (計算しない): (表示名, res/ の下のフォルダ名)
ps_refs = Tuple{String,String}[]
ps_report_tau = [30.0, 60.0, 120.0, 290.0, 500.0, 800.0]   # 表に出す τ [s] (ps_t_end − ps_t_heat 以下のもの)
ps_RUN = true
ps_SKIP_EXISTING = false
ps_FIGURE = true
# ----------------------------------------------------------------

using CSV, DataFrames, Statistics
import Plots

# ---------------------------- ケースの組み立て ----------------------------
ps_cases = Any[(name = "base", par = :base, val = NaN, pars = NamedTuple(), dir = nothing)]
for par in keys(ps_values)
    for v in ps_values[par]
        push!(ps_cases, (name = string(par) * "=" * string(Float64(v)), par = par, val = Float64(v), pars = NamedTuple{(par,)}((Float64(v),)), dir = nothing))
    end
end
ps_tag(c) = ps_tag_prefix * replace(c.name, "=" => "_", "+" => "")
ps_dir(c) = c.dir === nothing ? joinpath(@__DIR__, "..", "res", "transient_flat_ammonia" * ps_tag(c)) : joinpath(@__DIR__, "..", "res", c.dir)
println("実行するケース (", length(ps_cases), " 本):")
for c in ps_cases
    println("  ", c.name)
end

if ps_RUN
    for c in ps_cases
        if ps_SKIP_EXISTING && isfile(joinpath(ps_dir(c), "summary.csv")) && isfile(joinpath(ps_dir(c), "profiles.csv"))
            println("\n(飛ばす: 保存済み) ", c.name)
            continue
        end
        println("\n", "="^70, "\n  ケース ", c.name, " の計算を開始\n", "="^70)
        global SWEEP = merge((stage = ps_stage, t_end = ps_t_end, tag = ps_tag(c)), c.pars)   # t_end は :short のときだけ使われる
        try
            include(joinpath(@__DIR__, "transient_flat.jl"))
        catch err
            @warn "ケース $(c.name) の実行が途中で失敗 (保存済みの結果があれば比較に使う)" err
        end
    end
end

for (nm, folder) in ps_refs
    push!(ps_cases, (name = nm, par = :ref, val = NaN, pars = NamedTuple(), dir = folder))
end

# ---------------------------- 結果の読み込み ----------------------------
ps_lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
ps_res = Dict{String,Any}()
for c in ps_cases
    d = ps_dir(c)
    if !isfile(joinpath(d, "summary.csv")) || !isfile(joinpath(d, "profiles.csv"))
        @warn "結果がありません: $(d)"
        continue
    end
    summ = CSV.read(joinpath(d, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])
    idx = Dict{Tuple{String,Float64},Vector{Int}}()      # (comp, t) → profiles の行番号
    for i in 1:nrow(prof)
        push!(get!(idx, (String(prof.comp[i]), Float64(prof.t[i])), Int[]), i)
    end
    ps_res[c.name] = (summ = summ, prof = prof, idx = idx, tp = sort(unique(prof.t)))
end
haskey(ps_res, "base") || error("基準 (base) の結果がありません。ps_RUN = true で計算してください")
ps_ok = [c for c in ps_cases if haskey(ps_res, c.name)]

exp_file = joinpath(@__DIR__, "..", "data", "exp_5W_12min.csv")
isfile(exp_file) || error("exp_5W_12min.csv が見つかりません: ", exp_file)
ps_ex = CSV.read(exp_file, DataFrame)
ps_exv(cols, τ) = sum(ps_lin(ps_ex.t, Vector(ps_ex[!, cn]), ps_t_exp_onset + τ) for cn in cols) / length(cols)

K = 273.15
ps_tmax(r) = r.summ.t[end] - ps_t_heat
ps_col(r, col, τ) = ps_lin(r.summ.t, Vector(r.summ[!, col]), ps_t_heat + τ)
function ps_cells(r, comp, τ)
    t = r.tp[argmin(abs.(r.tp .- (ps_t_heat + τ)))]
    haskey(r.idx, (comp, Float64(t))) || return nothing
    return sort(r.prof[r.idx[(comp, Float64(t))], :], :cell)
end
function ps_wall(r, comp, τ)
    d = ps_cells(r, comp, τ)
    d === nothing && return NaN
    return mean(d.T_wall) - K
end
function ps_wallpos(r, comp, τ, pos)
    d = ps_cells(r, comp, τ)
    d === nothing && return NaN
    N = nrow(d)
    return ps_lin([(i - 0.5) / N for i in 1:N], Vector(d.T_wall), pos) - K
end
function ps_ec(r, τ)                   # EC のセンサ = 蒸発器ケースの側面の節点 S (comp = "eca")。なければ受熱面 P (summary の T_ec)
    v = ps_wall(r, "eca", τ)
    return isnan(v) ? ps_col(r, :T_ec, τ) - K : v
end
ps_sensors = [
    ("HB",   ["HB"],           (r, τ) -> ps_wall(r, "hb", τ)),
    ("EC",   ["EC1", "EC2"],   (r, τ) -> ps_ec(r, τ)),
    ("RES",  ["RES1", "RES2"], (r, τ) -> ps_wall(r, "cc", τ)),
    ("VL1",  ["VL1"],          (r, τ) -> ps_wallpos(r, "vl", τ, 1 / 6)),
    ("CON1", ["CON1"],         (r, τ) -> ps_wallpos(r, "cond", τ, 0.0)),
    ("CON5", ["CON5"],         (r, τ) -> ps_wallpos(r, "cond", τ, 1.0)),
    ("LL1",  ["LL1"],          (r, τ) -> ps_wallpos(r, "ll", τ, 0.10 / 1.15)),
    ("LL3",  ["LL3"],          (r, τ) -> ps_wallpos(r, "ll", τ, (1.15 - 0.04) / 1.15)),
]
ps_tau_rep = [τ for τ in ps_report_tau if τ <= minimum(ps_tmax(ps_res[c.name]) for c in ps_ok)]
rb = ps_res["base"]
fmt(x, d = 3) = lpad(isnan(x) ? "-" : string(round(x, digits = d)), 9)

# ---------------------------- 基準の実験との差 ----------------------------
println("\n基準 (base) の実験との差 (モデル − 実験) [K]:")
println(rpad("τ [s]", 7), join([lpad(s[1], 9) for s in ps_sensors]))
for τ in ps_tau_rep
    println(rpad(τ, 7), join([fmt(s[3](rb, τ) - ps_exv(s[2], τ)) for s in ps_sensors]))
end

# ---------------------------- 全ケースの実験との差 (モデル − 実験): パラメータを決めるときはこの表を見る ----------------------------
for τ in ps_tau_rep
    println("\nτ = ", τ, " s の実験との差 (モデル − 実験) [K]:")
    println(rpad("ケース", 24), join([lpad(s[1], 9) for s in ps_sensors]))
    for c in ps_ok
        r = ps_res[c.name]
        println(rpad(c.name, 24), join([fmt(s[3](r, τ) - ps_exv(s[2], τ)) for s in ps_sensors]))
    end
end

# ---------------------------- 感度: 基準との差 (ケース − 基準) ----------------------------
for τ in ps_tau_rep
    println("\nτ = ", τ, " s の基準との差 (ケース − 基準) [K]:")
    println(rpad("ケース", 24), join([lpad(s[1], 9) for s in ps_sensors]))
    for c in ps_ok
        c.name == "base" && continue
        r = ps_res[c.name]
        println(rpad(c.name, 24), join([fmt(s[3](r, τ) - s[3](rb, τ)) for s in ps_sensors]))
    end
end

# ---------------------------- 感度の指標: ln(値) あたりの変化 ----------------------------
println("\n感度指標 S = (最大値のケース − 最小値のケース) / ln(値の比) [K / ln 単位]。基準を含めた 3 点以上のうちの両端:")
for τ in ps_tau_rep
    println("τ = ", τ, " s")
    println(rpad("パラメータ", 12), join([lpad(s[1], 9) for s in ps_sensors]))
    for par in keys(ps_values)
        cs = [c for c in ps_ok if c.par == par]
        isempty(cs) && continue
        vals = vcat([Float64(getfield(ps_base, par))], [c.val for c in cs])
        names = vcat(["base"], [c.name for c in cs])
        imin, imax = argmin(vals), argmax(vals)
        println(rpad(string(par), 12), join([fmt((s[3](ps_res[names[imax]], τ) - s[3](ps_res[names[imin]], τ)) / log(vals[imax] / vals[imin])) for s in ps_sensors]))
    end
end

# ---------------------------- 立ち上がりの指標 ----------------------------
println("\n立ち上がりの指標 (τ = 10〜60 s):")
println(rpad("ケース", 24), lpad("n=0.5 の時刻", 14), lpad("Q_evap 最大", 13), lpad("その時刻", 10), lpad("ΔP(30 s)", 11))
ps_metrics = DataFrame(case = String[], par = String[], val = Float64[], t_n50 = Float64[], Q_evap_max = Float64[], t_Q_evap_max = Float64[], dP_30 = Float64[])
for c in ps_ok
    r = ps_res[c.name]
    τs = [τ for τ in 10.0:0.25:60.0 if τ <= ps_tmax(r)]
    nn = [ps_col(r, :n_nuc, τ) for τ in τs]
    i50 = findfirst(>=(0.5), nn)
    t50 = i50 === nothing ? NaN : (i50 == 1 ? τs[1] : τs[i50-1] + (0.5 - nn[i50-1]) / (nn[i50] - nn[i50-1]) * (τs[i50] - τs[i50-1]))
    qq = [ps_col(r, :Q_evap, τ) for τ in τs]
    iq = argmax(qq)
    P0 = ps_lin(r.summ.t, Vector(r.summ.P_grv), ps_t_heat - 0.5)
    dP30 = 30.0 <= ps_tmax(r) ? (ps_col(r, :P_grv, 30.0) - P0) / 1e3 : NaN
    push!(ps_metrics, (c.name, string(c.par), c.val, t50, qq[iq], τs[iq], dP30))
    println(rpad(c.name, 24), fmt(t50, 2) * "     ", fmt(qq[iq], 2) * "    ", fmt(τs[iq], 2), fmt(dP30, 2))
end

# ---------------------------- 熱収支 (budget.csv) の比較 ----------------------------
println("\n熱収支 (budget.csv) [W]: Q_amb = 蒸気管 + 液管の外気への漏れ、Q_house = 断熱材・PEEK、Q_plate = 凝縮器プレート → コールドプレート")
for τ in [60.0, 120.0, 290.0, 500.0, 800.0]
    τ <= minimum(ps_tmax(ps_res[c.name]) for c in ps_ok) || continue
    println("τ = ", τ, " s")
    println(rpad("ケース", 24), lpad("Q_apply", 9), lpad("Q_plate", 9), lpad("Q_amb", 9), lpad("Q_house", 9), lpad("dE/dt", 9), lpad("resid", 9))
    for c in ps_ok
        f = joinpath(ps_dir(c), "budget.csv")
        isfile(f) || (println(rpad(c.name, 24), "  budget.csv なし"); continue)
        b = CSV.read(f, DataFrame)
        g(col) = ps_lin(Vector(b.tau), Vector(b[!, col]), τ)
        println(rpad(c.name, 24), fmt(g(:Q_apply), 3), fmt(g(:Q_plate), 3), fmt(g(:Q_amb_vl) + g(:Q_amb_ll), 3), fmt(g(:Q_house), 3), fmt(g(:dE_total), 3), fmt(g(:resid), 3))
    end
end

# ---------------------------- 保存 ----------------------------
ps_out = joinpath(@__DIR__, "..", "res", ps_out_name)
mkpath(ps_out)
ps_tau_curve = collect(10.0:2.0:ps_t_end - ps_t_heat)
rows = DataFrame(case = String[], par = String[], val = Float64[], tau = Float64[], sensor = String[], model = Float64[], exp = Float64[])
for c in ps_ok
    r = ps_res[c.name]
    for τ in ps_tau_curve
        τ <= ps_tmax(r) || continue
        for s in ps_sensors
            push!(rows, (c.name, string(c.par), c.val, τ, s[1], s[3](r, τ), ps_exv(s[2], τ)))
        end
    end
end
CSV.write(joinpath(ps_out, "param_sweep.csv"), rows)
CSV.write(joinpath(ps_out, "param_metrics.csv"), ps_metrics)
println("\n保存しました: ", abspath(ps_out), " (param_sweep.csv, param_metrics.csv)")

if ps_FIGURE
    try
        cols = [:black, :royalblue, :orangered, :seagreen, :purple]
        for par in keys(ps_values)
            cs = [c for c in ps_ok if c.par == par || c.name == "base" || c.par == :ref]
            length(cs) < 2 && continue
            panels = []
            for sn in ("HB", "EC", "RES", "CON1")
                p = Plots.plot(xlabel = "tau [s]", ylabel = sn * " [degC]", legend = (sn == "HB" ? :bottomright : false), framestyle = :box)
                for (j, c) in enumerate(cs)
                    d = rows[(rows.case .== c.name) .& (rows.sensor .== sn), :]
                    Plots.plot!(p, d.tau, d.model, label = c.name, color = cols[mod1(j, length(cols))], lw = 2)
                end
                d = rows[(rows.case .== "base") .& (rows.sensor .== sn), :]
                Plots.plot!(p, d.tau, d.exp, label = "experiment", color = :gray, ls = :dash, lw = 2)
                push!(panels, p)
            end
            fig = Plots.plot(panels..., layout = (2, 2), size = (1100, 800), plot_title = string(par))
            Plots.savefig(fig, joinpath(ps_out, "param_sweep_" * string(par) * ".png"))
        end
        println("図を保存しました: ", abspath(ps_out), " (param_sweep_<パラメータ>.png)")
    catch err
        @warn "図の保存に失敗 (表と csv は出力済み)" err
    end
end