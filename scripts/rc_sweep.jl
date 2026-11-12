using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# rc_sweep.jl … r_c_wick を 1.5 / 1.85 / 2.2 µm に振って (現在の設定のまま) 計算し、HB・EVA1・VL1 を実験と比べる
#   使い方: julia rc_sweep.jl      (scripts/ に置く。transient_flat.jl と同じフォルダ。実験の csv は data/ に置く。引数なし)
#   しくみ: param_sweep.jl と同じ。Main.SWEEP を定義して transient_flat.jl を 1 回ずつ include する (stage = :short、t_end = rc_t_end、r_c_wick だけ置き換え)。
#           結果は res/transient_flat_ammonia<rc_tag_prefix><値>/ に保存される。
#   出すもの: res/<rc_out_name>/ に r_c_wick ごとの図 (rc_<値>um.png)、全ケース重ね図 (rc_all.png)、csv。ターミナルに表 (実験との差、n=0.5 の時刻)。
#   図: 横軸 τ = 加熱開始からの時間、縦軸 温度 [°C]。実線 = モデル、破線 = 実験。HB / EVA1 / VL1。
#   モデル側の位置: HB = hb の平均、EVA1 = 蒸発器ケース側面の節点 (eca)、VL1 = 蒸気管の 1/6 の位置の壁温 (param_sweep.jl と同じ)。
#   RUN = false にすると計算せず、保存済みの結果だけで図と表を作る。SKIP_EXISTING = true なら保存済みのケースは飛ばす。
# =============================================================================

# ---------------------------- 設定 ----------------------------
rc_values      = [1.5e-6, 1.85e-6, 2.2e-6]   # [m]
rc_t_end       = 120.0      # 計算の終わり [s] (加熱開始は 10 s)。長くするほど遅い
rc_t_heat      = 10.0       # transient_flat.jl の t_heat と同じ
rc_t_exp_onset = 303.8      # 実験 (exp_5W_12min.csv) の加熱開始時刻 [s]
rc_exp_file    = "exp_5W_12min.csv"
rc_tag_prefix  = "_rcS_"    # 結果フォルダ名の接尾辞の頭 (ほかの計算を上書きしないため)
rc_out_name    = "rc_sweep"
rc_xlim        = (0.0, 60.0)   # 図の横軸 τ [s]
rc_ylim        = :auto         # 縦軸 [°C]。例: (20, 30)
rc_report_tau  = [20.0, 30.0, 60.0, 110.0]   # 表に出す τ [s] (rc_t_end − rc_t_heat 以下のもの)
rc_size        = (900, 600)
rc_dpi         = 300
rc_font_label, rc_font_tick, rc_font_legend = 14, 12, 12
rc_lw          = 2.5
rc_RUN = true
rc_SKIP_EXISTING = false
# ----------------------------------------------------------------

using CSV, DataFrames, Statistics
import Plots

rc_name(v) = string(round(v * 1e6, digits = 2)) * "um"
rc_tag(v) = rc_tag_prefix * rc_name(v)
rc_dir(v) = joinpath(@__DIR__, "..", "res", "transient_flat_ammonia" * rc_tag(v))

if rc_RUN
    for v in rc_values
        if rc_SKIP_EXISTING && isfile(joinpath(rc_dir(v), "summary.csv")) && isfile(joinpath(rc_dir(v), "profiles.csv"))
            println("\n(飛ばす: 保存済み) r_c_wick = ", rc_name(v)); continue
        end
        println("\n", "="^70, "\n  r_c_wick = ", rc_name(v), " の計算を開始\n", "="^70)
        global SWEEP = (t_end = rc_t_end, tag = rc_tag(v), r_c_wick = v)
        try
            include(joinpath(@__DIR__, "transient_flat.jl"))
        catch err
            @warn "r_c_wick = $(rc_name(v)) の実行が途中で失敗 (保存済みの結果があれば使う)" err
        end
    end
end

# ---------------------------- 結果の読み込み ----------------------------
rc_lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
rc_res = Dict{Float64,Any}()
for v in rc_values
    d = rc_dir(v)
    if !isfile(joinpath(d, "summary.csv")) || !isfile(joinpath(d, "profiles.csv"))
        @warn "結果がありません: $(d)"; continue
    end
    summ = CSV.read(joinpath(d, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])
    idx = Dict{Tuple{String,Float64},Vector{Int}}()
    for i in 1:nrow(prof)
        push!(get!(idx, (String(prof.comp[i]), Float64(prof.t[i])), Int[]), i)
    end
    rc_res[v] = (summ = summ, prof = prof, idx = idx, tp = sort(unique(prof.t)))
end
isempty(rc_res) && error("結果がありません。rc_RUN = true で計算してください")
rc_ok = [v for v in rc_values if haskey(rc_res, v)]

exp_path = joinpath(@__DIR__, "..", "data", rc_exp_file)
isfile(exp_path) || error("実験データが見つかりません: ", exp_path)
rc_ex = CSV.read(exp_path, DataFrame)
rc_exv(cols, τ) = sum(rc_lin(rc_ex.t, Vector(rc_ex[!, cn]), rc_t_exp_onset + τ) for cn in cols) / length(cols)

K = 273.15
rc_tmax(r) = r.summ.t[end] - rc_t_heat
rc_col(r, col, τ) = rc_lin(r.summ.t, Vector(r.summ[!, col]), rc_t_heat + τ)
function rc_cells(r, comp, τ)
    t = r.tp[argmin(abs.(r.tp .- (rc_t_heat + τ)))]
    haskey(r.idx, (comp, Float64(t))) || return nothing
    return sort(r.prof[r.idx[(comp, Float64(t))], :], :cell)
end
function rc_wall(r, comp, τ)
    d = rc_cells(r, comp, τ)
    return d === nothing ? NaN : mean(d.T_wall) - K
end
function rc_wallpos(r, comp, τ, pos)
    d = rc_cells(r, comp, τ)
    d === nothing && return NaN
    N = nrow(d)
    return rc_lin([(i - 0.5) / N for i in 1:N], Vector(d.T_wall), pos) - K
end
function rc_eva(r, τ)
    v = rc_wall(r, "eca", τ)
    return isnan(v) ? rc_col(r, :T_ec, τ) - K : v
end
rc_sensors = [
    ("HB",   ["HB"],  (r, τ) -> rc_wall(r, "hb", τ),            :black),
    ("EVA1", ["EC1"], (r, τ) -> rc_eva(r, τ),                   :royalblue),
    ("VL1",  ["VL1"], (r, τ) -> rc_wallpos(r, "vl", τ, 1 / 6),  :crimson),
]

# ---------------------------- 表 ----------------------------
fmt(x, d = 2) = lpad(isnan(x) ? "-" : string(round(x, digits = d)), 9)
function rc_t50(r)
    τs = [τ for τ in 5.0:0.1:min(60.0, rc_tmax(r))]
    nn = [rc_col(r, :n_nuc, τ) for τ in τs]
    i = findfirst(>=(0.5), nn)
    return i === nothing ? NaN : (i == 1 ? τs[1] : τs[i-1] + (0.5 - nn[i-1]) / (nn[i] - nn[i-1]) * (τs[i] - τs[i-1]))
end
tau_rep = [τ for τ in rc_report_tau if τ <= minimum(rc_tmax(rc_res[v]) for v in rc_ok)]
println("\nn = 0.5 の時刻 [s] (実験の立ち上がり ≈ 20 s):")
for v in rc_ok
    println("  r_c_wick = ", rpad(rc_name(v), 8), fmt(rc_t50(rc_res[v])))
end
for τ in tau_rep
    println("\nτ = ", τ, " s の実験との差 (モデル − 実験) [K]:")
    println(rpad("r_c_wick", 12), join([lpad(s[1], 9) for s in rc_sensors]))
    for v in rc_ok
        println(rpad(rc_name(v), 12), join([fmt(s[3](rc_res[v], τ) - rc_exv(s[2], τ)) for s in rc_sensors]))
    end
end

# ---------------------------- csv と図 ----------------------------
out = joinpath(@__DIR__, "..", "res", rc_out_name); mkpath(out)
rows = DataFrame(r_c_wick_um = Float64[], tau = Float64[], sensor = String[], model = Float64[], exp = Float64[])
for v in rc_ok, τ in 0.0:0.5:rc_tmax(rc_res[v]), s in rc_sensors
    push!(rows, (v * 1e6, τ, s[1], s[3](rc_res[v], τ), rc_exv(s[2], τ)))
end
CSV.write(joinpath(out, "rc_sweep.csv"), rows)

function rc_base_plot(ttl)
    p = Plots.plot(; xlabel = "Time from heating onset [s]", ylabel = "Temperature [°C]", xlims = rc_xlim, title = ttl, size = rc_size, dpi = rc_dpi,
                   legend = :topleft, legend_columns = 2, framestyle = :box, guidefontsize = rc_font_label, tickfontsize = rc_font_tick, legendfontsize = rc_font_legend)
    rc_ylim === :auto || Plots.plot!(p, ylims = rc_ylim)
    return p
end
for v in rc_ok
    p = rc_base_plot("r_c,wick = " * string(round(v * 1e6, digits = 2)) * " µm")
    d = rows[rows.r_c_wick_um .== v * 1e6, :]
    for s in rc_sensors
        e = d[d.sensor .== s[1], :]
        Plots.plot!(p, e.tau, e.model; label = s[1] * " (model)", color = s[4], lw = rc_lw)
        Plots.plot!(p, e.tau, e.exp; label = s[1] * " (exp.)", color = s[4], lw = rc_lw, ls = :dash)
    end
    Plots.savefig(p, joinpath(out, "rc_" * rc_name(v) * ".png"))
end
# 全ケース重ね (センサごとに 1 枚): 実験 = 灰色の破線
for s in rc_sensors
    p = rc_base_plot(s[1] * ": effect of r_c,wick")
    shades = [:royalblue, :seagreen, :orangered, :purple]
    for (j, v) in enumerate(rc_ok)
        e = rows[(rows.r_c_wick_um .== v * 1e6) .& (rows.sensor .== s[1]), :]
        Plots.plot!(p, e.tau, e.model; label = string(round(v * 1e6, digits = 2)) * " µm", color = shades[mod1(j, 4)], lw = rc_lw)
    end
    e = rows[(rows.r_c_wick_um .== rc_ok[1] * 1e6) .& (rows.sensor .== s[1]), :]
    Plots.plot!(p, e.tau, e.exp; label = "experiment", color = :gray, lw = rc_lw, ls = :dash)
    Plots.savefig(p, joinpath(out, "rc_all_" * s[1] * ".png"))
end
println("\n保存しました: ", abspath(out), " (rc_<値>um.png, rc_all_<センサ>.png, rc_sweep.csv)")