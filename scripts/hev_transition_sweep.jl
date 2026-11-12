using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# hev_transition_sweep.jl … 蒸発の熱伝達係数 h_ev を定数として何本か動かし、実験と比べる (:long = τ ≈ 806 s まで)。
#   【未実行】 Claude は Julia を実行できないので、構文・動作は未確認。最初のエラーを貼ってください。
#   (h_ev を時間変化させる「遷移」のケースは削除した。ファイル名は経緯のまま。src 側の h_ev_at() は h_ev_end = NaN なら定数で、動作は従来どおり)
#
#   使い方: julia hev_transition_sweep.jl   (scripts/ に置く。transient_flat.jl と同じフォルダ。実験の csv は data/ に置く)
#   しくみ: h_sweep.jl と同じ。Main.SWEEP を定義して transient_flat.jl を 1 ケースずつ include する。
#           SWEEP のキー: h_ev (定数)。
#   ケース: 定数 (ev_const_values)。増減できる。
#   出すもの: 表 (実験との差、時間窓ごとの平均差と RMS)、res/<ev_out_name>/ に csv と図 (svg と png)。
#     hev_overview.*  … 全区間 (h_ev(τ)、EC−RES、HB、EC、RES、CON1)
#     hev_zoom.*      … 立ち上がりの拡大 (τ ≤ ev_zoom_max)
#     hev_error.*     … モデル − 実験 (同じ 6 枚。0 に近いほど良い)
#   注: ev_stage = :long は 816 s (τ ≈ 806 s) まで。数時間の時間変化を見るには :long2h (6910 s、計算に 4 時間ほど) に変える。
#       ev_SKIP_EXISTING = true なら、同じ名前の結果が保存済みのケースは計算しない (図だけ作り直したいとき用)。
#       ev_RUN = false にすると、計算せず、保存済みの結果だけで表と図を作る。
#   ケース名は h_ev から自動で決まる (別の設定は別フォルダ = 上書きしない)。
# =============================================================================

# ---------------------------- 設定 ----------------------------
ev_stage       = :long        # :long (816 s) | :long2h (6910 s)
ev_t_end       = ev_stage === :long2h ? 6910.0 : 816.0   # 計算の終わり [s] (transient_flat.jl の stage に合わせる)
ev_t_heat      = 10.0         # transient_flat.jl の t_heat と同じ (モデル時刻 = τ + 10 s)
ev_t_exp_onset = 303.8        # 実験の加熱開始時刻 [s] (τ = 実験時刻 − 303.8 s)
ev_tag_prefix  = "_hevT_hcon_"     # 結果フォルダ名の接尾辞の頭 (res/transient_flat_ammonia<接尾辞>/)
ev_out_name    = "hev_transition_hcon"   # 表・csv・図の保存先 (res/<名前>/)
ev_RUN           = true
ev_SKIP_EXISTING = true
ev_FIGURE        = true
ev_zoom_max    = 150.0        # 拡大図の τ の上限 [s]
ev_windows     = [(0.0, 50.0), (50.0, 150.0), (150.0, Inf)]   # 平均差・RMS を出す τ の窓 [s]
ev_report_tau  = [20.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0, 290.0, 500.0, 800.0]   # 表に出す τ [s] (計算範囲内のもの)

# 【計算するケース】定数 h_ev
#   h_con = 16500 × f_con_area 0.75 (要望書 2026-11-02) に変えたので、h_sweep.jl の保存済み (h_con = 20000 × 1.0) とは比べられない。全ケースをここで計算し直す。
ev_const_values = [100000.0, 50000.0, 25000.0]
ev_const_colors = [:royalblue, :seagreen, :darkorange]
# 【計算しないケース】h_sweep.jl で保存済みの結果を比較にだけ使う: (凡例, res/ の下のフォルダ名, h_ev, 色, 線種)
#   h_sweep.jl の tag = "_psHev_" + (base | h_ev_<値>)。フォルダが無ければ警告して飛ばす。h_sweep.jl の計算は τ ≈ 290 s までなので、それ以降は線が途切れる。
ev_refs = Tuple{String,String,Float64,Symbol,Symbol}[]   # 旧 h_con (20000×1.0) の h_sweep.jl の結果は使わない。使うなら (凝縮の h_con が違う旧結果として) 手で足す
# ----------------------------------------------------------------

using CSV, DataFrames, Statistics
import Plots

ev_fmtv(x) = string(round(Int, x))   # 100000.0 → "100000"

# ---------------------------- ケースの組み立て ----------------------------
# name = 結果フォルダの接尾辞に使う名前、label = 図の凡例、pars = SWEEP に入れる値
ev_cases = Any[]
for (j, hv) in enumerate(ev_const_values)
    push!(ev_cases, (name = "const" * ev_fmtv(hv), label = "h_ev = " * ev_fmtv(hv) * " (const)",
                     pars = (h_ev = hv,), color = ev_const_colors[mod1(j, length(ev_const_colors))],
                     lw = 1.6, ls = :solid, h_nominal = (hv, hv), trans = false, ref = false, dir = nothing))
end
for (lab, folder, hv, col, ls) in ev_refs
    push!(ev_cases, (name = "ref_" * folder, label = lab, pars = NamedTuple(), color = col, lw = 1.6, ls = ls,
                     h_nominal = (hv, hv), trans = false, ref = true, dir = folder))
end

ev_tag(c) = ev_tag_prefix * c.name
ev_dir(c) = c.dir !== nothing ? joinpath(@__DIR__, "..", "res", c.dir) : joinpath(@__DIR__, "..", "res", (ev_stage === :long2h ? "transient_flat_ammonia_2h" : "transient_flat_ammonia") * ev_tag(c))
println("ケース (", length(ev_cases), " 本), stage = ", ev_stage, ":")
for c in ev_cases
    println("  ", c.ref ? "(保存済みを読むだけ) " : "(計算) ", c.name, "   → ", ev_dir(c))
end

if ev_RUN
    for c in ev_cases
        c.ref && continue
        if ev_SKIP_EXISTING && isfile(joinpath(ev_dir(c), "summary.csv")) && isfile(joinpath(ev_dir(c), "profiles.csv"))
            println("\n(飛ばす: 保存済み) ", c.name)
            continue
        end
        println("\n", "="^70, "\n  ケース ", c.name, " の計算を開始\n", "="^70)
        global SWEEP = merge((t_end = ev_t_end, tag = ev_tag(c), stage = ev_stage), c.pars)
        try
            include(joinpath(@__DIR__, "transient_flat.jl"))
        catch err
            @warn "ケース $(c.name) の実行が途中で失敗 (保存済みの結果があれば比較に使う)" err
        end
    end
end

# ---------------------------- 結果の読み込み ----------------------------
ev_lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
const EV_K = 273.15
ev_res = Dict{String,Any}()
for c in ev_cases
    d = ev_dir(c)
    if !isfile(joinpath(d, "summary.csv")) || !isfile(joinpath(d, "profiles.csv"))
        @warn "結果がありません: $(d)"
        continue
    end
    summ = CSV.read(d * "/summary.csv", DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])
    idx = Dict{Tuple{String,Float64},Vector{Int}}()      # (comp, t) → profiles の行番号
    for i in 1:nrow(prof)
        push!(get!(idx, (String(prof.comp[i]), Float64(prof.t[i])), Int[]), i)
    end
    ev_res[c.name] = (summ = summ, prof = prof, idx = idx, tp = sort(unique(prof.t)))
end
ev_ok = [c for c in ev_cases if haskey(ev_res, c.name)]
isempty(ev_ok) && error("結果が 1 つもありません。ev_RUN = true で計算してください")

# 実験データ (列名は「101_HB」→「HB」、「RES_2」→「RES2」に正規化)
function ev_load_exp(stage)
    names_try = stage === :long2h ? ["exp_5W_startup.csv", "exp_5w_startup.csv", "exp_5W_2h.csv"] :
                                    ["exp_5W_12min.csv", "exp_5W_startup.csv", "exp_5w_startup.csv"]
    for nm in names_try
        f = joinpath(@__DIR__, "..", "data", nm)
        isfile(f) || continue
        df = CSV.read(f, DataFrame; silencewarnings = true)
        newn = [replace(replace(strip(String(n)), r"^\d+_" => ""), "RES_2" => "RES2") for n in names(df)]
        rename!(df, DataFrames.make_unique(Symbol.(newn)))
        need = [:t, :HB, :EC1, :EC2, :RES1, :RES2, :CON1]
        if !all(n -> n in propertynames(df), need)
            @warn "$(nm) に必要な列がありません (見つからない列: $(setdiff(need, propertynames(df))))"
            continue
        end
        dropmissing!(df, :t)
        out = DataFrame()
        for n in need
            out[!, n] = Float64.(coalesce.(df[!, n], NaN))
        end
        println("実験データ: ", f, " (", nrow(out), " 行)")
        return out
    end
    error("実験データが data/ に見つかりません (試した名前: $(names_try))")
end
ev_ex = ev_load_exp(ev_stage)
ev_exv(cols, τ) = sum(ev_lin(ev_ex.t, ev_ex[!, cn], ev_t_exp_onset + τ) for cn in cols) / length(cols)

# ---------------------------- モデルの値の取り出し (h_sweep.jl と同じ) ----------------------------
function ev_cells(r, comp, τ)
    t = r.tp[argmin(abs.(r.tp .- (ev_t_heat + τ)))]
    haskey(r.idx, (comp, Float64(t))) || return nothing
    return sort(r.prof[r.idx[(comp, Float64(t))], :], :cell)
end
function ev_wall(r, comp, τ)
    d = ev_cells(r, comp, τ)
    d === nothing && return NaN
    return mean(d.T_wall) - EV_K
end
function ev_wallpos(r, comp, τ, pos)
    d = ev_cells(r, comp, τ)
    d === nothing && return NaN
    N = nrow(d)
    return ev_lin([(i - 0.5) / N for i in 1:N], Vector(d.T_wall), pos) - EV_K
end
function ev_ec(r, τ)    # EC のセンサ = 蒸発器ケース側面の節点 S (comp = "eca")。なければ受熱面 P (summary の T_ec)
    v = ev_wall(r, "eca", τ)
    return isnan(v) ? ev_lin(r.summ.t, Vector(r.summ.T_ec), ev_t_heat + τ) - EV_K : v
end
ev_sensors = [      # (名前, 実験の列, モデルの関数)
    ("HB",   [:HB],            (r, τ) -> ev_wall(r, "hb", τ)),
    ("EC",   [:EC1, :EC2],     (r, τ) -> ev_ec(r, τ)),
    ("RES",  [:RES1, :RES2],   (r, τ) -> ev_wall(r, "cc", τ)),
    ("CON1", [:CON1],          (r, τ) -> ev_wallpos(r, "cond", τ, 0.0)),
]

# 共通の τ 格子 (全ケースの計算範囲の共通部分)
# 各ケースの計算範囲 (h_sweep.jl の保存済みは短い)。格子は最も長いケースに合わせ、短いケースはその先を NaN にする (線が途切れる)
ev_tend = Dict(c.name => ev_res[c.name].summ.t[end] - ev_t_heat for c in ev_ok)
ev_tmax = maximum(values(ev_tend))
ev_tmax = min(ev_tmax, ev_ex.t[end] - ev_t_exp_onset)
ev_step = ev_tmax > 2000 ? 5.0 : 2.0
ev_taus = collect(0.0:ev_step:ev_tmax)
println("比較する τ の範囲: 0 〜 ", round(ev_tmax, digits = 1), " s (刻み ", ev_step, " s)")

# 系列: ev_ser[ケース名][センサ] = モデルの温度 [°C]、"EC-RES" は EC − RES [K]、"h_ev" は h_ev(τ)
ev_ser = Dict{String,Dict{String,Vector{Float64}}}()
for c in ev_ok
    r = ev_res[c.name]
    s = Dict{String,Vector{Float64}}()
    for (sn, _, f) in ev_sensors
        s[sn] = [τ <= ev_tend[c.name] ? f(r, τ) : NaN for τ in ev_taus]
    end
    s["EC-RES"] = s["EC"] .- s["RES"]
    if "h_ev" in names(r.summ)
        s["h_ev"] = [τ <= ev_tend[c.name] ? ev_lin(r.summ.t, Vector{Float64}(r.summ.h_ev), ev_t_heat + τ) : NaN for τ in ev_taus]
    else    # 古い savefuncs.jl の結果 (h_sweep.jl の保存済みなど): 定数なので公称値で代用
        c.ref || @warn "$(c.name): summary.csv に h_ev 列がありません。公称値で代用します"
        s["h_ev"] = [τ <= ev_tend[c.name] ? c.h_nominal[1] : NaN for τ in ev_taus]
    end
    ev_ser[c.name] = s
end
ev_exs = Dict{String,Vector{Float64}}()
for (sn, cols, _) in ev_sensors
    ev_exs[sn] = [ev_exv(cols, τ) for τ in ev_taus]
end
ev_exs["EC-RES"] = ev_exs["EC"] .- ev_exs["RES"]
ev_quants = ["HB", "EC", "RES", "CON1", "EC-RES"]

# ---------------------------- 表 1: τ ごとの実験との差 (モデル − 実験) ----------------------------
ev_fmt(x, d = 3) = lpad(isnan(x) ? "-" : string(round(x, digits = d)), 9)
ev_taus_rep = [τ for τ in ev_report_tau if τ <= ev_tmax]
ev_iτ(τ) = argmin(abs.(ev_taus .- τ))
println("\n", "="^78, "\n実験との差 (モデル − 実験) [K]。EC-RES は (EC − RES) の差\n", "="^78)
for q in ev_quants
    println("\n", q, ":")
    println(rpad("case", 46), join([lpad(string(round(Int, τ)), 9) for τ in ev_taus_rep]))
    for c in ev_ok
        println(rpad(c.label, 46), join([ev_fmt(ev_ser[c.name][q][ev_iτ(τ)] - ev_exs[q][ev_iτ(τ)]) for τ in ev_taus_rep]))
    end
end
println("\n実験の値 (参考, EC-RES [K], HB [°C]):")
println(rpad("", 46), join([lpad(string(round(Int, τ)), 9) for τ in ev_taus_rep]))
println(rpad("EC-RES (exp)", 46), join([ev_fmt(ev_exs["EC-RES"][ev_iτ(τ)]) for τ in ev_taus_rep]))
println(rpad("HB (exp)", 46), join([ev_fmt(ev_exs["HB"][ev_iτ(τ)], 2) for τ in ev_taus_rep]))

# ---------------------------- 表 2: 時間窓ごとの平均差 (RMS) ----------------------------
println("\n", "="^78, "\n時間窓ごとの平均差 (モデル − 実験) と RMS [K]  形式: 平均差(RMS)\n", "="^78)
ev_metrics = DataFrame(case = String[], quantity = String[], tau_from = Float64[], tau_to = Float64[], bias = Float64[], rms = Float64[])
for q in ev_quants
    println("\n", q, ":")
    wl = [(a, min(b, ev_tmax)) for (a, b) in ev_windows if a < ev_tmax]
    println(rpad("case", 46), join([lpad("τ " * string(round(Int, a)) * "–" * string(round(Int, b)), 20) for (a, b) in wl]))
    for c in ev_ok
        cells = String[]
        for (a, b) in wl
            m = findall(τ -> a <= τ <= b, ev_taus)
            dd = ev_ser[c.name][q][m] .- ev_exs[q][m]
            dd = dd[.!isnan.(dd)]
            if isempty(dd)
                push!(cells, lpad("-", 20))
                continue
            end
            bias, rms = mean(dd), sqrt(mean(dd .^ 2))
            push!(ev_metrics, (c.name, q, a, b, bias, rms))
            push!(cells, lpad(string(round(bias, digits = 3)) * "(" * string(round(rms, digits = 3)) * ")", 20))
        end
        println(rpad(c.label, 46), join(cells))
    end
end

# ---------------------------- csv ----------------------------
ev_out = joinpath(@__DIR__, "..", "res", ev_out_name)
mkpath(ev_out)
ev_long = DataFrame(case = String[], tau = Float64[], quantity = String[], model = Float64[], exp = Float64[], diff = Float64[])
for c in ev_ok, q in ev_quants, (i, τ) in enumerate(ev_taus)
    push!(ev_long, (c.name, τ, q, ev_ser[c.name][q][i], ev_exs[q][i], ev_ser[c.name][q][i] - ev_exs[q][i]))
end
CSV.write(joinpath(ev_out, "hev_series.csv"), ev_long)
ev_hdf = DataFrame(tau = ev_taus)
for c in ev_ok
    ev_hdf[!, Symbol(c.name)] = ev_ser[c.name]["h_ev"]
end
CSV.write(joinpath(ev_out, "hev_h_ev_tau.csv"), ev_hdf)
CSV.write(joinpath(ev_out, "hev_window_metrics.csv"), ev_metrics)
println("\ncsv を保存しました: ", abspath(ev_out))

# ---------------------------- 図 ----------------------------
ev_tr = [c for c in ev_ok if c.trans]
ev_ylabels = Dict("HB" => "HB [°C]", "EC" => "EC [°C]", "RES" => "RES [°C]", "CON1" => "CON1 [°C]", "EC-RES" => "EC − RES [K]")
function ev_panel(q, mode, xmax; leg = false)
    m = findall(τ -> τ <= xmax, ev_taus)
    τp = ev_taus[m]
    yl = q == "h_ev" ? "h_ev [W/m²K]" : (mode === :err ? "model − exp  " * ev_ylabels[q] : ev_ylabels[q])
    pl = Plots.plot(xlabel = "τ [s]", ylabel = yl, framestyle = :box, xlims = (0, xmax), legend = (leg ? :topleft : false), legendfontsize = 7, foreground_color_legend = nothing, background_color_legend = :white,
                    yscale = q == "h_ev" ? :log10 : :identity, title = q == "h_ev" ? "prescribed h_ev(τ)" : q, titlefontsize = 10, guidefontsize = 9)
    if q != "h_ev" && mode === :abs
        Plots.plot!(pl, τp, ev_exs[q][m], label = "experiment", color = :black, ls = :solid, lw = 1.0, alpha = 0.55)
    end
    mode === :err && Plots.hline!(pl, [0.0], color = :gray, ls = :dot, label = "")
    for c in ev_ok
        y = q == "h_ev" ? ev_ser[c.name]["h_ev"][m] :
            (mode === :err ? ev_ser[c.name][q][m] .- ev_exs[q][m] : ev_ser[c.name][q][m])
        Plots.plot!(pl, τp, y, label = c.label, color = c.color, lw = c.lw, ls = c.ls)
    end
    return pl
end
function ev_figure(mode, xmax, fname, ttl)
    qs = ["h_ev", "EC-RES", "HB", "EC", "RES", "CON1"]
    panels = [ev_panel(q, mode, xmax; leg = (q == "EC-RES")) for q in qs]   # 凡例は EC−RES (2 枚目) にだけ出す (全ケース + 実験)
    fig = Plots.plot(panels..., layout = (3, 2), size = (1150, 1000), plot_title = ttl, plot_titlefontsize = 11,
                     left_margin = 5 * Plots.mm, bottom_margin = 3 * Plots.mm)
    Plots.savefig(fig, joinpath(ev_out, fname * ".svg"))
    Plots.savefig(fig, joinpath(ev_out, fname * ".png"))
    println("Saved: ", joinpath(ev_out, fname * ".svg"), " (.png)")
end
if ev_FIGURE
    ev_figure(:abs, ev_tmax, "hev_overview", "h_ev (const) sweep: model vs experiment (τ = t_exp − 303.8 s)")
    ev_figure(:abs, min(ev_zoom_max, ev_tmax), "hev_zoom", "h_ev (const) sweep: zoom on the onset")
    ev_figure(:err, ev_tmax, "hev_error", "model − experiment (closer to 0 is better)")
end
println("\n完了。結果: ", abspath(ev_out))