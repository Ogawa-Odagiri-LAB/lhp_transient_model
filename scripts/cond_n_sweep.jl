using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# cond_N_sweep.jl … 凝縮器の管の分割数 N (50 / 75 / 100) を変えて :long (τ ≈ 806 s) を計算し、長期の温度上昇・凝縮前線・CON1 が N で変わるかを見る。
#   【未実行】 Claude は Julia を実行できないので、構文・動作は未確認。最初のエラーを貼ってください。
#   使い方: julia cond_N_sweep.jl   (scripts/ に置く。transient_flat.jl と同じフォルダ。実験の csv は data/ に置く)
#   しくみ: hev_transition_sweep.jl と同じ。Main.SWEEP に cond_N を入れて transient_flat.jl を 1 ケースずつ include する。
#   条件: h_ev = 50000 (据え置き)、新しい幾何 (A_ec_wick = 364.5 mm²、C_ec_P = 8.374、溝長 27 mm)、h_con = 16500 × f_con_area 0.75。N 以外は全ケースで同じ。
#         チラー側の管 (ColdPlateGeometry の N = 50) とプレート格子 (30×6) は変えない。CON1 は全ケースで入口から 10 mm の位置を補間して読む (transient_flat.jl の con1_pos)。
#   見るもの: (1) CON1・Tsat−CON1・RES−CON1・HB の実験との差 (時間窓ごとの平均差と RMS)
#             (2) τ150 → 800 s の変化量 (HB・CON1・Tsat−CON1): 長期の上昇が N で変わるか
#             (3) 蒸気前線の位置 (蒸気セルの長さ) と二相長: N で収束するか
#             (4) N=50 を基準にした差 (N=75, 100 − N=50)。CON1 が 0.05 K 以内なら N は問題にならない
#   出すもの: 表、res/<cn_out_name>/ に csv と図 (cn_overview.svg/.png)
#   注: cn_stage = :long は 816 s まで。N=75, 100 は計算時間が N=50 (約 2200 s) の 1.5〜2.5 倍の見込み (未測定)。
#       cn_SKIP_EXISTING = true なら保存済みのケースは計算しない。cn_RUN = false なら保存済みの結果だけで表と図を作る。
# =============================================================================

# ---------------------------- 設定 ----------------------------
cn_stage       = :long        # :long (816 s) | :long2h (6910 s、4 時間以上)
cn_t_end       = cn_stage === :long2h ? 6910.0 : 816.0
cn_t_heat      = 10.0         # transient_flat.jl の t_heat と同じ (モデル時刻 = τ + 10 s)
cn_t_exp_onset = 303.8        # 実験の加熱開始時刻 [s]
cn_L_cond      = 0.56         # 凝縮器の管の長さ [m] (transient_flat.jl の CondGeometry の L)
cn_x_con1      = 0.010        # CON1 の位置 (入口から) [m]
cn_tag_prefix  = "_condN_"    # 結果フォルダ名の接尾辞の頭 (res/transient_flat_ammonia<接尾辞>/)
cn_out_name    = "cond_N_sweep"
cn_N_values    = [50, 75, 100]
cn_colors      = [:royalblue, :seagreen, :darkorange]
cn_RUN           = true
cn_SKIP_EXISTING = true
cn_FIGURE        = true
cn_windows     = [(0.0, 50.0), (50.0, 150.0), (150.0, Inf)]
cn_report_tau  = [20.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 300.0, 500.0, 800.0]
# ----------------------------------------------------------------

using CSV, DataFrames, Statistics
import Plots

cn_cases = [(name = "N" * string(N), N = N, color = cn_colors[mod1(j, length(cn_colors))]) for (j, N) in enumerate(cn_N_values)]
cn_tag(c) = cn_tag_prefix * c.name
cn_dir(c) = joinpath(@__DIR__, "..", "res", (cn_stage === :long2h ? "transient_flat_ammonia_2h" : "transient_flat_ammonia") * cn_tag(c))
println("ケース (", length(cn_cases), " 本), stage = ", cn_stage, ":")
for c in cn_cases
    println("  N = ", c.N, "   → ", cn_dir(c))
end

if cn_RUN
    for c in cn_cases
        if cn_SKIP_EXISTING && isfile(joinpath(cn_dir(c), "summary.csv")) && isfile(joinpath(cn_dir(c), "profiles.csv"))
            println("\n(飛ばす: 保存済み) ", c.name)
            continue
        end
        println("\n", "="^70, "\n  ケース ", c.name, " の計算を開始\n", "="^70)
        global SWEEP = (t_end = cn_t_end, tag = cn_tag(c), stage = cn_stage, cond_N = c.N)
        try
            include(joinpath(@__DIR__, "transient_flat.jl"))
        catch err
            @warn "ケース $(c.name) の実行が途中で失敗 (保存済みの結果があれば比較に使う)" err
        end
    end
end

# ---------------------------- 結果の読み込み ----------------------------
cn_lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
const CN_K = 273.15
cn_res = Dict{String,Any}()
for c in cn_cases
    d = cn_dir(c)
    if !isfile(joinpath(d, "summary.csv")) || !isfile(joinpath(d, "profiles.csv"))
        @warn "結果がありません: $(d)"
        continue
    end
    summ = CSV.read(joinpath(d, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :x, :T_wall])
    prof.comp = String.(prof.comp)      # InlineStrings → String (groupby のキーの検索のため)
    gd = groupby(prof, [:comp, :t])
    cn_res[c.name] = (summ = summ, gd = gd, tp = sort(unique(prof.t)), N = c.N)
end
cn_ok = [c for c in cn_cases if haskey(cn_res, c.name)]
isempty(cn_ok) && error("結果が 1 つもありません。cn_RUN = true で計算してください")

# 実験データ (列名は「101_HB」→「HB」、「RES_2」→「RES2」に正規化)
function cn_load_exp(stage)
    names_try = stage === :long2h ? ["exp_5W_startup.csv", "exp_5w_startup.csv", "exp_5W_2h.csv"] :
                                    ["exp_5W_12min.csv", "exp_5W_startup.csv", "exp_5w_startup.csv"]
    for nm in names_try
        f = joinpath(@__DIR__, "..", "data", nm)
        isfile(f) || continue
        df = CSV.read(f, DataFrame; silencewarnings = true)
        newn = [replace(replace(strip(String(n)), r"^\d+_" => ""), "RES_2" => "RES2") for n in names(df)]
        rename!(df, DataFrames.make_unique(Symbol.(newn)))
        need = [:t, :HB, :RES1, :RES2, :CON1]
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
cn_ex = cn_load_exp(cn_stage)
cn_exv(cols, τ) = sum(cn_lin(cn_ex.t, cn_ex[!, cn], cn_t_exp_onset + τ) for cn in cols) / length(cols)

# ---------------------------- モデルの値の取り出し ----------------------------
function cn_cells(r, comp, τ)
    t = r.tp[argmin(abs.(r.tp .- (cn_t_heat + τ)))]
    haskey(r.gd, (comp, Float64(t))) || return nothing
    return sort(DataFrame(r.gd[(comp, Float64(t))]), :cell)
end
function cn_wall(r, comp, τ)
    d = cn_cells(r, comp, τ)
    d === nothing && return NaN
    return mean(d.T_wall) - CN_K
end
function cn_con1(r, τ)      # 入口から 10 mm の壁温 (セル中心の線形補間)
    d = cn_cells(r, "cond", τ)
    d === nothing && return NaN
    N = nrow(d)
    return cn_lin([(i - 0.5) / N for i in 1:N], Vector(d.T_wall), cn_x_con1 / cn_L_cond) - CN_K
end
cn_tsat(r, τ) = cn_lin(r.summ.t, Vector(r.summ.T_cc), cn_t_heat + τ) - CN_K      # CC の飽和温度 (T_cc) [°C]
function cn_front(r, τ)     # (蒸気セルの長さ [mm] = x>0.5 のセル数 × dl, 二相長 [mm] = 0.001<x<0.999 のセル数 × dl)
    d = cn_cells(r, "cond", τ)
    d === nothing && return (NaN, NaN)
    dl = cn_L_cond / nrow(d) * 1e3
    return (count(>(0.5), d.x) * dl, count(x -> 0.001 < x < 0.999, d.x) * dl)
end

cn_tmax = min(minimum(cn_res[c.name].summ.t[end] - cn_t_heat for c in cn_ok), cn_ex.t[end] - cn_t_exp_onset)
cn_taus = collect(0.0:2.0:cn_tmax)
println("比較する τ の範囲: 0 〜 ", round(cn_tmax, digits = 1), " s")

sensors = [("HB", [:HB], (r, τ) -> cn_wall(r, "hb", τ)),
           ("RES", [:RES1, :RES2], (r, τ) -> cn_wall(r, "cc", τ)),
           ("CON1", [:CON1], (r, τ) -> cn_con1(r, τ))]
cn_ser = Dict{String,Dict{String,Vector{Float64}}}()
for c in cn_ok
    r = cn_res[c.name]
    s = Dict{String,Vector{Float64}}()
    for (nm, _, f) in sensors
        s[nm] = [f(r, τ) for τ in cn_taus]
    end
    s["Tsat-CON1"] = [cn_tsat(r, τ) - s["CON1"][i] for (i, τ) in enumerate(cn_taus)]
    s["RES-CON1"] = s["RES"] .- s["CON1"]
    s["len_vap"] = [cn_front(r, τ)[1] for τ in cn_taus]
    s["len_2ph"] = [cn_front(r, τ)[2] for τ in cn_taus]
    cn_ser[c.name] = s
end
cn_exs = Dict{String,Vector{Float64}}()
for (nm, cols, _) in sensors
    cn_exs[nm] = [cn_exv(cols, τ) for τ in cn_taus]
end
cn_exs["RES-CON1"] = cn_exs["RES"] .- cn_exs["CON1"]

# ---------------------------- 表 ----------------------------
fmt(x, d = 3) = lpad(isnan(x) ? "-" : string(round(x, digits = d)), 9)
idx_of(τ) = argmin(abs.(cn_taus .- τ))

println("\n", "="^100, "\n(1) 実験との差 (モデル − 実験) の時間窓ごとの平均差 (RMS) [K]\n", "="^100)
for nm in ("HB", "CON1", "RES-CON1")
    println("\n", nm, ":")
    println(rpad("N", 8), join([lpad(string(w[1]) * "–" * (isinf(w[2]) ? string(round(cn_taus[end])) : string(w[2])), 18) for w in cn_windows]))
    for c in cn_ok
        row = String[]
        for w in cn_windows
            m = [i for (i, τ) in enumerate(cn_taus) if w[1] <= τ <= w[2]]
            e = cn_ser[c.name][nm][m] .- cn_exs[nm][m]
            push!(row, lpad(string(round(mean(e), digits = 3)) * "(" * string(round(sqrt(mean(e .^ 2)), digits = 3)) * ")", 18))
        end
        println(rpad("N=" * string(c.N), 8), join(row))
    end
end

println("\n", "="^100, "\n(2) τ150 → 800 s の変化量 [K] (長期の上昇): モデル / 実験\n", "="^100)
i150, i800 = idx_of(150.0), idx_of(min(800.0, cn_taus[end]))
println(rpad("量", 12), join([lpad("N=" * string(c.N), 10) for c in cn_ok]), lpad("実験", 10))
for nm in ("HB", "CON1", "Tsat-CON1", "RES-CON1")
    row = [lpad(string(round(cn_ser[c.name][nm][i800] - cn_ser[c.name][nm][i150], digits = 3)), 10) for c in cn_ok]
    ex = haskey(cn_exs, nm) ? lpad(string(round(cn_exs[nm][i800] - cn_exs[nm][i150], digits = 3)), 10) : lpad("-", 10)
    println(rpad(nm, 12), join(row), ex)
end

println("\n", "="^100, "\n(3) 蒸気前線: 蒸気セルの長さ / 二相長 [mm] (x>0.5 のセル、0.001<x<0.999 のセル)\n", "="^100)
println(rpad("τ [s]", 8), join([lpad("N=" * string(c.N) * " 蒸気/二相", 20) for c in cn_ok]))
for τ in cn_report_tau
    τ <= cn_taus[end] || continue
    i = idx_of(τ)
    println(rpad(string(τ), 8), join([lpad(string(round(cn_ser[c.name]["len_vap"][i], digits = 1)) * " / " * string(round(cn_ser[c.name]["len_2ph"][i], digits = 1)), 20) for c in cn_ok]))
end

println("\n", "="^100, "\n(4) N=50 を基準にした差 (N − N=50) [K]: CON1 / Tsat−CON1 / HB\n", "="^100)
if any(c -> c.N == 50, cn_ok)
    base = first(c for c in cn_ok if c.N == 50)
    println(rpad("τ [s]", 8), join([lpad("N=" * string(c.N) * " CON1,Tsat-CON1,HB", 30) for c in cn_ok if c.N != 50]))
    for τ in cn_report_tau
        τ <= cn_taus[end] || continue
        i = idx_of(τ)
        cols = [lpad(join([string(round(cn_ser[c.name][nm][i] - cn_ser[base.name][nm][i], digits = 3)) for nm in ("CON1", "Tsat-CON1", "HB")], ", "), 30) for c in cn_ok if c.N != 50]
        println(rpad(string(τ), 8), join(cols))
    end
else
    println("(N=50 の結果がないので省略)")
end

# ---------------------------- csv と図 ----------------------------
out_dir = joinpath(@__DIR__, "..", "res", cn_out_name)
mkpath(out_dir)
df = DataFrame(tau = cn_taus)
for c in cn_ok, (k, v) in cn_ser[c.name]
    df[!, Symbol(c.name * "_" * k)] = v
end
for (k, v) in cn_exs
    df[!, Symbol("exp_" * k)] = v
end
CSV.write(joinpath(out_dir, "cond_N_series.csv"), df)
println("\ncsv を保存しました: ", out_dir)

if cn_FIGURE
    function panel(title, key; ylabel = "", exp_key = nothing)
        p = Plots.plot(; title = title, xlabel = "τ [s]", ylabel = ylabel, legend = :best, titlefontsize = 9, size = (560, 360))
        exp_key !== nothing && Plots.plot!(p, cn_taus, cn_exs[exp_key]; color = :black, lw = 1.6, label = "実験")
        for c in cn_ok
            Plots.plot!(p, cn_taus, cn_ser[c.name][key]; color = c.color, lw = 1.4, label = "N = " * string(c.N))
        end
        return p
    end
    fig = Plots.plot(panel("CON1 (入口から 10 mm) [°C]", "CON1"; ylabel = "°C", exp_key = "CON1"),
                     panel("HB [°C]", "HB"; ylabel = "°C", exp_key = "HB"),
                     panel("Tsat − CON1 [K]", "Tsat-CON1"; ylabel = "K"),
                     panel("RES − CON1 [K]", "RES-CON1"; ylabel = "K", exp_key = "RES-CON1"),
                     panel("蒸気セルの長さ (x>0.5) [mm]", "len_vap"; ylabel = "mm"),
                     panel("二相長 [mm]", "len_2ph"; ylabel = "mm");
                     layout = (3, 2), size = (1150, 1000))
    Plots.savefig(fig, joinpath(out_dir, "cn_overview.svg"))
    Plots.savefig(fig, joinpath(out_dir, "cn_overview.png"))
    println("Saved: ", joinpath(out_dir, "cn_overview.svg"), " (.png)")
end
println("\n完了。結果: ", out_dir)