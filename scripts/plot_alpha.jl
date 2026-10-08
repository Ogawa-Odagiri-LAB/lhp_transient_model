using Plots, CSV, DataFrames

const RES_DIR = "./res/vars_transient/"
const FIG_DIR = RES_DIR * "figs/"
# α を持つ流体コンポーネント（ループの流れ順）。wick は α が NaN なので除外
const FLUID_COMPONENTS = ["grv", "vl", "cond", "ll", "cc"]

ensure_fig_dir() = ispath(FIG_DIR) || mkpath(FIG_DIR)

# "<prefix>_<時刻>s.csv" を探して時刻順に並べる
# （"var_1.5s.csv" → 1.5。整数部と小数部をまとめて取り出す）
function list_var_files(dir_path=RES_DIR; prefix="var")
    ispath(dir_path) || error("ディレクトリ $dir_path が存在しません。")
    pattern = Regex("^" * prefix * "_(\\d+(?:\\.\\d+)?)s\\.csv\$")
    files = filter(f -> occursin(pattern, f), readdir(dir_path))
    isempty(files) && error("$(prefix)_*s.csv が $dir_path に見つかりません。")
    pairs = [(file=f, time=parse(Float64, match(pattern, f).captures[1])) for f in files]
    sort!(pairs, by = x -> x.time)
    return pairs
end

read_var(t; dir_path=RES_DIR, prefix="var") =
    CSV.read(joinpath(dir_path, "$(prefix)_$(t)s.csv"), DataFrame)

# 流体ノード数（wick を除く）
n_fluid_nodes(df) = count(in(FLUID_COMPONENTS), df.Component)

# 時系列データを読み込む。
# 別の分割数で計算した古い結果が同じフォルダに残っていると配列の大きさが
# そろわないため、最も新しく保存されたファイルと同じノード数のものだけを使う
function load_latest_run(dir_path=RES_DIR; prefix="var")
    pairs = list_var_files(dir_path; prefix)
    paths = [joinpath(dir_path, p.file) for p in pairs]
    dfs   = [CSV.read(path, DataFrame) for path in paths]

    newest = argmax(mtime.(paths))
    n_ref  = n_fluid_nodes(dfs[newest])
    keep   = [n_fluid_nodes(df) == n_ref for df in dfs]

    n_skip = count(!, keep)
    if n_skip > 0
        skipped = [p.file for (p, k) in zip(pairs, keep) if !k]
        msg = "ノード数が最新の結果（$(n_ref) ノード）と異なるファイルを $(n_skip) 個スキップしました。" *
              "古い計算結果が残っている可能性があります。"
        @warn msg skipped
    end
    times = [p.time for (p, k) in zip(pairs, keep) if k]
    return times, dfs[keep]
end

# ------------------------------------------------------------
# 1. 各コンポーネントの相の内訳を表示（どこが単相になっているかの確認用）
# ------------------------------------------------------------
function summarize_α(t; dir_path=RES_DIR, prefix="var")
    df = read_var(t; dir_path, prefix)
    println("t = $t s の α の内訳")
    for comp in FLUID_COMPONENTS
        α = df[df.Component .== comp, :α]
        isempty(α) && continue
        n_vap = count(>=(1.0), α)   # 蒸気単相（f_state で α ≥ 1 は 1 に丸められる）
        n_liq = count(<=(0.0), α)   # 液単相
        n_tp  = length(α) - n_vap - n_liq
        println(rpad(comp, 5),
                ": 二相 ", lpad(n_tp, 3),
                " / 蒸気単相(α≥1) ", lpad(n_vap, 3),
                " / 液単相(α≤0) ", lpad(n_liq, 3),
                "   min = ", round(minimum(α), sigdigits=4),
                ", max = ", round(maximum(α), sigdigits=4))
    end
end

# ------------------------------------------------------------
# 2. ある時刻の α の空間分布
# ------------------------------------------------------------
function plot_α_profile(t; dir_path=RES_DIR, prefix="var")
    df = read_var(t; dir_path, prefix)

    plot(
        title = "Void Fraction Distribution (t = $(t) s)",
        xlabel = "Position from Component Inlet [m]",
        ylabel = "α [-]",
        legend = :outerright,
        grid = true,
        size = (800, 500)
    )
    hline!([0.0, 1.0], color = :gray, linestyle = :dash, label = "")

    for comp in FLUID_COMPONENTS
        df_sub = df[df.Component .== comp, :]
        isempty(df_sub) && continue
        plot!(
            df_sub.Position,
            df_sub.α,
            label = comp,
            marker = :circle,
            markersize = 3,
            linewidth = 2
        )
    end

    # [0, 1] の範囲外になっている点を赤い×で強調
    df_out = df[in.(df.Component, Ref(FLUID_COMPONENTS)) .& ((df.α .< 0.0) .| (df.α .> 1.0)), :]
    if !isempty(df_out)
        scatter!(df_out.Position, df_out.α,
                 label = "out of [0,1]", marker = :x, markersize = 6, color = :red)
        @warn "t = $t s: α が [0, 1] の範囲外のノードが $(nrow(df_out)) 個あります"
    end

    ensure_fig_dir()
    savefig(joinpath(FIG_DIR, "alpha_profile_$(t)s.png"))
    return current()
end

# ------------------------------------------------------------
# 3. 代表ノードの α の時間変化
#    node には番号のほか :end（最終ノード）、:mid（中央ノード）も指定可
# ------------------------------------------------------------
const DEFAULT_NODES = [
    ("grv", 1), ("vl", 1), ("vl", :end),
    ("cond", 1), ("cond", :mid), ("cond", :end),
    ("ll", 1), ("ll", :end), ("cc", 1),
]

function plot_α_transient(; dir_path=RES_DIR, prefix="var", nodes=DEFAULT_NODES)
    times, dfs = load_latest_run(dir_path; prefix)

    plot(
        title = "Void Fraction History",
        xlabel = "Time [s]",
        ylabel = "α [-]",
        legend = :outerright,
        grid = true,
        size = (800, 500)
    )
    hline!([0.0, 1.0], color = :gray, linestyle = :dash, label = "")

    for (comp, node) in nodes
        vals = Float64[]
        n_used = 0
        for df in dfs
            sub = df[df.Component .== comp, :]
            n_max = maximum(sub.Node)
            n_used = node === :end ? n_max :
                     node === :mid ? (n_max + 1) ÷ 2 : node
            push!(vals, sub[sub.Node .== n_used, :α][1])
        end
        plot!(times, vals, label = "$(comp)[$(n_used)]",
              marker = :circle, markersize = 3, linewidth = 2)
    end

    display(current())
    ensure_fig_dir()
    savefig(joinpath(FIG_DIR, "alpha_transient.png"))
    return current()
end

# ------------------------------------------------------------
# 4. ループ全体の α を「ノード番号 × 時刻」のヒートマップで表示
#    （単相/二相の境界がどう動くかを一目で確認できる）
# ------------------------------------------------------------
function plot_α_heatmap(; dir_path=RES_DIR, prefix="var")
    times, dfs = load_latest_run(dir_path; prefix)

    rows = Vector{Vector{Float64}}()
    counts = Int[]
    for (k, df) in enumerate(dfs)
        subs = [df[df.Component .== c, :] for c in FLUID_COMPONENTS]
        push!(rows, reduce(vcat, [s.α for s in subs]))
        if k == 1
            counts = [nrow(s) for s in subs]
        end
    end
    A = reduce(vcat, permutedims.(rows))   # (時刻数) × (ノード数)

    # コンポーネント境界と目盛り位置
    edges = cumsum(counts)
    starts = [1; edges[1:end-1] .+ 1]
    mids = (starts .+ edges) ./ 2

    heatmap(
        1:size(A, 2), times, A;
        clims = (0.0, 1.0),       # 範囲外の値は端の色で表示される
        c = :viridis,
        title = "Void Fraction Map (grv → vl → cond → ll → cc)",
        xlabel = "Component",
        ylabel = "Time [s]",
        colorbar_title = "α",
        xticks = (mids, FLUID_COMPONENTS),
        size = (900, 500)
    )
    vline!(edges[1:end-1] .+ 0.5, color = :white, linestyle = :dash, label = "")

    display(current())
    ensure_fig_dir()
    savefig(joinpath(FIG_DIR, "alpha_heatmap.png"))
    return current()
end

# ------------------------------------------------------------
# 実行例
# ------------------------------------------------------------
summarize_α(1.0)
plot_α_profile(1.0)
plot_α_transient()
plot_α_heatmap()