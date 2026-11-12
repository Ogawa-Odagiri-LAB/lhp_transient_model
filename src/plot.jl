# =============================================================================
# plot.jl (プロット。単体で動く: LHPTransientModel に依存しない)
#
# save_results が書き出した CSV (summary.csv, profiles.csv, wick.csv, geometry.csv) を読んで、
#   ・コンポーネント別の時系列 (成分ごとに色・凡例)
#   ・任意の時刻の、LHP 全体の温度分布・乾き度 (vapor quality) の分布
# を描く。計算をやり直さずに、いつでも、何度でも、結果を眺められる。
#
# 【使い方】 plot_results.jl (スクリプト) を参照。REPL で、個別に呼ぶ場合:
#   include("src/plot.jl")
#   df_summary, df_prof, df_wick = load_results("res/transient_ammonia")
#   lengths = load_lengths("res/transient_ammonia")      # geometry.csv が無い古い結果は、手で指定:
#   # lengths = (grv = 0.035, vl = 2.0, cond = 1.0, ll = 2.0, cc = 0.07)
#
#   plt = plot_overview(df_summary, df_prof)                          # 時系列の一覧 (6 パネル)
#   plt = plot_distribution(df_prof, df_wick, lengths, 50.0; summary = df_summary)   # t=50 s の分布
#   plt = plot_distribution(df_prof, df_wick, lengths, [10.0, 50.0, 100.0]; summary = df_summary)
#   plt = plot_component_series(df_prof, :T; agg = :max)              # 任意の量・集約の時系列
#   savefig(plt, "x.png")
#
# 計算の直後に、REPL で、そのまま見る場合 (CSV を経由しない):
#   plot_overview(results_summary(sol, p, tables), results_profiles(sol, p))
#
# 【位置座標 s】 グルーブ → 蒸気管 → 凝縮器 → 液管 → CC の順に、セル中心の位置を
# 並べる (ウィックは半径方向なので、別のパネル)。s の原点はグルーブの入口。
# 各成分の長さ lengths = (grv, vl, cond, ll, cc) [m] は、geometry.csv (save_results が
# 書く) から読むか、手で指定する。セル数は、profiles.csv から数える。
#
# 【集約 agg】 成分内のセルが複数あるとき、時系列に使う代表値:
#   :mean (既定), :min, :max, :first (最上流のセル), :last (最下流のセル)
#   または cell=整数 で特定のセル。
#
# 【ラベル】 GR バックエンドの既定のフォントは、日本語が表示できないことがあるので、
# 既定は英語。日本語にするには:  set_plot_language!(:ja; fontfamily="Meiryo")
# =============================================================================

using Plots, DataFrames, CSV
using Statistics: mean

# ---------------------------------------------------------------------------
# ラベルと色
# ---------------------------------------------------------------------------

const PLOT_LABELS = Dict{String,String}(
    "grv" => "groove", "vl" => "vapor line", "cond" => "condenser", "ll" => "liquid line",
    "cc" => "CC", "evap" => "evaporator case", "wick" => "wick",
)
const PLOT_COMP_COLOR = Dict{String,Symbol}(
    "grv" => :purple, "vl" => :red, "cond" => :blue, "ll" => :green,
    "cc" => :orange, "evap" => :black, "wick" => :brown,
)
const PLOT_FLUID_COMPS = ("grv", "vl", "cond", "ll", "cc")
const PLOT_WALL_COMPS = ("vl", "cond", "ll", "cc")      # グルーブには壁がない (蒸発器ケースは別)

"""
    set_plot_language!(lang=:en; fontfamily=nothing)

ラベルの言語を切り替える。:ja で日本語。日本語の表示には、日本語のフォントが必要
(例: fontfamily="Meiryo" または "Yu Gothic")。
"""
function set_plot_language!(lang::Symbol=:en; fontfamily=nothing)
    if lang === :ja
        merge!(PLOT_LABELS, Dict("grv" => "グルーブ", "vl" => "蒸気管", "cond" => "凝縮器",
                                 "ll" => "液管", "cc" => "CC", "evap" => "蒸発器ケース", "wick" => "ウィック"))
    elseif lang === :en
        merge!(PLOT_LABELS, Dict("grv" => "groove", "vl" => "vapor line", "cond" => "condenser",
                                 "ll" => "liquid line", "cc" => "CC", "evap" => "evaporator case", "wick" => "wick"))
    else
        error("lang は :en または :ja")
    end
    fontfamily === nothing || Plots.default(fontfamily = fontfamily)
    return nothing
end

# ---------------------------------------------------------------------------
# 結果の読み込み (CSV)
# ---------------------------------------------------------------------------

"""
    load_results(dir) -> (summary, profiles, wick)

save_results が書き出した CSV を読む。
"""
function load_results(dir::AbstractString)
    summary = CSV.read(joinpath(dir, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(dir, "profiles.csv"), DataFrame)
    wick = CSV.read(joinpath(dir, "wick.csv"), DataFrame)
    return summary, prof, wick
end

"""
    load_lengths(dir) -> NamedTuple (grv, vl, cond, ll, cc) [m]

geometry.csv (save_results が書く) から、各成分の長さを読む。
"""
function load_lengths(dir::AbstractString)
    df = CSV.read(joinpath(dir, "geometry.csv"), DataFrame)
    d = Dict(String(r.comp) => Float64(r.L) for r in eachrow(df))
    return (grv = d["grv"], vl = d["vl"], cond = d["cond"], ll = d["ll"], cc = d["cc"])
end

"""
    available_times(prof) -> Vector{Float64}

保存されている時刻 (昇順)。
"""
available_times(prof::DataFrame) = sort(unique(Float64.(prof.t)))

# ---------------------------------------------------------------------------
# 成分別の時系列
# ---------------------------------------------------------------------------

_aggfun(agg::Symbol) = agg === :mean ? mean :
                       agg === :min ? minimum :
                       agg === :max ? maximum :
                       agg === :first ? first :
                       agg === :last ? last :
                       error("agg は :mean, :min, :max, :first, :last のいずれか")

"""
    component_series(prof, comp, var; agg=:mean, cell=nothing) -> (t, v)

成分 comp (例 "vl") の量 var (例 :T) の時系列。セルが複数あるときは、agg で代表値を取る。
cell に整数を渡すと、そのセルの値 (agg は無視)。
"""
function component_series(prof::DataFrame, comp::AbstractString, var::Symbol;
    agg::Symbol=:mean, cell::Union{Nothing,Integer}=nothing)

    sub = prof[prof.comp .== comp, :]
    isempty(sub) && return (Float64[], Float64[])
    if cell !== nothing
        sub = sub[sub.cell .== cell, :]
        isempty(sub) && return (Float64[], Float64[])
        order = sortperm(sub.t)
        return (collect(sub.t[order]), collect(sub[order, var]))
    end
    g = combine(groupby(sub, :t), var => _aggfun(agg) => :v)
    sort!(g, :t)
    return (collect(g.t), collect(g.v))
end

"""
    plot_component_series(prof, var; comps=PLOT_FLUID_COMPS, agg=:mean, cell=nothing,
                          ylabel=string(var), title="", rel_to_cc=false) -> Plot

量 var (:T, :P, :x, :u, :rho, :Em, :T_wall) の時系列を、成分ごとに色分けして 1 枚に描く。
rel_to_cc=true なら、各成分の値から CC の値 (同時刻) を引く (圧力差の表示用)。
"""
function plot_component_series(prof::DataFrame, var::Symbol;
    comps=PLOT_FLUID_COMPS, agg::Symbol=:mean, cell::Union{Nothing,Integer}=nothing,
    ylabel::AbstractString=string(var), title::AbstractString="", rel_to_cc::Bool=false)

    plt = plot(; xlabel = "t [s]", ylabel = ylabel, title = title, legend = :best)
    t_cc, v_cc = rel_to_cc ? component_series(prof, "cc", var; agg = agg, cell = cell) : (Float64[], Float64[])
    for comp in comps
        t, v = component_series(prof, comp, var; agg = agg, cell = cell)
        isempty(t) && continue
        if rel_to_cc
            length(v) == length(v_cc) || continue
            comp == "cc" && continue
            v = v .- v_cc
        end
        plot!(plt, t, v; label = get(PLOT_LABELS, comp, comp), color = get(PLOT_COMP_COLOR, comp, :gray), lw = 1.8)
    end
    return plt
end

"""
    plot_overview(summary, prof; agg=:mean, T_cond=nothing, T_amb=nothing, h_fg=nothing,
                  size=(1250, 950)) -> Plot

時系列の一覧 (6 パネル、成分ごとに色・凡例)。
  1. 流体温度 T (と、蒸発界面温度 T_out の破線)
  2. 壁温度 (蒸発器ケースは summary.T_ec。T_cond, T_amb を渡すと、基準線を引く)
  3. 乾き度 x
  4. 圧力差 P − P_cc (P_cc 自身の値は、図のタイトル欄に、最終時刻の値を表示)
  5. 流速 u (セルの下流側の面の速度。符号つき平均)
  6. 印加熱 Q_apply と、蒸発による熱 ṁ h_fg
agg で、成分内の代表値を選ぶ (:mean, :min, :max, :first, :last)。

summary に Q_evap の列が無い (古い結果) ときは、h_fg (T_out [K] → h_fg [J/kg] の関数)
を渡すと、Q_evap = mdot × h_fg(T_out) を計算する。例:
  h_fg = T -> PropsSI("H","T",T,"Q",1.0,"ammonia") - PropsSI("H","T",T,"Q",0.0,"ammonia")
"""
function plot_overview(summary::DataFrame, prof::DataFrame;
    agg::Symbol=:mean, T_cond=nothing, T_amb=nothing, h_fg=nothing, size=(1250, 950))

    summ = summary
    if !hasproperty(summ, :Q_evap) && h_fg !== nothing && hasproperty(summ, :T_out) && hasproperty(summ, :mdot)
        summ = copy(summary)
        summ[!, :Q_evap] = summ.mdot .* [h_fg(T) for T in summ.T_out]
    end

    # 1. 流体温度
    p1 = plot_component_series(prof, :T; comps = PLOT_FLUID_COMPS, agg = agg,
                               ylabel = "T [K]", title = "fluid temperature ($(agg))")
    hasproperty(summ, :T_out) &&
        plot!(p1, summ.t, summ.T_out; label = "T_out (interface)", color = :black, ls = :dash, lw = 1.2)

    # 2. 壁温度
    p2 = plot_component_series(prof, :T_wall; comps = PLOT_WALL_COMPS, agg = agg,
                               ylabel = "T [K]", title = "wall temperature ($(agg))")
    hasproperty(summ, :T_ec) &&
        plot!(p2, summ.t, summ.T_ec; label = PLOT_LABELS["evap"], color = PLOT_COMP_COLOR["evap"], lw = 1.8)
    T_cond === nothing || hline!(p2, [T_cond]; label = "T_cond (plate)", color = :gray, ls = :dot)
    T_amb === nothing || hline!(p2, [T_amb]; label = "T_amb", color = :gray, ls = :dashdot)

    # 3. 乾き度
    p3 = plot_component_series(prof, :x; comps = PLOT_FLUID_COMPS, agg = agg,
                               ylabel = "x [-]", title = "vapor quality ($(agg))")

    # 4. 圧力差
    t_pcc, v_pcc = component_series(prof, "cc", :P)
    ttl = isempty(v_pcc) ? "pressure difference" : "P - P_cc [Pa]  (P_cc(end) = $(round(v_pcc[end], digits = 0)) Pa)"
    p4 = plot_component_series(prof, :P; comps = ("grv", "vl", "cond", "ll"), agg = agg,
                               ylabel = "P - P_cc [Pa]", title = ttl, rel_to_cc = true)

    # 5. 流速
    p5 = plot_component_series(prof, :u; comps = PLOT_FLUID_COMPS, agg = agg,
                               ylabel = "u [m/s]", title = "velocity ($(agg))")

    # 6. 印加熱と蒸発による熱
    p6 = plot(; xlabel = "t [s]", ylabel = "Q [W]", title = "heat input / evaporation", legend = :best)
    hasproperty(summ, :Q_apply) &&
        plot!(p6, summ.t, summ.Q_apply; label = "Q_apply", color = :black, lw = 1.8)
    hasproperty(summ, :Q_evap) &&
        plot!(p6, summ.t, summ.Q_evap; label = "mdot * h_fg", color = :red, lw = 1.8)

    return plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = size, margin = 4Plots.mm)
end

# ---------------------------------------------------------------------------
# ループに沿った位置座標
# ---------------------------------------------------------------------------

"""
    loop_lengths(x) -> NamedTuple (grv, vl, cond, ll, cc)

各成分の長さ [m]。x は、(grv=..., vl=..., ...) の NamedTuple (数値)、または
LHPGeometry (geo.grv.L など)。
"""
function loop_lengths(x)
    if x isa NamedTuple && x.grv isa Real
        return (grv = Float64(x.grv), vl = Float64(x.vl), cond = Float64(x.cond),
                ll = Float64(x.ll), cc = Float64(x.cc))
    end
    return (grv = Float64(hasproperty(x.grv, :L_grv) ? x.grv.L_grv : x.grv.L), vl = Float64(x.vl.L), cond = Float64(x.cond.L),
            ll = Float64(x.ll.L), cc = Float64(x.cc.L))
end

"""
    loop_layout(prof, lengths) -> Vector{NamedTuple}

流体の成分 (グルーブ → 蒸気管 → 凝縮器 → 液管 → CC) の、ループに沿った位置。
各要素: (comp, N, dl, s_start, s_end)。セル数 N は prof から数える。s の原点は、グルーブの入口。
"""
function loop_layout(prof::DataFrame, lengths)
    Ls = loop_lengths(lengths)
    out = NamedTuple[]
    s0 = 0.0
    for name in PLOT_FLUID_COMPS
        L = getfield(Ls, Symbol(name))
        cells = prof.cell[prof.comp .== name]
        isempty(cells) && continue
        N = maximum(cells)
        push!(out, (comp = name, N = N, dl = L / N, s_start = s0, s_end = s0 + L))
        s0 += L
    end
    return out
end

# ---------------------------------------------------------------------------
# 任意の時刻の分布
# ---------------------------------------------------------------------------

# 保存されている時刻のうち、t に最も近いもの
function _nearest_times(avail::AbstractVector{<:Real}, ts)
    return [avail[argmin(abs.(avail .- t))] for t in ts]
end

# 時刻 tk の、成分ごとの (s, T, x, T_wall, u, s_u) を取り出す。
# 速度 u は、各セルの下流側の面の速度なので、面の位置 s_u = s_start + i Δl に置く。
function _distribution_at(prof::DataFrame, layout, tk::Real)
    prof_t = prof[prof.t .== tk, :]
    data = Dict{String,NamedTuple}()
    for c in layout
        sub = prof_t[prof_t.comp .== c.comp, :]
        isempty(sub) && continue
        sub = sort(sub, :cell)
        s = c.s_start .+ (sub.cell .- 0.5) .* c.dl
        s_u = c.s_start .+ sub.cell .* c.dl
        data[c.comp] = (s = collect(s), T = collect(sub.T), x = collect(sub.x), T_wall = collect(sub.T_wall),
                        u = collect(sub.u), s_u = collect(s_u))
    end
    return data
end

# 上端に、成分名のラベルと境界線を入れる (y_top: ラベルの y 位置)
function _decorate_components!(plt, layout, y_top)
    bounds = [c.s_end for c in layout[1:end-1]]
    isempty(bounds) || vline!(plt, bounds; color = :gray, ls = :dot, label = "")
    for c in layout
        annotate!(plt, 0.5 * (c.s_start + c.s_end), y_top,
                  Plots.text(get(PLOT_LABELS, c.comp, c.comp), 8, :gray))
    end
    return plt
end

"""
    plot_distribution(prof, wick, lengths, times; summary=nothing, color_by=:auto,
                      show_wall=true, show_velocity=true, size=(1050, 1250)) -> Plot

任意の時刻 (保存された時刻のうち、最も近いもの) の、LHP 全体の分布を描く。
  上: 温度 T (流体は実線、壁は破線。蒸発器ケース T_ec と、蒸発界面温度 T_out の点)
  2: 乾き度 x (vapor quality)
  3: 流速 u (各セルの下流側の面の速度。面の位置に描く。CC は、ループに沿った速度では
     なく、ウィックへの流出速度なので、含めない。show_velocity=false で、このパネルを省く)
  下: ウィックの温度 (0 = CC の流体、1..N = ウィックのセル、N+1 = 蒸発界面)
lengths: 各成分の長さ (grv, vl, cond, ll, cc) [m] の NamedTuple (load_lengths の戻り値、または手で指定)。
times は、スカラーまたはベクトル (重ね描き)。
color_by: :component (成分ごとに色。時刻が 1 つのとき既定) / :time (時刻ごとに色。複数時刻の既定)。
summary を渡すと、蒸発器ケースの温度と、蒸発界面の温度を描く。
"""
function plot_distribution(prof::DataFrame, wick::DataFrame, lengths, times;
    summary::Union{Nothing,DataFrame}=nothing, color_by::Symbol=:auto,
    show_wall::Bool=true, show_velocity::Bool=true, size=(1050, 1250))

    ts = times isa Real ? [Float64(times)] : collect(Float64, times)
    avail = sort(unique(prof.t))
    tk = _nearest_times(avail, ts)
    for (t, t1) in zip(ts, tk)
        abs(t - t1) > 1e-9 && println("  t = ", t, " s → 最も近い保存時刻 ", t1, " s を使います")
    end
    loop_pos = loop_layout(prof, lengths)
    nt = length(tk)
    if color_by === :auto
        color_by = nt == 1 ? :component : :time
    elseif color_by === :component && nt > 1
        println("  複数時刻のときは color_by=:time にします")
        color_by = :time
    end
    pal = nt > 1 ? palette(:viridis, nt) : nothing
    datas = [_distribution_at(prof, loop_pos, t) for t in tk]

    pT = plot(; xlabel = "position along loop s [m]", ylabel = "T [K]",
              title = "temperature distribution", legend = :outertopright)
    pX = plot(; xlabel = "position along loop s [m]", ylabel = "x [-]",
              title = "vapor quality distribution", legend = :outertopright)
    pU = plot(; xlabel = "position along loop s [m]", ylabel = "u [m/s]",
              title = "velocity distribution (downstream face of each cell)", legend = :outertopright)
    u_min = Inf; u_max = -Inf

    T_min = Inf; T_max = -Inf
    for (m, data) in enumerate(datas)
        tlabel = "t = $(tk[m]) s"
        if color_by === :component
            for c in loop_pos
                haskey(data, c.comp) || continue
                d = data[c.comp]
                col = get(PLOT_COMP_COLOR, c.comp, :gray)
                lab = get(PLOT_LABELS, c.comp, c.comp)
                plot!(pT, d.s, d.T; color = col, label = lab, lw = 2, marker = :circle, ms = 2.5)
                if show_wall && c.comp != "grv" && !all(isnan, d.T_wall)
                    plot!(pT, d.s, d.T_wall; color = col, ls = :dash, lw = 1.5, label = lab * " wall")
                end
                plot!(pX, d.s, d.x; color = col, label = lab, lw = 2, marker = :circle, ms = 2.5)
                if c.comp != "cc"
                    plot!(pU, d.s_u, d.u; color = col, label = lab, lw = 2, marker = :circle, ms = 2.5)
                end
            end
        else
            col = pal[m]
            s_all = Float64[]; T_all = Float64[]; x_all = Float64[]
            sw = Float64[]; Tw = Float64[]
            su_all = Float64[]; u_all = Float64[]
            for c in loop_pos
                haskey(data, c.comp) || continue
                d = data[c.comp]
                append!(s_all, d.s); append!(T_all, d.T); append!(x_all, d.x)
                if c.comp != "grv"
                    append!(sw, d.s); append!(Tw, d.T_wall)
                end
                if c.comp != "cc"
                    append!(su_all, d.s_u); append!(u_all, d.u)
                end
            end
            plot!(pT, s_all, T_all; color = col, label = tlabel, lw = 2)
            show_wall && plot!(pT, sw, Tw; color = col, ls = :dash, lw = 1.3, label = "")
            plot!(pX, s_all, x_all; color = col, label = tlabel, lw = 2)
            plot!(pU, su_all, u_all; color = col, label = tlabel, lw = 2)
        end
        for (cname, d) in data
            T_min = min(T_min, minimum(d.T)); T_max = max(T_max, maximum(d.T))
            if cname != "cc"
                u_min = min(u_min, minimum(d.u)); u_max = max(u_max, maximum(d.u))
            end
            if show_wall && !all(isnan, d.T_wall)
                T_min = min(T_min, minimum(filter(!isnan, d.T_wall)))
                T_max = max(T_max, maximum(filter(!isnan, d.T_wall)))
            end
        end

        # 蒸発器ケースの温度と、蒸発界面の温度 (グルーブの位置)
        if summary !== nothing
            k = argmin(abs.(summary.t .- tk[m]))
            s_g = 0.5 * (loop_pos[1].s_start + loop_pos[1].s_end)
            colm = color_by === :component ? PLOT_COMP_COLOR["evap"] : pal[m]
            if hasproperty(summary, :T_ec)
                scatter!(pT, [s_g], [summary.T_ec[k]]; marker = :diamond, ms = 7, color = colm,
                         label = m == 1 ? PLOT_LABELS["evap"] : "")
                T_min = min(T_min, summary.T_ec[k]); T_max = max(T_max, summary.T_ec[k])
            end
            if hasproperty(summary, :T_out)
                scatter!(pT, [s_g], [summary.T_out[k]]; marker = :star5, ms = 7, color = colm,
                         label = m == 1 ? "interface T_sat(P_grv)" : "")
            end
        end
    end

    # 範囲と、成分名のラベル・境界線
    if isfinite(T_min) && isfinite(T_max)
        rng = max(T_max - T_min, 1.0)
        plot!(pT; ylims = (T_min - 0.05 * rng, T_max + 0.14 * rng))
        _decorate_components!(pT, loop_pos, T_max + 0.09 * rng)
    end
    plot!(pX; ylims = (-0.05, 1.2))
    _decorate_components!(pX, loop_pos, 1.12)
    if isfinite(u_min) && isfinite(u_max)
        urng = max(u_max - u_min, 1e-3)
        plot!(pU; ylims = (u_min - 0.05 * urng, u_max + 0.14 * urng))
        hline!(pU, [0.0]; color = :gray, ls = :dot, label = "")
        _decorate_components!(pU, loop_pos, u_max + 0.09 * urng)
    end

    # ウィック: 0 = CC の流体 (平板型は secondary wick)、1..N = ウィックのセル、N+1 = 蒸発界面
    pW = plot(; xlabel = "wick node (0: CC fluid / secondary wick, 1..N: cells from CC side, N+1: interface)",
              ylabel = "T [K]", title = "wick temperature", legend = :outertopright)
    for (m, t1) in enumerate(tk)
        sub = sort(wick[wick.t .== t1, :], :cell)
        isempty(sub) && continue
        col = color_by === :component ? :brown : pal[m]
        xs = collect(Float64, sub.cell)
        Ts = collect(Float64, sub.T)
        Tcc = NaN
        cc = prof[(prof.t .== t1) .& (prof.comp .== "cc"), :]
        isempty(cc) || (Tcc = cc.T[1])
        # 平板型は、cell = 0 に secondary wick が入っている。円筒型は、CC の流体温度を 0 とする
        if !isnan(Tcc) && !(0 in sub.cell)
            xs = vcat(0.0, xs); Ts = vcat(Tcc, Ts)
        end
        if summary !== nothing && hasproperty(summary, :T_out)
            k = argmin(abs.(summary.t .- t1))
            xs = vcat(xs, maximum(sub.cell) + 1.0); Ts = vcat(Ts, summary.T_out[k])
        end
        plot!(pW, xs, Ts; color = col, marker = :circle, ms = 4, lw = 1.8, label = "t = $(t1) s")
    end

    if show_velocity
        return plot(pT, pX, pU, pW; layout = @layout([a; b; c; d{0.16h}]), size = size, margin = 4Plots.mm)
    end
    return plot(pT, pX, pW; layout = @layout([a; b; c{0.22h}]), size = (size[1], round(Int, size[2] * 0.76)), margin = 4Plots.mm)
end