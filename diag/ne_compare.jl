# ne_compare.jl — CC の気相の過剰圧力 ΔP_ne を液側へ伝える場合 / 伝えない場合の差分 (起動 τ = 0〜90 s)
#
# 【何を切るか】 CC の圧力 P (気相) は、界面の液温 Tl[1] の飽和圧より ΔP_ne = P − P_sat(Tl[1]) だけ高くなる (液が流入して蒸気を押すとき)。
#   ふつうは液側 (primary wick の出口、secondary wick、液管の出口) も同じ P を見る (自由界面なので液圧 = 蒸気圧)。
#   「伝えない」ケースでは、液側にだけ P_sat(Tl[1]) を渡す (ne_liq_decouple = true)。気相側 (CC の質量・エネルギー収支) は変えない。
#   → 液側の流れ (液管の押し出し・ウィックへの流入) に ΔP_ne がどれだけ効いているかが分かる。物理的に正しい設定ではなく、切り分け用。
# 【実行】 julia ne_compare.jl   (設定はこのファイルの冒頭。1 本 約 3 分)
#   「伝える」ケースは保存済みの単独実行 (res/transient_flat_ammonia_v9a) を使う。同じ設定 (r_c_wick=1.85 µm、h_con=20000 など) で :short を回した後に実行すること。
# 【出力】 表 (センサ・圧力・流量の差) と図 res/ne_compare/ne_compare.png、csv
nc_RUN      = true                                   # false: 保存済みの「伝えない」結果だけで表を作る
nc_t_end    = 100.0                                  # 計算の終わり [s]
nc_tag      = "_neoff"
nc_base_dir = "transient_flat_ammonia_v9a"           # 「伝える」ケースの結果フォルダ (res/ の下)
nc_t_heat   = 10.0
nc_t_exp    = 303.8
nc_tau      = [15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 90.0]

using CSV, DataFrames, Statistics
import Plots

if nc_RUN
    global SWEEP = (t_end = nc_t_end, tag = nc_tag, ne_off = 1.0)
    include(joinpath(@__DIR__, "transient_flat.jl"))
end

K = 273.15
lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
function load(dir)
    d = joinpath(@__DIR__, "..", "res", dir)
    isfile(joinpath(d, "summary.csv")) || error("結果がありません: " * d)
    summ = CSV.read(joinpath(d, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])
    idx = Dict{Tuple{String,Float64},Vector{Int}}()
    for i in 1:nrow(prof)
        push!(get!(idx, (String(prof.comp[i]), Float64(prof.t[i])), Int[]), i)
    end
    return (summ = summ, prof = prof, idx = idx, tp = sort(unique(prof.t)))
end
function cells(r, comp, τ)
    t = r.tp[argmin(abs.(r.tp .- (nc_t_heat + τ)))]
    haskey(r.idx, (comp, Float64(t))) || return nothing
    return sort(r.prof[r.idx[(comp, Float64(t))], :], :cell)
end
wallmean(r, comp, τ) = (d = cells(r, comp, τ); d === nothing ? NaN : mean(d.T_wall) - K)
function wallpos(r, comp, τ, pos)
    d = cells(r, comp, τ); d === nothing && return NaN
    N = nrow(d)
    return lin([(i - 0.5) / N for i in 1:N], Vector(d.T_wall), pos) - K
end
scol(r, col, τ) = (col in propertynames(r.summ)) ? lin(r.summ.t, Vector(r.summ[!, col]), nc_t_heat + τ) : NaN

on  = load(nc_base_dir)
off = load("transient_flat_ammonia" * nc_tag)

sens = [("HB", (r, τ) -> wallmean(r, "hb", τ)), ("EC", (r, τ) -> (v = wallmean(r, "eca", τ); isnan(v) ? scol(r, :T_ec, τ) - K : v)),
        ("RES", (r, τ) -> wallmean(r, "cc", τ)), ("VL1", (r, τ) -> wallpos(r, "vl", τ, 1/6)), ("VL3", (r, τ) -> wallpos(r, "vl", τ, 5/6)),
        ("CON1", (r, τ) -> wallpos(r, "cond", τ, 0.0)), ("LL1", (r, τ) -> wallpos(r, "ll", τ, 0.10/1.15)), ("LL2", (r, τ) -> wallpos(r, "ll", τ, 0.5/1.15))]
println("\n=== 差 (伝えない − 伝える) ===  壁温 [K]")
println(rpad("τ [s]", 8), join([lpad(s[1], 8) for s in sens]))
for τ in nc_tau
    println(rpad(string(τ), 8), join([lpad(string(round(s[2](off, τ) - s[2](on, τ), digits = 3)), 8) for s in sens]))
end

qs = [(:P_cc, "P_cc [kPa]", 1e-3), (:dP_ne, "ΔP_ne [kPa]", 1e-3), (:mdot, "ṁ_evap [mg/s]", 1e6), (:u_ll_end, "u_ll_end [mm/s]", 1e3),
      (:Q_hb_ec, "Q_hb_ec [W]", 1.0), (:n_nuc, "n", 1.0), (:T_ec, "T_ec [K]", 1.0)]
println("\n圧力・流量など:  伝える | 伝えない | 差")
for (col, lab, sc) in qs
    (col in propertynames(on.summ)) || continue
    println("\n", lab)
    println(rpad("τ [s]", 8), join([lpad(string(Int(τ)), 9) for τ in nc_tau]))
    a = [scol(on, col, τ) * sc for τ in nc_tau]; b = [scol(off, col, τ) * sc for τ in nc_tau]
    println(rpad("伝える", 8), join([lpad(string(round(v, digits = 3)), 9) for v in a]))
    println(rpad("伝えない", 8), join([lpad(string(round(v, digits = 3)), 9) for v in b]))
    println(rpad("差", 8), join([lpad(string(round(y - x, digits = 3)), 9) for (x, y) in zip(a, b)]))
end

# n = 0.5 の時刻
function t_half(r)
    t = r.summ.t; n = r.summ.n_nuc
    k = findfirst(>=(0.5), n)
    k === nothing && return NaN
    return k == 1 ? t[1] - nc_t_heat : (t[k-1] + (0.5 - n[k-1]) / (n[k] - n[k-1]) * (t[k] - t[k-1])) - nc_t_heat
end
println("\nn = 0.5 の時刻 τ [s]:  伝える ", round(t_half(on), digits = 2), "   伝えない ", round(t_half(off), digits = 2))

out = joinpath(@__DIR__, "..", "res", "ne_compare"); mkpath(out)
τc = collect(10.0:0.5:nc_t_end - nc_t_heat)
rows = DataFrame(tau = τc)
for (col, lab, sc) in qs
    (col in propertynames(on.summ)) || continue
    rows[!, Symbol(string(col) * "_on")]  = [scol(on, col, τ) * sc for τ in τc]
    rows[!, Symbol(string(col) * "_off")] = [scol(off, col, τ) * sc for τ in τc]
end
for (nm, f) in sens
    rows[!, Symbol(nm * "_on")]  = [f(on, τ) for τ in τc]
    rows[!, Symbol(nm * "_off")] = [f(off, τ) for τ in τc]
end
CSV.write(joinpath(out, "ne_compare.csv"), rows)
try
    pl = []
    for (col, lab, sc) in qs[1:4]
        (col in propertynames(on.summ)) || continue
        p1 = Plots.plot(τc, rows[!, Symbol(string(col) * "_on")], label = "ΔP_ne を伝える", lw = 2, xlabel = "τ [s]", ylabel = lab, xlims = (10, 70))
        Plots.plot!(p1, τc, rows[!, Symbol(string(col) * "_off")], label = "伝えない", lw = 2, ls = :dash)
        push!(pl, p1)
    end
    for nm in ("VL1", "LL1", "CON1", "HB")
        p1 = Plots.plot(τc, rows[!, Symbol(nm * "_on")], label = "伝える", lw = 2, xlabel = "τ [s]", ylabel = nm * " 壁温 [°C]", xlims = (10, 90))
        Plots.plot!(p1, τc, rows[!, Symbol(nm * "_off")], label = "伝えない", lw = 2, ls = :dash)
        push!(pl, p1)
    end
    Plots.savefig(Plots.plot(pl..., layout = (4, 2), size = (1100, 1300)), joinpath(out, "ne_compare.png"))
    println("図を保存しました: ", joinpath(out, "ne_compare.png"))
catch err
    @warn "図の作成に失敗 (表と csv は保存済み)" err
end
println("保存しました: ", out)