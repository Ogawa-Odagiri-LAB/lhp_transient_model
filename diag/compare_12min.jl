using Pkg; Pkg.activate(".")
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
using LHPTransientModel      # pchip_interp を使う
using CSV, DataFrames

# =============================================================================
# compare_12min.jl … 保存済みの計算結果 (stage = :long) と、実験 exp_5W_12min.csv を比べる (計算のやり直しは不要)。
#   使い方: julia compare_12min.jl      (この .jl は diag/ に置く。実験の csv は data/ に置く)
#   読むもの: res/<result_name>/summary.csv, profiles.csv
#   出すもの: 表 (モデル / 実験 / 差) と、res/<result_name>/compare_12min.csv
# 注: RES は CC の壁 (profiles.csv の comp = "cc" の T_wall の平均) と比べる。新しい CC では T_cc が界面側の薄い第 1 層の温度なので、流体の代用にはしない。
# =============================================================================
result_name = "transient_flat_ammonia_v9a"   # 800 s の結果は "transient_flat_ammonia" + run_tag、2 h は "transient_flat_ammonia_2h" + run_tag
long2h = occursin("_2h", result_name)
result_dir = joinpath(@__DIR__, "..", "res", result_name)
exp_name = long2h ? "exp_5W_2h.csv" : "exp_5W_12min.csv"
exp_file = joinpath(@__DIR__, "..", "data", exp_name)
isfile(exp_file) || error(exp_name, " が見つかりません: ", exp_file)
t_heat, t_exp_onset = 10.0, 303.8        # モデルの加熱開始時刻 (stage = :long) と、実験の対応する時刻
K = 273.15

summ = CSV.read(joinpath(result_dir, "summary.csv"), DataFrame)
prof = CSV.read(joinpath(result_dir, "profiles.csv"), DataFrame)
ex   = CSV.read(exp_file, DataFrame)

cellpos(N) = [(i - 0.5) / N for i in 1:N]
s_LL1, s_LL2 = 0.10 / 1.15, (1.15 - 0.50) / 1.15
tmodel = sort(unique(prof.t))
function wall(comp, t, pos)
    d = prof[(prof.comp .== comp) .& (prof.t .== t), :]
    d = sort(d, :cell)
    return pchip_interp(cellpos(nrow(d)), Vector(d.T_wall), pos) - K
end
nearest(t) = tmodel[argmin(abs.(tmodel .- t))]
# 蒸発器ケースを 2 節点にした結果 (comp = "eca" の行) では、EC1・EC2 のセンサは側面の節点 S と比べる (summary の T_ec は受熱面の節点 P)
has_eca = any((prof.comp .== "eca") .& (prof.T_wall .> 100))
eca_T(t) = only(prof[(prof.comp .== "eca") .& (prof.t .== t), :T_wall]) - K
srow(t) = summ[argmin(abs.(summ.t .- t)), :]
lin1(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
ecol(nm) = Vector(ex[!, nm])
exv(cols, τ) = sum(lin1(ex.t, ecol(c), t_exp_onset + τ) for c in cols) / length(cols)

sens = [
    ("HB",   ["HB"],           t -> srow(t).T_hb - K),
    ("EC",   ["EC1", "EC2"],   t -> has_eca ? eca_T(t) : srow(t).T_ec - K),
    ("RES",  ["RES1", "RES2"], t -> sum(prof[(prof.comp .== "cc") .& (prof.t .== t), :T_wall]) / count((prof.comp .== "cc") .& (prof.t .== t)) - K),
    ("VL1",  ["VL1"],  t -> wall("vl", t, 1/6)), ("VL2", ["VL2"], t -> wall("vl", t, 1/2)), ("VL3", ["VL3"], t -> wall("vl", t, 5/6)),
    ("CON1", ["CON1"], t -> wall("cond", t, 0.0)), ("CON2", ["CON2"], t -> wall("cond", t, 0.25)),
    ("CON3", ["CON3"], t -> wall("cond", t, 0.5)), ("CON4", ["CON4"], t -> wall("cond", t, 0.75)),
    ("CON5", ["CON5"], t -> wall("cond", t, 1.0)),
    ("LL1",  ["LL1"],  t -> wall("ll", t, s_LL1)), ("LL2", ["LL2"], t -> wall("ll", t, s_LL2)),
    ("LL3", ["LL3"], t -> wall("ll", t, (1.15 - 0.04) / 1.15)),   # 診断用: CC 手前 4 cm (熱電対が流体の温度を表すかは未確認。モデルの末端の壁と並べて見る)
]
τs = [τ for τ in (long2h ? (0, 90, 300, 600, 800, 1200, 1800, 2400, 3000, 3600, 4200, 4800, 5400, 6000, 6600, 6800) : (0, 30, 60, 90, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 800)) if t_heat + τ <= tmodel[end] && t_exp_onset + τ <= ex.t[end]]

println("実験との比較 (stage = :long、加熱開始からの時間 τ [s]):  モデル / 実験 / 差 (モデル − 実験) [°C or K]")
println(rpad("", 9), join([lpad("+" * string(τ), 7) for τ in τs]))
rows = DataFrame(sensor = String[], tau = Float64[], model = Float64[], exp = Float64[], diff = Float64[])
for (nm, cols, f) in sens
    mo = [f(nearest(t_heat + τ)) for τ in τs]
    e  = [exv(cols, τ) for τ in τs]
    println(rpad(nm * " 計算", 9), join([lpad(string(round(v, digits = 2)), 7) for v in mo]))
    println(rpad(nm * " 実験", 9), join([lpad(string(round(v, digits = 2)), 7) for v in e]))
    println(rpad(nm * " 差",   9), join([lpad(string(round(v, digits = 2)), 7) for v in (mo .- e)]))
    for (τ, a, b) in zip(τs, mo, e)
        push!(rows, (nm, τ, a, b, a - b))
    end
end
CSV.write(joinpath(result_dir, "compare_12min.csv"), rows)

# 実験の導通 (印加熱 / (Tsat − CON1) は圧力が無いので省略)。代わりに HB−CON2 を、時間で並べる (経路全体の温度差の目安)
println("\nHB − CON2 [K] (モデル / 実験):")
for τ in τs
    t = nearest(t_heat + τ)
    println("  +", lpad(τ, 4), "  ", lpad(round(srow(t).T_hb - K - wall("cond", t, 0.25), digits = 2), 6), "  ",
            lpad(round(exv(["HB"], τ) - exv(["CON2"], τ), digits = 2), 6))
end
# CON1 − CON2 [K]: 凝縮が起きている領域 (CON1 付近) の壁とプレート温度 (CON2 はほぼプレート温度) の差。圧力が要らない。
println("\nCON1 − CON2 [K] (モデル / 実験) と、モデルの凝縮領域の長さ (x > 0.01 の最後のセル / N、管全長に対する比):")
function cond_end(t)
    d = sort(prof[(prof.comp .== "cond") .& (prof.t .== t), :], :cell)
    k = findlast(>(0.01), d.x)
    return k === nothing ? 0.0 : k / nrow(d)
end
for τ in τs
    t = nearest(t_heat + τ)
    println("  +", lpad(τ, 4), "  ", lpad(round(wall("cond", t, 0.0) - wall("cond", t, 0.25), digits = 2), 6), "  ",
            lpad(round(exv(["CON1"], τ) - exv(["CON2"], τ), digits = 2), 6), "    凝縮領域の長さ (モデル): ", round(cond_end(t), digits = 3))
end

println("\n保存: ", joinpath(result_dir, "compare_12min.csv"))