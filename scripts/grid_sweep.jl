using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# grid_sweep.jl … 凝縮器プレートの格子 (nx × ny。コールドプレートも同じ格子) を変えて :short を続けて計算し、凝縮器の入口付近の温度差が格子で変わるかを見る
#   使い方: julia grid_sweep.jl      (scripts/ に置く。transient_flat.jl・heat_budget_calc.jl と同じフォルダ。実験の csv は data/ に置く)
#   しくみ: f_shah_sweep.jl・tamb_sweep.jl と同じ。Main.SWEEP に plate_nx, plate_ny を入れて transient_flat.jl を 1 回ずつ include する。
#   見るもの: Tsat − CON1、CON1 − プレート中央、プレート中央 − コールドプレート中央、プレート最高 − 中央 (どれも τ = 60〜290 s)。実験の RES − CON1、
#             CON1 − 117 (ps 用のプレートの実験ファイルがあれば) と比べる。格子で動かなければ、プレートの節点の大きさは問題にならない。
#   格子の選び方: 節点の幅は 300/nx × 60/ny mm。管の 2 本の線は y = 15, 45 mm: ny = 6 (10 mm) と ny = 10 (6 mm) は節点の中心に管がのる (ny = 12 は 2.5 mm ずれる)。
#                 管のセルは 11.2 mm 間隔なので、節点の幅が 11.2 mm より小さいと、管の接しない節点が出る (熱はプレートの中の熱伝導で広がる)。
#   注: :short の t_end = gs_t_end (既定 300 s)。プレートの熱の広がりは数十 s で準定常になるので、起動の 100〜290 s で比べられる (2 時間のゆっくりした上昇はこの確認の対象外)。
#       計算時間は節点が増えると増える (50 × 10 で約 2 倍と見込むが未測定)。gs_SKIP_EXISTING = true で保存済みのケースを飛ばす。gs_RUN = false なら比較だけ。
# =============================================================================

# ---------------------------- 設定 ----------------------------
gs_t_end        = 300.0         # :short の終わり [s] (加熱開始は 10 s。τ = 290 s まで)
gs_t_heat       = 10.0          # transient_flat.jl の t_heat と同じ
gs_t_exp_onset  = 303.8         # 実験の加熱開始時刻 [s]
gs_tag_prefix   = "_psGrid_"
gs_out_name     = "grid_sweep"
gs_cases = [("30x6", 30, 6),    # 基準 (10 × 10 mm)
            ("30x10", 30, 10),  # 幅方向だけ細かく (10 × 6 mm)。管の線が節点の中心にのる
            ("50x10", 50, 10)]  # 6 × 6 mm
gs_exp_file      = "exp_5W_12min.csv"           # data/ (τ ≤ 800 s の実験。RES, CON1 など)
gs_exp_plate     = "exp_5W_plate.csv"           # data/ (タブ区切り: t, 117_cond, 118_coldplate。なければ飛ばす)
gs_report_tau   = [60.0, 120.0, 200.0, 290.0]
gs_RUN = true
gs_SKIP_EXISTING = false
# ----------------------------------------------------------------

using CSV, DataFrames, Statistics

gs_tag(c) = gs_tag_prefix * c[1]
gs_dir(c) = joinpath(@__DIR__, "..", "res", "transient_flat_ammonia" * gs_tag(c))

println("実行するケース (", length(gs_cases), " 本): ", join([c[1] for c in gs_cases], ", "))
if gs_RUN
    for c in gs_cases
        if gs_SKIP_EXISTING && isfile(joinpath(gs_dir(c), "summary.csv")) && isfile(joinpath(gs_dir(c), "profiles.csv"))
            println("\n(飛ばす: 保存済み) ", c[1])
            continue
        end
        println("\n", "="^70, "\n  格子 ", c[1], " の計算を開始\n", "="^70)
        global SWEEP = (stage = :short, t_end = gs_t_end, tag = gs_tag(c), plate_nx = c[2], plate_ny = c[3])
        try
            include(joinpath(@__DIR__, "transient_flat.jl"))
        catch err
            @warn "格子 $(c[1]) の実行が途中で失敗 (保存済みの結果があれば比較に使う)" err
        end
    end
end

# ---------------------------- 読み込み ----------------------------
K = 273.15
lin(xs, ys, x) = x <= xs[1] ? ys[1] : (x >= xs[end] ? ys[end] : (k = searchsortedlast(xs, x); ys[k] + (ys[k+1] - ys[k]) * (x - xs[k]) / (xs[k+1] - xs[k])))
res = Dict{String,Any}()
for c in gs_cases
    d = gs_dir(c)
    if !isfile(joinpath(d, "summary.csv")) || !isfile(joinpath(d, "profiles.csv"))
        @warn "結果がありません: $(d)"
        continue
    end
    summ = CSV.read(joinpath(d, "summary.csv"), DataFrame)
    prof = CSV.read(joinpath(d, "profiles.csv"), DataFrame, select = [:t, :comp, :cell, :T_wall])
    res[c[1]] = (summ = summ, prof = prof, nx = c[2], ny = c[3], tp = sort(unique(prof.t)))
end
isempty(res) && error("結果がありません。gs_RUN = true で計算してください")

# プレート中央の 4 節点 (loop_heat_diag と同じ定義)
function centre_nodes(nx, ny)
    ic = (cld(nx, 2), cld(nx, 2) + 1 - (nx % 2))
    jc = (cld(ny, 2), cld(ny, 2) + 1 - (ny % 2))
    return [(j - 1) * nx + i for i in ic for j in jc]
end

function metrics(r, τ)
    t = r.tp[argmin(abs.(r.tp .- (gs_t_heat + τ)))]
    snap = r.prof[r.prof.t .== t, :]
    wall(comp) = (d = sort(snap[snap.comp .== comp, :], :cell); Vector(d.T_wall))
    plate, cplate, cond, cc = wall("plate"), wall("cplate"), wall("cond"), wall("cc")
    cn = centre_nodes(r.nx, r.ny)
    Tsat = lin(r.summ.t, Vector(r.summ.T_cc), t) - K
    con1 = cond[1] - K
    pc = mean(plate[cn]) - K
    cpc = isempty(cplate) ? NaN : mean(cplate[cn]) - K
    return (Tsat = Tsat, CON1 = con1, Tpl_c = pc, Tcp_c = cpc, Tpl_max = maximum(plate) - K, RES = mean(cc) - K)
end

# 実験
ex = CSV.read(joinpath(@__DIR__, "..", "data", gs_exp_file), DataFrame)
exv(cols, τ) = mean(lin(ex.t, Vector(ex[!, cn]), gs_t_exp_onset + τ) for cn in cols)
exp_plate = nothing
let f = joinpath(@__DIR__, "..", "data", gs_exp_plate)
    isfile(f) && (global exp_plate = CSV.read(f, DataFrame; delim = '\t'))
end
exp_plate === nothing && println("(", gs_exp_plate, " がないので、実験の 117・118 は表に出しません)")
fmt(x, d = 3) = lpad(isnan(x) ? "-" : string(round(x, digits = d)), 9)

tmax = minimum(r.summ.t[end] - gs_t_heat for r in values(res))
taus = [τ for τ in gs_report_tau if τ <= tmax]

rows = DataFrame(case = String[], tau = Float64[], Tsat_CON1 = Float64[], CON1_Tplc = Float64[], Tplc_Tcpc = Float64[], Tplmax_Tplc = Float64[], RES_CON1 = Float64[])
for τ in taus
    println("\nτ = ", τ, " s:   (差はすべて [K])")
    println(rpad("格子", 10), lpad("Tsat-CON1", 10), lpad("CON1-Tpl中央", 13), lpad("Tpl中央-Tcp中央", 16), lpad("Tpl最高-中央", 13), lpad("RES-CON1", 10), lpad("CON1 [°C]", 10), lpad("Tpl中央", 9))
    for c in gs_cases
        haskey(res, c[1]) || continue
        m = metrics(res[c[1]], τ)
        println(rpad(c[1], 10), fmt(m.Tsat - m.CON1), "   ", fmt(m.CON1 - m.Tpl_c), "      ", fmt(m.Tpl_c - m.Tcp_c), "       ", fmt(m.Tpl_max - m.Tpl_c), fmt(m.RES - m.CON1), fmt(m.CON1, 2), fmt(m.Tpl_c, 2))
        push!(rows, (c[1], τ, m.Tsat - m.CON1, m.CON1 - m.Tpl_c, m.Tpl_c - m.Tcp_c, m.Tpl_max - m.Tpl_c, m.RES - m.CON1))
    end
    resv = exv(["RES1", "RES2"], τ); c1 = exv(["CON1"], τ)
    line = "実験       " * lpad("-", 10) * "   " * lpad("-", 9) * "      " * lpad("-", 9) * "       " * lpad("-", 9) * fmt(resv - c1) * fmt(c1, 2)
    if exp_plate !== nothing
        p117 = lin(exp_plate.t, Vector(exp_plate[!, "117_cond"]), gs_t_exp_onset + τ)
        p118 = lin(exp_plate.t, Vector(exp_plate[!, "118_coldplate"]), gs_t_exp_onset + τ)
        line *= fmt(p117, 2) * "   (CON1−117 = " * string(round(c1 - p117, digits = 3)) * ",  117−118 = " * string(round(p117 - p118, digits = 3)) * ")"
    end
    println(line)
end

# 他のセンサの実験との差 (格子の影響が HB・EC・RES・VL1 に広がっていないか)
println("\n実験との差 (モデル − 実験) [K]:   HB | EC | RES | VL1 | CON1")
function wallpos(r, comp, t, pos)
    snap = r.prof[(r.prof.t .== t) .& (r.prof.comp .== comp), :]
    snap = sort(snap, :cell); N = nrow(snap)
    return lin([(i - 0.5) / N for i in 1:N], Vector(snap.T_wall), pos) - K
end
for τ in taus
    println("τ = ", τ, " s")
    for c in gs_cases
        haskey(res, c[1]) || continue
        r = res[c[1]]
        t = r.tp[argmin(abs.(r.tp .- (gs_t_heat + τ)))]
        hb = wallpos(r, "hb", t, 0.5); ec = wallpos(r, any((r.prof.t .== t) .& (r.prof.comp .== "eca")) ? "eca" : "evap", t, 0.5); rs = wallpos(r, "cc", t, 0.5); vl1 = wallpos(r, "vl", t, 1 / 6); c1 = wallpos(r, "cond", t, 0.0)
        println("  ", rpad(c[1], 10), fmt(hb - exv(["HB"], τ)), fmt(ec - exv(["EC1", "EC2"], τ)), fmt(rs - exv(["RES1", "RES2"], τ)), fmt(vl1 - exv(["VL1"], τ)), fmt(c1 - exv(["CON1"], τ)))
    end
end

out = joinpath(@__DIR__, "..", "res", gs_out_name)
mkpath(out)
CSV.write(joinpath(out, "grid_metrics.csv"), rows)
println("\n保存しました: ", abspath(out), " (grid_metrics.csv)")