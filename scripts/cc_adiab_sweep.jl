using Pkg; Pkg.activate(".") # プロジェクト環境のアクティベート

# =============================================================================
# cc_adiab_sweep.jl … CC ケースを完全断熱 (CC ケース → 断熱材 節点 3 の熱伝導を 0) にして :long (816 s) を計算し、基準と実験に比べる 【未実行】
#   目的: 長期の温度上昇 (実験 HB +2.5 K/2 h に対しモデルは τ≈90 s で飽和) に対し、CC ケースから周囲への熱損失が効いているかの上限試験。
#         CC は還流の過冷度で熱を吸収する (感度 約 50 K/W) ので、CC に熱が残れば T_cc が上がり、ループ全体の温度が上がる、という仮説の確認。
#   使い方: julia scripts/cc_adiab_sweep.jl      (scripts/ に置く。transient_flat.jl・heat_budget_calc.jl と同じフォルダ。実験の csv は data/exp_5W_12min.csv)
#   しくみ: T_amb_sweep.jl と同じ。Main.SWEEP を定義して transient_flat.jl を include する。cc_adiabatic = true (src/types.jl, src/residuals_housing.jl の追加) で節点 3 の壁側の伝導 Gi を 0 にする。
#          節点 3 には、断熱材の横方向の伝導 (蒸発器ケース側から)、液管の末端 2 セルの熱 (n_ll_ci)、外表面の損失は残る (熱容量 0 の節点は DAE が特異になるので、節点自体は消さない)。
#   ケース: base (cc_adiabatic = false。ps_base_ref が保存済みなら計算せずそれを使う) / cc_adiabatic (true)。
#           ps_extra_cases に追加ケース (例: f_ec_cc = 0 = 蒸発器ケース→CC ケースの伝導も切る) を足せる。
#   結果は res/transient_flat_ammonia<ps_tag_prefix><名前>/ に保存 (budget.csv も出る)。表・csv は res/<ps_out_name>/。 
#   注: src/ を変更したので Julia を再起動してから実行する。1 本 :long で約 1400〜1600 s。
# =============================================================================

# ---------------------------- 設定 ----------------------------
ps_stage       = :long      # :long (816 s) / :long2h (6910 s。差が 0.2 K 以上のときに)
ps_t_end       = 816.0      # 表・csv の範囲に使う終わり [s] (:long は 816 s。:long2h なら 6910.0)
ps_t_heat      = 10.0       # transient_flat.jl の t_heat と同じ
ps_t_exp_onset = 303.8      # 実験の加熱開始時刻 [s]
ps_tag_prefix  = "_psCCad_"          # 結果フォルダ名の接尾辞の頭 (前の結果を上書きしないため。:long2h にするときは "_psCCad2h_" などに変える)
ps_out_name    = "cc_adiab_sweep"    # 表・csv の保存先 (res/<名前>/)
ps_base_ref    = "transient_flat_ammonia_psTamb_base"   # base に使う保存済みの結果 (res/ の下)。"" なら base も計算する。中身が現在のコードの既定値 (T_amb 21.5 °C) と同じであること
ps_extra_cases = Any[
    # (name = "cc_adiab_fecc0", pars = (cc_adiabatic = true, f_ec_cc = 0.0)),   # 追加ケースの例 (外すと計算する)
]
ps_report_tau  = [30.0, 60.0, 120.0, 290.0, 500.0, 800.0]   # 表に出す τ [s] (ps_t_end − ps_t_heat 以下のもの)
ps_RUN = true
ps_SKIP_EXISTING = false
# ----------------------------------------------------------------

using CSV, DataFrames, Statistics

# ---------------------------- ケースの組み立て ----------------------------
ps_cases = Any[]
if ps_base_ref != "" && isfile(joinpath(@__DIR__, "..", "res", ps_base_ref, "summary.csv"))
    push!(ps_cases, (name = "base", pars = NamedTuple(), dir = ps_base_ref))
else
    push!(ps_cases, (name = "base", pars = (cc_adiabatic = false,), dir = nothing))
end
push!(ps_cases, (name = "cc_adiabatic", pars = (cc_adiabatic = true,), dir = nothing))
for c in ps_extra_cases
    push!(ps_cases, (name = c.name, pars = c.pars, dir = nothing))
end
ps_tag(c) = ps_tag_prefix * c.name
ps_dir(c) = c.dir === nothing ? joinpath(@__DIR__, "..", "res", "transient_flat_ammonia" * ps_tag(c)) : joinpath(@__DIR__, "..", "res", c.dir)
println("ケース (", length(ps_cases), " 本):")
for c in ps_cases
    println("  ", c.name, c.dir === nothing ? "" : "   (保存済み: res/" * c.dir * ")")
end

if ps_RUN
    for c in ps_cases
        c.dir === nothing || continue      # 保存済みの結果を使うケースは計算しない
        if ps_SKIP_EXISTING && isfile(joinpath(ps_dir(c), "summary.csv")) && isfile(joinpath(ps_dir(c), "profiles.csv"))
            println("\n(飛ばす: 保存済み) ", c.name)
            continue
        end
        println("\n", "="^70, "\n  ケース ", c.name, " の計算を開始\n", "="^70)
        global SWEEP = merge((stage = ps_stage, t_end = ps_t_end, tag = ps_tag(c)), c.pars)
        try
            include(joinpath(@__DIR__, "transient_flat.jl"))
        catch err
            @warn "ケース $(c.name) の実行が途中で失敗 (保存済みの結果があれば比較に使う)" err
        end
    end
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
ps_has(r, col) = hasproperty(r.summ, col)
ps_col(r, col, τ) = ps_has(r, col) ? ps_lin(r.summ.t, Vector{Float64}(r.summ[!, col]), ps_t_heat + τ) : NaN
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

# ---------------------------- 表 ----------------------------
println("\n基準 (base) の実験との差 (モデル − 実験) [K]:")
println(rpad("τ [s]", 7), join([lpad(s[1], 9) for s in ps_sensors]))
for τ in ps_tau_rep
    println(rpad(τ, 7), join([fmt(s[3](rb, τ) - ps_exv(s[2], τ)) for s in ps_sensors]))
end
for τ in ps_tau_rep
    println("\nτ = ", τ, " s の実験との差 (モデル − 実験) [K]:")
    println(rpad("ケース", 24), join([lpad(s[1], 9) for s in ps_sensors]))
    for c in ps_ok
        println(rpad(c.name, 24), join([fmt(s[3](ps_res[c.name], τ) - ps_exv(s[2], τ)) for s in ps_sensors]))
    end
end
for τ in ps_tau_rep
    println("\nτ = ", τ, " s の基準との差 (ケース − 基準) [K]:")
    println(rpad("ケース", 24), join([lpad(s[1], 9) for s in ps_sensors]))
    for c in ps_ok
        c.name == "base" && continue
        println(rpad(c.name, 24), join([fmt(s[3](ps_res[c.name], τ) - s[3](rb, τ)) for s in ps_sensors]))
    end
end

# 昇温 (τ = 60 s を基準にした変化) … 長期の傾きを見る。実験の HB は 2 h で約 +2.5 K
println("\nτ = ", first(ps_tau_rep), " s からの昇温 (τ − ", first(ps_tau_rep), " s) [K]。実験と比べる:")
println(rpad("ケース", 24), join([lpad(s[1], 9) for s in ps_sensors]))
for τ in ps_tau_rep[2:end]
    println("τ = ", τ)
    println(rpad("実験", 24), join([fmt(ps_exv(s[2], τ) - ps_exv(s[2], first(ps_tau_rep))) for s in ps_sensors]))
    for c in ps_ok
        r = ps_res[c.name]
        println(rpad(c.name, 24), join([fmt(s[3](r, τ) - s[3](r, first(ps_tau_rep))) for s in ps_sensors]))
    end
end

# CC まわりの熱 [W] (summary.csv にある列だけ): Q_hi_cc = CC ケース → 節点 3、Q_ec_cc = 蒸発器ケース → CC ケース、Q_br = 熱橋、Q_cc = CC の流体が受ける熱
ps_qcols = [:Q_hi_cc, :Q_ec_cc, :Q_br, :Q_cc, :Q_hi_ev, :Q_hi_hb, :Q_peek]
println("\nCC まわりの熱 [W] (+ は CC ケースから出る/入る向きは列名どおり。cc_adiabatic では Q_hi_cc = 0):")
for τ in ps_tau_rep
    println("τ = ", τ, " s")
    println(rpad("ケース", 24), join([lpad(string(q), 10) for q in ps_qcols]))
    for c in ps_ok
        println(rpad(c.name, 24), join([lpad(let v = ps_col(ps_res[c.name], q, τ); isnan(v) ? "-" : string(round(v, sigdigits = 3)) end, 10) for q in ps_qcols]))
    end
end

# 熱収支 (budget.csv)
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
rows = DataFrame(case = String[], tau = Float64[], sensor = String[], model = Float64[], exp = Float64[])
for c in ps_ok
    r = ps_res[c.name]
    for τ in 10.0:2.0:ps_t_end - ps_t_heat
        τ <= ps_tmax(r) || continue
        for s in ps_sensors
            push!(rows, (c.name, τ, s[1], s[3](r, τ), ps_exv(s[2], τ)))
        end
    end
end
CSV.write(joinpath(ps_out, "cc_adiab_sweep.csv"), rows)
println("\n保存しました: ", abspath(ps_out), " (cc_adiab_sweep.csv)")