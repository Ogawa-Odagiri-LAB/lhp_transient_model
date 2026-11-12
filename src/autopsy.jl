# =============================================================================
# autopsy.jl (IDA が止まったときの、残差の滑らかさとヤコビアンの構造の検査。診断のみ)
#
# failure_residual_diag の F は、du を 0.4 s 間隔の保存点の差分で近似しており、信頼できない。
# ここでは du に依存しない検査をする:
#   (0) 前線まわりのセルの飽和液線からの距離 δ = ρ − ρ_sl(e) (tables を渡したとき)
#   (1) 滑らかさの検査: 最後の保存点の状態で、focus に合う変数 j を相対 h = 1e-9, 1e-7, 1e-5, 1e-3 だけ
#       動かし (中心差分)、残差の各行 i の偏微分 ∂F_i/∂u_j を求める。
#       Newton が使う小さい刻み (1e-9 と 1e-5) で偏微分が大きく変わる (i, j) は、
#       その幅の中に折れ・急な変化がある。h=1e-3 での変化は参考 (大きな刻みで折れをまたぐだけ)。
#   (2) ヤコビアンの構造の検査: 初期状態で検出した疎パターンに無い (i, j) のうち偏微分が大きいもの。
# full=true で全変数を調べる (数分)。既定は focus の変数だけ (数秒)。
# =============================================================================

const AUTOPSY_FOCUS = r"^(vl\.\w+\[(10|11|12|13|14)\]|ll\.\w+\[(18|19|20)\]|nuc\.|cc\.|grv\.)"

"""
    failure_autopsy(sol, p, geo, condition, u0; tables = nothing, focus = AUTOPSY_FOCUS, full = false, n_top = 25) -> Nothing

診断のみ。IDA が止まったとき、sol の最後の保存点で (0) δ (1) 滑らかさ (2) パターンの見落とし を調べる。
"""
function failure_autopsy(sol, p::AbstractLHPParams, geo, condition, u0; tables = nothing, focus::Regex = AUTOPSY_FOCUS, full::Bool = false, n_top::Int = 25)
    length(sol.t) < 2 && return nothing
    k = length(sol.t)
    u = collect(Float64, sol.u[k])
    du = (u .- collect(Float64, sol.u[k-1])) ./ (sol.t[k] - sol.t[k-1])
    t = sol.t[k]
    n = length(u)
    labels = variable_labels(p.axes, n)
    idx(name) = findfirst(==(name), labels)

    # ---- (0) 飽和液線からの距離 ----
    if tables !== nothing
        println("\n[飽和液線からの距離] t = ", round(t, digits = 4), " s,  δ = ρm − ρ_sl(Em) [kg/m³] (正: 液側 = 剛、負: 二相側)")
        for (blk, cells) in (("vl", 9:14), ("ll", 17:20), ("cond", (1, 2, 25, 50)))
            for c in cells
                a = idx("$blk.ρm[$c]"); b = idx("$blk.Em[$c]")
                (a === nothing || b === nothing) && continue
                δ = u[a] - tables.P_of_ρe.ρ_ref(u[b])
                println("    ", rpad("$blk[$c]", 10), "ρ = ", rpad(round(u[a], digits = 3), 10), "δ = ", round(δ, sigdigits = 3))
            end
        end
    end

    hs = (1.0e-9, 1.0e-7, 1.0e-5, 1.0e-3)
    F_p = zeros(n); F_m = zeros(n)
    safe!(F, du_, u_, p_, t_) = try
        f_dae_residual!(F, du_, u_, p_, t_)
    catch err
        err isa DomainError || rethrow()
        fill!(F, NaN)
    end
    D = [zeros(n) for _ in hs]
    nl_list = Tuple{Float64,Float64,Float64,Int,Int,NTuple{4,Float64}}[]   # (impact, nl_small, nl_big, i, j, partials)
    mi_list = Tuple{Float64,Int,Int}[]
    diffmask = Bool.(getdata(initialize_differential_vars(geo)))
    scales = variable_scales(geo, condition)
    pattern = detect_jacobian_pattern(p, collect(Float64, getdata(u0)), condition.t_init, scales, diffmask)
    cols = full ? (1:n) : [j for j in 1:n if occursin(focus, labels[j])]
    println("\n[滑らかさの検査] t = ", round(t, digits = 4), " s。調べる変数 ", length(cols), " 個 (", full ? "全部" : "focus", ")。相対 h = ", hs)
    for j in cols
        base = max(abs(u[j]), 1.0e-6)
        for (a, h) in enumerate(hs)
            δ = h * base; uj = u[j]
            u[j] = uj + δ; safe!(F_p, du, u, p, t)
            u[j] = uj - δ; safe!(F_m, du, u, p, t)
            u[j] = uj
            @. D[a] = (F_p - F_m) / (2δ)
        end
        step = 1.0e-5 * base
        for i in 1:n
            d1, d3, d4 = D[1][i], D[3][i], D[4][i]
            (isfinite(d1) && isfinite(d3) && isfinite(d4)) || continue
            mag = max(abs(d1), abs(d3))
            mag * step < 1.0e-9 && continue
            nls = abs(d3 - d1) / (abs(d3) + abs(d1) + 1.0e-300)
            nlb = abs(d4 - d1) / (abs(d4) + abs(d1) + 1.0e-300)
            nls > 0.1 && push!(nl_list, (mag * step, nls, nlb, i, j, (D[1][i], D[2][i], D[3][i], D[4][i])))
            if !pattern[i, j] && i != j
                push!(mi_list, (maximum(abs.(getindex.(D, i))) * step, i, j))
            end
        end
    end
    sort!(nl_list; by = x -> -x[1])
    println("  Newton の刻み (h=1e-9 と 1e-5) で偏微分が 10% 以上変わる (i, j) の数: ", length(nl_list), "  (上位 ", n_top, " を F の変化の大きさ順に)")
    println("    行 i (残差)           列 j (変数)           F の変化(相対1e-5)  変化率(1e-9↔1e-5)  変化率(1e-9↔1e-3)   ∂F/∂u (h=1e-9 | 1e-7 | 1e-5 | 1e-3)")
    for (imp, nls, nlb, i, j, ds) in nl_list[1:min(n_top, length(nl_list))]
        println("    ", rpad(labels[i], 22), rpad(labels[j], 22), lpad(string(round(imp, sigdigits = 3)), 12), lpad(string(round(nls, digits = 2)), 16), lpad(string(round(nlb, digits = 2)), 18), "   ",
                join([string(round(d, sigdigits = 3)) for d in ds], " | "))
    end
    sort!(mi_list; by = x -> -x[1])
    println("\n[パターンの検査] 疎パターンに無いが偏微分が大きい (i, j): ", length(mi_list), " 個 (上位 ", min(n_top, length(mi_list)), ")")
    for (imp, i, j) in mi_list[1:min(n_top, length(mi_list))]
        println("    行 ", rpad(labels[i], 22), " ← 変数 ", rpad(labels[j], 22), "  F の変化 = ", round(imp, sigdigits = 3))
    end
    return nothing
end