# =============================================================================
# test_diagnose.jl
#
# IDA の刻みを制限している原因を調べる診断。
#   スイープ (test_solver_sweep.jl) の結果:
#     ・BDF 次数 2 で、解を変えずに、ステップが 45 万 → 3.4 万、棄却が 6.3 万 → 96
#     ・しかし、Newton 法の収束失敗がステップの 42% で起き、ヤコビアンの再計算で
#       時間の大半を使う (1 ステップ 約 23 ms)
#     ・許容誤差を緩めると (rtol 1e-4 → 1e-3)、収束失敗が増え、遅くなる
#       → 刻みを制限しているのは、精度ではなく、非線形方程式の収束
#     ・P_grv の ±650 Pa の揺らぎが、全設定で残る → 方程式が持つ振動
#
# このスクリプトは、
#   A. 圧力・速度・乾き度の時系列 (5 ms 間隔) と、P_grv の揺らぎの周期
#   B. 同じ状態で、後退オイラーの非線形方程式を Newton 法 (毎反復で新しい
#      ヤコビアン) で解き、収束の様子と、収束を妨げる変数 (変数名) を表示
# を行う。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials, Statistics
using DiffEqBase: DAEProblem, DAEFunction, solve

# ---------------------------- 設定 ----------------------------
const T_PROBE = 3.0        # 診断する時刻 [s] (起動の過渡が落ち着いた後)
const T_WINDOW = 0.5       # 時系列を取る区間の長さ [s]
const SAVE_DT = 0.005      # 時系列の間隔 [s]
const RTOL = 1.0e-3        # IDA の許容誤差 (診断用の前進計算。次数 2 で十分)
const ATOL_FRAC = 1.0e-9
const ATOL_FRAC_U = 1.0e-7
const H_LIST = (1.0e-5, 1.0e-4, 1.0e-3)   # Newton 法を試す刻み [s]
# --------------------------------------------------------------

const SRC = joinpath(@__DIR__, "..", "src")
for f in ("physics", "smoothing", "types", "proptables", "propaccess", "correlations",
          "heattransfer", "geometry_cc", "residuals_transport", "residuals_wall",
          "residuals_wick", "residuals_cc", "residuals_evapcase", "residuals_grv",
          "residuals_main", "initialcondition", "sparsity", "solverfuncs")
    include(joinpath(SRC, f * ".jl"))
end

fluid = "ammonia"
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)
wick = WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0)
geo = LHPGeometry(wick,
                  GrvGeometry(N=1),
                  TLGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  CondGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  TLGeometry(N=5, L=1.0, D_i=0.00218, D_o=0.00318),
                  CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037))

# 保存時刻を指定すると、IDA が保存時刻ごとに停止・再開して、刻みが歪む疑いがある
# (test_solver_sweep.jl 参照)。saveat は終了時刻だけにして、IDA の内部ステップを
# 全て保存し (save_everystep=true)、そこから時系列を作る。次数2の刻みは約 150 μs で、
# 5 ms 間隔より細かい。
t_end = T_PROBE + T_WINDOW
condition = LHPTransientCondition(Q_apply=20.0, t_init=0.0, t_end=t_end, fluid=fluid, saveat=t_end)

println("前進計算 (次数2) で t=", t_end, " s まで進めます (内部ステップを全て保存)...")
t_run = @elapsed (sol, p, u0) = solve_transient_lhp!(geo, condition, tables;
    linear_solver=:KLU, reltol=RTOL, atol_frac=ATOL_FRAC, atol_frac_u=ATOL_FRAC_U,
    max_order=2, verbose=false, save_everystep=true)
println("  完了 (", round(t_run, digits=1), " s), retcode=", sol.retcode, ", 保存ステップ数 ", length(sol.t))
n = length(sol.u[1])
labels = variable_labels(p.axes, n)

# =============================================================================
println()
println("="^76)
println("A. 時系列 (t=", T_PROBE, "〜", t_end, " s, ", SAVE_DT * 1e3, " ms 間隔)")
println("="^76)
# 各ステップの値 (窓 [T_PROBE, t_end] 内)
idx = findall(t -> T_PROBE - 1e-12 <= t <= t_end + 1e-12, sol.t)
println("  窓内の内部ステップ数 ", length(idx), " (平均刻み ", round(T_WINDOW / max(length(idx) - 1, 1) * 1e6, digits=1), " μs)")
tt = sol.t[idx]
getcomp(f) = [f(components_at(sol, p, k)) for k in idx]
Pg_raw = getcomp(c -> c.grv.P[1])
function interp1(ts, ys, t)       # 線形補間
    k = searchsortedlast(ts, t)
    k >= length(ts) && return ys[end]
    k < 1 && return ys[1]
    w = (t - ts[k]) / (ts[k+1] - ts[k])
    return (1 - w) * ys[k] + w * ys[k+1]
end
println("   t[s]    P_grv[Pa]  P_vl1[Pa]  P_vl5[Pa]  P_cond1[Pa]  u_grv[m/s]  u_vl1[m/s]  x_grv    x_vl1    x_cond1   ṁ[kg/s]")
t_print = T_PROBE
for k in idx
    if sol.t[k] >= t_print - 1e-12
        c = components_at(sol, p, k)
        ṁ = c.cc.u[1] * tables.ρ_l_sat(c.cc.T[1]) * p.p_cc.A
        println("  ", lpad(round(sol.t[k], digits=3), 6), "  ", lpad(round(c.grv.P[1], digits=0), 10), " ",
                lpad(round(c.vl.P[1], digits=0), 10), " ", lpad(round(c.vl.P[end], digits=0), 10), " ",
                lpad(round(c.cond.P[1], digits=0), 11), "  ", lpad(round(c.grv.u[1], digits=4), 10), "  ",
                lpad(round(c.vl.u[1], digits=4), 10), "  ", round(c.grv.x[1], digits=5), "  ",
                round(c.vl.x[1], digits=5), "  ", round(c.cond.x[1], digits=5), "  ", round(ṁ, sigdigits=3))
        t_print += 0.02          # 20 ms ごとに表示
    end
end
# 均等な 5 ms 間隔にリサンプルして、離散フーリエ変換に使う
ts = collect(range(T_PROBE, t_end, step=SAVE_DT))
Pg = [interp1(tt, Pg_raw, t) for t in ts]
# P_grv の離散フーリエ変換 (直流を除いて、周波数 0.5〜100 Hz の振幅)
Pg_ac = Pg .- mean(Pg)
Nt = length(Pg_ac); T_total = ts[end] - ts[1]
freqs = collect(range(1.0, 1.0 / (2 * SAVE_DT), length=200))
amp = [abs(sum(Pg_ac[j] * cis(-2π * f * (ts[j] - ts[1])) for j in 1:Nt)) * 2 / Nt for f in freqs]
imax = argmax(amp)
println("  P_grv の揺らぎ: 標準偏差 ", round(std(Pg), digits=1), " Pa, 最大−最小 ", round(maximum(Pg) - minimum(Pg), digits=1), " Pa")
println("  支配的な周波数 ≈ ", round(freqs[imax], digits=1), " Hz (周期 ", round(1000 / freqs[imax], digits=1),
        " ms), その振幅 ", round(amp[imax], digits=1), " Pa")
println("  (参考: 蒸気管 (L=1 m, 音速≈400 m/s) の 1/4 波長共振 ≈ 100 Hz、液管/凝縮器 (音速≈1200 m/s,")
println("   セル長 0.2 m) の格子の最短波長 ≈ 3 kHz。5 ms 間隔では、100 Hz 以上は折り返しで見える)")

# =============================================================================
println()
println("="^76)
println("B. Newton 法の収束の様子 (後退オイラー、毎反復で新しい密ヤコビアン)")
println("="^76)
kp = argmin(abs.(sol.t .- T_PROBE))      # T_PROBE に最も近い内部ステップ
u_n = copy(sol.u[kp])
tn = sol.t[kp]
scales = variable_scales(geo, condition)
abstol = absolute_tolerances(geo, condition; atol_frac=ATOL_FRAC, atol_frac_u=ATOL_FRAC_U)
res = zeros(n); res1 = zeros(n)

function wrms(v, u)
    ewt_inv = RTOL .* abs.(u) .+ abstol
    return sqrt(sum((v ./ ewt_inv) .^ 2) / length(v))
end
function dense_jacobian(u, du, γ, t)
    f_dae_residual!(res, du, u, p, t)
    F0 = copy(res)
    J = zeros(n, n)
    for j in 1:n
        d = sqrt(eps(Float64)) * max(abs(u[j]), 1.0e-3 * scales[j])
        u1 = copy(u); du1 = copy(du)
        u1[j] += d; du1[j] += γ * d
        f_dae_residual!(res1, du1, u1, p, t)
        J[:, j] = (res1 .- F0) ./ d
    end
    return F0, J
end

println("  t_n = ", tn, " s の状態から、刻み h で 1 ステップ進める (予測子: u_n そのまま)")
for h in H_LIST
    println("  --- h = ", h, " s ---")
    u = copy(u_n)
    converged = false
    for it in 1:20
        du = (u .- u_n) ./ h
        F0, J = dense_jacobian(u, du, 1.0 / h, tn + h)
        δ = -(J \ F0)
        u .+= δ
        nrm = wrms(δ, u)
        # 更新量への寄与が最大の変数と、残差が最大の方程式 (重み付き)
        ewt_inv = RTOL .* abs.(u) .+ abstol
        iδ = argmax(abs.(δ) ./ ewt_inv)
        iF = argmax(abs.(F0) ./ (ewt_inv .* (1 / h)))
        println("    反復 ", lpad(it, 2), ": WRMS(更新量) = ", lpad(round(nrm, sigdigits=3), 10),
                "   最大寄与: ", rpad(labels[iδ], 12), " (δ=", round(δ[iδ], sigdigits=3), ")",
                "   最大残差: ", labels[iF])
        if nrm < 0.1
            println("    → 収束 (WRMS<0.1, ", it, " 回)")
            converged = true
            break
        end
    end
    converged || println("    → 20 回で収束せず")
end
println()
println("  (毎反復でヤコビアンを更新しても、反復が多い・滞る場合は、その変数の方程式の非線形性が")
println("   収束を妨げている。3〜4 回で収束するなら、IDA が古いヤコビアンを使い回すために、")
println("   収束失敗が多いと考えられる。)")