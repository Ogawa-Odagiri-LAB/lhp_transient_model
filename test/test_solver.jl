# =============================================================================
# test_solver.jl
#
# IDA による過渡計算の動作確認 (実走)。冒頭の設定を変えて、短時間で比較できる。
#   1. 完走するか (retcode, 計算時間, IDA の統計)
#   2. t=0 の解が、初期値と一致するか
#   3. 物理的な妥当性 (時間変化の表: 圧力・温度・蒸発量・乾き度)
#   4. 質量保存
#   5. エネルギー収支 (蓄積の変化 ΔE が、入力 ∫Q を上回らないこと)
#   6. :KLU と :Dense の比較 (差が最大の変数を、変数名と時刻で特定する)
#
# 実行時間の目安: テーブル構築が初回のみ数分。IDA は、設定とエネルギー式の
# 離散化によっては、数分〜1時間以上かかる (前回: Dense で 20 s の計算に 3661 s)。
# まず T_END を短く (例: 3 s) して、設定の効果を比べるとよい。
# =============================================================================

using CoolProp, Interpolations, Serialization, ComponentArrays, Parameters, StaticArrays,
      SparseArrays, LinearAlgebra, Random, Sundials
using DiffEqBase: DAEProblem, DAEFunction, solve

# ---------------------------- 設定 ----------------------------
const T_END    = 20.0       # 計算する時間 [s]。蒸気が凝縮器に届くのは t≈7〜11 s
const SAVEAT   = 0.5        # 保存間隔 [s]
const RELTOL   = 1.0e-6     # 相対許容誤差 (緩める: 1e-4, 1e-3)
const ATOL_FRAC = 1.0e-9    # 絶対許容誤差 = ATOL_FRAC × 変数の代表的な大きさ
const ATOL_FRAC_U = 1.0e-9  # 流体の速度 u だけ別の比率 (音波の振動を追う必要を減らす: 1e-7 など)
const MAX_ORDER = 5         # BDF の最大次数 (2 にすると、振動に強く、刻みが大きく取れることがある)
const RUN_KLU  = true       # :KLU でも走らせて、:Dense と比較する
# --------------------------------------------------------------

const SRC = joinpath(@__DIR__, "..", "src")   # テストは ./test/、ソースは ./src/
for f in ("physics", "smoothing", "types", "proptables", "propaccess", "correlations",
          "heattransfer", "geometry_cc", "residuals_transport", "residuals_wall",
          "residuals_wick", "residuals_cc", "residuals_evapcase", "residuals_grv",
          "residuals_main", "initialcondition", "sparsity", "solverfuncs")
    include(joinpath(SRC, f * ".jl"))
end

fluid = "ammonia"
println("テーブル構築中 (初回は数分かかる場合があります。キャッシュされます)...")
tables = init_fluid_lookups(fluid; T_min=230.0, T_max=360.0, P_min=1.0e5, P_max=4.5e6,
                            Ne=400, Δδ=1.0, N_T_sat=200)

# 小さめの分割。ウィックの固相の実効熱伝導率 k_eff は、入力せず、SUS304 の熱伝導率、
# 液の熱伝導率、空隙率 φ から、マクスウェルの複合則 (式53) で計算される
# (φ=0.65 で約 3.6 W/mK)。
wick = WickGeometry(N=3, φ=0.65, ρ_wick=8900.0, Cp_wick=444.0)
geo = LHPGeometry(wick,
                  GrvGeometry(N=1),
                  TLGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  CondGeometry(N=5, L=1.0, D_i=0.00435, D_o=0.00635),
                  TLGeometry(N=5, L=1.0, D_i=0.00218, D_o=0.00318),
                  CCGeometry(N=1, L=0.07, D_i=0.035, D_o=0.037))
condition = LHPTransientCondition(Q_apply=5.0, t_init=0.0, t_end=T_END, saveat=SAVEAT, fluid=fluid)

println("設定: T_END=", T_END, " s, RELTOL=", RELTOL, ", ATOL_FRAC=", ATOL_FRAC,
        ", ATOL_FRAC_U=", ATOL_FRAC_U, ", MAX_ORDER=", MAX_ORDER)

println("="^72)
println("テスト1: :Dense で完走するか")
println("="^72)
t_dense = @elapsed (sol, p, u0) = solve_transient_lhp!(geo, condition, tables;
    linear_solver=:Dense, reltol=RELTOL, atol_frac=ATOL_FRAC, atol_frac_u=ATOL_FRAC_U, max_order=MAX_ORDER)
println("  計算時間: ", round(t_dense, digits=1), " s")
println("  retcode = ", sol.retcode, ", 保存点 ", length(sol.t), " 個")
success = occursin("Success", string(sol.retcode))
println("  成功: ", success)

println()
println("="^72)
println("テスト2: t=0 の解が、初期値と一致するか")
println("="^72)
println("  max|sol.u[1] − u0| = ", maximum(abs, sol.u[1] .- getdata(u0)), " (0 のはず)")

println()
println("="^72)
println("テスト3: 物理的な妥当性 (時間変化)")
println("="^72)
function check_state(k)
    c = components_at(sol, p, k)
    comps = (:grv, :vl, :cond, :ll, :cc)
    ρ_min = minimum(minimum(getproperty(getproperty(c, comp), :ρm)) for comp in comps)
    x_min = minimum(minimum(getproperty(getproperty(c, comp), :x)) for comp in comps)
    x_max = maximum(maximum(getproperty(getproperty(c, comp), :x)) for comp in comps)
    u_max = maximum(maximum(abs, getproperty(getproperty(c, comp), :u)) for comp in comps)
    return (ρ_min = ρ_min, x_min = x_min, x_max = x_max, u_max = u_max, finite = all(isfinite, sol.u[k]))
end
println("  t[s]  | T_ec[K] | P_cc[Pa] | T_cc[K] | P_grv[Pa] | P_cond1[Pa] | ṁ[kg/s] | max|u|[m/s] | x∈[min,max] | 有限")
nshow = min(length(sol.t), 24)
for k in unique(round.(Int, range(1, length(sol.t), length=nshow)))
    c = components_at(sol, p, k)
    ρ_l = tables.ρ_l_sat(c.cc.T[1])
    ṁ = c.cc.u[1] * ρ_l * p.p_cc.A
    s = check_state(k)
    println("  ", lpad(round(sol.t[k], digits=2), 6), " | ", round(c.T_wall.evap[1], digits=3),
            " | ", round(c.cc.P[1], digits=0), " | ", round(c.cc.T[1], digits=3),
            " | ", round(c.grv.P[1], digits=0), " | ", round(c.cond.P[1], digits=0),
            " | ", round(ṁ, sigdigits=3), " | ", round(s.u_max, sigdigits=3),
            " | [", round(s.x_min, digits=3), ",", round(s.x_max, digits=3), "] | ", s.finite)
end
c_end = components_at(sol, p, length(sol.t))
println("  蒸発器ケースの温度上昇 ΔT_ec = ", c_end.T_wall.evap[1] - u0.T_wall.evap[1], " K (正のはず)")

println()
println("="^72)
println("テスト4: 質量保存")
println("="^72)
function total_mass(c)
    m = 0.0
    for (comp, pc) in ((:grv, p.p_grv), (:vl, p.p_vl), (:cond, p.p_cond), (:ll, p.p_ll), (:cc, p.p_cc))
        m += sum(getproperty(getproperty(c, comp), :ρm)) * pc.V
    end
    return m
end
m0 = total_mass(components_at(sol, p, 1)); mend = total_mass(c_end)
println("  全流体の質量: t=0 で ", m0, " kg, t=", sol.t[end], " で ", mend, " kg, 相対変化 ", (mend - m0) / m0)

println()
println("="^72)
println("テスト5: エネルギー収支 (ΔE は ∫Q_apply dt − (放熱) を超えないはず)")
println("="^72)
# 壁の内部エネルギー m ∫_{T_ref}^{T} c(T') dT'。c は温度の関数なので、C(T)·T とは違う
# (T が変わると、T0 ΔC の誤差が出る: 凝縮器の壁が約9 K冷えるとき、約35 J)。
# 基準温度 T_ref は、差だけが意味を持つので、何でもよい。台形則 (20 分割) で積分する。
function wall_internal_energy(T, V_w; T_ref=273.15, n=20)
    Ts = range(T_ref, T, length=n + 1)
    f = [ρ_sus304(t) * c_sus304(t) for t in Ts]
    Iint = sum(0.5 * (f[k] + f[k+1]) * (Ts[k+1] - Ts[k]) for k in 1:n)
    return V_w * Iint
end

function stored_energy(c)
    E = 0.0
    for (comp, pc) in ((:grv, p.p_grv), (:vl, p.p_vl), (:cond, p.p_cond), (:ll, p.p_ll), (:cc, p.p_cc))
        cm = getproperty(c, comp)
        E += sum(cm.ρm .* cm.Em) * pc.V
    end
    E += sum(c.wick.E)
    for (wall, wp) in ((:vl, p.wp_vl), (:cond, p.wp_cond), (:ll, p.wp_ll))
        E += sum(wall_internal_energy(T, wp.V_w) for T in getproperty(c.T_wall, wall))
    end
    E += wall_internal_energy(c.T_wall.cc[1], p.cw.V_w)
    E += wall_internal_energy(c.T_wall.evap[1], p.ec.V_w)
    return E
end
E0 = stored_energy(components_at(sol, p, 1)); Eend = stored_energy(c_end)
Q_vals = [applied_heat(t, condition) for t in sol.t]
trapz(f) = sum(0.5 * (f[k] + f[k+1]) * (sol.t[k+1] - sol.t[k]) for k in 1:length(sol.t)-1)
Q_in = trapz(Q_vals)

# 系から出る熱 [W]: 蒸気管・液管の壁→周囲 (自然対流)、凝縮器の壁→凝縮器プレート。
# CCケース・蒸発器ケースは断熱。Q_amb_* は「壁が受け取る向き」が正なので、符号を反転する。
# 残差 (residuals_wall.jl) と同じ関数を使うので、残差が仮定する放熱と一致する。
function heat_out(c)
    Q = 0.0
    for (wall, wp) in ((:vl, p.wp_vl), (:ll, p.wp_ll))
        for T in getproperty(c.T_wall, wall)
            Q -= Q_amb_pipe(T, wp, condition.T_amb)
        end
    end
    for T in c.T_wall.cond
        Q -= Q_amb_condenser(T, p.wp_cond, condition.T_cond; h_con=condition.h_con)
    end
    return Q
end
Q_out_vals = [heat_out(components_at(sol, p, k)) for k in 1:length(sol.t)]
Q_out = trapz(Q_out_vals)
println("  蓄積したエネルギーの変化 ΔE            = ", Eend - E0, " J")
println("  印加熱の積分 ∫Q_apply dt               = ", Q_in, " J")
println("  放熱の積分 ∫Q_out dt (周囲+凝縮器プレート) = ", Q_out, " J")
println("  収支の残差 ΔE − (∫Q_in − ∫Q_out)       = ", (Eend - E0) - (Q_in - Q_out), " J",
        "  (入力に対して ", round(100 * ((Eend - E0) - (Q_in - Q_out)) / Q_in, digits=2), " %)")
println("  (壁の蓄熱は m∫c dT で計算 (C(T)·T の誤差を避けた)。残差は、台形則の積分誤差 (保存間隔 ", SAVEAT, " s) と、",
        "非保存的な部分 (運動エネルギーの扱いなど) の誤差だけのはず)")

if RUN_KLU
    println()
    println("="^72)
    println("テスト6: :KLU が :Dense と一致するか")
    println("="^72)
    try
        t_klu = @elapsed (sol_k, _, _) = solve_transient_lhp!(geo, condition, tables;
            linear_solver=:KLU, reltol=RELTOL, atol_frac=ATOL_FRAC, atol_frac_u=ATOL_FRAC_U, max_order=MAX_ORDER, verbose=true)
        println("  KLU: retcode = ", sol_k.retcode, ", 計算時間 ", round(t_klu, digits=1), " s (Dense: ", round(t_dense, digits=1), " s)")
        if length(sol_k.t) == length(sol.t)
            scl = variable_scales(geo, condition)
            labels = variable_labels(p.axes, length(scl))
            Uk = hcat(sol_k.u...); Ud = hcat(sol.u...)
            rel = abs.(Uk .- Ud) ./ scl                      # 代表的な大きさに対する相対差
            (rmax, idx) = findmax(rel)
            i, k = Tuple(idx)
            println("  最終解の差: 代表的な大きさに対する相対差の最大 = ", maximum(rel[:, end]))
            println("  全時刻での最大 = ", rmax, "  (変数 ", labels[i], ", t = ", sol.t[k], " s)")
            println("    その値: KLU ", Uk[i, k], ", Dense ", Ud[i, k])
            println("  時刻ごとの、相対差の最大 (振動の位相がずれると、ある時刻から急に大きくなる):")
            for k2 in unique(round.(Int, range(1, length(sol.t), length=min(12, length(sol.t)))))
                (r, ii) = findmax(rel[:, k2])
                println("    t=", lpad(round(sol.t[k2], digits=2), 6), ": ", round(r, sigdigits=3), "  (", labels[ii], ")")
            end
        else
            println("  保存点の数が違う: KLU ", length(sol_k.t), ", Dense ", length(sol.t))
        end
    catch err
        println("  KLU は実行できませんでした: ", sprint(showerror, err)[1:min(300, end)])
    end
end