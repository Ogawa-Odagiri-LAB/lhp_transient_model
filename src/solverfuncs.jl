# =============================================================================
# solverfuncs.jl (IDA による時間積分)
#
# 旧実装 (時間ステップごとに、後退オイラーの非線形方程式を TrustRegion で解く)
# から、完全陰的な DAE ソルバー IDA (Sundials.jl) に置き換えた。
# 残差 F(t, u, du) = 0 (residuals_main.jl の f_dae_residual!) を、IDA が、
# 可変ステップ・可変次数の BDF で時間積分する。時間ステップや Newton 反復の
# 制御は IDA が行うので、時間刻み dt の指定は不要 (許容誤差で制御される)。
#
# 使い方:
#   tables = init_fluid_lookups("ammonia"; T_min=230.0, T_max=360.0, P_max=4.5e6)
#   sol, p, u0 = solve_transient_lhp!(geo, condition, tables)
#
# 【既定の設定と、その根拠】 (test_solver.jl, test_solver_sweep.jl の実測。アンモニア、
#  15 s の起動の計算、125 変数。基準 = 次数5、rtol 1e-6、KLU: 122 s、47.8 万ステップ)
#   linear_solver = :KLU   Dense と一致 (最終解の差 7.5e-6)。次数2 で、Dense より約20倍速い
#                          (収束失敗: Dense 14 278 回 (ステップの42%)、KLU 65 回)
#   max_order = 2          次数5 の刻みは、精度でなく安定性で決まる (rtol を 1e-4 → 1e-6 に
#                          しても、ステップ数が 56.7 万 → 47.7 万と、ほぼ変わらない)。
#                          次数2: 12 622 ステップ、6 s (基準の約 1/20)、基準との差は
#                          P_cc で 0.1 Pa、蒸発量で 3e-4。次数3 は、次数5 より悪かった
#                          (240 万ステップ、2555 s。原因不明) ので使わない
#   reltol = 1e-4          次数2 で、基準との差が上の通り
#   atol_frac_u = 1e-7     流体の速度 u の絶対許容誤差 (→ 1e-6 m/s)。流れの精度に必要な値
#   max_nonlinear_iters=4  上限を 8, 15 に増やしても、効果がなかった (ステップ数が増えた)
#   saveat                 刻みへの影響はなかった (saveat=0.5 と 0.05 で、ステップ・棄却が完全に同じ)
#
# 【線形ソルバー】
#   :KLU (既定)    … 疎行列 (sparsity.jl)。パターンを数値的に検出し、彩色した差分
#                    ヤコビアンを自前で作って渡す (評価回数が色の数まで減る)。
#   :Dense         … 密行列。ヤコビアンは Sundials 内部の差分で作る (n 回の残差評価)。
#                    次数2 では、収束失敗が多く、KLU より大幅に遅い (理由は未特定)。
#
# 【許容誤差】 reltol と、変数ごとの絶対許容誤差 atol_frac × (変数の代表的な大きさ)。
# 液体の圧力は、わずかな密度変化 (1 kg/m³ → 約 1e6 Pa) で大きく変わるので、
# 圧力・密度の絶対許容誤差を厳しめにしてある。
#
# 【初期化】 initial_condition が整合な (u0, du0) を作るので、IDA 側の初期条件の
# 補正 (IDACalcIC) は、ほとんど何も動かさない。
# =============================================================================

using Sundials
using DiffEqBase: DAEProblem, DAEFunction, solve

"""
    absolute_tolerances(geo, condition; atol_frac=1e-9, atol_frac_u=atol_frac) -> Vector{Float64}

変数ごとの絶対許容誤差 = (変数ごとの比率) × (その変数の代表的な大きさ)。
状態ベクトルと同じ並び。代表的な大きさは types.jl の *_scale
(ρm:50 kg/m³, u:10 m/s, Em:4e5 J/kg, P:1e6 Pa, T:300 K, x:1)。

流体の速度 u の比率だけ atol_frac_u で別に指定できる。液管の流速は 7.5e-3 m/s 程度
(20 W) で、必要な精度は、その 1% (7.5e-5 m/s) 程度。一方、圧力の揺らぎ ±50 Pa
が作る音波の速度振幅は δP/(ρc) = 50/(600×1200) ≈ 7e-5 m/s で、既定の
atol_frac=1e-9 (→ 1e-8 m/s) では、許容誤差の約100倍になり、音波の振動を精度よく
追うために、刻みが数十 μs まで小さくなる。u の比率を 1e-7 (→ 1e-6 m/s) 程度に
緩めると、流れの精度を保ったまま、音波を追う必要が減る。
"""
function absolute_tolerances(geo, condition; atol_frac::Float64=1.0e-9,
    atol_frac_u::Float64=atol_frac)
    frac = initialize_state(geo)
    getdata(frac) .= atol_frac
    for comp in condition.fluid_comp
        getproperty(frac, comp).u .= atol_frac_u
    end
    return copy(getdata(frac)) .* variable_scales(geo, condition)
end

"""
    DOMAIN_ERROR_COUNT

ガード付きの残差 (f_dae_residual_guarded!) が、DomainError を捕まえた回数。
solve_transient_lhp! が、計算の最初に 0 に戻し、最後に、1 回以上あれば表示する。
"""
const DOMAIN_ERROR_COUNT = Ref(0)

"""
    f_dae_residual_guarded!(res, du, u, p, t)

f_dae_residual! を、DomainError (sqrt や小数べきの底が負、など) から守る。

IDA の Newton 反復の途中の試行点は、収束しない (刻みが大きすぎる、ヤコビアンが古い)
ときに、1〜2 回の反復で、物理的にありえない状態 (ρm<0、物性テーブルの範囲外の T など)
へ飛ぶことがある。その試行点で、相関式が DomainError を投げると、IDA は、その刻みを
縮めて再試行する前に、計算全体が中断してしまう (実測: Q 5→10 W の低熱負荷で、
t=170 s の凝縮の前線のセルで、約 200 万回目の残差評価で発生。解そのものは、直前まで
正常で、ṁ>0、圧力も一定だった)。

DomainError を捕まえて残差を NaN にして返すと、IDA は、その反復が収束しなかったと
見なして、刻みを縮めて再試行する (SciML の推奨: 解が定義域を出る試行点は、NaN で棄却させる)。
DomainError 以外の例外 (バグ) は、そのまま投げる。
"""
function f_dae_residual_guarded!(res, du, u, p, t)
    # u/du に NaN/Inf が入っていたら、物性テーブルの補間に渡す前に NaN 残差を返す
    # (補間に NaN を渡すと範囲外読み出しで Julia ごとクラッシュする)。
    if !(all(isfinite, u) && all(isfinite, du))
        DOMAIN_ERROR_COUNT[] += 1
        fill!(res, NaN)
        return nothing
    end
    try
        f_dae_residual!(res, du, u, p, t)
    catch e
        e isa DomainError || rethrow()
        DOMAIN_ERROR_COUNT[] += 1
        fill!(res, NaN)
    end
    return nothing
end

"""
    build_dae_problem(p, u0, du0, diffmask, tspan; linear_solver=:Dense, ...) -> DAEProblem

DAEProblem を作る。linear_solver=:KLU のときは、疎行列のパターンを検出し、
彩色した差分ヤコビアン (ColoredFDJacobian) を jac として渡す。
"""
function build_dae_problem(p::AbstractLHPParams, u0::AbstractVector, du0::AbstractVector,
    diffmask::AbstractVector{Bool}, tspan, scales; linear_solver::Symbol=:Dense,
    nprobe::Int=3, verbose::Bool=true, guard_domain_errors::Bool=true)

    fres = guard_domain_errors ? f_dae_residual_guarded! : f_dae_residual!

    if linear_solver === :KLU
        verbose && println("疎行列のパターンを検出中 (変数 $(length(u0)) 個、ランダムな状態 $(nprobe) 個)...")
        pattern = detect_jacobian_pattern(p, u0, tspan[1], scales, diffmask; nprobe=nprobe)
        color = greedy_coloring(pattern)
        check_coloring(pattern, color) || error("彩色が妥当ではありません (同じ色の列が行を共有)")
        jf = ColoredFDJacobian(p, pattern, color, scales)
        if verbose
            println("  非ゼロ: $(count(pattern)) / $(length(u0)^2) ($(round(100 * count(pattern) / length(u0)^2, digits=2)) %), ",
                    "色の数: $(maximum(color)) (密行列の差分なら $(length(u0)) 回の評価)")
        end
        f = DAEFunction(fres; jac=jf, jac_prototype=jf.prototype)
        return DAEProblem(f, du0, u0, tspan, p; differential_vars=diffmask)
    elseif linear_solver === :Dense
        return DAEProblem(fres, du0, u0, tspan, p; differential_vars=diffmask)
    else
        error("linear_solver は :Dense または :KLU")
    end
end

"""
    solve_transient_lhp!(geo, condition, tables;
                         linear_solver=:KLU, reltol=1e-4, atol_frac=1e-9, atol_frac_u=1e-7,
                         dtmax=Inf, maxiters=10^7, max_order=2, verbose=true, kwargs...) -> (sol, p, u0)

過渡計算を実行する。
  戻り値  sol: IDA の解 (sol.t, sol.u)。sol.u[k] は、状態ベクトル (Vector) で、
                ComponentArray(sol.u[k], p.axes) で、コンポーネントごとに取り出せる
          p  : LHPParams (後処理や、残差の評価に使う)
          u0 : 初期値 (ComponentArray)
  condition.saveat ごとに解を保存する。ただし、保存時刻ごとに IDA が停止・再開して、
  刻みが歪む疑いがあるので、saveat は粗く (0.5 s 以上) 取るのが無難。内部ステップを
  全て保存したいときは、kwargs に save_everystep=true を渡す (saveat は使われない)。
  Q_apply のステップの区切り点 (applied_heat_breakpoints) を tstops に渡す。
  atol_frac_u: 流体の速度 u の絶対許容誤差の比率 (既定 1e-7)。absolute_tolerances 参照。
  max_nonlinear_iters: 1 ステップあたりの Newton 反復の上限 (IDA の既定 4)。収束しないと
             収束失敗として、ヤコビアンを再計算して刻みを縮める (113 回の残差評価を伴う)。
             増やすと、失敗が減ることがある。
  guard_domain_errors: true (既定) なら、残差の DomainError を捕まえて NaN を返し、IDA に、
             その刻みを棄却させる (f_dae_residual_guarded! 参照)。false なら、例外で計算が中断する。
  max_order: BDF の最大次数 (1〜5)。3〜5次は、減衰の弱い振動 (音響波) に対して安定
             領域が狭く、刻みが強く制限されることがある。2 にすると A 安定になり、
             振動に強くなる (精度は下がる)。
  kwargs は、IDA の solve にそのまま渡される (例: initializealg=...)。
"""
function solve_transient_lhp!(geo::AnyLHPGeometry, condition::LHPTransientCondition,
    tables::FluidTables; linear_solver::Symbol=:KLU, reltol::Float64=1.0e-4,
    atol_frac::Float64=1.0e-9, atol_frac_u::Float64=1.0e-7, dtmax::Float64=Inf,
    maxiters::Int=10^7, max_order::Int=2, max_nonlinear_iters::Int=4,
    verbose::Bool=true, guard_domain_errors::Bool=true, kwargs...)

    DOMAIN_ERROR_COUNT[] = 0
    # --- 整合初期条件 ---
    u0, du0, p = initial_condition(geo, condition, tables; verbose=verbose)
    diffmask = Bool.(getdata(initialize_differential_vars(geo)))
    scales = variable_scales(geo, condition)
    abstol = absolute_tolerances(geo, condition; atol_frac=atol_frac, atol_frac_u=atol_frac_u)

    tspan = (condition.t_init, condition.t_end)
    prob = build_dae_problem(p, copy(getdata(u0)), du0, diffmask, tspan, scales;
                             linear_solver=linear_solver, verbose=verbose,
                             guard_domain_errors=guard_domain_errors)

    # 既定値のキーワードは渡さない (Sundials.jl のバージョンで IDA のキーワード名が
    # 違っても、既定の動作が影響を受けないように)
    ida_kwargs = Pair{Symbol,Any}[:linear_solver => linear_solver]
    max_order != 5 && push!(ida_kwargs, :max_order => max_order)
    max_nonlinear_iters != 4 && push!(ida_kwargs, :max_nonlinear_iters => max_nonlinear_iters)
    solver = IDA(; ida_kwargs...)
    tstops = applied_heat_breakpoints(condition)
    verbose && println("IDA で時間積分を開始: t = ", tspan[1], " 〜 ", tspan[2], " s, ",
                       "線形ソルバー ", linear_solver, ", 最大次数 ", max_order,
                       ", Newton 反復の上限 ", max_nonlinear_iters,
                       ", reltol=", reltol, ", atol_frac=", atol_frac, ", atol_frac_u=", atol_frac_u,
                       isempty(tstops) ? "" : ", tstops=$(tstops)")

    # save_everystep=true (IDA の内部ステップを全て保存) を指定したときは、saveat を
    # 渡さない。保存時刻の指定は、IDA が保存時刻ごとに停止・再開して、刻みが歪む疑い
    # がある (test_solver_sweep.jl 参照)。
    save_kwargs = haskey(kwargs, :save_everystep) ? NamedTuple() : (saveat = condition.saveat,)
    sol = solve(prob, solver; abstol=abstol, reltol=reltol,
                tstops=tstops, dtmax=dtmax, maxiters=maxiters,
                save_kwargs..., kwargs...)

    if DOMAIN_ERROR_COUNT[] > 0
        println("注意: 残差の評価で DomainError を ", DOMAIN_ERROR_COUNT[], " 回捕まえ、その試行点を棄却させました ",
                "(IDA の Newton 反復の試行点が、定義域の外へ飛んだもの。刻みを縮めて再試行されます)")
    end
    if verbose
        println("IDA 終了: retcode = ", sol.retcode, ", 保存点 ", length(sol.t),
                " 個, 最終時刻 t = ", sol.t[end], " s")
        # 時間が何に使われたかの診断: 残差の評価回数 nf、ヤコビアンの評価回数 njacs、
        # 受理された/棄却されたステップ数 naccept/nreject など
        println("IDA の統計:")
        for f in fieldnames(typeof(sol.stats))
            v = getfield(sol.stats, f)
            (v isa Number && v != 0) && println("    ", rpad(string(f), 18), v)
        end
    end
    return sol, p, u0
end

"""
    components_at(sol, p, k) -> ComponentArray

k 番目の保存点の解を、コンポーネントごとに取り出せる ComponentArray にする
(例: components_at(sol, p, 5).cc.P[1])。
"""
components_at(sol, p::AbstractLHPParams, k::Int) = ComponentArray(copy(sol.u[k]), p.axes)