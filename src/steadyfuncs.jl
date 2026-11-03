# =============================================================================
# steadyfuncs.jl (定常状態の計算: 擬似時間発展)
#
# 現在の過渡モデル (残差 F(t,u,du)=0) を、そのまま使って定常状態を求める。
# 一定の印加熱量 Q のもとで、IDA で時間積分を続け、状態の時間変化が十分小さく
# なったところで止める (擬似時間発展)。
#
# 【なぜ F(u, 0)=0 を直接ニュートン法で解かないか】
#   閉じた系なので、作動流体の総質量は保存される。F(u,0)=0 は総質量を決める式を
#   含まない (質量はどの値でも定常になりうる) ので、ヤコビアンが特異に近くなる。
#   時間発展なら、初期条件で決めた総質量 (= 充填量) が保存されたまま、その総質量
#   での定常状態に収束する。LHP の起動は複数の定常解を持ちうる (起動前の静止状態と
#   循環状態など) が、時間発展なら、初期条件から物理的に到達できる側に落ち着く。
#
# 【使い方】 scripts/steady.jl 参照。
#   sol, p, u_ss, info = solve_steady_lhp!(geo, condition, tables)
#   save_steady_state(path, u_ss, p, condition)
#   u_init = load_steady_state(path, geo)    # → solve_transient_lhp!(...; u_init=u_init)
#
# 時間発展は t_chunk [s] ごとの区間に分け、各区間の終わりで、du0 を残差から逆算して
# (initialcondition.jl の compute_du0)、時間変化率を調べる。収束していれば終了、
# そうでなければ、その状態から次の区間を続ける。
# =============================================================================

using DiffEqBase: remake
import DiffEqBase

"""
    steady_rates(u, p, t, diffmask) -> (dT_wall, dρ_comp, du, loc)

状態 u での時間変化率。
  dT_wall : max|dT_wall/dt| [K/s]    (壁温度の変化率)
  dρ_comp : max_c |Σ_i dρm_i/dt| / Σ_i ρm_i [1/s]
            コンポーネント c (蒸気管・凝縮器・液管・グルーブ・CC) ごとの、質量の変化率
            (総和の変化率 / 総和)。コンポーネントの中の気液界面の移動 (1セルの密度が
            大きく変わる) は、そのコンポーネントの総質量の変化としてしか効かないので、
            1セルだけの激しい変化 (凝縮器入口の前線など) に引きずられない。
            セルごとの max|dρ/dt|/ρ は、loc.ρ_cell に表示用として出す。
  du      : 微分変数の du (compute_du0 で、残差が du について線形なことから逆算)
  loc     : 最大値の場所 (T: 壁温度、ρ: コンポーネント、ρ_cell: セルごとの最大)
"""
function steady_rates(u::AbstractVector, p::AbstractLHPParams, t, diffmask::AbstractVector{Bool})
    du, _ = compute_du0(u, p, t, diffmask)
    d = ComponentArray(du, p.axes)
    x = ComponentArray(copy(u), p.axes)
    dT_wall = 0.0
    locT = "-"
    for w in keys(d.T_wall)
        v = getproperty(d.T_wall, w)
        for i in eachindex(v)
            if abs(v[i]) > dT_wall
                dT_wall = abs(v[i]); locT = "T_wall.$(w)[$(i)]"
            end
        end
    end
    dρ_comp = 0.0
    locρ = "-"
    dρ_cell = 0.0
    locρc = "-"
    for comp in p.cond.fluid_comp
        ρ = getproperty(x, comp).ρm
        dρ = getproperty(d, comp).ρm
        r = abs(sum(dρ)) / sum(ρ)
        if r > dρ_comp
            dρ_comp = r; locρ = string(comp)
        end
        for i in eachindex(ρ)
            if ρ[i] > 1.0e-6 && abs(dρ[i]) / ρ[i] > dρ_cell
                dρ_cell = abs(dρ[i]) / ρ[i]; locρc = "$(comp).ρm[$(i)]"
            end
        end
    end
    return dT_wall, dρ_comp, du, (T = locT, ρ = locρ, ρ_cell = locρc, dρ_cell = dρ_cell)
end

"""
    solve_steady_lhp!(geo, condition, tables;
                      t_chunk=50.0, t_max=20000.0, tol_dT=1e-3, tol_dρ=1e-5,
                      linear_solver=:KLU, reltol=1e-3, atol_frac=1e-9, atol_frac_u=1e-7,
                      max_order=1, max_nonlinear_iters=4, guard_domain_errors=true,
                      u_init=nothing, verbose=true) -> (sol, p, u_ss, info)

一定の印加熱量 (condition.Q_apply_init を使う。ステップはつけない) で、定常状態まで
擬似時間発展させる。
  終了条件: max|dT_wall/dt| < tol_dT [K/s] かつ、コンポーネントごとの質量変化率
            max_c |Σ dρm/dt|/Σ ρm < tol_dρ [1/s] (steady_rates 参照)
            (または時刻が t_max に達する、IDA が失敗する)
  戻り値:   sol   … 最後の区間の IDA の解 (保存点は、区間の終わりだけ)
            p     … パラメータ
            u_ss  … 定常状態の状態ベクトル (Vector)
            info  … (converged, t, dT_wall, dρ_rel, history)  history は区間ごとの記録
  max_order = 1 (後退オイラー): 定常状態だけが欲しいので、時間精度は要らない。1次の
            後退オイラーは L 安定で、減衰の弱い音響振動を、追わずに潰して大きな刻みで進める
            (次数2 は、振動を追うために刻みが小さくなる。過渡計算の 4.8 s/s の遅さの一因)。
            定常解 F(u,0)=0 は、積分法に依らない。reltol も 1e-3 に緩めてある。
  u_init:   指定すると、その状態から始める (前回の続き、別の熱負荷の定常解など)。
"""
function solve_steady_lhp!(geo::AnyLHPGeometry, condition::LHPTransientCondition,
    tables::FluidTables; t_chunk::Float64=50.0, t_max::Float64=20000.0,
    tol_dT::Float64=1.0e-3, tol_dρ::Float64=1.0e-5,
    linear_solver::Symbol=:KLU, reltol::Float64=1.0e-3,
    atol_frac::Float64=1.0e-9, atol_frac_u::Float64=1.0e-7,
    max_order::Int=1, max_nonlinear_iters::Int=4, guard_domain_errors::Bool=true,
    u_init::Union{Nothing,AbstractVector}=nothing, verbose::Bool=true)

    DOMAIN_ERROR_COUNT[] = 0
    # 一定の熱負荷 (ステップなし)・時刻 0 から。T_cond, T_amb, ΔT_startup などは condition のまま
    cs = LHPTransientCondition(condition; Q_apply = condition.Q_apply_init,
                               Q_apply_init = condition.Q_apply_init,
                               t_init = 0.0, t_end = t_chunk, t_step = 0.0)

    u0, du0, p = initial_condition(geo, cs, tables; verbose=verbose, u_init=u_init)
    diffmask = Bool.(getdata(initialize_differential_vars(geo)))
    scales = variable_scales(geo, cs)
    abstol = absolute_tolerances(geo, cs; atol_frac=atol_frac, atol_frac_u=atol_frac_u)

    u = copy(getdata(u0))
    prob = build_dae_problem(p, copy(u), du0, diffmask, (0.0, t_chunk), scales;
                             linear_solver=linear_solver, verbose=verbose,
                             guard_domain_errors=guard_domain_errors)
    ida_kwargs = Pair{Symbol,Any}[:linear_solver => linear_solver]
    max_order != 5 && push!(ida_kwargs, :max_order => max_order)
    max_nonlinear_iters != 4 && push!(ida_kwargs, :max_nonlinear_iters => max_nonlinear_iters)
    solver = IDA(; ida_kwargs...)

    verbose && println("定常計算 (擬似時間発展): Q = ", cs.Q_apply, " W, T_cond = ", cs.T_cond,
                       " K, T_amb = ", cs.T_amb, " K, 区間 ", t_chunk, " s, 最大 ", t_max,
                       " s, 終了条件 dT_wall < ", tol_dT, " K/s かつ 質量変化率 < ", tol_dρ, " 1/s")

    history = NamedTuple[]
    sol = nothing
    t = 0.0
    converged = false
    du = du0
    while t < t_max
        t1 = t + t_chunk
        prob_k = remake(prob; u0=copy(u), du0=copy(du), tspan=(t, t1))
        # initializealg = NoInit: 区間の初期値は、前の区間の IDA の解 (代数拘束の残差は IDA の
        # 許容誤差の大きさ) なので、既定の CheckInit (残差を abstol 以下に要求する) には通らない。
        # 検査をやめて、そのまま IDA に渡す (最初のステップの Newton 反復が残差を詰める)。
        t_w = @elapsed sol = solve(prob_k, solver; abstol=abstol, reltol=reltol, saveat=[t1],
                    maxiters=10^7, initializealg=DiffEqBase.SciMLBase.NoInit())
        if !(string(sol.retcode) in ("Success", "Default")) || isempty(sol.u)
            @warn "IDA が失敗しました (retcode = $(sol.retcode), t = $(t) 〜 $(t1) s)。直前の状態を返します。" *
                  (isempty(history) ? "" : " 直前の区間の最大変化率の位置: $(history[end].loc)")
            break
        end
        u = copy(sol.u[end])
        t = sol.t[end]
        dT_wall, dρ_rel, du_all, loc = steady_rates(u, p, t, diffmask)
        du = du_all
        c = ComponentArray(copy(u), p.axes)
        push!(history, (t=t, wall_s=t_w, dT_wall=dT_wall, dρ_rel=dρ_rel, loc=string(loc), P_cc=c.cc.P[1],
                        T_evap=c.T_wall.evap[1], Q_out=heat_out(c, p)))
        if verbose
            println("  t = ", lpad(round(t, digits=1), 8), " s:  max|dT_wall/dt| = ",
                    rpad(round(dT_wall, sigdigits=3), 9), " K/s,  質量変化率 = ",
                    rpad(round(dρ_rel, sigdigits=3), 9), " 1/s,  P_cc = ",
                    round(c.cc.P[1], digits=0), " Pa,  T_evap = ", round(c.T_wall.evap[1], digits=3),
                    " K,  Q_out = ", round(history[end].Q_out, digits=3), " W")
            println("      凝縮器入口: x = ", round(c.cond.x[1], sigdigits=3), ", u = ",
                    round(c.cond.u[1], sigdigits=3), " m/s;  凝縮器出口: x = ", round(c.cond.x[end], sigdigits=3),
                    ", u = ", round(c.cond.u[end], sigdigits=3), " m/s;  蒸気管出口 u = ",
                    round(c.vl.u[end], sigdigits=3), " m/s")
            st = sol.stats
            println("      最大: ", loc.T, ", 質量変化 ", loc.ρ, ", セル ", loc.ρ_cell, " (",
                    round(loc.dρ_cell, sigdigits=3), " 1/s)   計算 ", round(t_w, digits=1), " s,  受理 ",
                    st.naccept, ", 棄却 ", st.nreject, ", 収束失敗 ", st.nnonlinconvfail, ", nf ", st.nf)
        end
        if dT_wall < tol_dT && dρ_rel < tol_dρ
            converged = true
            break
        end
    end

    if DOMAIN_ERROR_COUNT[] > 0
        println("注意: 残差の評価で DomainError を ", DOMAIN_ERROR_COUNT[], " 回捕まえ、その試行点を棄却させました")
    end
    verbose && println(converged ? "定常状態に収束しました (t = $(t) s)" :
                       "収束しませんでした (t = $(t) s)。t_max を増やすか、history を確認してください。")
    info = (converged=converged, t=t, dT_wall=isempty(history) ? NaN : history[end].dT_wall,
            dρ_rel=isempty(history) ? NaN : history[end].dρ_rel, history=history)
    return sol, p, u, info
end

# -----------------------------------------------------------------------------
# 定常状態の保存・読み込み (Serialization。素の配列と数値だけを保存する)
# -----------------------------------------------------------------------------

"""
    save_steady_state(path, u, p, condition)

定常状態の状態ベクトル u を保存する。幾何の取り違えを検出するため、変数の並び
(コンポーネント名と長さ) と、計算条件の主な値も一緒に保存する。
"""
function save_steady_state(path::AbstractString, u::AbstractVector, p::AbstractLHPParams,
    condition::LHPTransientCondition)
    d = ComponentArray(copy(u), p.axes)
    layout = [(string(k), [(string(v), length(getproperty(getproperty(d, k), v)))
                           for v in keys(getproperty(d, k))]) for k in keys(d)]
    meta = (Q = condition.Q_apply_init, T_cond = condition.T_cond, T_amb = condition.T_amb,
            T_init = condition.T_init, fluid = condition.fluid)
    mkpath(dirname(abspath(path)))
    Serialization.serialize(path, (u = collect(Float64, u), layout = layout, meta = meta))
    println("定常状態を保存しました: ", abspath(path))
    return nothing
end

"""
    load_steady_state(path, geo) -> Vector{Float64}

save_steady_state で保存した状態ベクトルを読み込む。変数の並びが geo のものと違えば
エラー (幾何・分割数が保存時と違う)。
"""
function load_steady_state(path::AbstractString, geo::AnyLHPGeometry)
    isfile(path) || error("定常状態のファイルが見つかりません: $(path)")
    s = Serialization.deserialize(path)
    d = initialize_state(geo)
    length(getdata(d)) == length(s.u) ||
        error("保存された状態ベクトルの長さ ($(length(s.u))) が、現在の幾何 ($(length(getdata(d)))) と違います")
    layout = [(string(k), [(string(v), length(getproperty(getproperty(d, k), v)))
                           for v in keys(getproperty(d, k))]) for k in keys(d)]
    layout == s.layout ||
        error("保存された定常状態の変数の並びが、現在の幾何と違います (分割数や型を確認してください)")
    println("定常状態を読み込みました: ", abspath(path), "  (Q = ", s.meta.Q, " W, T_cond = ",
            s.meta.T_cond, " K, T_amb = ", s.meta.T_amb, " K)")
    return s.u
end