# =============================================================================
# initialcondition.jl (整合初期条件)
#
# IDA は、初期値 (u0, du0) が、残差 F(t0, u0, du0) = 0 を満たす (整合) ことを
# 要求する。定常状態は求めず、静止 (u=0) の平衡に近い状態を物理的に決めて、
# 代数拘束 (P, T, x の状態方程式、CC出口の流速) が満たされるように代数変数を
# 決め、微分変数の du0 を、残差から逆算する。
#
# 【初期状態】
#   CC         … T_init の気液二相、液体積分率 f_l_cc_init。圧力 P0 は、表の
#                P(ρ,e) をそのまま使う (表と整合した値。厳密な飽和圧力と数十Pa違う)
#   液 (液管、凝縮器 [cond_init_phase=:liquid]) … 過冷却度 ΔT_sub_init の圧縮液
#   蒸気 (グルーブ、蒸気管、[凝縮器 :vapor]) … 過熱度 ΔT_sup_init の蒸気
#   ウィック   … 蒸発界面の温度 T_sat(P_grv) (E は式46)。T_init からの mK のずれが、
#                大きなコンダクタンス G_N (約64 W/K) で増幅され、ṁ≠0 の凝縮から
#                始まってしまうのを避ける
#   壁         … 接する流体と同じ温度 (蒸発器ケースは T_sat(P_grv)、CCケースは T_init)
#   速度       … u=0 (CC出口の u だけは、代数拘束 u=ṁ/(ρ_l A) から求める)
#   圧力       … 全域で P0 (液・蒸気とも、表の P(ρ,e) = P0 となる ρ を二分法で求める)。
#                傾斜・標高差があるときは、静水圧 P0 − ρ g (z − z_cc) をつける。
#
# 液・蒸気を、飽和線 (T_init) の上ではなく、少しずらして置く理由: 飽和線の上は、
# 表の補間誤差のため、液が二相側にはみ出して液の剛性が使えなくなる
# (δ=ρ−ρ_sl(e) が負になる)。過冷却・過熱にして、境界から離す。
#
# 【du0 の求め方】 残差 F は du について厳密に線形 (du は、状態だけに依存する
# 係数を持つ形でしか現れない)。du=0 の残差 F0 と、微分変数の du を1つずつ1にした
# ときの残差から、du についての係数行列 J を作り、微分変数の行について
#     J du_diff = −F0_diff
# を解く。線形関数なので、差分の誤差はなく、検算で F(u0, du0) ≈ 0 を確認する。
# 代数変数の行は、du を含まず、u0 を整合に作ってあるので F0≈0 になっているはず
# (consistency_report で確認する)。
# =============================================================================

using ComponentArrays: getaxes, getdata

# -----------------------------------------------------------------------------
# 1. 表の P(ρ,e) = P_target となる ρ を求める
# -----------------------------------------------------------------------------

"""
    solve_ρ_for_P(tables, e, P_target, ρ_lo, ρ_hi; nscan=400) -> ρ

(ρ,e) の表で、e を固定し、P(ρ,e) = P_target となる ρ を求める。
ρ_lo〜ρ_hi を nscan 点で走査し、P が P_target を (小→大に) またぐ最初の区間を
見つけて二分法で絞る。液・蒸気とも、ρ とともに P が単調に増える範囲を
走査する想定。またぐ区間がなければエラー (ρ の範囲を確認すること)。
"""
function solve_ρ_for_P(tables::FluidTables, e, P_target, ρ_lo, ρ_hi; nscan::Int=400)
    f(ρ) = tables.P_of_ρe(ρ, e) - P_target
    ρs = range(ρ_lo, ρ_hi, length=nscan)
    fprev = f(ρs[1])
    for k in 2:nscan
        fk = f(ρs[k])
        if fprev < 0 <= fk
            a, b = ρs[k-1], ρs[k]
            for _ in 1:200
                m = 0.5 * (a + b)
                f(m) < 0 ? (a = m) : (b = m)
                (b - a) < 1.0e-14 * b && break
            end
            return 0.5 * (a + b)
        end
        fprev = fk
    end
    error("P(ρ,e)=$(P_target) Pa となる ρ が、[$(ρ_lo), $(ρ_hi)] kg/m3 に見つかりません " *
          "(e=$(e) J/kg)。過冷却度・過熱度、温度範囲を確認してください。")
end

# -----------------------------------------------------------------------------
# 2. 流体コンポーネントの初期状態
# -----------------------------------------------------------------------------

"""
    fill_fluid_component!(cmp, p_c, kind, e, ρ_hs, P0, z_ref, g, tables, T_ref)

流体コンポーネント cmp (ρm,u,Em,P,T,x の ComponentVector) を、一様な比内部
エネルギー e、圧力 P0 (+静水圧) の液 (kind=:liquid) または蒸気 (:vapor) で
埋める。u=0。P, T, x は、表の flash (flash_ρe) から取る
(残差の状態拘束と厳密に整合する)。

  ρ_hs   : 静水圧の見積もりに使う密度 (液は飽和液、蒸気は飽和蒸気の密度)
  z_ref  : 圧力 P0 の基準標高 (CC中心)
  T_ref  : 飽和蒸気密度の基準温度 (蒸気の ρ の探索範囲の目安)
"""
function fill_fluid_component!(cmp, p_c, kind::Symbol, e, ρ_hs, P0, z_ref, g,
    tables::FluidTables, T_ref)
    ρ_vs = tables.ρ_v_sat(T_ref)
    for i in 1:p_c.N
        z = cell_elevation(p_c, i)
        P_t = P0 - ρ_hs * g * (z - z_ref)
        ρ = if kind === :liquid
            ρ_ref = tables.ρ_ref(e)     # 飽和液線上の密度 (δ=0)
            solve_ρ_for_P(tables, e, P_t, ρ_ref - 1.0, ρ_ref + 0.9 * last(tables.δ_range))
        else
            solve_ρ_for_P(tables, e, P_t, 0.3 * ρ_vs, 3.0 * ρ_vs)
        end
        st = flash_ρe(ρ, e, tables)
        cmp.ρm[i] = ρ
        cmp.u[i]  = 0.0
        cmp.Em[i] = e + g * z          # Em = em + u²/2 + g z  (u=0)
        cmp.P[i]  = st.P
        cmp.T[i]  = st.T
        cmp.x[i]  = st.x
    end
    return nothing
end

"""
    fill_initial_state!(u0, p::AbstractLHPParams, tables)

u0 (initialize_state の戻り値) を、物理的に決めた初期状態で埋める
(CC出口の流速 u を除く。それは、残差から決める: set_cc_outlet_velocity!)。
"""
function fill_initial_state!(u0, p::AbstractLHPParams, tables::FluidTables)
    c = p.cond
    g = c.g
    T0 = c.T_init
    (c.ΔT_sub_init > 0 && c.ΔT_sup_init > 0) ||
        error("ΔT_sub_init, ΔT_sup_init は正の値にしてください (飽和線の上を避けるため)")
    c.cond_init_phase in (:liquid, :vapor) ||
        error("cond_init_phase は :liquid または :vapor")
    T_liq = T0 - c.ΔT_sub_init
    T_vap = T0 + c.ΔT_sup_init

    # --- CC: T0 の気液二相。圧力 P0 は、表の P(ρ,e) (表と整合) ---
    ρ_l = tables.ρ_l_sat(T0)
    ρ_v = tables.ρ_v_sat(T0)
    f_l = c.f_l_cc_init
    0 < f_l < 1 || error("f_l_cc_init は 0〜1 の間にしてください")
    ρ_cc = f_l * ρ_l + (1 - f_l) * ρ_v
    x_cc = (1 - f_l) * ρ_v / ρ_cc
    Psat = tables.Psat_T(T0)
    e_l = tables.h_l_sat(T0) - Psat / ρ_l
    e_v = tables.h_v_sat(T0) - Psat / ρ_v
    e_cc = e_l + x_cc * (e_v - e_l)
    z_cc = cell_elevation(p.p_cc, 1)
    st_cc = flash_ρe(ρ_cc, e_cc, tables)
    P0 = st_cc.P
    u0.cc.ρm[1] = ρ_cc
    u0.cc.u[1]  = 0.0
    u0.cc.Em[1] = e_cc + g * z_cc
    u0.cc.P[1]  = st_cc.P
    u0.cc.T[1]  = st_cc.T
    u0.cc.x[1]  = st_cc.x

    # --- 液・蒸気の比内部エネルギー (圧縮液・過熱蒸気、CoolProp から) ---
    e_liq = CoolProp.PropsSI("Umass", "T", T_liq, "P", P0, c.fluid)
    e_vap = CoolProp.PropsSI("Umass", "T", T_vap, "P", P0, c.fluid)
    ρ_liq_hs = ρ_l
    ρ_vap_hs = ρ_v

    # --- 各コンポーネント ---
    fill_fluid_component!(u0.grv,  p.p_grv,  :vapor,  e_vap, ρ_vap_hs, P0, z_cc, g, tables, T0)
    fill_fluid_component!(u0.vl,   p.p_vl,   :vapor,  e_vap, ρ_vap_hs, P0, z_cc, g, tables, T0)
    if c.cond_init_phase === :liquid
        fill_fluid_component!(u0.cond, p.p_cond, :liquid, e_liq, ρ_liq_hs, P0, z_cc, g, tables, T0)
    else
        fill_fluid_component!(u0.cond, p.p_cond, :vapor, e_vap, ρ_vap_hs, P0, z_cc, g, tables, T0)
    end
    fill_fluid_component!(u0.ll,   p.p_ll,   :liquid, e_liq, ρ_liq_hs, P0, z_cc, g, tables, T0)

    # --- ウィック: 蒸発界面の温度 T_out = T_sat(P_grv) から始める ---
    # ウィック最外セルと界面の間のコンダクタンス G_N は大きい (約64 W/K) ので、
    # 初期のウィックの温度が T_sat(P_grv) と 1 mK ずれるだけで、0.06 W の熱が流れ、
    # 蒸発量 ṁ が (印加熱 20 W の蒸発量の約0.6% の) 凝縮から始まってしまう
    # (表の補間誤差の mK 程度のずれで起こる)。ウィックと蒸発器ケースを、
    # 界面と同じ温度 T_sat(P_grv) から始めると、ṁ=0 の平衡から始められる。
    # E は式46 (residual_wick! の状態拘束と同じ式)。
    wp = p.wick_p
    T_w = Tsat(groove_pressure(u0.grv), tables)
    e_w = liquid_props(T_w, P0, tables).e
    for i in 1:wp.N
        u0.wick.T[i] = T_w
        u0.wick.E[i] = wp.V[i] * ((1 - wp.por) * wp.ρ_w * wp.C_w * T_w +
                                  wp.por * p.ρ_l_ref * e_w)
    end
    # 平板型: secondary wick (常に液で満たされている) も同じ温度から始める
    if :swick in keys(u0)
        sp = p.swick_p
        u0.swick.T[1] = T_w
        u0.swick.E[1] = sp.V * ((1 - sp.por) * sp.ρ_w * sp.C_w * T_w +
                                sp.por * p.ρ_l_ref * e_w)
    end

    # --- 壁: 接する流体と同じ温度 (蒸発器ケースは界面と同じ T_w、CCケースは T_init) ---
    T_cond_wall = c.cond_init_phase === :liquid ? T_liq : T_vap
    u0.T_wall.evap .= T_w
    u0.T_wall.cc   .= T0
    u0.T_wall.vl   .= T_vap
    u0.T_wall.cond .= T_cond_wall
    u0.T_wall.ll   .= T_liq
    return P0
end

# -----------------------------------------------------------------------------
# 3. CC出口の流速 (代数拘束から決める)
# -----------------------------------------------------------------------------

"""
    set_cc_outlet_velocity!(u0, p, t0)

CC出口の面の流速 u.cc.u[1] は、運動量式ではなく代数拘束
    u = ṁ_out/(ρ_l A)       (residual_cc! の resid.u[1] = u − ṁ_out/(ρ_l A))
で決まる。ṁ_out はウィックの蒸発量で、状態から計算される。残差の u 行を
u=0 で評価すると resid.u[1] = −ṁ_out/(ρ_l A) なので、u = −resid.u[1] を代入すれば、
この拘束が満たされる (ṁ_out は u.cc.u に依存しないので、1回で厳密に決まる)。
"""
function set_cc_outlet_velocity!(u0, p::AbstractLHPParams, t0)
    n = length(getdata(u0))
    res = zeros(n)
    f_dae_residual!(res, zeros(n), getdata(u0), p, t0)
    r = ComponentArray(res, p.axes)
    u0.cc.u[1] -= r.cc.u[1]
    return nothing
end

# -----------------------------------------------------------------------------
# 4. du0 の計算
# -----------------------------------------------------------------------------

"""
    compute_du0(u0vec, p, t0, diffmask) -> (du0, res0)

F(t0, u0, du) が du について線形であることを使って、微分変数の du0 を求める。
  res0 = F(t0, u0, 0)
  J[:,k] = F(t0, u0, e_k) − res0   (k: 微分変数)     … 線形なので誤差なし
  J du_diff = −res0_diff            (微分変数の行について)
res0 は、du=0 での残差 (代数変数の行は、整合に作ってあれば ≈0)。
"""
function compute_du0(u0vec::AbstractVector, p::AbstractLHPParams, t0, diffmask::AbstractVector{Bool})
    n = length(u0vec)
    D = findall(diffmask)
    du = zeros(n)
    res0 = zeros(n)
    f_dae_residual!(res0, du, u0vec, p, t0)

    J = zeros(length(D), length(D))
    res = zeros(n)
    for (k, j) in enumerate(D)
        du[j] = 1.0
        f_dae_residual!(res, du, u0vec, p, t0)
        du[j] = 0.0
        J[:, k] = res[D] .- res0[D]
    end

    du0 = zeros(n)
    du0[D] = -(J \ res0[D])
    return du0, res0
end

# -----------------------------------------------------------------------------
# 5. 整合初期条件の生成 (エントリポイント)
# -----------------------------------------------------------------------------

"""
    initial_condition(geo, condition, tables; verbose=true) -> (u0, du0, p)

整合初期条件を生成する。戻り値:
  u0  … ComponentArray (初期値)
  du0 … Vector{Float64} (u0 と同じ並びの、整合な du0)
  p   … LHPParams (f_dae_residual! に渡すパラメータ)

手順: (1) 物理的に初期状態を決め、代数変数を表の flash と整合させる
      (2) CC出口の流速を、代数拘束から決める
      (3) du0 を、残差が du について線形であることを使って求める
      (4) F(t0, u0, du0) を評価して、整合性を確認する (consistency_report)
"""
function initial_condition(geo::AnyLHPGeometry, condition::LHPTransientCondition,
    tables::FluidTables; verbose::Bool=true)

    u0 = initialize_state(geo)
    p = build_params(geo, condition, tables, getaxes(u0))
    t0 = condition.t_init

    P0 = fill_initial_state!(u0, p, tables)
    set_cc_outlet_velocity!(u0, p, t0)

    diffmask = Bool.(getdata(initialize_differential_vars(geo)))
    du0, res0 = compute_du0(getdata(u0), p, t0, diffmask)

    # 検算: F(t0, u0, du0)
    res = zeros(length(res0))
    f_dae_residual!(res, du0, getdata(u0), p, t0)
    if verbose
        println("初期圧力 P0 = ", P0, " Pa")
        consistency_report(res0, res, p.axes, diffmask)
    end
    return u0, du0, p
end

"""
    consistency_report(res0, res, axes, diffmask)

整合性の検算を表示する。res0 は du=0 での残差、res は du=du0 での残差。
  - 代数変数の行: res0 と res の両方が ≈0 であること (du を含まないので同じ値)
  - 微分変数の行: res が ≈0 (res0 は、du=0 のときの不整合の大きさで、
                  それが du0 で打ち消される)
"""
function consistency_report(res0, res, axes, diffmask)
    n = length(res)
    alg = .!diffmask
    println("整合性の検算:")
    println("  代数変数の行: max|F| = ", maximum(abs, res[alg]),
            "  (≈0 のはず。P[Pa], T[K], x[-], u[m/s] の各拘束の残差)")
    println("  微分変数の行: max|F(u0, 0)| = ", maximum(abs, res0[diffmask]),
            " → max|F(u0, du0)| = ", maximum(abs, res[diffmask]))
    r0 = ComponentArray(res0, axes)
    r  = ComponentArray(res, axes)
    println("  ブロックごとの max|F(u0, du0)| (微分変数のもの) と、max|F(u0, 0)|:")
    for comp in keys(r)
        sub = getproperty(r, comp)
        sub0 = getproperty(r0, comp)
        for var in keys(sub)
            println("    ", rpad("$(comp).$(var)", 14), " ",
                    lpad(round(maximum(abs, getproperty(sub, var)), sigdigits=3), 10),
                    "   (du=0: ", round(maximum(abs, getproperty(sub0, var)), sigdigits=3), ")")
        end
    end
end