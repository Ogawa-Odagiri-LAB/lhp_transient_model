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
#   壁         … 接する流体と同じ温度 (蒸発器ケースは T_sat(P_grv)、CCケースは T_init)。
#                T_wall_init_* を指定した壁は、その温度から始める (実験の初期の壁温度)
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
function fill_fluid_component!(cmp, p_c, kind, e, ρ_hs, P0, z_ref, g,
    tables::FluidTables, T_ref)
    ρ_vs = tables.ρ_v_sat(T_ref)
    for i in 1:p_c.N
        # kind, e はセルごとのベクトルでもよい (液と蒸気が混在する管: initial_from_measurements.jl)
        kind_i = kind isa AbstractVector ? kind[i] : kind
        e_i    = e isa AbstractVector ? e[i] : e
        ρ_hs_i = ρ_hs isa NamedTuple ? (kind_i === :liquid ? ρ_hs.liquid : ρ_hs.vapor) : ρ_hs
        z = cell_elevation(p_c, i)
        P_t = P0 - ρ_hs_i * g * (z - z_ref)
        ρ = if kind_i === :liquid
            ρ_ref = tables.ρ_ref(e_i)     # 飽和液線上の密度 (δ=0)
            solve_ρ_for_P(tables, e_i, P_t, ρ_ref - 1.0, ρ_ref + 0.9 * last(tables.δ_range))
        else
            solve_ρ_for_P(tables, e_i, P_t, 0.3 * ρ_vs, 3.0 * ρ_vs)
        end
        st = flash_ρe(ρ, e_i, tables)
        cmp.ρm[i] = ρ
        cmp.u[i]  = 0.0
        cmp.Em[i] = e_i + g * z          # Em = em + u²/2 + g z  (u=0)
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
    if :ccl in keys(u0)
        # 平板型の新しい CC (気相 1 セル + 液 16 層): 全層 T0、名目の空き容積 X = V (1 − f_l)、気相の質量 mv = ρ_v,sat(T0) X。
        # 圧力は、CC の圧力の閉じた式 (cc_solve_P) から。要約 (cc.ρm, Em, T, x) は残差と同じ式 (cc_summary_values)。u は set_cc_outlet_velocity! で決める。
        pc = p.p_cc
        u0.ccl.Tl .= T0
        u0.ccl.X[1]  = pc.V * (1 - f_l)
        u0.ccl.mv[1] = ρ_v * u0.ccl.X[1]
        P0 = cc_solve_P(u0.ccl.Tl, u0.ccl.X[1], u0.ccl.mv[1], pc, tables)
        sm = cc_summary_values(P0, u0.ccl.Tl, u0.ccl.X[1], u0.ccl.mv[1], pc, g, tables)
        u0.cc.ρm[1] = sm.ρm
        u0.cc.u[1]  = 0.0
        u0.cc.Em[1] = sm.Em
        u0.cc.P[1]  = P0
        u0.cc.T[1]  = sm.T
        u0.cc.x[1]  = sm.x
    else
        st_cc = flash_ρe(ρ_cc, e_cc, tables)
        P0 = st_cc.P
        u0.cc.ρm[1] = ρ_cc
        u0.cc.u[1]  = 0.0
        u0.cc.Em[1] = e_cc + g * z_cc
        u0.cc.P[1]  = st_cc.P
        u0.cc.T[1]  = st_cc.T
        u0.cc.x[1]  = st_cc.x
    end

    # --- 液・蒸気の温度と比内部エネルギー (圧縮液・過熱蒸気、CoolProp から) ---
    # 壁の初期温度 T_wall_init_* を指定した管は、流体の温度も壁に合わせる (壁と流体の温度差が
    # あると、初期に、壁-流体の熱流束が桁違いに大きくなる。特に蒸気管の壁が飽和温度付近だと、
    # 凝縮側の熱伝達率が効いて、蒸気が一気に凝縮する)。ただし、飽和線から離す制約を守る:
    #   蒸気 ≥ T_init + ΔT_sup_init、液 ≤ T_init − ΔT_sub_init
    # (制約でずれる分だけ、壁と流体に温度差が残る)。
    _sc(x) = x isa AbstractVector ? sum(x) / length(x) : x    # セルごとのベクトルは、平均 (位相ベクトル未指定のときの代表値)
    T_vap_vl   = c.T_wall_init_vl   === nothing ? T_vap : max(_sc(c.T_wall_init_vl), T_vap)
    T_liq_ll   = c.T_wall_init_ll   === nothing ? T_liq : min(_sc(c.T_wall_init_ll), T_liq)
    T_cond_fl  = if c.cond_init_phase === :liquid
        c.T_wall_init_cond === nothing ? T_liq : min(_sc(c.T_wall_init_cond), T_liq)
    else
        c.T_wall_init_cond === nothing ? T_vap : max(_sc(c.T_wall_init_cond), T_vap)
    end
    e_vap    = CoolProp.PropsSI("Umass", "T", T_vap,    "P", P0, c.fluid)   # グルーブ
    e_vap_vl = CoolProp.PropsSI("Umass", "T", T_vap_vl, "P", P0, c.fluid)
    e_liq_ll = CoolProp.PropsSI("Umass", "T", T_liq_ll, "P", P0, c.fluid)
    e_cond   = CoolProp.PropsSI("Umass", "T", T_cond_fl, "P", P0, c.fluid)
    ρ_liq_hs = ρ_l
    ρ_vap_hs = ρ_v

    # --- 各コンポーネント ---
    # グルーブ: 既定は蒸気 (過熱 ΔT_sup_init)。grv_init=:liquid なら満水 (過冷却 ΔT_sub_init の液)
    c.grv_init in (:vapor, :liquid) || error("grv_init は :vapor または :liquid")
    if c.grv_init === :liquid
        e_liq_grv = CoolProp.PropsSI("Umass", "T", T_liq, "P", P0, c.fluid)
        fill_fluid_component!(u0.grv, p.p_grv, :liquid, e_liq_grv, ρ_liq_hs, P0, z_cc, g, tables, T0)
    else
        fill_fluid_component!(u0.grv, p.p_grv, :vapor, e_vap, ρ_vap_hs, P0, z_cc, g, tables, T0)
    end
    # グルーブに少量の液を入れる (grv_init=:vapor かつ f_l_grv_init > 0): 飽和温度 T0 の気液二相 (CC と同じ置き方)。
    # 核沸騰が始まると、蒸気の発生と同時に、液が出口から押し出される (f_full_grv を超えた分)。
    if c.grv_init === :vapor && c.f_l_grv_init > 0
        c.f_l_grv_init < 1 || error("f_l_grv_init は 0〜1")
        fg = c.f_l_grv_init
        ρ_g = fg * ρ_l + (1 - fg) * ρ_v
        x_g = (1 - fg) * ρ_v / ρ_g
        e_g = e_l + x_g * (e_v - e_l)
        for i in 1:p.p_grv.N
            st_g = flash_ρe(ρ_g, e_g, tables)
            u0.grv.ρm[i] = ρ_g;  u0.grv.u[i] = 0.0
            u0.grv.Em[i] = e_g + g * cell_elevation(p.p_grv, i)
            u0.grv.P[i]  = st_g.P;  u0.grv.T[i] = st_g.T;  u0.grv.x[i] = st_g.x
        end
        println("  グルーブの初期液体積分率 f_l = ", fg, " (x = ", round(x_g, digits = 5), ",  T = ",
                round(u0.grv.T[1] - 273.15, digits = 3), " °C,  P = ", round(u0.grv.P[1] / 1e3, digits = 2), " kPa)")
    end
    # 位相ベクトル (init_phase_*) が与えられた管は、セルごとに液/蒸気を置く。流体の温度は壁に合わせ、
    # 飽和線から離す (液 ≤ T_init − ΔT_sub_init、蒸気 ≥ T_init + ΔT_sup_init)。
    ρ_hs_mix = (liquid = ρ_liq_hs, vapor = ρ_vap_hs)
    function fill_mixed!(cmp, p_c, phase, T_wall_init, name)
        length(phase) == p_c.N || error("init_phase_$(name) の長さ $(length(phase)) ≠ N=$(p_c.N)")
        Tw = T_wall_init === nothing ? fill(T0, p_c.N) :
             (T_wall_init isa AbstractVector ? T_wall_init : fill(T_wall_init, p_c.N))
        # 蒸気の e: 表の飽和蒸気の e_v (CC と同じ式) を基準に、CoolProp の過熱分を足す。
        # CoolProp の Umass をそのまま使うと、表の飽和蒸気線とのずれ (約 2 kJ/kg ≈ 0.8 K 相当) で、
        # 過熱 0.1 K の蒸気が表では二相 (x=0.998、T が 0.7 K 高い) になってしまう。
        ρ_v0 = tables.ρ_v_sat(T0)
        e_v_tab = tables.h_v_sat(T0) - tables.Psat_T(T0) / ρ_v0
        e_v_cp  = CoolProp.PropsSI("Umass", "T", T0, "Q", 1, c.fluid)
        # 蒸気: 表の補間誤差のため、飽和線のすぐ上 (過熱 0.1 K) の蒸気は、表では二相 (x<1) になり、
        # flash の T が 0.7 K 以上ずれる。表で単相蒸気 (x ≥ 0.9999) になる最小の過熱度まで、
        # 目標温度に上乗せして探す。
        function vapor_e(T_target)
            e_last = 0.0
            # まず、飽和蒸気そのもの (CC の x = 1 と同じ表の e_v) とそのすぐ上 (+0, +100, +250 J/kg ≈ 過熱 0〜0.15 K) を試す。
            # 以前は T_target (飽和 + 0.1 K) から探し始めて、表の誤差のため過熱 0.8 K の蒸気 (T = 飽和 + 0.87 K) になった。
            # 0.8 K 相当 (約 2 kJ/kg) より小さい過熱で単相になる最小の δe を、50 J/kg 刻みで探す。
            for δe in 0.0:50.0:3000.0
                e = e_v_tab + δe
                ρ = try
                    solve_ρ_for_P(tables, e, P0, 0.3 * ρ_v0, 3.0 * ρ_v0)
                catch
                    continue
                end
                if flash_ρe(ρ, e, tables).x >= 0.9999
                    println("  蒸気セルの初期状態 (", name, "): 飽和蒸気の e に +", δe, " J/kg で単相蒸気")
                    return e
                end
            end
            for extra in (0.0, 0.1, 0.2, 0.4, 0.8, 1.5, 3.0)
                e = e_v_tab + (CoolProp.PropsSI("Umass", "T", T_target + extra, "P", P0, c.fluid) - e_v_cp)
                e_last = e
                ρ = try
                    solve_ρ_for_P(tables, e, P0, 0.3 * ρ_v0, 3.0 * ρ_v0)
                catch
                    continue
                end
                flash_ρe(ρ, e, tables).x >= 0.9999 && return e
            end
            return e_last
        end
        # 飽和に近い液 (:sat): 過冷却度 ΔT_sub_ll_sat の液。表の P(ρ,e) で単相の液 (x ≤ 1e-6) にならなければ、過冷却度を大きくする。
        function liquid_e_sat(Tw_i)
            for dT in sort(unique(filter(<=(c.ΔT_sub_init), [c.ΔT_sub_ll_sat, 0.04, 0.06, c.ΔT_sub_init])))
                e = CoolProp.PropsSI("Umass", "T", min(Tw_i, T0 - dT), "P", P0, c.fluid)
                ρ_r = tables.ρ_ref(e)
                ok = try
                    ρ = solve_ρ_for_P(tables, e, P0, ρ_r - 1.0, ρ_r + 0.9 * last(tables.δ_range))
                    flash_ρe(ρ, e, tables).x <= 1.0e-6
                catch
                    false
                end
                ok && return e
            end
            return CoolProp.PropsSI("Umass", "T", min(Tw_i, T_liq), "P", P0, c.fluid)
        end
        e_cells = [phase[i] === :vapor ? vapor_e(max(Tw[i], T_vap)) :
                   phase[i] === :sat   ? liquid_e_sat(Tw[i]) :
                   CoolProp.PropsSI("Umass", "T", min(Tw[i], T_liq), "P", P0, c.fluid)
                   for i in 1:p_c.N]
        kinds = [phase[i] === :vapor ? :vapor : :liquid for i in 1:p_c.N]
        fill_fluid_component!(cmp, p_c, kinds, e_cells, ρ_hs_mix, P0, z_cc, g, tables, T0)
        ns = count(==(:sat), phase)
        ns > 0 && println("  飽和に近い液のセル (", name, "): ", ns, " 個 (セル ", findfirst(==(:sat), phase), "〜", findlast(==(:sat), phase),
                          ")、過冷却度 ", round(minimum(T0 - cmp.T[i] for i in eachindex(phase) if phase[i] === :sat), digits = 3), "〜",
                          round(maximum(T0 - cmp.T[i] for i in eachindex(phase) if phase[i] === :sat), digits = 3), " K、x 最大 ",
                          maximum(cmp.x[i] for i in eachindex(phase) if phase[i] === :sat))
    end
    if c.init_phase_vl !== nothing
        fill_mixed!(u0.vl, p.p_vl, c.init_phase_vl, c.T_wall_init_vl, "vl")
    else
        fill_fluid_component!(u0.vl, p.p_vl, :vapor, e_vap_vl, ρ_vap_hs, P0, z_cc, g, tables, T0)
    end
    if c.init_phase_cond !== nothing
        fill_mixed!(u0.cond, p.p_cond, c.init_phase_cond, c.T_wall_init_cond, "cond")
    elseif c.cond_init_phase === :liquid
        fill_fluid_component!(u0.cond, p.p_cond, :liquid, e_cond, ρ_liq_hs, P0, z_cc, g, tables, T0)
    else
        fill_fluid_component!(u0.cond, p.p_cond, :vapor, e_cond, ρ_vap_hs, P0, z_cc, g, tables, T0)
    end
    if c.init_phase_ll !== nothing
        fill_mixed!(u0.ll, p.p_ll, c.init_phase_ll, c.T_wall_init_ll, "ll")
    else
        fill_fluid_component!(u0.ll, p.p_ll, :liquid, e_liq_ll, ρ_liq_hs, P0, z_cc, g, tables, T0)
    end
    # 液管の CC 側の末端を、飽和の気液二相にする (蒸気の体積分率 f_v_ll_init)。CC と同じ置き方 (T0 の飽和、P は表の flash)。
    if c.f_v_ll_init > 0
        0 < c.f_v_ll_init < 1 || error("f_v_ll_init は 0〜1")
        f_l_ll = 1 - c.f_v_ll_init
        ρ_ll2 = f_l_ll * ρ_l + c.f_v_ll_init * ρ_v
        x_ll2 = c.f_v_ll_init * ρ_v / ρ_ll2
        e_ll2 = e_l + x_ll2 * (e_v - e_l)
        N_ll = p.p_ll.N
        for i in max(1, N_ll - c.n_ll_2ph + 1):N_ll
            st_l2 = flash_ρe(ρ_ll2, e_ll2, tables)
            u0.ll.ρm[i] = ρ_ll2;  u0.ll.u[i] = 0.0
            u0.ll.Em[i] = e_ll2 + g * cell_elevation(p.p_ll, i)
            u0.ll.P[i]  = st_l2.P;  u0.ll.T[i] = st_l2.T;  u0.ll.x[i] = st_l2.x
        end
        println("  液管の末端 ", c.n_ll_2ph, " セルを二相に: f_v = ", c.f_v_ll_init, " (x = ", round(x_ll2, digits = 5),
                ",  T = ", round(u0.ll.T[N_ll] - 273.15, digits = 3), " °C,  P = ", round(u0.ll.P[N_ll] / 1e3, digits = 2), " kPa)")
    end

    # 蒸気管の初期蒸気セルを、飽和の気液二相 (乾き度 x_vap_vl_init < 1) にする。
    # 表の限界で、単相蒸気は飽和から約 0.7 K 上でしか始められない (実験の VL1 は飽和に近い 22.6 °C)。
    # この過熱蒸気が前線の二相セルに触れて、起動 1 s 後に前線セルが一気に凝縮し、+1.7 kPa の圧力の跳ねが出た。
    if c.x_vap_vl_init < 1 && c.init_phase_vl !== nothing
        0 < c.x_vap_vl_init || error("x_vap_vl_init は 0〜1")
        x_v2 = c.x_vap_vl_init
        ρ_v2 = 1 / (x_v2 / ρ_v + (1 - x_v2) / ρ_l)
        e_v2 = e_l + x_v2 * (e_v - e_l)
        # 表の P(ρ,e) の誤差で、(ρ_v2, e_v2) の flash の P が P0 から数 kPa ずれる (x=0.9995 で +1.8 kPa)。
        # ずれた圧力で始めると、起動後 5 s の過渡で蒸気が凝縮し、液管の二相も潰れる。e を固定して、表の P が P0 になる ρ を探す。
        ρ_v2 = try
            solve_ρ_for_P(tables, e_v2, P0, 0.3 * ρ_v, 3.0 * ρ_v)
        catch
            ρ_v2
        end
        st_v2 = flash_ρe(ρ_v2, e_v2, tables)
        for i in eachindex(c.init_phase_vl)
            c.init_phase_vl[i] === :vapor || continue
            u0.vl.ρm[i] = ρ_v2;  u0.vl.u[i] = 0.0
            u0.vl.Em[i] = e_v2 + g * cell_elevation(p.p_vl, i)
            u0.vl.P[i]  = st_v2.P;  u0.vl.T[i] = st_v2.T;  u0.vl.x[i] = st_v2.x
        end
        println("  蒸気管の蒸気セルを二相に: x = ", x_v2, " (flash: x = ", round(st_v2.x, digits = 5),
                ",  T = ", round(st_v2.T - 273.15, digits = 3), " °C (T_init = ", round(T0 - 273.15, digits = 3),
                "),  P = ", round(st_v2.P / 1e3, digits = 2), " kPa,  ρ = ", round(ρ_v2, digits = 3), ")")
    end

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
        for (i, V_i) in enumerate((sp.V_A, sp.V_B))     # 節点 A、B
            u0.swick.T[i] = T_w
            u0.swick.E[i] = V_i * ((1 - sp.por) * sp.ρ_w * sp.C_w * T_w +
                                   sp.por * p.ρ_l_ref * e_w)
        end
    end

    # --- 壁: 接する流体と同じ温度 (蒸発器ケースは界面と同じ T_w、CCケースは T_init) ---
    u0.T_wall.evap .= T_w
    u0.T_wall.cc   .= T0
    :hb in keys(u0.T_wall) && (u0.T_wall.hb .= c.T_wall_init_hb === nothing ? T_w : c.T_wall_init_hb)
    u0.T_wall.vl   .= T_vap_vl
    u0.T_wall.cond .= T_cond_fl
    u0.T_wall.ll   .= T_liq_ll
    :cplate in keys(u0) && (u0.cplate.T .= c.T_cond)    # コールドプレート (2次元) も T_cond (冷却液の温度) から始める
    :plate in keys(u0) && (u0.plate.T .= c.T_cond)    # 凝縮器プレートはコールドプレートの温度から始める (実測は使わない)
    # 壁の初期温度の指定 (LHPTransientCondition の T_wall_init_*) があれば、上書きする
    for (name, T_set) in ((:evap, c.T_wall_init_evap), (:cc, c.T_wall_init_cc),
                          (:vl, c.T_wall_init_vl), (:cond, c.T_wall_init_cond),
                          (:ll, c.T_wall_init_ll))
        T_set === nothing || (getproperty(u0.T_wall, name) .= T_set)
    end
    :cca in keys(u0) && c.cc_2node && (u0.cca.T .= u0.T_wall.cc[1])    # 2 節点の CC ケース: 節点 A も B と同じ温度から始める
    :eca in keys(u0) && c.ec_2node && (u0.eca.T .= u0.T_wall.evap[1])  # 2 節点の蒸発器ケース: 側面の節点 S も P と同じ温度から始める
    :hi in keys(u0) && init_housing!(u0, c)    # 断熱材・PEEK は、壁と周囲の間の定常の温度から始める
    # 位相ベクトルを与えた管の蒸気セル: 壁を、表の flash が返した流体の温度に合わせる。
    # 飽和線の近く (過熱 0.1 K) の蒸気は、表 P(ρ,e) の補間誤差のため、flash の T が目標より
    # 数百 mK ずれる (実測: 目標 22.7°C に対し 23.4°C)。壁と流体に温度差があると、蒸気側の
    # 熱伝達 (凝縮) が大きく、起動直後の刻みが潰れるので、壁を流体に合わせる。
    for (cmp, ph, wname) in ((u0.vl, c.init_phase_vl, :vl), (u0.cond, c.init_phase_cond, :cond))
        ph === nothing && continue
        Tw = getproperty(u0.T_wall, wname)
        for i in eachindex(ph)
            ph[i] === :vapor && (Tw[i] = cmp.T[i])
        end
        nv = count(==(:vapor), ph)
        nv > 0 && println("  蒸気セルの初期温度 (", wname, "): 最大 ", round(maximum(cmp.T[i] for i in eachindex(ph) if ph[i] === :vapor) - T0, digits = 3),
                          " K (飽和温度 T_init からの差。表の誤差で 0.8 K 以上になることがある)")
    end
    # 液管の飽和に近い液セル (:sat) と末端の二相セル: 流体の方が壁より最大 0.065 K 低い (または二相は壁が高い)。
    # 二相セルは壁との熱伝達が大きく、飽和液セルも h≈6000 W/m²K 相当で 0.4 K/s 以上の速さで加熱されて
    # 体積膨張し、末端の圧力が +1.4 kPa 跳ねて二相 (蒸気) が起動後 5 s で潰れる (du=0 残差 ll.Em[18..20] が 1e6 級)。
    # 壁を流体の温度に合わせて、起動直後の過渡をなくす (差は 0.07 K 以下で、センサの不確かさの範囲内)。
    if c.init_phase_ll !== nothing
        for i in eachindex(c.init_phase_ll)
            c.init_phase_ll[i] === :sat && (u0.T_wall.ll[i] = u0.ll.T[i])
        end
    end
    if c.f_v_ll_init > 0
        for i in max(1, p.p_ll.N - c.n_ll_2ph + 1):p.p_ll.N
            u0.T_wall.ll[i] = u0.ll.T[i]
        end
    end
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
u_init を渡すと、手順 (1)(2) の代わりに、その状態ベクトルを u0 とする
(定常計算の結果から過渡計算を始めるときなど)。

手順: (1) 物理的に初期状態を決め、代数変数を表の flash と整合させる
      (2) CC出口の流速を、代数拘束から決める
      (3) du0 を、残差が du について線形であることを使って求める
      (4) F(t0, u0, du0) を評価して、整合性を確認する (consistency_report)
"""
function initial_condition(geo::AnyLHPGeometry, condition::LHPTransientCondition,
    tables::FluidTables; verbose::Bool=true, u_init::Union{Nothing,AbstractVector}=nothing)

    u0 = initialize_state(geo)
    p = build_params(geo, condition, tables, getaxes(u0))
    t0 = condition.t_init

    if u_init === nothing
        P0 = fill_initial_state!(u0, p, tables)
        set_cc_outlet_velocity!(u0, p, t0)
    else
        # 定常計算 (steadyfuncs.jl) などで得た状態ベクトルを、そのまま初期値にする。
        # du0 は、下で、この状態から残差を使って逆算する (整合初期条件)。
        length(u_init) == length(getdata(u0)) ||
            error("u_init の長さ ($(length(u_init))) が、この幾何の状態ベクトルの長さ " *
                  "($(length(getdata(u0)))) と違います。保存時と同じ幾何・分割数か確認してください。")
        getdata(u0) .= u_init
        P0 = u0.cc.P[1]
    end

    diffmask = Bool.(getdata(initialize_differential_vars(geo)))
    du0, res0 = compute_du0(getdata(u0), p, t0, diffmask)

    # 検算: F(t0, u0, du0)
    res = zeros(length(res0))
    f_dae_residual!(res, du0, getdata(u0), p, t0)
    if verbose
        println("初期圧力 P0 = ", P0, " Pa")
        consistency_report(res0, res, p.axes, diffmask)
        # 診断: du=0 の残差が大きい変数の上位と、蒸気管のセルごとの状態
        labels = variable_labels(p.axes, length(res0))
        println("  du=0 の残差の大きい変数 (上位 10):")
        for j in sortperm(abs.(res0), rev = true)[1:min(10, length(res0))]
            println("    ", rpad(labels[j], 18), " ", round(res0[j], sigdigits = 4))
        end
        if hasproperty(u0, :vl)
            println("  蒸気管のセル: i, T_fluid[°C], x, ρ[kg/m3], T_wall[°C], F_wall(du=0)[W], F_Em(du=0)")
            r0 = ComponentArray(res0, p.axes)
            for i in 1:p.p_vl.N
                println("    ", lpad(i, 2), "  ", round(u0.vl.T[i] - 273.15, digits = 3), "  ", round(u0.vl.x[i], digits = 3),
                        "  ", round(u0.vl.ρm[i], digits = 2), "  ", round(u0.T_wall.vl[i] - 273.15, digits = 3),
                        "  ", round(r0.T_wall.vl[i], sigdigits = 3), "  ", round(r0.vl.Em[i], sigdigits = 3))
            end
        end
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