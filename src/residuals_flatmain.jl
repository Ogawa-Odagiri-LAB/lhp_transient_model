# =============================================================================
# residuals_flatmain.jl (平板型の全体の残差 f_dae_residual!)
#
# 円筒型の residuals_main.jl と同じ構成で、平板型蒸発器 (secondary wick 付き) 用の
# メソッド。残差関数の名前は同じ (f_dae_residual!) で、p の型 (FlatLHPParams) で
# 切り替わる。initialcondition.jl, sparsity.jl, solverfuncs.jl, postprocess.jl,
# savefuncs.jl は、AbstractLHPParams に対して書いてあるので、そのまま使える。
#
# 【円筒型との違い】
#   流体:   grv → vl → cond → ll → cc →(ṁ)→ secondary wick →(ṁ)→ primary wick →(蒸発)→ grv
#           (液管はCCに直接入る。CCから出る液 ṁ は secondary wick を通って primary wick へ)
#   共有する量 (1回だけ計算して複数の残差に配る):
#     Q_nc   = flatcc_natural_convection の Q_fluid : secondary wick → CC流体
#              (CC流体の残差と secondary wick の残差で共有。CCケース → 流体 の対流は無い)
#     Q_cs   = Q_case_swick : CCケース → secondary wick
#              (CCケースの残差と secondary wick の残差で共有)
#     Q_leak = residual_wick! の戻り値 : primary wick → secondary wick
#              (primary wick の第1層の残差と secondary wick の残差で共有)
#     Q_ec_cc = Q_ec_wall(T_ec, ec_h, T_ccw, cw_h) : 蒸発器ケース → CCケース
#   その他 (iface, Q_ec_vl, Q_ec_grv, h_conv, ε_u) は円筒型と同じ。
#   CC流体の Q_leak (ウィック内表面 → CC) は無い (primary wick は CC に直接接しない) ので 0。
# =============================================================================

using ComponentArrays: ComponentArray, getaxes

"""
    FlatLHPParams

LHPParams (円筒型) と同じフィールド名に、secondary wick の量を加えたもの。
"""
struct FlatLHPParams{A,T} <: AbstractLHPParams
    geo::FlatLHPGeometry
    cond::LHPTransientCondition
    tables::T
    axes::A
    p_grv::Any
    p_vl::Any
    p_cond::Any
    p_ll::Any
    p_cc::Any
    wp_vl::Any
    wp_cond::Any
    wp_ll::Any
    cw::Any                  # CCケース (半円筒の壁)
    ec::Any                  # 蒸発器ケース (蒸気管との熱伝導用)
    ec_h::Any                # 蒸発器ケース → CCケースの経路 (蒸発器ケース側)
    cw_h::Any                # 同 (CCケース側)
    wick_p::Any              # primary wick
    swick_p::Any             # secondary wick (V, k_eff, ... と CCケースとの熱伝導の量)
    ρ_l_ref::Float64
end

"""
    build_params(geo::FlatLHPGeometry, condition, tables, axes) -> FlatLHPParams

k_s (SUS304) と k_l (液) は、基準温度 T_init での値で、primary / secondary wick の
k_eff (マクスウェル) に使う (円筒型と同じ扱い)。
"""
function build_params(geo::FlatLHPGeometry, condition::LHPTransientCondition,
    tables::FluidTables, axes)
    k_s = k_sus304(condition.T_init)
    k_l = tables.k_l_sat(condition.T_init)
    ec_h, cw_h = flat_ec_cc_params(geo.grv, geo.cc)
    return FlatLHPParams(
        geo, condition, tables, axes,
        grv_params(geo.grv), pipe_params(geo.vl), pipe_params(geo.cond),
        pipe_params(geo.ll), flatcc_params(geo.cc),
        wall_params(geo.vl), wall_params(geo.cond), wall_params(geo.ll),
        flatcc_wall_params(geo.cc), flat_evap_case_params(geo.grv), ec_h, cw_h,
        flatwick_params(geo.wick, geo.cc, geo.swick, k_s, k_l),
        swick_params(geo.swick, geo.cc, k_s, k_l),
        tables.ρ_l_sat(condition.T_init),
    )
end

"""
    f_dae_residual!(res, du, u, p::FlatLHPParams, t)

平板型の全体の残差 (円筒型の f_dae_residual! と同じ形。上の説明を参照)。
"""
function f_dae_residual!(res, du, u, p::FlatLHPParams, t)
    r  = ComponentArray(res, p.axes)
    d  = ComponentArray(du,  p.axes)
    x  = ComponentArray(u,   p.axes)
    geo, c, tables = p.geo, p.cond, p.tables
    g = c.g

    # (1) 共有する量
    T_ec  = x.T_wall.evap[1]
    T_ccw = x.T_wall.cc[1]
    T_cc  = x.cc.T[1]
    P_cc  = x.cc.P[1]
    T_sw  = x.swick.T[1]
    P_grv = groove_pressure(x.grv)

    # 蒸発界面: primary wick の最外層 (蒸発器側) と界面。ṁ は secondary wick を通ってくる
    iface = wick_interface(x.wick.T[end], P_grv, T_ec, c.h_ev, geo.grv.A_ec_wick,
                           p.wick_p.G[end], P_cc, tables;
                           T_cc = T_cc, ΔT_startup = c.ΔT_startup, ε_startup = c.ε_startup)
    S_m_total, S_E_total = wick_groove_source(iface)

    Q_ec_vl = Q_ec_wall(T_ec, p.ec, x.T_wall.vl[1], p.wp_vl)
    Q_ec_cc = Q_ec_wall(T_ec, p.ec_h, T_ccw, p.cw_h)
    h_grv    = groove_h_ec(x.grv, p.p_grv, iface.ṁ, g, tables)
    Q_ec_grv = groove_evap_heat(x.grv, T_ec, h_grv, p.p_grv)

    up_vl   = transport_ghost_up(x.grv,  p.p_grv,  p.p_vl)
    up_cond = transport_ghost_up(x.vl,   p.p_vl,   p.p_cond)
    up_ll   = transport_ghost_up(x.cond, p.p_cond, p.p_ll)
    h_vl   = zeros(p.p_vl.N);   compute_h_conv_transport!(h_vl,   x.vl,   up_vl,   x.T_wall.vl,   p.p_vl,   g, tables; ε_Re = 200.0)
    h_cond = zeros(p.p_cond.N); compute_h_conv_transport!(h_cond, x.cond, up_cond, x.T_wall.cond, p.p_cond, g, tables; ε_Re = 200.0)
    h_ll   = zeros(p.p_ll.N);   compute_h_conv_transport!(h_ll,   x.ll,   up_ll,   x.T_wall.ll,   p.p_ll,   g, tables; ε_Re = 200.0)

    # CC: secondary wick → 流体の自然対流、CCケース → secondary wick の熱伝導、液管からの流入
    nc    = flatcc_natural_convection(x.cc, T_sw, p.p_cc, g, tables)
    Q_cs  = Q_case_swick(T_ccw, T_sw, p.swick_p)
    inlet = cc_inlet(x.ll, p.p_ll, x.cc; ε_u = c.ε_u)

    # (2) 流体の残差
    down_grv = transport_ghost_down(x.vl, p.p_vl, p.p_grv)
    residual_grv!(r.grv, d.grv, x.grv, down_grv, Q_ec_grv, S_m_total, S_E_total,
                  p.p_grv, g, tables; ε_u = c.ε_u)

    down_vl = transport_ghost_down(x.cond, p.p_cond, p.p_vl)
    residual_transport_line!(r.vl, d.vl, x.vl, up_vl, down_vl, x.T_wall.vl, h_vl,
                             p.p_vl, g, tables; ε_u = c.ε_u)

    down_cond = transport_ghost_down(x.ll, p.p_ll, p.p_cond)
    residual_transport_line!(r.cond, d.cond, x.cond, up_cond, down_cond, x.T_wall.cond, h_cond,
                             p.p_cond, g, tables; ε_u = c.ε_u)

    down_ll = transport_ghost_down(x.cc, p.p_cc, p.p_ll)
    residual_transport_line!(r.ll, d.ll, x.ll, up_ll, down_ll, x.T_wall.ll, h_ll,
                             p.p_ll, g, tables; ε_u = c.ε_u)

    # primary wick: ノード0 = secondary wick (温度 T_sw)。Q_leak = primary wick → secondary wick
    Q_leak = residual_wick!(r.wick, d.wick, x.wick, p.wick_p, iface, T_sw, P_cc, tables;
                            ρ_l_ref = p.ρ_l_ref)

    # secondary wick: CCから ṁ を受け、primary wick へ渡す
    residual_swick!(r.swick, d.swick, x.swick, p.swick_p, iface.ṁ, T_cc, x.wick.T[1], P_cc,
                    Q_cs, nc.Q_fluid, Q_leak, tables; ρ_l_ref = p.ρ_l_ref)

    # CC (貯水池): 入口 = 液管 (直接)、出口 = secondary wick へ ṁ。Q_leak (CC直結) は 0
    residual_cc!(r.cc, d.cc, x.cc, inlet, iface.ṁ, nc.Q_fluid, 0.0, p.p_cc, g, tables)

    # (3) 壁の残差
    T_wall_ll_last = x.T_wall.ll[end]

    residual_wall_line!(r.T_wall.vl, d.T_wall.vl, x.T_wall.vl, x.vl.T, h_vl, p.wp_vl,
                        wall_ghost(T_ec, p.ec), wall_ghost(x.T_wall.cond[1], p.wp_cond),
                        c.T_amb, :pipe)
    residual_wall_line!(r.T_wall.cond, d.T_wall.cond, x.T_wall.cond, x.cond.T, h_cond, p.wp_cond,
                        wall_ghost(x.T_wall.vl[end], p.wp_vl), wall_ghost(x.T_wall.ll[1], p.wp_ll),
                        c.T_cond, :condenser; h_con = c.h_con)
    residual_wall_line!(r.T_wall.ll, d.T_wall.ll, x.T_wall.ll, x.ll.T, h_ll, p.wp_ll,
                        wall_ghost(x.T_wall.cond[end], p.wp_cond), wall_ghost(T_ccw, p.cw),
                        c.T_amb, :pipe)

    # CCケース: 蒸発器ケース・液管の壁との熱伝導、secondary wick への熱伝導 Q_cs
    Q_cc_ll = Q_cc_ll_gain(T_ccw, p.cw, T_wall_ll_last, p.wp_ll)
    residual_cc_wall!(r.T_wall.cc, d.T_wall.cc, x.T_wall.cc, Q_cs, Q_ec_cc, Q_cc_ll, p.cw)

    r.T_wall.evap[1] = residual_evapcase(d.T_wall.evap[1], T_ec, p.ec, applied_heat(t, c),
                                         Q_ec_vl, Q_ec_cc, iface.Q_ec, sum(Q_ec_grv))
    return nothing
end