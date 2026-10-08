# =============================================================================
# residuals_main.jl (全体の残差 f_dae_residual!)
#
# IDA (Sundials.jl) に渡す残差関数 F(t, u, du) = 0 を、各コンポーネントの残差
# (residuals_*.jl) を束ねて組み立てる。
#
#   DAEProblem(f_dae_residual!, du0, u0, tspan, p; differential_vars=...)
#   f_dae_residual!(res, du, u, p, t)   … res, du, u は Vector{Float64}
#
# 【ループの構成】 蒸発器 (ウィック・蒸発器ケース) → グルーブ → 蒸気管 → 凝縮器
#   → 液管 → CC → ウィック … と循環する。
#   流体:   grv → vl → cond → ll → cc →(蒸発量 ṁ)→ wick →(蒸発)→ grv
#   壁:     evap ⇄ vl ⇄ cond ⇄ ll ⇄ cc ⇄ evap  (軸方向の熱伝導)
#
# 【共有する量 (1回だけ計算して、複数の残差に配る)】
#   別々に計算し直すと、失う熱と受け取る熱がずれて、エネルギーが保存しなくなる。
#     iface       … wick_interface: Q_ec-wick, ṁ, T_out, Q_extra
#                    (蒸発器ケース・ウィック・グルーブ・CC で使う)
#     Q_ec_vl     … 蒸発器ケース→蒸気管の壁 (蒸発器ケース・蒸気管の壁)
#     Q_ec_cc     … 蒸発器ケース→CCケース (蒸発器ケース・CCの壁)
#     Q_ec_grv    … 蒸発器ケース→グルーブ蒸気 (蒸発器ケース・グルーブ)
#     h_conv      … 流体-壁の熱伝達係数 (vl, cond, ll: 流体の q と壁の Q_fluid)
#     Q_fluid(cc) … CCケース→CC流体 (CCの流体・CCの壁)
#     Q_leak      … ウィック内表面→CC流体 (ウィック・CC)
#     ε_u         … 風上ブレンドの幅 (液管・CC の境界の質量を一致させる)
#
# 【コンポーネント間の接続 (ゴースト値)】
#   流体: 隣のコンポーネントの最終/最初のセルの値。コンポーネントの境界で
#         断面積が違うので、面の速度を、質量流量 ρuA が保たれるよう、この成分の
#         断面積に換算する (transport_ghost_up/down)。
#   壁:   上流・下流の壁セルの温度と幾何 (wall_ghost)。
#         vl の上流は蒸発器ケース、ll の下流は CCケース。
#   ただし、蒸気管の流体の上流はグルーブ、壁の上流は蒸発器ケースで、
#   流体と壁で上流の相手が違う。
#
# 【凝縮器の壁】 周囲 T_amb ではなく、凝縮器プレート T_cond と、接触熱伝達係数
# h_con (式51) で結合する。
# =============================================================================

using ComponentArrays: ComponentArray, getaxes

# -----------------------------------------------------------------------------
# 1. パラメータ (事前計算)
# -----------------------------------------------------------------------------

"""
    LHPParams

残差関数が毎回使う、事前計算した量。幾何パラメータ (pipe_params など) は、
最初に1回だけ作って、ここに入れておく。
"""
struct LHPParams{A,T}
    geo::LHPGeometry
    cond::LHPTransientCondition
    tables::T
    axes::A                  # ComponentArray の軸 (Vector を ComponentArray に包むために使う)
    # 流体の幾何
    p_grv::Any
    p_vl::Any
    p_cond::Any
    p_ll::Any
    p_cc::Any
    # 壁の幾何
    wp_vl::Any
    wp_cond::Any
    wp_ll::Any
    cw::Any                  # CCケース
    ec::Any                  # 蒸発器ケース
    # ウィック
    wick_p::Any
    ρ_l_ref::Float64         # 式46の液の密度 (定数)
end

"""
    build_params(geo, condition, tables, axes) -> LHPParams

幾何パラメータを事前計算する。axes は ComponentArrays.getaxes(u0)。
ρ_l_ref (ウィックの液の基準密度) は、初期温度 T_init での飽和液密度。
"""
function build_params(geo::LHPGeometry, condition::LHPTransientCondition,
    tables::FluidTables, axes)
    return LHPParams(
        geo, condition, tables, axes,
        grv_params(geo.grv), pipe_params(geo.vl), pipe_params(geo.cond),
        pipe_params(geo.ll), cc_params(geo.cc),
        wall_params(geo.vl), wall_params(geo.cond), wall_params(geo.ll),
        cc_wall_params(geo.cc), evap_case_params(geo.grv),
        wick_params(geo.wick, k_sus304(condition.T_init), tables.k_l_sat(condition.T_init)),
        tables.ρ_l_sat(condition.T_init),
    )
end

# -----------------------------------------------------------------------------
# 2. 全体の残差
# -----------------------------------------------------------------------------

"""
    f_dae_residual!(res, du, u, p::LHPParams, t)

IDA に渡す全体の残差。res, du, u は Vector{Float64} で、p.axes で
ComponentArray に包んで、各コンポーネントの残差関数に渡す。
"""
function f_dae_residual!(res, du, u, p::LHPParams, t)
    r  = ComponentArray(res, p.axes)
    d  = ComponentArray(du,  p.axes)
    x  = ComponentArray(u,   p.axes)
    geo, c, tables = p.geo, p.cond, p.tables
    g = c.g

    # -------------------------------------------------------------------------
    # (1) 共有する量を、状態から1回だけ計算する
    # -------------------------------------------------------------------------
    T_ec  = x.T_wall.evap[1]         # 蒸発器ケース (長さ1の配列。types.jl 参照)
    T_ccw = x.T_wall.cc[1]
    T_cc  = x.cc.T[1]
    P_cc  = x.cc.P[1]
    P_grv = groove_pressure(x.grv)

    # ウィック: 蒸発界面 (蒸発器ケースから界面への熱 Q_ec、蒸発量 ṁ など)
    iface = wick_interface(x.wick.T[end], P_grv, T_ec, c.h_ev, geo.grv.A_ec_wick,
                           p.wick_p.G[end], P_cc, tables;
                           T_cc = T_cc, ΔT_startup = c.ΔT_startup, ε_startup = c.ε_startup)
    S_m_total, S_E_total = wick_groove_source(iface)

    # 蒸発器ケースから、隣の壁・グルーブ蒸気への熱
    Q_ec_vl = Q_ec_wall(T_ec, p.ec, x.T_wall.vl[1], p.wp_vl)
    Q_ec_cc = Q_ec_wall(T_ec, p.ec, T_ccw, p.cw)
    h_grv    = groove_h_ec(x.grv, p.p_grv, iface.ṁ, g, tables)
    Q_ec_grv = groove_evap_heat(x.grv, T_ec, h_grv, p.p_grv)

    # 輸送管 (vl, cond, ll) の流体-壁の熱伝達係数 (流体の q と壁の Q_fluid で共有)
    up_vl   = transport_ghost_up(x.grv,  p.p_grv,  p.p_vl)
    up_cond = transport_ghost_up(x.vl,   p.p_vl,   p.p_cond)
    up_ll   = transport_ghost_up(x.cond, p.p_cond, p.p_ll)
    h_vl   = zeros(p.p_vl.N);   compute_h_conv_transport!(h_vl,   x.vl,   up_vl,   x.T_wall.vl,   p.p_vl,   g, tables; ε_Re = 200.0)
    h_cond = zeros(p.p_cond.N); compute_h_conv_transport!(h_cond, x.cond, up_cond, x.T_wall.cond, p.p_cond, g, tables; ε_Re = 200.0)
    h_ll   = zeros(p.p_ll.N);   compute_h_conv_transport!(h_ll,   x.ll,   up_ll,   x.T_wall.ll,   p.p_ll,   g, tables; ε_Re = 200.0)

    # CC: 自然対流 (CCケース→流体) と液管からの流入
    nc = cc_natural_convection(x.cc, T_ccw, p.p_cc, g, tables)
    inlet = cc_inlet(x.ll, p.p_ll, x.cc; ε_u = c.ε_u)

    # -------------------------------------------------------------------------
    # (2) 流体の残差
    # -------------------------------------------------------------------------
    # グルーブ (下流 = 蒸気管の最初のセル)
    down_grv = transport_ghost_down(x.vl, p.p_vl, p.p_grv)
    residual_grv!(r.grv, d.grv, x.grv, down_grv, Q_ec_grv, S_m_total, S_E_total,
                  p.p_grv, g, tables; ε_u = c.ε_u)

    # 蒸気管 (上流 = グルーブ、下流 = 凝縮器)
    down_vl = transport_ghost_down(x.cond, p.p_cond, p.p_vl)
    residual_transport_line!(r.vl, d.vl, x.vl, up_vl, down_vl, x.T_wall.vl, h_vl,
                             p.p_vl, g, tables; ε_u = c.ε_u)

    # 凝縮器 (上流 = 蒸気管、下流 = 液管)
    down_cond = transport_ghost_down(x.ll, p.p_ll, p.p_cond)
    residual_transport_line!(r.cond, d.cond, x.cond, up_cond, down_cond, x.T_wall.cond, h_cond,
                             p.p_cond, g, tables; ε_u = c.ε_u)

    # 液管 (上流 = 凝縮器、下流 = CC)
    down_ll = transport_ghost_down(x.cc, p.p_cc, p.p_ll)
    residual_transport_line!(r.ll, d.ll, x.ll, up_ll, down_ll, x.T_wall.ll, h_ll,
                             p.p_ll, g, tables; ε_u = c.ε_u)

    # ウィック (CC の流出 ṁ を受けて蒸発する)。Q_leak は CC の熱源になる
    Q_leak = residual_wick!(r.wick, d.wick, x.wick, p.wick_p, iface, T_cc, P_cc, tables;
                            ρ_l_ref = p.ρ_l_ref)

    # CC (貯水池): 入口 = 液管、出口 = ウィックへ ṁ
    residual_cc!(r.cc, d.cc, x.cc, inlet, iface.ṁ, nc.Q_fluid, Q_leak, p.p_cc, g, tables)

    # -------------------------------------------------------------------------
    # (3) 壁の残差
    # -------------------------------------------------------------------------
    T_wall_ll_last = x.T_wall.ll[end]

    # 蒸気管の壁 (上流 = 蒸発器ケース、下流 = 凝縮器の最初の壁)
    residual_wall_line!(r.T_wall.vl, d.T_wall.vl, x.T_wall.vl, x.vl.T, h_vl, p.wp_vl,
                        wall_ghost(T_ec, p.ec), wall_ghost(x.T_wall.cond[1], p.wp_cond),
                        c.T_amb, :pipe)
    # 凝縮器の壁 (上流 = 蒸気管の最終セル、下流 = 液管の最初の壁)。環境は凝縮器プレート
    residual_wall_line!(r.T_wall.cond, d.T_wall.cond, x.T_wall.cond, x.cond.T, h_cond, p.wp_cond,
                        wall_ghost(x.T_wall.vl[end], p.wp_vl), wall_ghost(x.T_wall.ll[1], p.wp_ll),
                        c.T_cond, :condenser; h_con = c.h_con)
    # 液管の壁 (上流 = 凝縮器の最終セル、下流 = CCケース)
    residual_wall_line!(r.T_wall.ll, d.T_wall.ll, x.T_wall.ll, x.ll.T, h_ll, p.wp_ll,
                        wall_ghost(x.T_wall.cond[end], p.wp_cond), wall_ghost(T_ccw, p.cw),
                        c.T_amb, :pipe)

    # CCケース: 蒸発器ケース・液管の壁との熱伝導、CC流体との自然対流
    Q_cc_ll = Q_cc_ll_gain(T_ccw, p.cw, T_wall_ll_last, p.wp_ll)
    residual_cc_wall!(r.T_wall.cc, d.T_wall.cc, x.T_wall.cc, nc.Q_fluid, Q_ec_cc, Q_cc_ll, p.cw)

    # 蒸発器ケース (スカラー): 印加熱 Q_apply(t) を受け、隣の壁・ウィック・グルーブへ渡す
    r.T_wall.evap[1] = residual_evapcase(d.T_wall.evap[1], T_ec, p.ec, applied_heat(t, c),
                                         Q_ec_vl, Q_ec_cc, iface.Q_ec, sum(Q_ec_grv))
    return nothing
end