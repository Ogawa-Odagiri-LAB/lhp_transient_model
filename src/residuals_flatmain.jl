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
#     Q_cs   = Q_case_swick2 : CCケース → secondary wick (2 節点 A、B の合計)
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
    cc2::Any                 # CC ケースの 2 節点の経路 (cc_2node のとき。flat_cc_junction_params)
    wick_p::Any              # primary wick
    swick_p::Any             # secondary wick (V, k_eff, ... と CCケースとの熱伝導の量)
    ρ_l_ref::Float64
    plate_nodes::Vector{Int}  # 凝縮器の管の各セルが接するプレートの節点
    cplate_nodes::Vector{Int} # チラー配管の各セルが接するコールドプレートの節点 (コールドプレートがなければ空)
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
        pipe_params(geo.ll),
        flatcc_params(geo.cc; K_liq = liquid_bulk_modulus(condition.T_init, tables),
                      M_mol = CoolProp.PropsSI("molar_mass", condition.fluid)),
        wall_params(geo.vl), wall_params(geo.cond), wall_params(geo.ll),
        flatcc_wall_params(geo.cc), flat_evap_case_params(geo.grv), ec_h, cw_h,
        flat_cc_junction_params(geo.grv, geo.cc),
        flatwick_params(geo.wick, geo.cc, geo.swick, k_s, k_l),
        swick_params(geo.swick, geo.cc, k_s, k_l),
        tables.ρ_l_sat(condition.T_init),
        plate_tube_nodes(geo.plate, geo.cond.L, geo.cond.N),
        n_cplate(geo) > 0 ? cplate_pipe_nodes(geo.plate, geo.cplate) : Int[],
    )
end

"""
    h_ev_at(c, t) -> h_ev(t) [W/m²K]

c.h_ev_end が NaN なら定数 c.h_ev (従来どおり)。そうでなければ、c.h_ev (早い時刻) から c.h_ev_end (遅い時刻) へ、
t = c.h_ev_t_mid を中心に tanh で滑らかに移る。補間は対数空間 (h の比が滑らかに変わる)。試験用 (what-if)。
"""
function h_ev_at(c, t)
    isnan(c.h_ev_end) && return c.h_ev
    s = 0.5 * (1 - tanh((t - c.h_ev_t_mid) / c.h_ev_width))    # 早い時刻 1 → 遅い時刻 0
    return exp(log(c.h_ev_end) + (log(c.h_ev) - log(c.h_ev_end)) * s)
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
    # CC (気相 1 セル + 液 N_CC_LAYER 層、residuals_flatcc.jl): 液管の流入・ウィックへの流出・secondary wick の液は最下層 (層 N) につながる
    swp_cc = c.cc_wall_cond ? p.swick_p : nothing      # 壁面凝縮 (cc_wall_cond) を入れるときだけ secondary wick の物性を渡す
    T_cc_liq = x.ccl.Tl[end]          # 液相の最下層の温度
    P_cc  = x.cc.P[1]                 # CC の圧力 (気相の圧力。液の圧力も同じ = 静水圧は無視)
    # 【数値実験】c.ne_liq_decouple = true: 液側 (ウィック・二次ウィック・液管出口) には、気相の過剰圧力 ΔP_ne = P − P_sat(Tl[1]) を伝えず、界面の飽和圧 P_sat(Tl[1]) を渡す。ΔP_ne の影響の切り分け用 (物理的根拠なし)
    P_cc_liq = c.ne_liq_decouple ? tables.Psat_T(x.ccl.Tl[1]) : P_cc
    bt    = cc_bottom_state(T_cc_liq, P_cc, p.p_cc, g, tables)
    T_sw  = x.swick.T[1]              # secondary wick の節点 A (平面側 + 両端面。primary wick・液管の流れに接する)
    T_swB = x.swick.T[2]              # 節点 B (半円環 = 曲面。蒸気が凝縮する壁)
    P_grv = groove_pressure(x.grv)

    # 蒸発界面: primary wick の最外層 (蒸発器側) と界面。ṁ は secondary wick を通ってくる
    # グルーブの濡れ (液がある) の重み。grv_boiling=false なら 0 (従来どおり)
    w_wet_c = c.grv_boiling ?
        groove_wet_weight(x.grv, p.p_grv, g, tables; f_th = c.f_wet_th, ε_wet = c.ε_wet) :
        zeros(p.p_grv.N)
    w_wet = sum(w_wet_c) / length(w_wet_c)
    # 界面蒸発を止める重み (:dry は液がある間、:saturation はサブクールの液がある間だけ)
    # 核沸騰の活性度 n (状態 nuc.n)。grv_evap_gate = :nucleation のとき、ゲートとグルーブの沸騰の h の両方がこれを使う
    nuc_on = c.grv_boiling && c.grv_evap_gate === :nucleation
    n_nuc  = nuc_on ? [x.nuc.n[i] for i in 1:p.p_grv.N] : nothing
    w_stop, _ = evap_stop_weight(x.grv, p.p_grv, c, tables, w_wet_c; T_ec = T_ec, n_nuc = n_nuc)
    # r_c_wick > 0: グルーブの沸騰の h は別の状態 nc (キャビティ半径 r_c_grv) で切り替える
    sep_nuc = nuc_on && c.r_c_wick > 0
    n_cav   = sep_nuc ? [x.nuc.nc[i] for i in 1:p.p_grv.N] : n_nuc

    iface = wick_interface(x.wick.T[end], P_grv, T_ec, h_ev_at(c, t), geo.grv.A_ec_wick,
                           p.wick_p.G[end], P_cc_liq, tables;
                           T_cc = T_cc_liq, ΔT_startup = c.ΔT_startup, ε_startup = c.ε_startup,
                           w_wet = w_stop, free_interface = c.grv_boiling, q_evap_max = c.q_evap_max)
    S_m_total, S_E_total = wick_groove_source(iface)

    Q_ec_vl = Q_ec_wall(T_ec, p.ec, x.T_wall.vl[1], p.wp_vl)
    T_cca   = c.cc_2node ? x.cca.T[1] : T_ccw
    # ec_2node = true: CC ケースへ逃げる熱は側面の節点 S から出る (熱の流れは P → S → CC ケース)。S の温度で伝導を計算する
    T_ec_src = c.ec_2node ? x.eca.T[1] : T_ec
    ecc     = flat_ec_cc_flows(c, p, T_ec_src, T_cca, T_ccw)
    Q_foot  = c.G_foot_cc * (x.T_wall.hb[1] - T_cca)    # 蒸発器ケースの足: HB → CC ケース (A 節点。cc_2node = false なら CC ケース)
    Q_ec_cc = ecc.Q_ec_in      # 蒸発器ケース (ec_2node: 側面 S) が失う熱 (1 節点: CC ケース B へ、2 節点: 節点 A へ)
    h_grv    = c.grv_boiling ?
        groove_h_ec_boiling(x.grv, p.p_grv, iface.ṁ, T_ec, w_wet_c, g, tables;
                            r_c = c.r_c_grv, ε_onset = c.ε_onset,
                            f_th_h = c.f_wet_th_h, ε_h = c.ε_wet_h, n_nuc = n_cav) :
        groove_h_ec(x.grv, p.p_grv, iface.ṁ, g, tables)
    Q_ec_grv = groove_evap_heat(x.grv, T_ec, h_grv, p.p_grv)

    down_grv = transport_ghost_down(x.vl, p.p_vl, p.p_grv)
    up_vl   = transport_ghost_up(x.grv,  p.p_grv,  p.p_vl)
    # グルーブの成層した出口: 蒸気管へ渡す流束から液の成分を除き、グルーブには差分を戻す
    strat = nothing
    if c.grv_stratified
        strat = groove_stratified_face(x.grv, p.p_grv, down_grv, g, tables;
                                       ε_u = c.ε_u, f_full = c.f_full_grv, ε_full = c.ε_full_grv)
        up_vl = merge(up_vl, (ρ = strat.ρ_eff, Em = strat.Em_eff))
    end
    up_cond = transport_ghost_up(x.vl,   p.p_vl,   p.p_cond)
    up_ll   = transport_ghost_up(x.cond, p.p_cond, p.p_ll)
    h_vl   = zeros(p.p_vl.N);   compute_h_conv_transport!(h_vl,   x.vl,   up_vl,   x.T_wall.vl,   p.p_vl,   g, tables; ε_Re = 200.0, ε_phase = c.ε_phase_h, xps_mode = c.h_xps_mode, f_cond = c.f_h_shah)
    h_cond = zeros(p.p_cond.N); compute_h_conv_transport!(h_cond, x.cond, up_cond, x.T_wall.cond, p.p_cond, g, tables; ε_Re = 200.0, ε_phase = c.ε_phase_h, xps_mode = c.h_xps_mode, f_cond = c.f_h_shah)
    h_ll   = zeros(p.p_ll.N);   compute_h_conv_transport!(h_ll,   x.ll,   up_ll,   x.T_wall.ll,   p.p_ll,   g, tables; ε_Re = 200.0, ε_phase = c.ε_phase_h, xps_mode = c.h_xps_mode, f_cond = (c.f_h_shah_ll < 0 ? c.f_h_shah : c.f_h_shah_ll))
    # 管壁の半径方向の熱伝導抵抗 (内面〜管壁の中心) を、流体−壁の熱伝達 h に直列に入れる (wall_radial_R = true)。壁の節点は管壁の中心の温度。
    # 流体側の q と壁側の Q_fluid が同じ h を共有するので、熱は厳密に保存される。false なら従来 (壁は半径方向に等温)。
    if c.wall_radial_R
        for i in eachindex(h_vl);   h_vl[i]   = h_series(h_vl[i],   wall_R_in(x.T_wall.vl[i],   p.wp_vl));   end
        for i in eachindex(h_cond); h_cond[i] = h_series(h_cond[i], wall_R_in(x.T_wall.cond[i], p.wp_cond)); end
        for i in eachindex(h_ll);   h_ll[i]   = h_series(h_ll[i],   wall_R_in(x.T_wall.ll[i],   p.wp_ll));   end
    end

    # CC: secondary wick → 流体の自然対流、CCケース → secondary wick の熱伝導、液管からの流入
    nc    = flatcc_natural_convection(x.ccl, P_cc, T_sw, T_swB, p.p_cc, g, tables;
                                      boiling = c.cc_boiling, r_c = c.r_c_cc, ε_onset = c.ε_onset_cc)
    qcs   = Q_case_swick2(T_ccw, T_sw, T_swB, p.swick_p)     # CCケース → swick の節点 A、B
    Q_cs  = qcs.Q
    inlet = cc_inlet(x.ll, p.p_ll, bt.ρ, bt.φ; ε_u = c.ε_u)     # 液管 → 最下層

    # 核沸騰の活性度 n の式: τ_nuc dn/dt = n_target(T_ec, T_sat, n) − n (:nucleation 以外では n = 一定)
    for i in 1:p.p_grv.N
        if nuc_on
            Ts_i = Tsat(x.grv.P[i], tables)
            sp_i = sat_props(Ts_i, tables)
            r_c_if = sep_nuc ? c.r_c_wick : c.r_c_grv         # 界面蒸発のゲート: 細孔半径 (分けたとき)
            ε_if   = sep_nuc ? c.ε_onset_wick : c.ε_onset
            r.nuc.n[i] = c.τ_nuc * d.nuc.n[i] -
                         (nucleation_target(T_ec, Ts_i, sp_i, r_c_if, ε_if, c.f_hold_onset, x.nuc.n[i];
                                         x = x.grv.x[i], x_hold = (sep_nuc ? c.x_hold_wick : c.x_hold_nuc), ε_x = c.ε_x_nuc,
                                         n_lat = c.n_lat_nuc, ε_n = c.ε_n_nuc,
                                         T_onset = (sep_nuc && c.onset_ref === :wick) ? x.wick.T[end] : nothing) - x.nuc.n[i])
            if sep_nuc    # グルーブの気泡核生成 (キャビティ半径 r_c_grv)
                r.nuc.nc[i] = c.τ_nuc * d.nuc.nc[i] -
                              (nucleation_target(T_ec, Ts_i, sp_i, c.r_c_grv, c.ε_onset, c.f_hold_onset, x.nuc.nc[i];
                                              x = x.grv.x[i], x_hold = c.x_hold_nuc, ε_x = c.ε_x_nuc,
                                              n_lat = c.n_lat_nuc, ε_n = c.ε_n_nuc) - x.nuc.nc[i])
            else
                r.nuc.nc[i] = d.nuc.nc[i] + x.nuc.nc[i]    # 使わない変数: nc → 0
            end
        else
            r.nuc.n[i]  = d.nuc.n[i]  + x.nuc.n[i]     # 使わない変数: n → 0 (ヤコビアンを特異にしない)
            r.nuc.nc[i] = d.nuc.nc[i] + x.nuc.nc[i]
        end
    end

    # (2) 流体の残差
    residual_grv!(r.grv, d.grv, x.grv, down_grv, Q_ec_grv, S_m_total, S_E_total,
                  p.p_grv, g, tables; ε_u = c.ε_u, strat = strat)

    down_vl = transport_ghost_down(x.cond, p.p_cond, p.p_vl)
    residual_transport_line!(r.vl, d.vl, x.vl, up_vl, down_vl, x.T_wall.vl, h_vl,
                             p.p_vl, g, tables; ε_u = c.ε_u)

    down_cond = transport_ghost_down(x.ll, p.p_ll, p.p_cond)
    residual_transport_line!(r.cond, d.cond, x.cond, up_cond, down_cond, x.T_wall.cond, h_cond,
                             p.p_cond, g, tables; ε_u = c.ε_u)

    # 液管出口のゴースト値: 最下層の液 (ρ, Em)、CC の圧力、出口の流速 (cc.u)。cc_inlet と同じ ρ・φ なので、質量・エネルギーが厳密に一致する
    down_ll = (ρ = bt.ρ, Em = bt.Em, P = P_cc_liq, u = x.cc.u[1] * p.p_cc.A / p.p_ll.A, dl = p.p_cc.dl)
    residual_transport_line!(r.ll, d.ll, x.ll, up_ll, down_ll, x.T_wall.ll, h_ll,
                             p.p_ll, g, tables; ε_u = c.ε_u)

    # primary wick: ノード0 = secondary wick (温度 T_sw)。Q_leak = primary wick → secondary wick
    Q_leak = residual_wick!(r.wick, d.wick, x.wick, p.wick_p, iface, T_sw, P_cc_liq, tables;
                            ρ_l_ref = p.ρ_l_ref)

    # secondary wick: CC (最下層) から ṁ を受け、primary wick へ渡す。CC の流体が受け取る熱 nc.Q_fluid は、swick が失う熱 (残差の −Q_nc)
    # 壁面凝縮の潜熱 ccr.Q_wall は secondary wick へ入る (swick が流体へ渡す熱から引く)。ccr の計算は CC の液の残差と同じ
    ccr = residual_ccl!(r.ccl, d.ccl, x.ccl, P_cc, inlet, iface.ṁ, nc, c.Q_cc_fluid, p.p_cc, g, tables; swp = swp_cc, T_sw = (c.cc_wall_dT ? T_swB : NaN))
    residual_swick2!(r.swick, d.swick, x.swick, p.swick_p, iface.ṁ, T_cc_liq, x.wick.T[1], P_cc_liq,
                     qcs.Q_A, qcs.Q_B, nc.Q_A, nc.Q_B, ccr.Q_wall, Q_leak, tables; ρ_l_ref = p.ρ_l_ref)

    # CC (気相 1 セル + 液 16 層): 入口 = 液管 (最下層へ直接)、出口 = 最下層から secondary wick へ ṁ。
    # ΔP_ne は状態ではなく、気相の圧力 P と界面の液温 Tl[1] の飽和圧の差として結果に出る (界面の凝縮 (97)〜(99) 式)
    residual_cc_summary!(r.cc, x.cc, x.ccl, iface.ṁ, p.p_cc, g, tables)

    # (3) 壁の残差
    T_wall_ll_last = x.T_wall.ll[end]

    residual_wall_line!(r.T_wall.vl, d.T_wall.vl, x.T_wall.vl, x.vl.T, h_vl, p.wp_vl,
                        wall_ghost(T_ec, p.ec), wall_ghost(x.T_wall.cond[1], p.wp_cond),
                        c.T_amb, :pipe; amb = pipe_amb_spec(c, :vl))
    residual_wall_line!(r.T_wall.cond, d.T_wall.cond, x.T_wall.cond, x.cond.T, h_cond, p.wp_cond,
                        wall_ghost(x.T_wall.vl[end], p.wp_vl), wall_ghost(x.T_wall.ll[1], p.wp_ll),
                        [x.plate.T[k] for k in p.plate_nodes], :condenser; h_con = c.h_con, f_con = c.f_con_area, wall_R = c.wall_radial_R)
    if isempty(p.cplate_nodes)
        residual_plate!(r.plate.T, d.plate.T, x.plate.T, x.T_wall.cond, p.wp_cond, c, geo.plate, p.plate_nodes)
    else
        residual_plate!(r.plate.T, d.plate.T, x.plate.T, x.T_wall.cond, p.wp_cond, c, geo.plate, p.plate_nodes; Tcp = x.cplate.T)
        residual_cplate!(r.cplate.T, d.cplate.T, x.cplate.T, x.plate.T, c, geo.plate, geo.cplate, p.cplate_nodes)
    end
    # 熱橋: 熱源 (HB または蒸発器ケース) → 液管の CC 側の末端
    T_br_src = c.br_src === :hb ? x.T_wall.hb[1] : T_ec
    Q_br_ll = thermal_bridge_ll(c, T_br_src, x.T_wall.ll)
    Q_br = sum(Q_br_ll)
    # 液管の CC 側の末端 n_ll_ci セルは CC 側の断熱材の中: 周囲温度 = 断熱材の節点 3 の温度。その熱は節点 3 へ入る
    T_hi3 = c.hi_on ? x.hi.T[3] : c.T_amb
    T_env_ll = ll_env_temps(c, length(x.T_wall.ll), T_hi3)
    Q_ll_to_hi3 = 0.0
    if c.hi_on && c.n_ll_ci > 0
        amb_ll = pipe_amb_spec(c, :ll)
        N_ll = length(x.T_wall.ll)
        for k in 0:min(c.n_ll_ci, N_ll) - 1
            Q_ll_to_hi3 -= Q_amb_pipe(x.T_wall.ll[N_ll - k], p.wp_ll, T_hi3; amb = amb_ll)
        end
    end
    residual_wall_line!(r.T_wall.ll, d.T_wall.ll, x.T_wall.ll, x.ll.T, h_ll, p.wp_ll,
                        wall_ghost(x.T_wall.cond[end], p.wp_cond), wall_ghost(T_ccw, p.cw),
                        T_env_ll, :pipe; amb = pipe_amb_spec(c, :ll), Q_extra = Q_br_ll)

    # 周囲との遅い熱経路: 断熱材 3 節点 (HB 側・蒸発器ケース側・CC 側) と PEEK。壁が失う熱 (失う向きが正)
    Q_hb_loss, Q_ev_loss, Q_cc_loss = residual_housing!(r.hi.T, d.hi.T, x.hi.T, r.peek.T, d.peek.T, x.peek.T,
                                                        x.T_wall.hb[1], T_ec_src, T_ccw, c; Q_ll_in = Q_ll_to_hi3)    # ec_2node: 断熱材に接するのは側面 S
    # CCケース: 蒸発器ケース・液管の壁との熱伝導、secondary wick への熱伝導 Q_cs
    Q_cc_ll = Q_cc_ll_gain(T_ccw, p.cw, T_wall_ll_last, p.wp_ll)
    if c.cc_2node
        c.C_cc_extra > 0 || error("cc_2node = true には C_cc_extra > 0 (節点 A の熱容量) が必要")
        r.cca.T[1] = c.C_cc_extra * d.cca.T[1] - (Q_ec_cc + Q_foot - ecc.Q_ab)
        residual_cc_wall!(r.T_wall.cc, d.T_wall.cc, x.T_wall.cc, Q_cs, ecc.Q_ab, Q_cc_ll, p.cw; C_extra = 0.0, Q_amb = -Q_cc_loss)
    else
        r.cca.T[1] = d.cca.T[1] + x.cca.T[1]    # 使わない変数: → 0
        residual_cc_wall!(r.T_wall.cc, d.T_wall.cc, x.T_wall.cc, Q_cs, Q_ec_cc + Q_foot, Q_cc_ll, p.cw; C_extra = c.C_cc_extra, Q_amb = -Q_cc_loss)
    end

    # ヒーターブロック HB: ヒーターの熱 Q_apply が入り、接触熱伝達で蒸発器ケースへ伝わる
    T_hb  = x.T_wall.hb[1]
    Q_hb_ec = c.h_hb_ec * c.A_hb_ec * (T_hb - T_ec)
    r.T_wall.hb[1] = c.C_hb * d.T_wall.hb[1] - (applied_heat(t, c) - Q_hb_ec - Q_foot - (c.br_src === :hb ? Q_br : 0.0) - Q_hb_loss)
    # 蒸発器ケース: ec_2node = true のとき、P (T_wall.evap) の熱容量は C_ec_P (一定)、側面の節点 S (eca) は P だけと G_ec_ps で結ぶ
    if c.ec_2node
        (c.C_ec_P > 0 && c.C_ec_S > 0 && c.G_ec_ps > 0) || error("ec_2node = true には C_ec_P, C_ec_S, G_ec_ps > 0 が必要")
        Q_ps = c.G_ec_ps * (T_ec - x.eca.T[1])
        r.eca.T[1] = c.C_ec_S * d.eca.T[1] - (Q_ps - Q_ec_cc - Q_ev_loss)    # S: P から Q_ps が入り、CC ケースへ Q_ec_cc、断熱材へ Q_ev_loss が出る
        C_extra_ec = c.C_ec_P - evap_case_heat_capacity(T_ec, p.ec)    # residual_evapcase の C = 幾何 + C_extra = C_ec_P
    else
        Q_ps = 0.0
        r.eca.T[1] = d.eca.T[1] + x.eca.T[1]    # 使わない変数: → 0
        C_extra_ec = c.C_ec_extra
    end
    r.T_wall.evap[1] = residual_evapcase(d.T_wall.evap[1], T_ec, p.ec, Q_hb_ec,
                                         Q_ec_vl + (c.br_src === :hb ? 0.0 : Q_br) + (c.ec_2node ? 0.0 : Q_ev_loss) + Q_ps, c.ec_2node ? 0.0 : Q_ec_cc, iface.Q_ec, sum(Q_ec_grv); C_extra = C_extra_ec)
    return nothing
end