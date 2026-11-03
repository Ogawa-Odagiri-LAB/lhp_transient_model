# =============================================================================
# initial_from_measurements.jl (実測の温度・封入量からの初期条件)
#
# 実験の初期のセンサ温度と封入量から、LHPTransientCondition の初期条件の引数を作る。
#
#   measured_initial(geo, tables, sensors; charge, θ_tol, ...) -> NamedTuple
#
# 戻り値の NamedTuple を、LHPTransientCondition(...; init...) に展開して渡す:
#     init = measured_initial(geo, tables, sensors; charge = 7.6e-3)
#     condition = LHPTransientCondition(; Q_apply = 5.0, ..., init...)
#
# 【手順】
#  1. 基準温度: CC の壁 (RES1, RES2 の平均) を CC 流体の温度 T_init とし、P0 = P_sat(T_init)。
#     蒸発器ケース = EC1, EC2 の平均、凝縮器板 = cond。
#  2. 壁温度: 各管のセンサ位置の温度を PCHIP で補間し、セルの中心の壁温度とする (端は平坦に外挿)。
#       蒸気管 VL1,2,3 : 蒸発器側の根元から 1/6, 1/2, 5/6
#       凝縮器 CON1..5 : 入口から出口まで 5 点の等間隔 (0, 1/4, 1/2, 3/4, 1)
#       液管   LL1,2,3 : 凝縮器から約 10 cm、CC から約 50 cm、CC から約 4 cm (全長 1.15 m)
#  3. 位相: Θ = T_wall − T_sat(P0)。
#       Θ < −θ_tol (既定 0.15 K) なら液 (蒸気なら壁で凝縮して壁が T_sat まで温まるはずなので)。
#       それ以外は蒸気。ただし、液管は常に液 (CC へ入る液)。
#     センサの誤差 (約 0.05〜0.1 K) 以内の Θ は飽和 (蒸気) とみなす。θ_tol は、蒸気管の液の前線の位置を
#     動かす主なパラメータ。結果は表示するので、確認すること。
#  4. 封入量: M = Σ(液・蒸気の管) + ウィック (空隙が液) + グルーブ (満水) + CC。CC の二相の液体積分率
#     f_l を、M が封入量に一致するように決める (CC 以外が決まれば f_l は一意)。
# =============================================================================

"""
    pchip_interp(xs, ys, x) -> y

単調な PCHIP (Fritsch–Carlson)。xs は昇順。範囲外は端の値 (平坦な外挿)。
"""
function pchip_interp(xs::AbstractVector, ys::AbstractVector, x::Real)
    n = length(xs)
    x <= xs[1] && return float(ys[1])
    x >= xs[n] && return float(ys[n])
    h = diff(xs)
    δ = diff(ys) ./ h
    d = zeros(n)
    d[1] = δ[1]; d[n] = δ[n-1]
    for k in 2:n-1
        if δ[k-1] * δ[k] > 0
            w1 = 2h[k] + h[k-1]; w2 = h[k] + 2h[k-1]
            d[k] = (w1 + w2) / (w1 / δ[k-1] + w2 / δ[k])
        end
    end
    k = searchsortedlast(xs, x)
    k = min(k, n - 1)
    t = (x - xs[k]) / h[k]
    h00 = 2t^3 - 3t^2 + 1; h10 = t^3 - 2t^2 + t; h01 = -2t^3 + 3t^2; h11 = t^3 - t^2
    return h00 * ys[k] + h10 * h[k] * d[k] + h01 * ys[k+1] + h11 * h[k] * d[k+1]
end

"セル中心の位置 (長さに対する比) 0.5/N, 1.5/N, ... "
cell_positions(N::Int) = [(i - 0.5) / N for i in 1:N]

"""
    measured_initial(geo, tables, s; charge, θ_tol=0.15, ΔT_sub_init=0.1, ΔT_sup_init=0.1,
                     verbose=true) -> NamedTuple

s: センサ温度 [°C] の NamedTuple (HB, EC1, EC2, VL1..3, CON1..5, LL1..3, RES1, RES2, cond)。
charge: 封入量 [kg]。
f_l_grv: grv_init=:vapor のときのグルーブの初期液体積分率 (LHPTransientCondition の f_l_grv_init と同じにする)。
θ_sat_ll: 液管の壁温が T_sat − θ_sat_ll [K] 以上のセルを、飽和に近い液 (:sat、過冷却度 ΔT_sub_ll_sat) にする (0 なら使わない。センサの誤差 0.05〜0.1 K 程度)。
f_v_ll, n_ll_2ph: 液管の CC 側の末端 n_ll_2ph セルを、蒸気の体積分率 f_v_ll の飽和二相にする (LHPTransientCondition の f_v_ll_init, n_ll_2ph と同じにする)。
n_vap_vl: 蒸気管の蒸気セル数 (蒸発器側から。nothing なら θ_tol で判定)。
grv_init: グルーブの初期状態 (:liquid = 満水、:vapor = 蒸気)。封入量の収支 (CC に割り当てる質量) に効く。LHPTransientCondition の grv_init と同じにする。
"""
function measured_initial(geo::FlatLHPGeometry, tables::FluidTables, s;
    charge::Float64, θ_tol::Float64 = 0.15, ΔT_sub_init::Float64 = 0.1, ΔT_sup_init::Float64 = 0.1,
    grv_init::Symbol = :liquid, n_vap_vl::Union{Nothing,Int} = nothing, f_l_grv::Float64 = 0.0, θ_sat_ll::Float64 = 0.0, f_v_ll::Float64 = 0.0, n_ll_2ph::Int = 1, ΔT_sub_ll_sat::Float64 = 0.02, fluid::String = "ammonia", verbose::Bool = true)

    K = 273.15
    T_cc = (s.RES1 + s.RES2) / 2 + K
    T_ec = (s.EC1 + s.EC2) / 2 + K
    T_cond_plate = s.cond + K

    # 体積・物性を得るため、仮の条件でパラメータを作る
    cond0 = LHPTransientCondition(Q_apply = 0.0, t_init = 0.0, t_end = 1.0, T_init = T_cc, fluid = fluid)
    p0 = build_params(geo, cond0, tables, getaxes(initialize_state(geo)))

    Ts = T_cc    # CC は気液二相なので、飽和温度 = CC の流体温度 (P0 = P_sat(T_cc))

    # --- 壁温度 (PCHIP) ---
    # 液管のセンサ位置 (凝縮器側の根元から、全長に対する比): LL1 = 凝縮器から約 10 cm、LL2 = CC から約 50 cm、
    # LL3 = CC から約 4 cm (全長 1.15 m)
    s_vl, s_ll = [1/6, 1/2, 5/6], [0.10 / 1.15, (1.15 - 0.50) / 1.15, (1.15 - 0.04) / 1.15]
    s_cd = [0.0, 0.25, 0.5, 0.75, 1.0]
    Tw_vl   = [pchip_interp(s_vl, [s.VL1, s.VL2, s.VL3], x) + K for x in cell_positions(geo.vl.N)]
    Tw_cond = [pchip_interp(s_cd, [s.CON1, s.CON2, s.CON3, s.CON4, s.CON5], x) + K for x in cell_positions(geo.cond.N)]
    Tw_ll   = [pchip_interp(s_ll, [s.LL1, s.LL2, s.LL3], x) + K for x in cell_positions(geo.ll.N)]

    # --- 位相 ---
    phase_of(Tw) = [(T - Ts) < -θ_tol ? :liquid : :vapor for T in Tw]
    ph_vl   = phase_of(Tw_vl)
    # n_vap_vl を与えたら、蒸気管の蒸気セルは蒸発器側から n_vap_vl 個 (残りは液)。θ_tol による判定は使わない。
    if n_vap_vl !== nothing
        0 <= n_vap_vl <= geo.vl.N || error("n_vap_vl は 0〜", geo.vl.N)
        ph_vl = [i <= n_vap_vl ? :vapor : :liquid for i in 1:geo.vl.N]
    end
    ph_cond = phase_of(Tw_cond)
    ph_ll   = fill(:liquid, geo.ll.N)
    if θ_sat_ll > 0
        ph_ll = [Tw_ll[i] >= Ts - θ_sat_ll ? :sat : :liquid for i in 1:geo.ll.N]
    end

    # 蒸気のセルの壁は、飽和温度 + ΔT_sup_init 以上に持ち上げる。壁が飽和温度より低い蒸気のセルは、
    # 壁での凝縮の熱伝達が桁違いに大きく (t=0 で 200 W 級の熱流)、起動直後に一気に凝縮して、
    # 時間刻みが 数十 μs に潰れる。センサの誤差 (0.05〜0.1 K) 程度の補正である。
    T_vap_min = Ts + ΔT_sup_init
    n_lift, max_lift = 0, 0.0
    for i in eachindex(Tw_vl)
        if ph_vl[i] === :vapor && Tw_vl[i] < T_vap_min
            max_lift = max(max_lift, T_vap_min - Tw_vl[i]); n_lift += 1; Tw_vl[i] = T_vap_min
        end
    end
    for i in eachindex(Tw_cond)
        if ph_cond[i] === :vapor && Tw_cond[i] < T_vap_min
            max_lift = max(max_lift, T_vap_min - Tw_cond[i]); n_lift += 1; Tw_cond[i] = T_vap_min
        end
    end

    # --- 封入量の収支 ---
    T_liq = T_cc - ΔT_sub_init
    ρ_l_cc, ρ_v_cc = tables.ρ_l_sat(T_cc), tables.ρ_v_sat(T_cc)
    function tube_mass(pc, phase, Tw)
        m = 0.0
        for i in 1:pc.N
            if phase[i] === :sat
                m += tables.ρ_l_sat(min(Tw[i], T_cc - ΔT_sub_ll_sat)) * pc.V
            elseif phase[i] === :liquid
                m += tables.ρ_l_sat(min(Tw[i], T_liq)) * pc.V   # pc.V: 1セルの体積
            else
                m += tables.ρ_v_sat(max(Tw[i], T_cc + ΔT_sup_init)) * pc.V
            end
        end
        return m
    end
    m_vl, m_cond, m_ll = tube_mass(p0.p_vl, ph_vl, Tw_vl), tube_mass(p0.p_cond, ph_cond, Tw_cond),
                         tube_mass(p0.p_ll, ph_ll, Tw_ll)
    # 液管の末端を二相にした分の質量の補正 (液 → 蒸気の体積分率 f_v_ll の混合)
    if f_v_ll > 0
        for i in max(1, geo.ll.N - n_ll_2ph + 1):geo.ll.N
            ρ_prev = ph_ll[i] === :sat ? tables.ρ_l_sat(min(Tw_ll[i], T_cc - ΔT_sub_ll_sat)) :
                                         tables.ρ_l_sat(min(Tw_ll[i], T_liq))
            m_ll += ((1 - f_v_ll) * ρ_l_cc + f_v_ll * ρ_v_cc - ρ_prev) * p0.p_ll.V
        end
    end
    ρ_l_ref = p0.ρ_l_ref
    m_pw = sum(p0.wick_p.V) * p0.wick_p.por * ρ_l_ref
    m_sw = p0.swick_p.V * p0.swick_p.por * ρ_l_ref
    # グルーブ: grv_init = :liquid なら満水、:vapor なら蒸気 (x = 1) で、質量は蒸気の分だけ
    ρ_grv = grv_init === :liquid ? tables.ρ_l_sat(T_liq) :
            f_l_grv * ρ_l_cc + (1 - f_l_grv) * ρ_v_cc      # :vapor で f_l_grv > 0 なら二相 (飽和)
    m_grv = ρ_grv * p0.p_grv.V * p0.p_grv.N
    m_rest = m_vl + m_cond + m_ll + m_pw + m_sw + m_grv
    m_cc = charge - m_rest
    V_cc = p0.p_cc.V * p0.p_cc.N
    f_l_raw = (m_cc / V_cc - ρ_v_cc) / (ρ_l_cc - ρ_v_cc)
    f_l = clamp(f_l_raw, 0.02, 0.98)

    if verbose
        println("=== 実測からの初期条件 ===")
        println("T_cc (RES平均) = ", round(T_cc - K, digits = 3), " °C,  T_ec (EC平均) = ", round(T_ec - K, digits = 3),
                " °C,  凝縮器板 = ", round(T_cond_plate - K, digits = 3), " °C")
        println("P0 = P_sat(T_cc) = ", round(tables.Psat_T(T_cc) / 1e6, digits = 4), " MPa")
        nl(ph) = count(x -> x === :liquid || x === :sat, ph)
        println("位相 (θ_tol = ", θ_tol, " K):  蒸気管 液 ", nl(ph_vl), "/", geo.vl.N, " セル",
                "(前線: 蒸発器側から ", round(count(==(:vapor), ph_vl) / geo.vl.N, digits = 2), "),  凝縮器 液 ",
                nl(ph_cond), "/", geo.cond.N, ",  液管 液 ", nl(ph_ll), "/", geo.ll.N,
                θ_sat_ll > 0 ? " (うち飽和に近い液 " * string(count(==(:sat), ph_ll)) * " セル)" : "")
        println("蒸気セルの壁を T_sat+", ΔT_sup_init, " K まで持ち上げたセル数: ", n_lift, " (最大 ", round(max_lift, digits = 3), " K)")
        println("封入量の収支 [g]:  蒸気管 ", round(m_vl * 1e3, digits = 3), "  凝縮器 ", round(m_cond * 1e3, digits = 3),
                "  液管 ", round(m_ll * 1e3, digits = 3), "  1次ウィック ", round(m_pw * 1e3, digits = 3),
                "  2次ウィック ", round(m_sw * 1e3, digits = 3), "  グルーブ ", round(m_grv * 1e3, digits = 3))
        println("  → CC に割り当てる質量 ", round(m_cc * 1e3, digits = 3), " g  (CC 体積 ", round(V_cc * 1e6, digits = 2),
                " mL) → 液体積分率 f_l = ", round(f_l_raw, digits = 3),
                f_l != f_l_raw ? "  [範囲外のため " * string(f_l) * " に制限。封入量・体積・位相を見直すこと]" : "")
        println("熱容量 [J/K]:  蒸発器ケース(幾何の V_w) ", round(evap_case_heat_capacity(T_ec, p0.ec), digits = 2),
                "  CC ケース ", round(cc_wall_heat_capacity(T_cc, p0.cw), digits = 2),
                "  (実測の見積り: 蒸発器ケース 6.2 + 足 1.0、CC ケース 6.4)")
        println("  感度: 蒸気管の前線が 1 セル動くと CC の質量が ",
                round(p0.p_vl.V * (ρ_l_cc - ρ_v_cc) * 1e3, digits = 3), " g 動く")
    end

    return (T_init = T_cc, T_cond = T_cond_plate,
            ΔT_sub_init = ΔT_sub_init, ΔT_sup_init = ΔT_sup_init, ΔT_sub_ll_sat = ΔT_sub_ll_sat,
            T_wall_init_evap = T_ec, T_wall_init_cc = T_cc,
            T_wall_init_vl = Tw_vl, T_wall_init_cond = Tw_cond, T_wall_init_ll = Tw_ll,
            init_phase_vl = ph_vl, init_phase_cond = ph_cond, init_phase_ll = ph_ll,
            f_l_cc_init = f_l, T_wall_init_hb = s.HB + K)
end

"""
    total_charge(u0, p) -> M [kg]

初期状態の全質量 (流体 Σρm V + ウィックの空隙の液)。measured_initial の封入量と比べる検算用。
"""
function total_charge(u0, p::FlatLHPParams)
    M = 0.0
    for (comp, pc) in ((:grv, p.p_grv), (:vl, p.p_vl), (:cond, p.p_cond), (:ll, p.p_ll), (:cc, p.p_cc))
        M += sum(getproperty(u0, comp).ρm) * pc.V
    end
    M += sum(p.wick_p.V) * p.wick_p.por * p.ρ_l_ref + p.swick_p.V * p.swick_p.por * p.ρ_l_ref
    return M
end