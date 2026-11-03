# =============================================================================
# savefuncs.jl (結果の保存)
#
# solve_transient_lhp! の結果 (sol, p) を、CSV に保存する。
# 旧実装は、時間ステップごとに状態を保存していたが、IDA は時間積分を内部で行うので、
# condition.saveat ごとに保存された解 (sol.t, sol.u) を、まとめて書き出す。
#
#   save_results(sol, p, tables, dir)
#     dir/summary.csv … 時系列の要約 (1 保存点 = 1 行)
#     dir/profiles.csv … 空間分布 (縦持ち: 時刻 × 成分 × セル)
#     dir/wick.csv     … ウィックの温度・エネルギーの分布 (時刻 × セル)
#     dir/geometry.csv … 各成分のセル数 N と長さ L (プロットで、位置座標を作るため)
#
# 保存した CSV は、計算をやり直さずに、plot.jl (単体で動く) で、いつでも眺められる。
#
# profiles.csv は縦持ち (long format)。成分ごとにセル数が違うので、横持ちにすると
# 列が不揃いになる。DataFrames で、成分・時刻を絞り込んで使う:
#   using CSV, DataFrames
#   prof = CSV.read("profiles.csv", DataFrame)
#   vl = prof[(prof.comp .== "vl") .& (prof.t .== 10.0), :]
# =============================================================================

using CSV, DataFrames

"""
    results_summary(sol, p, tables) -> DataFrame

時系列の要約。1 保存点 = 1 行。
  t [s], Q_apply [W], T_ec [K] (蒸発器ケース), P_cc [Pa], T_cc [K], P_grv [Pa] (グルーブの平均圧力),
  T_out [K] (蒸発界面の温度 = T_sat(P_grv)), mdot [kg/s] (CC からウィックへの流量 = 蒸発量),
  Q_evap [W] (蒸発による熱 = mdot × h_fg(T_out)。印加熱 Q_apply と比べると、熱が蒸発に使われているか分かる)
"""
function results_summary(sol, p::AbstractLHPParams, tables::FluidTables)
    df = DataFrame(t = Float64[], Q_apply = Float64[], T_ec = Float64[], P_cc = Float64[],
                   T_cc = Float64[], P_grv = Float64[], T_out = Float64[], mdot = Float64[],
                   Q_evap = Float64[],
                   grv_rho = Float64[], grv_x = Float64[], grv_fl = Float64[],
                   w_wet = Float64[], h_grv = Float64[], Q_ec_grv = Float64[],
                   g_s = Float64[], dT_sub = Float64[], T_out_if = Float64[],
                   h_cc = Float64[], w_on_cc = Float64[], Q_cc = Float64[], x_cc = Float64[], u_ll_end = Float64[], T_hb = Float64[],
                   Q_hb_ec = Float64[], Q_ec_vl = Float64[], Q_ec_cc = Float64[], Q_ec_wick = Float64[], Q_br = Float64[], n_nuc = Float64[], dP_ne = Float64[],
                   Q_hi_hb = Float64[], Q_hi_ev = Float64[], Q_hi_cc = Float64[], Q_peek = Float64[],
                   Q_wall = Float64[], Q_cs = Float64[], Q_leak = Float64[], T_sw = Float64[], T_swB = Float64[],
                   T_l1 = Float64[], T_lN = Float64[], f_v_cc = Float64[], m_v_cc = Float64[], mdot_cond = Float64[], h_i = Float64[], h_ev = Float64[],
                   # --- CC まわりの熱経路 (診断。物理は変えない)。符号は各列の説明の向き ---
                   Q_foot = Float64[], T_ecS = Float64[], T_ccA = Float64[], T_ccB = Float64[], T_llend = Float64[], T_llfl = Float64[],
                   Q_ab = Float64[], Q_cs_A = Float64[], Q_cs_B = Float64[], Q_cc_ll = Float64[],
                   Q_nc_A = Float64[], Q_nc_B = Float64[], Q_nc_l = Float64[], Q_nc_v = Float64[], Q_l1 = Float64[],
                   Q_AB_sw = Float64[], Q_hflow = Float64[], Q_ret = Float64[], mdot_in = Float64[],
                   E_swA = Float64[], E_swB = Float64[],
                   bal_caseA = Float64[], bal_caseB = Float64[], bal_swA = Float64[], bal_swB = Float64[], bal_fluid = Float64[])
    for k in eachindex(sol.t)
        c = components_at(sol, p, k)
        # グルーブの診断: ρm, x, 液体積分率 f_l, 濡れ重み, 蒸発器ケース→グルーブの h と Q (grv_boiling のとき)
        gq = p.cond.grv_boiling
        T_ec_k = c.T_wall.evap[1]
        # 核沸騰の活性度 n (grv_evap_gate = :nucleation のときだけ。平板型の状態 nuc.n)
        n_k = (gq && p.cond.grv_evap_gate === :nucleation && hasproperty(c, :nuc)) ?
              [c.nuc.n[i] for i in 1:p.p_grv.N] : nothing
        n_h = (n_k !== nothing && p.cond.r_c_wick > 0) ? [c.nuc.nc[i] for i in 1:p.p_grv.N] : n_k    # グルーブの沸騰の h 用
        ρ_l_g = tables.ρ_l_sat(c.grv.T[1])
        fl_g = (1 - c.grv.x[1]) * c.grv.ρm[1] / ρ_l_g
        ṁ_k = c.cc.u[1] * tables.ρ_l_sat(c.cc.T[1]) * p.p_cc.A
        if gq
            w_k = groove_wet_weight(c.grv, p.p_grv, p.cond.g, tables;
                                    f_th = p.cond.f_wet_th, ε_wet = p.cond.ε_wet)
            h_k = groove_h_ec_boiling(c.grv, p.p_grv, ṁ_k, T_ec_k, w_k, p.cond.g, tables;
                                      r_c = p.cond.r_c_grv, ε_onset = p.cond.ε_onset,
                                      f_th_h = p.cond.f_wet_th_h, ε_h = p.cond.ε_wet_h, n_nuc = n_h)
        else
            w_k = zeros(p.p_grv.N)
            h_k = groove_h_ec(c.grv, p.p_grv, ṁ_k, p.cond.g, tables)
        end
        # 界面のゲート g_s と、グルーブ第1セルのサブクール度 (T_sat(P) − T)、界面の温度 T_out_if
        w_stop_k = gq ? evap_stop_weight(c.grv, p.p_grv, p.cond, tables, w_k; T_ec = T_ec_k, n_nuc = n_k)[1] : 0.0
        newcc = hasproperty(c, :ccl)    # 平板型の新しい CC (気相 1 + 液 16 層): 液はウィック・液管につながる最下層の温度、気相の飽和温度は cc.T
        ccd = newcc ? cc_state_diag(c, p, tables) : nothing
        iface_k = wick_interface(c.wick.T[end], groove_pressure(c.grv), T_ec_k, h_ev_at(p.cond, sol.t[k]),
                                 p.geo.grv.A_ec_wick, p.wick_p.G[end], c.cc.P[1], tables;
                                 T_cc = newcc ? c.ccl.Tl[end] : c.cc.T[1], ΔT_startup = p.cond.ΔT_startup,
                                 ε_startup = p.cond.ε_startup, w_wet = w_stop_k,
                                 free_interface = gq, q_evap_max = p.cond.q_evap_max)
        dT_sub_k = Tsat(c.grv.P[1], tables) - c.grv.T[1]
        # 平板型のみ: CC の自然対流・核沸騰の診断 (円筒型は NaN)
        nc_k = hasproperty(c, :swick) ?
            flatcc_natural_convection(c.ccl, c.cc.P[1], c.swick.T[1], c.swick.T[2], p.p_cc, p.cond.g, tables;
                                      boiling = p.cond.cc_boiling, r_c = p.cond.r_c_cc,
                                      ε_onset = p.cond.ε_onset_cc) :
            (h_l = NaN, w_on = NaN, Q_fluid = NaN)
        Q_k = sum(groove_evap_heat(c.grv, T_ec_k, h_k, p.p_grv))
        # 蒸発器ケースの熱収支の各項 (蒸発器ケースが失う向きが正。残差と同じ関数・同じ引数)
        flat = hasproperty(p, :ec_h)
        Q_hb_ec_k = flat ? p.cond.h_hb_ec * p.cond.A_hb_ec * (c.T_wall.hb[1] - T_ec_k) : NaN
        Q_ec_vl_k = flat ? Q_ec_wall(T_ec_k, p.ec, c.T_wall.vl[1], p.wp_vl) : NaN
        # ec_2node のとき CC ケースへ熱が出るのは側面 S (c.eca) であり、P (T_ec_k) ではない (残差 residuals_flatmain.jl と同じ引数にする)
        T_ec_src_k = (flat && p.cond.ec_2node) ? c.eca.T[1] : T_ec_k
        T_cca_k = (:cca in keys(c) && p.cond.cc_2node) ? c.cca.T[1] : c.T_wall.cc[1]
        ecc_k = flat ? flat_ec_cc_flows(p.cond, p, T_ec_src_k, T_cca_k, c.T_wall.cc[1]) : nothing
        Q_ec_cc_k = flat ? ecc_k.Q_ec_in : NaN
        Q_br_k = flat ? sum(thermal_bridge_ll(p.cond, p.cond.br_src === :hb ? c.T_wall.hb[1] : T_ec_k, c.T_wall.ll)) : NaN
        # secondary wick の熱収支の診断 (平板型): 壁面凝縮で swick に入る熱 Q_wall、CC ケース → swick の Q_cs、primary wick → swick の Q_leak
        swk = flat && hasproperty(c, :swick)
        Q_wall_k = swk ? (newcc ? ccd.Q_wall : 0.0) : NaN    # 壁面凝縮の潜熱 (cc_wall_cond = false なら 0)
        Q_cs_k   = swk ? Q_case_swick2(c.T_wall.cc[1], c.swick.T[1], c.swick.T[2], p.swick_p).Q : NaN
        Q_leak_k = swk ? p.wick_p.G[1] * (c.wick.T[1] - c.swick.T[1]) : NaN
        # CC まわりの熱経路 (平板型・新 CC のみ。残差と同じ関数・同じ引数)
        cch = swk && newcc
        T_ccB_k = flat ? c.T_wall.cc[1] : NaN
        qcs_k = cch ? Q_case_swick2(T_ccB_k, c.swick.T[1], c.swick.T[2], p.swick_p) : nothing
        Q_ab_k = (flat && p.cond.cc_2node) ? ecc_k.Q_ab : NaN
        T_wll_k = flat ? c.T_wall.ll[end] : NaN
        Q_ccll_k = flat ? Q_cc_ll_gain(T_ccB_k, p.cw, T_wll_k, p.wp_ll) : NaN
        # secondary wick の A→B 方向の内部熱 Q_AB (B → A が正)、液の流れが運ぶ ṁ(hf_in − hf_out)
        T_swA_k = cch ? c.swick.T[1] : NaN
        T_swB_k = cch ? c.swick.T[2] : NaN
        Q_ABsw_k = cch ? p.swick_p.G_AB * (T_swB_k - T_swA_k) : NaN
        T_ccliq_k = newcc ? c.ccl.Tl[end] : NaN
        P_wl_k = newcc ? (p.cond.ne_liq_decouple ? tables.Psat_T(c.ccl.Tl[1]) : c.cc.P[1]) : NaN
        ṁ_if_k = iface_k.ṁ
        wm_k = cch ? smooth_step(ṁ_if_k, 1e-8) : NaN
        hf_in_k  = cch ? wm_k * liquid_props(T_ccliq_k, P_wl_k, tables).h + (1 - wm_k) * liquid_props(T_swA_k, P_wl_k, tables).h : NaN
        hf_out_k = cch ? wm_k * liquid_props(T_swA_k, P_wl_k, tables).h + (1 - wm_k) * liquid_props(c.wick.T[1], P_wl_k, tables).h : NaN
        Q_hflow_k = cch ? ṁ_if_k * (hf_in_k - hf_out_k) : NaN
        # 液管から CC へ戻る液が、最下層の液から受ける顕熱 ṁ_in(φ_in − φ_N) と、その流量
        bt_k = newcc ? cc_bottom_state(T_ccliq_k, c.cc.P[1], p.p_cc, p.cond.g, tables) : nothing
        inl_k = newcc ? cc_inlet(c.ll, p.p_ll, bt_k.ρ, bt_k.φ; ε_u = p.cond.ε_u) : nothing
        Q_ret_k = newcc ? inl_k.ṁ * (inl_k.φ - bt_k.φ) : NaN
        # 層 1 に入る熱 (蒸気側の壁から + 界面凝縮の潜熱)
        Q_l1_k = newcc ? nc_k.Q_v + ccd.ṁ_cond * ccd.λ : NaN
        # 液管の最後のセルの流体温度
        T_llfl_k = c.ll.T[end]
        P_grv = groove_pressure(c.grv)
        ρ_l = tables.ρ_l_sat(c.cc.T[1])
        T_out = Tsat(P_grv, tables)
        mdot = c.cc.u[1] * ρ_l * p.p_cc.A
        h_fg = tables.h_v_sat(T_out) - tables.h_l_sat(T_out)
        push!(df, (t = sol.t[k], Q_apply = applied_heat(sol.t[k], p.cond),
                   T_ec = c.T_wall.evap[1], P_cc = c.cc.P[1], T_cc = c.cc.T[1],
                   P_grv = P_grv, T_out = T_out, mdot = mdot, Q_evap = mdot * h_fg,
                   grv_rho = c.grv.ρm[1], grv_x = c.grv.x[1], grv_fl = fl_g,
                   w_wet = w_k[1], h_grv = h_k[1], Q_ec_grv = Q_k,
                   g_s = iface_k.g_s, dT_sub = dT_sub_k, T_out_if = iface_k.T_out,
                   h_cc = nc_k.h_l, w_on_cc = nc_k.w_on, Q_cc = nc_k.Q_fluid, x_cc = c.cc.x[1],
                   u_ll_end = c.ll.u[end],
                   T_hb = :hb in keys(c.T_wall) ? c.T_wall.hb[1] : NaN,
                   Q_hb_ec = Q_hb_ec_k, Q_ec_vl = Q_ec_vl_k, Q_ec_cc = Q_ec_cc_k,
                   Q_ec_wick = iface_k.Q_ec, Q_br = Q_br_k, n_nuc = n_k === nothing ? NaN : n_k[1],
                   dP_ne = newcc ? ccd.dP_ne : NaN,     # ΔP_ne = P_cc − P_sat(層 1 の温度) [Pa]
                   Q_hi_hb = (:hi in keys(c) && p.cond.hi_on) ? hi_G_in(p.cond, p.cond.hi_A_hb) * (c.T_wall.hb[1] - c.hi.T[1]) : NaN,
                   Q_hi_ev = (:hi in keys(c) && p.cond.hi_on) ? hi_G_in(p.cond, p.cond.hi_A_ev) * (T_ec_k - c.hi.T[2]) : NaN,
                   Q_hi_cc = (:hi in keys(c) && p.cond.hi_on) ? (p.cond.cc_adiabatic ? 0.0 : hi_G_in(p.cond, p.cond.hi_A_cc) * (c.T_wall.cc[1] - c.hi.T[3])) : NaN,
                   Q_peek = (:hi in keys(c) && p.cond.peek_on) ? p.cond.peek_k * p.cond.peek_A / (peek_dz(p.cond)[1] / 2) * (c.T_wall.hb[1] - c.peek.T[1]) : NaN,
                   Q_wall = Q_wall_k, Q_cs = Q_cs_k, Q_leak = Q_leak_k, T_sw = swk ? c.swick.T[1] : NaN, T_swB = swk ? c.swick.T[2] : NaN,
                   T_l1 = newcc ? ccd.T1 : NaN, T_lN = newcc ? ccd.TN : NaN, f_v_cc = newcc ? ccd.f_v : NaN,
                   m_v_cc = newcc ? ccd.mv : NaN, mdot_cond = newcc ? ccd.ṁ_cond : NaN, h_i = newcc ? ccd.h_i : NaN,
                   h_ev = h_ev_at(p.cond, sol.t[k]),
                   Q_foot = flat ? p.cond.G_foot_cc * (c.T_wall.hb[1] - T_cca_k) : NaN,
                   T_ecS = (flat && p.cond.ec_2node) ? c.eca.T[1] : NaN, T_ccA = T_cca_k, T_ccB = T_ccB_k,
                   T_llend = T_wll_k, T_llfl = T_llfl_k,
                   Q_ab = Q_ab_k, Q_cs_A = cch ? qcs_k.Q_A : NaN, Q_cs_B = cch ? qcs_k.Q_B : NaN, Q_cc_ll = Q_ccll_k,
                   Q_nc_A = newcc ? nc_k.Q_A : NaN, Q_nc_B = newcc ? nc_k.Q_B : NaN,
                   Q_nc_l = newcc ? nc_k.Q_l : NaN, Q_nc_v = newcc ? nc_k.Q_v : NaN, Q_l1 = Q_l1_k,
                   Q_AB_sw = Q_ABsw_k, Q_hflow = Q_hflow_k, Q_ret = Q_ret_k, mdot_in = newcc ? inl_k.ṁ : NaN,
                   E_swA = cch ? c.swick.E[1] : NaN, E_swB = cch ? c.swick.E[2] : NaN,
                   bal_caseA = NaN, bal_caseB = NaN, bal_swA = NaN, bal_swB = NaN, bal_fluid = NaN))
    end
    # 閉じ方の確認 (各節点・CC 流体の熱収支の残り [W]。0 に近いほど経路が漏れなく入っている)。
    #   dT/dt は保存点の中心差分 (保存間隔が粗いと精度が落ちる)。
    #   bal_caseA : C_cc_extra dT_A/dt − (Q_ec_cc + Q_foot − Q_ab)                       (cc_2node のとき)
    #   bal_caseB : C_B dT_B/dt − (Q_ab + Q_cc_ll − Q_cs_A − Q_cs_B − Q_hi_cc)
    #   bal_swA   : dE_A/dt − (Q_hflow + Q_cs_A − Q_nc_A + Q_leak + Q_AB_sw)
    #   bal_swB   : dE_B/dt − (Q_cs_B − Q_nc_B + Q_wall − Q_AB_sw)
    #   bal_fluid : Q_nc_l + Q_nc_v − (Q_nc_A + Q_nc_B)  (CC 流体が受ける熱と swick が失う熱の一致。流体側の全収支は列 Q_nc_*, Q_l1, Q_ret から手で確認する)
    n = nrow(df)
    if n >= 3 && all(!isnan, df.T_ccB) && hasproperty(p, :cw)
        dTdt(v, i) = i == 1 ? (v[2] - v[1]) / (df.t[2] - df.t[1]) :
                     i == n ? (v[n] - v[n-1]) / (df.t[n] - df.t[n-1]) :
                     (v[i+1] - v[i-1]) / (df.t[i+1] - df.t[i-1])
        for i in 1:n
            if p.cond.cc_2node
                df.bal_caseA[i] = p.cond.C_cc_extra * dTdt(df.T_ccA, i) - (df.Q_ec_cc[i] + df.Q_foot[i] - df.Q_ab[i])
                df.bal_caseB[i] = cc_wall_heat_capacity(df.T_ccB[i], p.cw) * dTdt(df.T_ccB, i) -
                                  (df.Q_ab[i] + df.Q_cc_ll[i] - df.Q_cs_A[i] - df.Q_cs_B[i] - df.Q_hi_cc[i])
            end
            if !isnan(df.E_swA[i])
                df.bal_swA[i] = dTdt(df.E_swA, i) - (df.Q_hflow[i] + df.Q_cs_A[i] - df.Q_nc_A[i] + df.Q_leak[i] + df.Q_AB_sw[i])
                df.bal_swB[i] = dTdt(df.E_swB, i) - (df.Q_cs_B[i] - df.Q_nc_B[i] + df.Q_wall[i] - df.Q_AB_sw[i])
            end
            df.bal_fluid[i] = df.Q_nc_l[i] + df.Q_nc_v[i] - (df.Q_nc_A[i] + df.Q_nc_B[i])
        end
    end
    return df
end

"""
    results_profiles(sol, p) -> DataFrame

空間分布 (縦持ち)。1 行 = (保存時刻, 成分, セル)。
  comp: "grv", "vl", "cond", "ll", "cc" (流体の成分)、"evap" (蒸発器ケース)、"hb" (ヒーターブロック、平板型のみ。T_wall だけ入る)
  cell: 成分内のセル番号 (上流 = 1)
  rho [kg/m³], u [m/s] (セルの下流側の面の速度), Em [J/kg], P [Pa], T [K], x [-] (乾き度),
  T_wall [K] (壁温度。グルーブは壁を持たないので NaN)
"""
function results_profiles(sol, p::AbstractLHPParams)
    df = DataFrame(t = Float64[], comp = String[], cell = Int[], rho = Float64[], u = Float64[],
                   Em = Float64[], P = Float64[], T = Float64[], x = Float64[], T_wall = Float64[])
    for k in eachindex(sol.t)
        c = components_at(sol, p, k)
        for comp in (:grv, :vl, :cond, :ll, :cc)
            cm = getproperty(c, comp)
            n = length(cm.ρm)
            Tw = comp === :grv ? fill(NaN, n) : getproperty(c.T_wall, comp)
            for i in 1:n
                push!(df, (t = sol.t[k], comp = String(comp), cell = i, rho = cm.ρm[i], u = cm.u[i],
                           Em = cm.Em[i], P = cm.P[i], T = cm.T[i], x = cm.x[i], T_wall = Tw[i]))
            end
        end
        # CC の液の層 (平板型の新しい CC): comp = "ccl"、cell = 層 (1 = 気液界面側、N = 底)。T = 層の温度、rho = 飽和液密度、P = CC の圧力、x = 0
        if :ccl in keys(c)
            for (i, Tk) in enumerate(c.ccl.Tl)
                push!(df, (t = sol.t[k], comp = "ccl", cell = i, rho = p.tables.ρ_l_sat(Tk), u = NaN, Em = NaN,
                           P = c.cc.P[1], T = Tk, x = 0.0, T_wall = NaN))
            end
        end
        # 凝縮器プレート (平板型): comp = "plate"、T_wall にプレートの温度
        if :plate in keys(c)
            for (i, Tpl) in enumerate(c.plate.T)
                push!(df, (t = sol.t[k], comp = "plate", cell = i, rho = NaN, u = NaN, Em = NaN,
                           P = NaN, T = NaN, x = NaN, T_wall = Tpl))
            end
        end
        # コールドプレート (2次元): comp = "cplate"
        if :cplate in keys(c)
            for (i, Tcp) in enumerate(c.cplate.T)
                push!(df, (t = sol.t[k], comp = "cplate", cell = i, rho = NaN, u = NaN, Em = NaN,
                           P = NaN, T = NaN, x = NaN, T_wall = Tcp))
            end
        end
        # 断熱材 (hi: 3 節点) と PEEK (peek): 節点の温度を T_wall に入れる
        for wname in (:hi, :peek, :cca, :eca)
            wname in keys(c) || continue
            for (i, Tv) in enumerate(getproperty(c, wname).T)
                push!(df, (t = sol.t[k], comp = String(wname), cell = i, rho = NaN, u = NaN, Em = NaN,
                           P = NaN, T = NaN, x = NaN, T_wall = Tv))
            end
        end
        # 集中容量の壁: 蒸発器ケース (evap)、平板型はヒーターブロック (hb)。流体の量は NaN、T_wall に壁温度
        for wname in (:evap, :hb)
            wname in keys(c.T_wall) || continue
            push!(df, (t = sol.t[k], comp = String(wname), cell = 1, rho = NaN, u = NaN, Em = NaN,
                       P = NaN, T = NaN, x = NaN, T_wall = getproperty(c.T_wall, wname)[1]))
        end
    end
    return df
end

"""
    results_wick(sol, p) -> DataFrame

ウィックの分布。1 行 = (保存時刻, セル)。cell は内側 (CC側) = 1。平板型は、secondary wick を cell = 0 で含む。T [K], E [J] (式49)。
"""
function results_wick(sol, p::AbstractLHPParams)
    df = DataFrame(t = Float64[], cell = Int[], T = Float64[], E = Float64[])
    for k in eachindex(sol.t)
        c = components_at(sol, p, k)
        # 平板型: secondary wick は cell = 0 の行として保存する (primary wick の内側のノード0)
        :swick in keys(c) && begin
            push!(df, (t = sol.t[k], cell = 0, T = c.swick.T[1], E = c.swick.E[1]))      # 節点 A
            push!(df, (t = sol.t[k], cell = -1, T = c.swick.T[2], E = c.swick.E[2]))     # 節点 B (半円環)
        end
        for i in eachindex(c.wick.T)
            push!(df, (t = sol.t[k], cell = i, T = c.wick.T[i], E = c.wick.E[i]))
        end
    end
    return df
end

"""
    results_geometry(p) -> DataFrame

各成分のセル数 N と長さ L [m] (comp = "grv", "vl", "cond", "ll", "cc")。
プロットで、ループに沿った位置座標を作るのに使う (plot.jl の load_lengths が読む)。
"""
function results_geometry(p::AbstractLHPParams)
    g = p.geo
    return DataFrame(comp = ["grv", "vl", "cond", "ll", "cc"],
                     N = [g.grv.N, g.vl.N, g.cond.N, g.ll.N, g.cc.N],
                     L = [groove_length(g.grv), g.vl.L, g.cond.L, g.ll.L, g.cc.L])
end

"""
    save_results(sol, p, tables, dir) -> dir

summary.csv, profiles.csv, wick.csv, geometry.csv を dir に書き出す (dir が無ければ作る)。
"""
function save_results(sol, p::AbstractLHPParams, tables::FluidTables, dir::AbstractString)
    mkpath(dir)
    CSV.write(joinpath(dir, "summary.csv"), results_summary(sol, p, tables))
    CSV.write(joinpath(dir, "profiles.csv"), results_profiles(sol, p))
    CSV.write(joinpath(dir, "wick.csv"), results_wick(sol, p))
    CSV.write(joinpath(dir, "geometry.csv"), results_geometry(p))
    println("結果を保存しました: ", abspath(dir), " (summary.csv, profiles.csv, wick.csv, geometry.csv)")
    return dir
end