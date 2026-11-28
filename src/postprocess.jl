# =============================================================================
# postprocess.jl (後処理: エネルギー収支)
#
# solve_transient_lhp! の結果 (sol, p) から、系全体のエネルギー収支を確認する。
#
#   ΔE = E(t_end) − E(t_init)           … 蓄積したエネルギーの変化
#   Q_in  = ∫ Q_apply dt                 … 印加熱
#   Q_out = ∫ (周囲への放熱 + 凝縮器プレートへの放熱) dt
#   残差 = ΔE − (Q_in − Q_out)
#
# 流体のエネルギーは保存形 (式3) の離散化で、面のフラックスを隣り合うセルが共有する
# ので、厳密に保存する。残差は、保存間隔 (saveat) の粗さによる台形則の積分誤差と、
# IDA の許容誤差による。残差が入力に対して数 % を超えるなら、どこかでエネルギーが
# 生成・消滅している。
#
# 【蓄積エネルギー E の定義】
#   流体: Σ ρm Em V      (Em = e + u²/2 + g z。保存形の左辺そのもの)
#   ウィック: Σ E_i       (式49)
#   壁: m ∫_{T_ref}^{T} c(T') dT'   (c は温度の関数なので、C(T)·T とは違う。
#       C(T)·T だと、温度が変わると T0 ΔC の誤差が出る: 凝縮器の壁が約9 K 冷えるとき約35 J)
# 基準温度 T_ref は、差だけが意味を持つので、何でもよい。
# =============================================================================

"""
    wall_internal_energy(T, V_w; T_ref=273.15, n=20) -> E [J]

SUS304 の壁 (体積 V_w) の内部エネルギー m ∫ c dT (台形則、n 分割)。
"""
function wall_internal_energy(T, V_w; T_ref=273.15, n::Int=20)
    Ts = range(T_ref, T, length=n + 1)
    f = [ρ_sus304(t) * c_sus304(t) for t in Ts]
    Iint = sum(0.5 * (f[k] + f[k+1]) * (Ts[k+1] - Ts[k]) for k in 1:n)
    return V_w * Iint
end

"""
    stored_energy(c, p) -> E [J]

状態 c (ComponentArray) に蓄えられている、系全体のエネルギー。
"""
function stored_energy(c, p::AbstractLHPParams)
    E = 0.0
    for (comp, pc) in ((:grv, p.p_grv), (:vl, p.p_vl), (:cond, p.p_cond), (:ll, p.p_ll), (:cc, p.p_cc))
        cm = getproperty(c, comp)
        E += sum(cm.ρm .* cm.Em) * pc.V
    end
    E += sum(c.wick.E)
    :swick in keys(c) && (E += sum(c.swick.E))   # 平板型: secondary wick
    for (wall, wp) in ((:vl, p.wp_vl), (:cond, p.wp_cond), (:ll, p.wp_ll))
        E += sum(wall_internal_energy(T, wp.V_w) for T in getproperty(c.T_wall, wall))
    end
    E += wall_internal_energy(c.T_wall.cc[1], p.cw.V_w)
    if p.cond.ec_2node    # 蒸発器ケース 2 節点: P は一定の熱容量 C_ec_P、S は C_ec_S
        E += p.cond.C_ec_P * (c.T_wall.evap[1] - 273.15) + p.cond.C_ec_S * (c.eca.T[1] - 273.15)
    else
        E += wall_internal_energy(c.T_wall.evap[1], p.ec.V_w)
        hasproperty(p.cond, :C_ec_extra) && (E += p.cond.C_ec_extra * (c.T_wall.evap[1] - 273.15))   # 蒸発器ケースの追加熱容量
    end
    hasproperty(p.cond, :C_cc_extra) && (E += p.cond.C_cc_extra * ((p.cond.cc_2node ? c.cca.T[1] : c.T_wall.cc[1]) - 273.15))   # CC ケースの追加熱容量 (2 節点では節点 A)
    :hb in keys(c.T_wall) && (E += p.cond.C_hb * (c.T_wall.hb[1] - 273.15))   # 平板型: ヒーターブロック
    :plate in keys(c) && (E += sum(plate_cell_params(p.cond, p.geo.plate).C .* (c.plate.T .- 273.15)))   # 凝縮器プレート
    n_cplate(p.geo) > 0 && (E += sum(cplate_cell_params(p.cond, p.geo.plate, p.geo.cplate).C .* (c.cplate.T .- 273.15)))   # コールドプレート
    E += housing_energy(c, p.cond)    # 断熱材・PEEK
    return E
end

"""
    stored_energy_parts(c, p) -> NamedTuple [J]

stored_energy を部位ごとに分けたもの (合計は stored_energy と一致する)。
  HB / EC ケース / CC ケース / 蒸気管・凝縮器・液管の壁 / グルーブ・蒸気管・凝縮器・液管・CC の流体 / 1次ウィック / 2次ウィック
"""
function stored_energy_parts(c, p::AbstractLHPParams)
    fl(comp, pc) = sum(getproperty(c, comp).ρm .* getproperty(c, comp).Em) * pc.V
    wl(wall, wp) = sum(wall_internal_energy(T, wp.V_w) for T in getproperty(c.T_wall, wall))
    E_cc = wall_internal_energy(c.T_wall.cc[1], p.cw.V_w)
    hasproperty(p.cond, :C_cc_extra) && (E_cc += p.cond.C_cc_extra * ((p.cond.cc_2node ? c.cca.T[1] : c.T_wall.cc[1]) - 273.15))
    if p.cond.ec_2node
        E_ec = p.cond.C_ec_P * (c.T_wall.evap[1] - 273.15) + p.cond.C_ec_S * (c.eca.T[1] - 273.15)
    else
        E_ec = wall_internal_energy(c.T_wall.evap[1], p.ec.V_w)
        hasproperty(p.cond, :C_ec_extra) && (E_ec += p.cond.C_ec_extra * (c.T_wall.evap[1] - 273.15))
    end
    E_hb = :hb in keys(c.T_wall) ? p.cond.C_hb * (c.T_wall.hb[1] - 273.15) : 0.0
    return (HB = E_hb, EC_case = E_ec, CC_case = E_cc,
            wall_vl = wl(:vl, p.wp_vl), wall_cond = wl(:cond, p.wp_cond), wall_ll = wl(:ll, p.wp_ll),
            plate = :plate in keys(c) ? sum(plate_cell_params(p.cond, p.geo.plate).C .* (c.plate.T .- 273.15)) : 0.0,
            cplate = n_cplate(p.geo) > 0 ? sum(cplate_cell_params(p.cond, p.geo.plate, p.geo.cplate).C .* (c.cplate.T .- 273.15)) : 0.0,
            fl_grv = fl(:grv, p.p_grv), fl_vl = fl(:vl, p.p_vl), fl_cond = fl(:cond, p.p_cond),
            fl_ll = fl(:ll, p.p_ll), fl_cc = fl(:cc, p.p_cc),
            wick = sum(c.wick.E), swick = :swick in keys(c) ? sum(c.swick.E) : 0.0,
            housing = housing_energy(c, p.cond))
end

"""
    energy_partition(sol, p; t0, ts) -> Nothing

時刻 t0 から各時刻 ts [s] までに、印加熱がどこへ入ったか (部位ごとの蓄積エネルギーの増分 [J]) を表で表示する。
同時に、印加熱の積分と放熱の積分、その差 (= 蓄積の合計に等しいはず) を出す。
流体の増分には、顕熱と蒸発・凝縮の潜熱 (相変化に伴う内部エネルギー変化) が両方入っている。
EC の昇温が実験より速いとき、足りない熱容量 (熱が入っているはずの部位) を探すのに使う。
"""
function energy_partition(sol, p::AbstractLHPParams; t0::Real, ts::AbstractVector{<:Real}, tables = nothing)
    idx(t) = argmin(abs.(sol.t .- t))
    k0 = idx(t0)
    P0 = stored_energy_parts(components_at(sol, p, k0), p)
    names = keys(P0)
    # 印加熱・放熱の積分 (台形則、保存点)
    Qin_t  = [applied_heat(t, p.cond) for t in sol.t]
    Qout_t = [heat_out(components_at(sol, p, k), p) for k in 1:length(sol.t)]
    cum(f, k) = sum(0.5 * (f[j] + f[j+1]) * (sol.t[j+1] - sol.t[j]) for j in k0:k-1; init = 0.0)
    println("\n熱の行き先: 加熱開始 t0 = ", t0, " s からの蓄積エネルギーの増分 [J]")
    println(rpad("部位", 11), join([lpad("+" * string(round(Int, t - t0)) * " s", 9) for t in ts]))
    rows = Vector{Vector{Float64}}()
    for nm in names
        row = Float64[]
        for t in ts
            Pk = stored_energy_parts(components_at(sol, p, idx(t)), p)
            push!(row, getproperty(Pk, nm) - getproperty(P0, nm))
        end
        push!(rows, row)
        println(rpad(string(nm), 11), join([lpad(string(round(v, digits = 2)), 9) for v in row]))
    end
    tot = [sum(r[i] for r in rows) for i in eachindex(ts)]
    Qin  = [cum(Qin_t, idx(t)) for t in ts]
    Qout = [cum(Qout_t, idx(t)) for t in ts]
    # 流体の各行は、液が部位の間を移動するときの内部エネルギーの出入りを含む (熱ではない)。部位をまたぐ移動は合計で消えるので、流体の合計を別に出す
    ifl = [i for (i, nm) in enumerate(names) if startswith(string(nm), "fl_")]
    println(rpad("流体 合計", 11), join([lpad(string(round(sum(rows[i][j] for i in ifl), digits = 2)), 9) for j in eachindex(ts)]))
    println(rpad("合計 ΔE", 11), join([lpad(string(round(v, digits = 2)), 9) for v in tot]))
    println(rpad("∫Q_apply", 11), join([lpad(string(round(v, digits = 2)), 9) for v in Qin]))
    println(rpad("∫Q_out", 11), join([lpad(string(round(v, digits = 2)), 9) for v in Qout]))
    println(rpad("残差", 11), join([lpad(string(round(tot[i] - (Qin[i] - Qout[i]), digits = 2)), 9) for i in eachindex(ts)]))

    # --- 蒸気の質量 (Σ ρm x V) の t0 からの増分 [mg]: 蒸発した量の目安。蒸気 1 mg の潜熱は約 1.17 J、体積は約 0.14 mL (蒸気管の液を押し出す量) ---
    vap(c) = (grv = sum(c.grv.ρm .* c.grv.x) * p.p_grv.V, vl = sum(c.vl.ρm .* c.vl.x) * p.p_vl.V,
              cond = sum(c.cond.ρm .* c.cond.x) * p.p_cond.V, ll = sum(c.ll.ρm .* c.ll.x) * p.p_ll.V,
              cc = sum(c.cc.ρm .* c.cc.x) * p.p_cc.V)
    V0 = vap(components_at(sol, p, k0))
    println("\n蒸気の質量の増分 [mg] (t0 から。乾き度 x が質量基準のとき)")
    println(rpad("部位", 11), join([lpad("+" * string(round(Int, t - t0)) * " s", 9) for t in ts]))
    Vk = [vap(components_at(sol, p, idx(t))) for t in ts]
    for nm in keys(V0)
        println(rpad(string(nm), 11), join([lpad(string(round(1e6 * (getproperty(v, nm) - getproperty(V0, nm)), digits = 2)), 9) for v in Vk]))
    end
    dV_tot = [1e6 * sum(getproperty(v, nm) - getproperty(V0, nm) for nm in keys(V0)) for v in Vk]
    println(rpad("合計", 11), join([lpad(string(round(v, digits = 2)), 9) for v in dV_tot]))

    # --- 蒸気の収支: 蒸発した質量 ∫ṁ dt (CC → ウィックへの流量。蒸発器で液→蒸気になる量) と、凝縮した質量 = ∫ṁ dt − 蒸気の増分の合計 ---
    # 凝縮は、蒸気管の前線・壁、凝縮器、CC (圧縮による凝縮) のどこかで起きる。蒸気の増分の部位別の表と合わせて見る。
    if tables !== nothing
        mdot_k = [let c = components_at(sol, p, k)
                      c.cc.u[1] * tables.ρ_l_sat(c.cc.T[1]) * p.p_cc.A
                  end for k in 1:length(sol.t)]
        cumm(k) = sum(0.5 * (mdot_k[j] + mdot_k[j+1]) * (sol.t[j+1] - sol.t[j]) for j in k0:k-1; init = 0.0)
        m_ev = [1e6 * cumm(idx(t)) for t in ts]
        println(rpad("∫ṁ dt", 11), join([lpad(string(round(v, digits = 2)), 9) for v in m_ev]), "   (蒸発した質量 [mg]。1 mg の潜熱 ≈ 1.17 J)")
        println(rpad("凝縮した量", 11), join([lpad(string(round(m_ev[i] - dV_tot[i], digits = 2)), 9) for i in eachindex(ts)]),
                "   (= ∫ṁ dt − 蒸気の増分の合計 [mg]。蒸気管・凝縮器・CC で凝縮した質量)")
    end
    # 前線の位置と蒸気の体積: x > 0.5 のセル数 (蒸気管、凝縮器) と、蒸気管の蒸気の体積
    println("\n蒸気管・凝縮器の蒸気セル数 (x > 0.5):")
    println(rpad("部位", 11), join([lpad("+" * string(round(Int, t - t0)) * " s", 9) for t in [t0; ts]]))
    cks = [components_at(sol, p, idx(t)) for t in [t0; ts]]
    println(rpad("vl (/" * string(p.p_vl.N) * ")", 11), join([lpad(string(count(>(0.5), c.vl.x)), 9) for c in cks]))
    println(rpad("cond (/" * string(p.p_cond.N) * ")", 11), join([lpad(string(count(>(0.5), c.cond.x)), 9) for c in cks]))
    return nothing
end

"""
    ll_end_budget(sol, p, tables; ts, cell = :end) -> Nothing

液管の CC 側の末端セル (二相セルを置くセル、LL3 の下流) の、質量と熱の収支を、時刻 ts [s] ごとに表示する (診断用、物理は変えない)。

  x, ρ, m   … 末端セルの乾き度、密度、質量 [mg]
  dm/dt     … 保存点の差分 (±1 保存点) で求めた質量の変化率 [mg/s]
  ṁ_in      … 上流の面 (セル N−1 との間) から入る質量流量 [mg/s]  (u>0 で上流から入る。風上の密度)
  ṁ_out     … CC 側の面から出る質量流量 [mg/s]                   (u>0 で CC へ出る。u<0 なら CC から入る)
  Σ         … ṁ_in − ṁ_out − dm/dt (質量保存の確認。0 に近いはず)
  dx/dt     … 乾き度の変化率 [1/s]
  壁の熱     … 壁の温度の変化から、壁の熱収支を逆算 [mW]:
                 C_w dT_w/dt (壁の蓄熱) = Q_amb + Q_cond,L + Q_br + (Q_fluid + Q_cond,R)
               Q_amb: 周囲 (断熱材つきの自然対流) との熱、Q_cond,L: 上流の壁との軸方向熱伝導、Q_br: 熱橋。
               残り = Q_fluid + Q_cond,R (流体との熱 + CC 側の壁へ) を逆算する。正が壁が受ける向き
               (残り < 0 なら、壁が流体へ熱を渡している = 蒸気が凝縮するほうに働く。ただし CC 側の壁との伝導が混ざる)。
"""
function ll_end_budget(sol, p::AbstractLHPParams, tables; ts::AbstractVector{<:Real})
    N = p.p_ll.N
    A, V = p.p_ll.A, p.p_ll.V
    wp = p.wp_ll
    amb = pipe_amb_spec(p.cond, :ll)
    idx(t) = argmin(abs.(sol.t .- t))
    println("\n液管の末端セル (セル ", N, ") の質量と熱の収支 (診断):")
    println("  ", lpad("t", 7), lpad("x", 9), lpad("ρ", 8), lpad("m[mg]", 8), lpad("dm/dt", 8), lpad("ṁ_in", 8), lpad("ṁ_out", 8), lpad("Σ", 8),
            lpad("dx/dt", 10), lpad("T-Tsat[mK]", 11), lpad("C dTw/dt", 10), lpad("Q_amb", 8), lpad("Q_cond,L", 9), lpad("Q_br", 8), lpad("残り", 8))
    for t in ts
        k = idx(t)
        k1, k2 = max(1, k - 1), min(length(sol.t), k + 1)
        k1 == k2 && continue
        c  = components_at(sol, p, k)
        c1 = components_at(sol, p, k1)
        c2 = components_at(sol, p, k2)
        dt = sol.t[k2] - sol.t[k1]
        m(cc) = 1e6 * cc.ll.ρm[N] * V
        dmdt = (m(c2) - m(c1)) / dt
        dxdt = (c2.ll.x[N] - c1.ll.x[N]) / dt
        ρ_up_cc = :ccl in keys(c) ? tables.ρ_l_sat(c.ccl.Tl[end]) : c.cc.ρm[1]    # 新しい CC: 液管がつながる最下層の液の密度
        # 面 N−1|N (セル N−1 の下流面、速度 u[N−1]) と面 N|CC (速度 u[N])
        u_in  = c.ll.u[N-1]
        u_out = c.ll.u[N]
        ṁ_in  = 1e6 * A * u_in  * (u_in  >= 0 ? c.ll.ρm[N-1] : c.ll.ρm[N])
        ṁ_out = 1e6 * A * u_out * (u_out >= 0 ? c.ll.ρm[N]   : ρ_up_cc)
        Σ = ṁ_in - ṁ_out - dmdt
        dTsat = 1e3 * (c.ll.T[N] - Tsat(c.ll.P[N], tables))
        # 壁の熱収支の逆算
        Tw = c.T_wall.ll
        dTwdt = (c2.T_wall.ll[N] - c1.T_wall.ll[N]) / dt
        Cw = wall_heat_capacity(Tw[N], wp)
        Q_amb = Q_amb_pipe(Tw[N], wp, ll_env_temps(p.cond, N, (p.cond.hi_on && :hi in keys(c)) ? c.hi.T[3] : p.cond.T_amb)[N]; amb = amb)
        Q_L = wall_axial_conductance(Tw[N-1], wp, Tw[N], wp) * (Tw[N-1] - Tw[N])
        T_src = p.cond.br_src === :hb ? c.T_wall.hb[1] : c.T_wall.evap[1]
        Q_br = thermal_bridge_ll(p.cond, T_src, Tw)[N]
        rest = Cw * dTwdt - Q_amb - Q_L - Q_br
        f(v, d = 3) = lpad(string(round(v, digits = d)), 8)
        println("  ", lpad(string(round(sol.t[k], digits = 2)), 7), lpad(string(round(c.ll.x[N], digits = 5)), 9), f(c.ll.ρm[N], 1), f(m(c), 2), f(dmdt), f(ṁ_in), f(ṁ_out), f(Σ),
                lpad(string(round(dxdt, sigdigits = 3)), 10), lpad(string(round(dTsat, digits = 1)), 11),
                lpad(string(round(1e3 * Cw * dTwdt, digits = 2)), 10), f(1e3 * Q_amb, 2), lpad(string(round(1e3 * Q_L, digits = 2)), 9), f(1e3 * Q_br, 2), f(1e3 * rest, 2))
    end
    return nothing
end

"""
    cc_interface_diag(sol, p, tables; t0, ts, dT_vl_res = nothing, dT_vl_res_exp = nothing) -> Nothing

平板型の新しい CC (気相 1 セル + 液 16 層) の、気液界面の診断 (診断用、物理は変えない)。

  f_v, s        … 気相の体積分率と、気相のキャップの高さ s = R − h_liq [mm]
  A_if          … 液面 (弦の幅 × 長さ) の面積 [cm²]
  u_LL, Re_j    … 液管出口の面の速度 [mm/s] と、噴流のレイノルズ数 (ρ u d / μ、d = 液管の内径)
  q_in          … 液の流入量 A u_LL [mm³/s]
  T1−T0         … 層 1 (気液界面側) の温度の、加熱前 (τ = 0) からの上昇 [K]。T_v−T1 は、界面の運動論的な温度差 [mK] (ΔP_ne の原因ではなく、ほぼ 0)
  ΔP_ne         … P_cc − P_sat(T1) [kPa] (気相の過剰圧力。層 1 の温度が上がった分)
  ṁ_cond, Q_c   … 界面の凝縮量 [mg/s] と、その潜熱 ṁ λ [W]。E_in は、その積算 [J]
  T_lN, T_mid   … 最下層 (流入側) と中間層 (層 8) の温度の、τ = 0 からの変化 [K]
  VL1−RES       … 壁の温度差 (加熱前 τ = 0 からの変化) [K]。計算と実験。×dP_sat/dT で圧力に換算 [kPa]
"""
function cc_interface_diag(sol, p::AbstractLHPParams, tables; t0::Real, ts::AbstractVector{<:Real},
                           dT_vl_res = nothing, dT_vl_res_exp = nothing)
    if !(:ccl in keys(components_at(sol, p, 1)))
        println("\ncc_interface_diag: 新しい CC (状態 ccl) の結果ではありません。スキップします")
        return nothing
    end
    A_ll = p.p_ll.A
    d_ll = sqrt(4 * A_ll / π)
    ks = [k for k in eachindex(sol.t) if sol.t[k] >= t0]
    tk = sol.t[ks]
    rows = map(ks) do k
        c = components_at(sol, p, k)
        dg = cc_state_diag(c, p, tables)
        sp = sat_props(dg.T_v, tables)
        u = c.ll.u[end]
        (dg = dg, u = u, q_in = A_ll * u, Re = sp.ρ_l * abs(u) * d_ll / sp.μ_l, Q = dg.ṁ_cond * dg.λ,
         dPdT = sp.h_fg * sp.ρ_v * sp.ρ_l / (dg.T_v * (sp.ρ_l - sp.ρ_v)), T_mid = c.ccl.Tl[max(1, N_CC_LAYER ÷ 2)],
         s = p.p_cc.R - dg.h_liq, dT = dT_vl_res === nothing ? NaN : dT_vl_res(c))
    end
    Q = [r.Q for r in rows]
    E_in = vcat(0.0, cumsum((Q[2:end] .+ Q[1:end-1]) ./ 2 .* diff(tk)))
    r0 = rows[1]
    dT0 = r0.dT
    idx(τ) = argmin(abs.(tk .- (t0 + τ)))
    f1(x, d = 2) = lpad(string(round(x, digits = d)), 8)
    println("\nCC 気液界面の診断 (加熱開始からの時間 τ [s]。診断のみ、物理は変えない):")
    println("  液管の内径 d = ", round(1e3 * d_ll, digits = 2), " mm,  A_if(τ=0) = ", round(1e4 * r0.dg.A_if, digits = 2), " cm²,  V_v(τ=0) = ",
            round(1e6 * r0.dg.f_v * p.p_cc.V, digits = 3), " mL,  h_i(τ=0) = ", round(r0.dg.h_i, sigdigits = 3), " W/m²K (σ = ", round(r0.dg.σ, digits = 3), ")")
    println("  ", lpad("τ", 4), lpad("f_v", 8), lpad("s[mm]", 8), lpad("A_if", 8), lpad("u_LL", 8), lpad("Re_j", 8), lpad("q_in", 8),
            lpad("T1-T0", 8), lpad("Tv-T1[mK]", 10), lpad("ΔP_ne", 8), lpad("ṁ_c[mg/s]", 10), lpad("Q_c[W]", 8), lpad("Q_w[W]", 8), lpad("g_T", 6), lpad("E_in", 8),
            lpad("T_lN-T0", 9), lpad("T_mid-T0", 9), lpad("VL1-RES", 9), lpad("換算", 8), lpad("実験", 8), lpad("実験換算", 9))
    for τ in ts
        t0 + τ > sol.t[end] && continue
        i = idx(τ)
        r = rows[i]
        dg = r.dg
        dTm = isnan(r.dT) ? NaN : r.dT - dT0
        dTe = dT_vl_res_exp === nothing ? NaN : dT_vl_res_exp(τ)
        println("  ", lpad(τ, 4), f1(dg.f_v, 3), f1(1e3 * r.s), f1(1e4 * dg.A_if), f1(1e3 * r.u, 1), f1(r.Re, 0), f1(1e9 * r.q_in, 1),
                f1(dg.T1 - r0.dg.T1, 3), lpad(string(round(1e3 * (dg.T_v - dg.T1), digits = 3)), 10), f1(dg.dP_ne / 1e3, 2),
                lpad(string(round(1e6 * dg.ṁ_cond, digits = 3)), 10), f1(r.Q, 3), f1(dg.Q_wall, 3), f1(dg.g_T, 2)[3:end], f1(E_in[i], 3),
                f1(dg.TN - r0.dg.TN, 3), f1(r.T_mid - r0.T_mid, 3),
                f1(dTm, 3), f1(dTm * r.dPdT / 1e3, 1), f1(dTe, 3), f1(dTe * r.dPdT / 1e3, 1))
    end
    return nothing
end

"""
    loop_heat_diag(sol, p, tables; t0, ts) -> Nothing

ループ全体の熱の収支 [W] を、加熱開始からの時間 τ ごとに表示する (診断用、物理は変えない)。前線が通った後の飽和温度が実験より高い原因を調べる。

  Q_apply, Q_plate, Q_amb … 印加熱、凝縮器プレート → コールドプレートへ出る熱、蒸気管・液管から周囲へ出る熱
  dE/dt の内訳           … 蓄積エネルギーの増加率 [W]: HB、蒸発器ケース、CC ケース、管の壁 (蒸気管 + 凝縮器 + 液管)、
                             流体 (グルーブ + 蒸気管 + 凝縮器 + 液管 + CC)、ウィック (1次 + 2次)、プレート
  Σ                      … Q_apply − Q_plate − Q_amb − Σ dE/dt (収支の確認。0 に近いはず)
  Tsat−CON1, Tsat−Tpl    … P_cc の飽和温度と、凝縮器の最初のセルの壁の温度の差、プレートの最高温度との差 [K]
  蒸気セル               … 蒸気管 / 凝縮器の蒸気セル数 (x > 0.5)
"""
function loop_heat_diag(sol, p::AbstractLHPParams, tables; t0::Real, ts::AbstractVector{<:Real})
    cnd = p.cond
    K = 273.15
    idx(t) = argmin(abs.(sol.t .- t))
    f1(x, d = 2) = lpad(string(round(x, digits = d)), 8)
    println("\nループの熱の収支 [W] (加熱開始からの時間 τ [s]。診断のみ):")
    println("  ", lpad("τ", 4), lpad("Q_apply", 8), lpad("Q_plate", 8), lpad("Q_amb", 8), lpad("dE/dt", 8), lpad("Σ", 8),
            lpad("HB", 8), lpad("EC", 8), lpad("CCケース", 9), lpad("管の壁", 8), lpad("流体", 8), lpad("ウィック", 9), lpad("プレート", 9),
            lpad("Tsat-CON1", 10), lpad("Tsat-Tpl", 9), lpad("T_pl最高", 9), lpad("T_pl中心", 9), lpad("T_cp最高", 9), lpad("T_cp中心", 9), lpad("蒸気セル", 10))
    for τ in ts
        k = idx(t0 + τ)
        k1, k2 = max(1, k - 1), min(length(sol.t), k + 1)
        k1 == k2 && continue
        c  = components_at(sol, p, k)
        ca = components_at(sol, p, k1)
        cb = components_at(sol, p, k2)
        dt = sol.t[k2] - sol.t[k1]
        Ea, Eb = stored_energy_parts(ca, p), stored_energy_parts(cb, p)
        dE(f) = (f(Eb) - f(Ea)) / dt
        d_hb   = dE(e -> e.HB)
        d_ec   = dE(e -> e.EC_case)
        d_cc   = dE(e -> e.CC_case)
        d_wall = dE(e -> e.wall_vl + e.wall_cond + e.wall_ll)
        d_fl   = dE(e -> e.fl_grv + e.fl_vl + e.fl_cond + e.fl_ll + e.fl_cc)
        d_wick = dE(e -> e.wick + e.swick)
        d_pl   = dE(e -> e.plate + e.cplate + e.housing)
        dEdt = d_hb + d_ec + d_cc + d_wall + d_fl + d_wick + d_pl
        Qa = applied_heat(sol.t[k], cnd)
        Q_plate = n_cplate(p.geo) > 0 ? sum(cplate_pipe_G(cnd, p.geo.cplate) * (c.cplate.T[kk] - cnd.T_cond) for kk in p.cplate_nodes) :
                  (:plate in keys(c) ? sum(plate_cell_params(cnd, p.geo.plate).Gc * (T - cnd.T_cond) for T in c.plate.T) : NaN)   # 冷却液へ出る熱 (コールドプレート (2次元) があれば配管から、なければ接触面から)
        Tcp_max = n_cplate(p.geo) > 0 ? maximum(c.cplate.T) - K : NaN
        Tcp_c = NaN
        if n_cplate(p.geo) > 0
            pg = p.geo.plate
            ic = (cld(pg.nx, 2), cld(pg.nx, 2) + 1 - (pg.nx % 2))
            jc = (cld(pg.ny, 2), cld(pg.ny, 2) + 1 - (pg.ny % 2))
            Tcp_c = sum(c.cplate.T[(j - 1) * pg.nx + i] for i in ic for j in jc) / (length(ic) * length(jc)) - K   # 118 (コールドプレート下面の中心) に対応
        end
        Q_amb = 0.0
        for (wall, wp) in ((:vl, p.wp_vl), (:ll, p.wp_ll))
            amb = pipe_amb_spec(cnd, wall)
            Tw_pipe = getproperty(c.T_wall, wall)
            T_env = wall === :ll ? ll_env_temps(cnd, length(Tw_pipe), (cnd.hi_on && :hi in keys(c)) ? c.hi.T[3] : cnd.T_amb) : fill(cnd.T_amb, length(Tw_pipe))
            for (j, T) in enumerate(Tw_pipe)
                Q_amb -= Q_amb_pipe(T, wp, T_env[j]; amb = amb)
            end
        end
        Ts = Tsat(c.cc.P[1], tables) - K
        Tpl = :plate in keys(c) ? maximum(c.plate.T) - K : NaN
        # プレート中央 (x = L/2, y = W/2) の 4 節点の平均。実験の 117_cond (プレート上面の中央) と比べる
        Tpl_c = NaN
        if :plate in keys(c)
            pg = p.geo.plate
            ic = (cld(pg.nx, 2), cld(pg.nx, 2) + 1 - (pg.nx % 2))
            jc = (cld(pg.ny, 2), cld(pg.ny, 2) + 1 - (pg.ny % 2))
            Tpl_c = sum(c.plate.T[(j - 1) * pg.nx + i] for i in ic for j in jc) / (length(ic) * length(jc)) - K
        end
        println("  ", lpad(τ, 4), f1(Qa), f1(Q_plate), f1(Q_amb), f1(dEdt), f1(Qa - Q_plate - Q_amb - dEdt),
                f1(d_hb), f1(d_ec), lpad(string(round(d_cc, digits = 2)), 9), f1(d_wall), f1(d_fl), lpad(string(round(d_wick, digits = 2)), 9), lpad(string(round(d_pl, digits = 2)), 9),
                lpad(string(round(Ts - (c.T_wall.cond[1] - K), digits = 2)), 10), lpad(string(round(Ts - Tpl, digits = 2)), 9), lpad(string(round(Tpl, digits = 2)), 9), lpad(string(round(Tpl_c, digits = 2)), 9), lpad(string(round(Tcp_max, digits = 2)), 9), lpad(string(round(Tcp_c, digits = 2)), 9),
                lpad(string(count(>(0.5), c.vl.x)) * "/" * string(count(>(0.5), c.cond.x)), 10))
    end
    return nothing
end

"""
    ll_thermal_diag(sol, p, tables; t0, ts, s_LL1, s_LL2) -> Nothing

液管 (LL) の温度の前線を、液の前線と比べる診断 (診断用、物理は変えない)。LL1 の温度低下が早く、LL2 が遅い原因を調べる。

  s           … 液が進んだ距離 [m] (中央のセルの速度の積算)
  x_th        … 熱の前線の位置 [m] (凝縮器側から。初期より 0.3 K 以上冷えたセルのうち、CC に最も近いセルの中心)
  x_th/s      … 熱の前線の進みの遅れ (1 なら液と同じ速さ)。壁の熱容量が流体より大きいと 1 より小さくなる
  C_w/C_f     … 液管 1 セルの、壁と流体の熱容量の比。熱の前線の速さ ≈ 液の速さ × 1/(1 + C_w/C_f) (壁と流体が密に結合しているとき)
  T_f, T_w    … LL1、LL2 の位置の、流体と壁の温度 [°C]
  Q_amb, Q_br … LL2 のセルの、周囲から壁へ入る熱と熱橋の熱 [mW]
"""
function ll_thermal_diag(sol, p::AbstractLHPParams, tables; t0::Real, ts::AbstractVector{<:Real}, s_LL1::Real, s_LL2::Real)
    cnd = p.cond
    K = 273.15
    N = p.p_ll.N
    L = N * p.p_ll.dl
    wp = p.wp_ll
    amb = pipe_amb_spec(cnd, :ll)
    ks = [k for k in eachindex(sol.t) if sol.t[k] >= t0]
    tk = sol.t[ks]
    cs = [components_at(sol, p, k) for k in ks]
    um = [c.ll.u[cld(N, 2)] for c in cs]
    sk = vcat(0.0, cumsum((um[2:end] .+ um[1:end-1]) ./ 2 .* diff(tk)))
    c0 = cs[1]
    pos = [(i - 0.5) / N for i in 1:N]
    interp(v, s) = begin
        s <= pos[1] && return v[1]
        s >= pos[end] && return v[end]
        j = searchsortedlast(pos, s)
        v[j] + (v[j+1] - v[j]) * (s - pos[j]) / (pos[j+1] - pos[j])
    end
    f1(x, d = 2) = lpad(string(round(x, digits = d)), 8)
    i2 = clamp(round(Int, s_LL2 * N + 0.5), 1, N)
    sp0 = sat_props(c0.ll.T[i2], tables)
    C_f = sp0.ρ_l * sp0.cp_l * p.p_ll.V
    C_w = wall_heat_capacity(c0.T_wall.ll[i2], wp)
    println("\n液管の温度の前線と液の前線 (加熱開始からの時間 τ [s]。診断のみ):")
    println("  液管 1 セルの熱容量: 流体 ", round(C_f, digits = 3), " J/K,  壁 ", round(C_w, digits = 3), " J/K  (C_w/C_f = ", round(C_w / C_f, digits = 2), ")")
    println("  ", lpad("τ", 4), lpad("s[m]", 8), lpad("x_th[m]", 8), lpad("x_th/s", 8), lpad("u_mid", 8),
            lpad("LL1 T_f", 8), lpad("LL1 T_w", 8), lpad("LL2 T_f", 8), lpad("LL2 T_w", 8), lpad("Q_amb", 8), lpad("Q_br", 8))
    for τ in ts
        t0 + τ > sol.t[end] && continue
        i = argmin(abs.(tk .- (t0 + τ)))
        c = cs[i]
        ith = 0
        for j in 1:N
            c.ll.T[j] < c0.ll.T[j] - 0.3 && (ith = j)
        end
        x_th = ith == 0 ? 0.0 : pos[ith] * L
        T_src = cnd.br_src === :hb ? c.T_wall.hb[1] : c.T_wall.evap[1]
        Q_br = thermal_bridge_ll(cnd, T_src, c.T_wall.ll)[i2]
        Q_amb = Q_amb_pipe(c.T_wall.ll[i2], wp, ll_env_temps(cnd, length(c.T_wall.ll), (cnd.hi_on && :hi in keys(c)) ? c.hi.T[3] : cnd.T_amb)[i2]; amb = amb)
        ratio = sk[i] > 1.0e-4 ? x_th / sk[i] : NaN
        println("  ", lpad(τ, 4), f1(sk[i]), f1(x_th), f1(ratio), f1(1e3 * um[i], 1),
                f1(interp(c.ll.T, s_LL1) - K), f1(interp(c.T_wall.ll, s_LL1) - K), f1(interp(c.ll.T, s_LL2) - K), f1(interp(c.T_wall.ll, s_LL2) - K),
                f1(1e3 * Q_amb, 2), f1(1e3 * Q_br, 2))
    end
    return nothing
end

"""
    startup_onset_diag(sol, p, tables; t0, ts) -> Nothing

起動 (グルーブの核生成) の立ち上がりの診断 (診断用、物理は変えない)。VL1・VL2・LL1 が実験より早い原因を調べる。
  EC−Tsat … 蒸発器ケースとグルーブの飽和温度の差 [K]。n, nc … 界面蒸発・気泡核生成の活性度。u_end … 液管の CC 側の速度 [mm/s]。
  ΔP      … P_cc − P_cc(τ = 0) [kPa]。VL1 − Tsat(P_cc) … 蒸気管の壁温度が飽和温度に追従する遅れ [K] (壁の熱容量による)。
"""
function startup_onset_diag(sol, p::AbstractLHPParams, tables; t0::Real, ts::AbstractVector{<:Real}, vl1_pos::Real = 1/6)
    K = 273.15
    idx(t) = argmin(abs.(sol.t .- t))
    c0 = components_at(sol, p, idx(t0))
    P0 = c0.cc.P[1]
    N = p.p_vl.N
    pos = [(i - 0.5) / N for i in 1:N]
    interp(v, s) = begin
        s <= pos[1] && return v[1]
        s >= pos[end] && return v[end]
        j = searchsortedlast(pos, s)
        v[j] + (v[j+1] - v[j]) * (s - pos[j]) / (pos[j+1] - pos[j])
    end
    f1(x, d = 3) = lpad(string(round(x, digits = d)), 8)
    println("\n起動の立ち上がり (加熱開始からの時間 τ [s]。診断のみ):")
    println("  ", lpad("τ", 4), lpad("EC-Tsat", 8), lpad("n", 8), lpad("nc", 8), lpad("u_end", 8), lpad("ΔP[kPa]", 8), lpad("VL1壁", 8), lpad("Tsat(Pcc)", 10), lpad("VL1-Tsat", 9))
    for τ in ts
        t0 + τ > sol.t[end] && continue
        k = idx(t0 + τ)
        c = components_at(sol, p, k)
        Ts = Tsat(c.grv.P[1], tables)
        Tsc = Tsat(c.cc.P[1], tables)
        Tv = interp(c.T_wall.vl, vl1_pos)
        println("  ", lpad(τ, 4), f1(c.T_wall.evap[1] - Ts), f1(maximum(c.nuc.n)), f1(maximum(c.nuc.nc)), f1(1e3 * c.ll.u[end], 1),
                f1((c.cc.P[1] - P0) / 1e3, 2), f1(Tv - K, 2), f1(Tsc - K, 2), f1(Tv - Tsc, 3))
    end
    return nothing
end

"""
    heat_out(c, p) -> Q [W]

系から出ていく熱。蒸気管・液管の壁→周囲 (自然対流)、凝縮器の壁→凝縮器プレート。
CCケース、蒸発器ケースは断熱。Q_amb_* は「壁が受け取る向き」が正なので、符号を反転する。
残差 (residuals_wall.jl) と同じ関数を使うので、残差が仮定する放熱と一致する。
"""
function heat_out(c, p::AbstractLHPParams)
    cnd = p.cond
    Q = 0.0
    T_hi3 = (cnd.hi_on && :hi in keys(c)) ? c.hi.T[3] : cnd.T_amb
    for (wall, wp) in ((:vl, p.wp_vl), (:ll, p.wp_ll))
        amb = pipe_amb_spec(cnd, wall)
        Tw_pipe = getproperty(c.T_wall, wall)
        T_env = wall === :ll ? ll_env_temps(cnd, length(Tw_pipe), T_hi3) : fill(cnd.T_amb, length(Tw_pipe))
        for (j, T) in enumerate(Tw_pipe)
            Q -= Q_amb_pipe(T, wp, T_env[j]; amb = amb)
        end
    end
    if n_cplate(p.geo) > 0    # コールドプレート (2次元): チラーの冷却液へ出る熱
        Gp = cplate_pipe_G(cnd, p.geo.cplate)
        Q += sum(Gp * (c.cplate.T[k] - cnd.T_cond) for k in p.cplate_nodes)
    elseif :plate in keys(c)    # 平板型: 凝縮器プレート→コールドプレート (T_cond 一定) へ出る熱
        pp = plate_cell_params(cnd, p.geo.plate)
        Q += sum(pp.Gc * (T - cnd.T_cond) for T in c.plate.T)
    else
        for T in c.T_wall.cond
            Q -= Q_amb_condenser(T, p.wp_cond, cnd.T_cond; h_con = cnd.h_con, f_con = cnd.f_con_area)
        end
    end
    Q += housing_heat_out(c, cnd)    # 断熱材 → 周囲、PEEK → 下板
    return Q
end

"""
    energy_balance(sol, p; verbose=true) -> NamedTuple

系全体のエネルギー収支。保存点 (condition.saveat) の値を、台形則で積分する。
戻り値: (ΔE, Q_in, Q_out, residual, E, Q_out_t)。E, Q_out_t は、保存点ごとの時系列。
"""
function energy_balance(sol, p::AbstractLHPParams; verbose::Bool=true, history::Bool=false, t_hist = nothing)
    n = length(sol.t)
    E = [stored_energy(components_at(sol, p, k), p) for k in 1:n]
    Q_in_t = [applied_heat(t, p.cond) + p.cond.Q_cc_fluid for t in sol.t]    # + 診断用の CC 流体への入熱 (既定 0)
    Q_out_t = [heat_out(components_at(sol, p, k), p) for k in 1:n]
    trapz(f) = sum(0.5 * (f[k] + f[k+1]) * (sol.t[k+1] - sol.t[k]) for k in 1:n-1)
    ΔE = E[end] - E[1]
    Q_in = trapz(Q_in_t)
    Q_out = trapz(Q_out_t)
    residual = ΔE - (Q_in - Q_out)
    if verbose
        println("エネルギー収支 (t = ", sol.t[1], " 〜 ", sol.t[end], " s):")
        println("  蓄積したエネルギーの変化 ΔE            = ", round(ΔE, digits=2), " J")
        println("  印加熱の積分 ∫Q_apply dt               = ", round(Q_in, digits=2), " J")
        println("  放熱の積分 ∫Q_out dt (周囲+凝縮器プレート) = ", round(Q_out, digits=2), " J")
        println("  残差 ΔE − (∫Q_in − ∫Q_out)             = ", round(residual, digits=2), " J",
                "  (入力に対して ", round(100 * residual / max(Q_in, eps()), digits=2), " %)")
    end
    if verbose && history
        # 残差の時間履歴: R(t) = [E(t) − E(t_0)] − [∫Q_in − ∫Q_out] (t_0 = 最初の保存点)。どの時刻の区間で残差が生じるかを見る。
        # 区間ごとの増分 dR も出す。IDA の許容誤差 (reltol) による誤差なら、加熱の立ち上がり (t_step〜t_step+τ_step) など、変化の速い区間に出る。
        ts_h = t_hist === nothing ? [t for t in (1.0, 5.0, 10.0, 12.0, 13.0, 13.5, 14.0, 14.5, 15.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0)
                                        if sol.t[1] <= t <= sol.t[end]] : t_hist
        cumt(f, k) = sum(0.5 * (f[j] + f[j+1]) * (sol.t[j+1] - sol.t[j]) for j in 1:k-1; init = 0.0)
        println("\nエネルギー残差の時間履歴 R(t) [J]  (t は計算の絶対時刻 [s]、dR は前の行からの増分):")
        println("  ", lpad("t", 7), lpad("R", 9), lpad("dR", 9), lpad("E-E0", 10), lpad("∫Qin", 9), lpad("∫Qout", 9))
        R_prev = 0.0
        for t in ts_h
            k = argmin(abs.(sol.t .- t))
            Rk = (E[k] - E[1]) - (cumt(Q_in_t, k) - cumt(Q_out_t, k))
            println("  ", lpad(string(round(sol.t[k], digits = 2)), 7), lpad(string(round(Rk, digits = 3)), 9), lpad(string(round(Rk - R_prev, digits = 3)), 9),
                    lpad(string(round(E[k] - E[1], digits = 2)), 10), lpad(string(round(cumt(Q_in_t, k), digits = 2)), 9), lpad(string(round(cumt(Q_out_t, k), digits = 2)), 9))
            R_prev = Rk
        end
    end
    return (ΔE = ΔE, Q_in = Q_in, Q_out = Q_out, residual = residual, E = E, Q_out_t = Q_out_t)
end

"""
    h_failure_diag(sol, p, tables; n_last = 3, cells_vl, cells_cond) -> Nothing

IDA が途中で止まったとき (ConvergenceFailure) に、最後に保存された数点で、蒸気管 (VL) と凝縮器のセルの
熱伝達 h の内訳を表示する診断 (診断のみ、物理は変えない)。
  u, G      … セル速度 [mm/s] (面の平均。入口側セルは u.u[1] で代用), G = ρ|u| [kg/m²s]
  x, x_ps   … フラッシュの乾き度, 疑似乾き度
  Tsat−T    … 過冷却度 [K],  Tw−T … 壁−流体 [K]
  Re_LS     … Shah の Re_LS (領域 II の h_Nu はこれの −1/3 乗)
  h_I, h_Nu … Shah の項 [W/m²K] (f 倍する前),  h_1ph … 単相の h
  w_2ph     … 二相側の重み,  h_tot … 実際に使われる h (f 倍後),  h_tot(f=1) … 参考
"""
function h_failure_diag(sol, p::AbstractLHPParams, tables; n_last::Int = 3,
                        cells_vl = 10:20, cells_cond = vcat(1:3, 44:50))
    c0 = p.cond
    g = 9.80665
    isempty(sol.t) && return nothing
    f1(x, d = 2) = lpad(string(round(x, digits = d)), 9)
    fe(x) = lpad(string(round(x, sigdigits = 3)), 10)
    ε_ph = c0.ε_phase_h
    for k in max(1, length(sol.t) - n_last + 1):length(sol.t)
        c = components_at(sol, p, k)
        println("\n[h の診断] t = ", round(sol.t[k], digits = 4), " s  (f_h_shah = ", c0.f_h_shah, ", ε_phase_h = ", ε_ph, ")")
        for (nm, comp, pp, cells) in (("vl", c.vl, p.p_vl, cells_vl), ("cond", c.cond, p.p_cond, cells_cond))
            println("  ", nm, ":  ", lpad("i", 3), lpad("u[mm/s]", 9), lpad("G", 9), lpad("x", 9), lpad("x_ps", 9),
                    lpad("Tsat-T", 9), lpad("Tw-T", 9), lpad("Re_LS", 10), lpad("h_I", 10), lpad("h_Nu", 10), lpad("h_1ph", 10), lpad("w_2ph", 9), lpad("h_tot", 10), lpad("h_tot(f=1)", 12))
            for i in cells
                i > pp.N && continue
                u_c = i == 1 ? comp.u[1] : 0.5 * (comp.u[i-1] + comp.u[i])
                ρ, Em = comp.ρm[i], comp.Em[i]
                G = ρ * soft_abs(u_c, 1.0e-4)
                em = internal_energy(Em, u_c, cell_elevation(pp, i), g)
                P, T, x = flash_ρe(ρ, em, tables)
                xps = c0.h_xps_mode === :P ? pseudo_quality_P(em, ρ, P, tables) : pseudo_quality(em, ρ, P, T, tables)
                props = sat_props(T, tables)
                Pr_ = soft_floor(P / tables.Pcrit, 1.0e-4)
                xs = quality_for_correlation(x)
                Re_LS = G * (1 - xs) * pp.D / props.μ_l
                Re_LT = G * pp.D / props.μ_l
                h_LT = h_LT_shah(props.k_l, pp.D, Re_LT, props.Pr_l)
                h_I = h_I_shah(h_LT, props.μ_l, props.μ_v, n_shah(Pr_), xs, Pr_)
                h_Nu = h_Nu_shah(Re_LS, props.ρ_l, props.ρ_v, g, props.k_l, props.μ_l)
                μ1, k1, cp1, Pr1 = transport_props(ρ, em, tables)
                h_1ph = h_from_Nu(Nu_pipe_smooth(reynolds_number_safe(ρ, u_c, pp.D, μ1), Pr1,
                                 dittus_boelter_exponent_smooth(getproperty(c.T_wall, Symbol(nm))[i], T); ε_Re = 200.0), k1, pp.D)
                Tw = getproperty(c.T_wall, Symbol(nm))[i]
                w_hi = smooth_step(xps - 1.0, ε_ph); w_lo = smooth_step(-xps, ε_ph)
                w2 = (1 - w_hi) * (1 - w_lo)
                h_tot = h_transport_wall(ρ, u_c, Em, cell_elevation(pp, i), Tw, pp.D, G, g, tables;
                                         ε_Re = 200.0, ε_phase = ε_ph, xps_mode = c0.h_xps_mode, f_cond = c0.f_h_shah)
                h_tot1 = h_transport_wall(ρ, u_c, Em, cell_elevation(pp, i), Tw, pp.D, G, g, tables;
                                          ε_Re = 200.0, ε_phase = ε_ph, xps_mode = c0.h_xps_mode, f_cond = 1.0)
                println("       ", lpad(i, 3), f1(u_c * 1e3, 3), fe(G), f1(x, 4), f1(xps, 4),
                        f1(Tsat(P, tables) - T, 3), f1(Tw - T, 3), fe(Re_LS), fe(h_I), fe(h_Nu), fe(h_1ph), f1(w2, 4), fe(h_tot), fe(h_tot1))
            end
        end
    end
    return nothing
end

"""
    failure_state_diag(sol, p, tables; n_last = 6) -> Nothing

IDA が止まったとき、最後の保存点 n_last 個で、各成分の状態の範囲を表示する診断 (診断のみ)。
どの成分の何が動いているか (速度、圧力、温度、乾き度、核生成、CC の ΔP_ne) を見る。
"""
function failure_state_diag(sol, p::AbstractLHPParams, tables; n_last::Int = 6)
    K = 273.15
    isempty(sol.t) && return nothing
    f(x, d = 3) = lpad(string(round(x, digits = d)), 9)
    println("\n[状態の診断] 最後の保存点 (成分ごとの範囲)")
    for k in max(1, length(sol.t) - n_last + 1):length(sol.t)
        c = components_at(sol, p, k)
        println("t = ", round(sol.t[k], digits = 4), " s")
        for nm in (:vl, :cond, :ll, :grv, :cc)
            comp = getproperty(c, nm)
            println("   ", rpad(string(nm), 5), " max|u|=", f(1e3 * maximum(abs.(comp.u))), " mm/s   P=[", f(minimum(comp.P) / 1e3, 2), ",", f(maximum(comp.P) / 1e3, 2),
                    "] kPa   T=[", f(minimum(comp.T) - K), ",", f(maximum(comp.T) - K), "]   x=[", f(minimum(comp.x), 4), ",", f(maximum(comp.x), 4), "]")
        end
        println("   nuc.n max=", f(maximum(c.nuc.n), 4), "  nc max=", f(maximum(c.nuc.nc), 4), "   ccl: ΔP_ne=", f((c.cc.P[1] - tables.Psat_T(c.ccl.Tl[1])) / 1e3, 4), " kPa  Tl1=", f(c.ccl.Tl[1] - K), "  TlN=", f(c.ccl.Tl[end] - K), "  X=", f(1e9 * c.ccl.X[1], 3), " mm³  mv=", f(1e6 * c.ccl.mv[1], 3), " mg",
                "   wick.T=[", f(minimum(c.wick.T) - K), ",", f(maximum(c.wick.T) - K), "]   swick.T=[", f(minimum(c.swick.T) - K), ",", f(maximum(c.swick.T) - K), "]")
        println("   T_wall: hb=", f(c.T_wall.hb[1] - K), " evap=", f(c.T_wall.evap[1] - K), " cc=", f(c.T_wall.cc[1] - K),
                "   EC−Tsat(P_grv)=", f(c.T_wall.evap[1] - Tsat(c.grv.P[1], tables), 4), " K")
    end
    return nothing
end

"""
    failure_residual_diag(sol, p; n_top = 15) -> Nothing

IDA が止まったとき、最後の保存点で、残差 F(t, u, du) を評価して、残差の大きい変数と、変化の速い変数を表示する診断 (診断のみ)。
du は、最後の 2 点の差分 (u_k − u_{k−1})/Δt で近似する。ステップが刻めなくなった原因の変数を探す。
"""
function failure_residual_diag(sol, p::AbstractLHPParams; n_top::Int = 15)
    length(sol.t) < 2 && return nothing
    k = length(sol.t)
    u1, u0 = sol.u[k], sol.u[k-1]
    dt = sol.t[k] - sol.t[k-1]
    du = (u1 .- u0) ./ dt
    res = zeros(length(u1))
    f_dae_residual!(res, du, copy(u1), p, sol.t[k])
    labels = variable_labels(p.axes, length(u1))
    println("\n[残差の診断] t = ", round(sol.t[k], digits = 4), " s,  Δt (差分) = ", round(dt, digits = 4), " s")
    println("  残差 |F| の大きい変数:   名前 | F | u | du(差分)")
    for i in sortperm(abs.(res); rev = true)[1:min(n_top, length(res))]
        println("    ", rpad(labels[i], 16), lpad(string(round(res[i], sigdigits = 4)), 12), lpad(string(round(u1[i], sigdigits = 7)), 14), lpad(string(round(du[i], sigdigits = 4)), 12))
    end
    rel = abs.(du) ./ (abs.(u1) .+ 1.0e-6 .* (1 .+ abs.(u1)))
    println("  変化の速い変数 (|du|/|u|):   名前 | u | du(差分) | |du|/|u| [1/s]")
    for i in sortperm(rel; rev = true)[1:min(n_top, length(rel))]
        println("    ", rpad(labels[i], 16), lpad(string(round(u1[i], sigdigits = 7)), 14), lpad(string(round(du[i], sigdigits = 4)), 12), lpad(string(round(rel[i], sigdigits = 3)), 12))
    end
    return nothing
end
# =============================================================================
# 毛細管圧の余裕 (平板型、事後の診断。物理は変えない)
#
# モデルは「primary wick と secondary wick は常に液で満たされている」と仮定していて、毛細管限界・ドライアウトは
# 残差に入っていない。保存点ごとに、毛細管圧が必要な圧力差を上回っているかを事後に確認する。
# 上面加熱 (orientation = :heat_top) では、液を CC の液面から primary wick の蒸発面まで、重力に逆らって持ち上げる。
#
# 【高さ】 CC の液面を基準 (0)、上向きを正。
#   z_s    : secondary wick の、CC 内の蒸気に面する平面までの高さ
#            :heat_top → h_vap (液面 → 平面)、 :heat_bottom → 0 (平面は液に浸かっている)
#   z_evap : primary wick の蒸発面 (グルーブ側) までの高さ
#            :heat_top →  h_vap + T_swick + H + height、 :heat_bottom → −(h_liq + T_swick + H + height) (重力が助ける)
#   (T_swick: secondary wick 平面部の厚さ、H: primary wick の高さ、height: グルーブ (フィン) の高さ)
# 【必要な圧力差】
#   primary wick の蒸発面のメニスカス:   ΔP_req_p = (P_grv − P_cc) + ρ_l g z_evap
#     (P_grv − P_cc はループ全体の圧力損失。ウィック内の流動損失はモデルに無いので含めない)
#   secondary wick の蒸気側のメニスカス: ΔP_req_s = ρ_l g z_s
# 【毛細管圧】 ΔP_cap = 2σ/r (接触角 0)。r に最大細孔半径を入れると保守側。margin = ΔP_cap − ΔP_req < 0 なら満液の仮定が破れる。
# =============================================================================

"""
    capillary_margin(sol, p::FlatLHPParams; r_pore_p = 5.1e-6, r_pore_s = p.cond.r_c_cc, verbose = true) -> DataFrame

保存点ごとの毛細管圧の余裕 (上の説明)。r_pore_p, r_pore_s は primary / secondary wick の細孔半径 [m]
(既定: primary は実測の 5.1 µm、secondary は CC の核沸騰に使っている r_c_cc)。列:
  t, f_v (CC の気相の体積分率), h_liq, h_vap, z_s, z_evap [m], dP_loop (= P_grv − P_cc),
  dP_req_p, dP_cap_p, margin_p, dP_req_s, dP_cap_s, margin_s [Pa]
"""
function capillary_margin(sol, p::FlatLHPParams; r_pore_p = 5.1e-6, r_pore_s = p.cond.r_c_cc, verbose::Bool = true)
    geo, tables, g = p.geo, p.tables, p.cond.g
    top = p.p_cc.orientation === :heat_top
    stack = geo.cc.T_swick + geo.wick.H + geo.grv.height
    df = DataFrame(t = Float64[], f_v = Float64[], h_liq = Float64[], h_vap = Float64[], z_s = Float64[],
                   z_evap = Float64[], dP_loop = Float64[], dP_req_p = Float64[], dP_cap_p = Float64[],
                   margin_p = Float64[], dP_req_s = Float64[], dP_cap_s = Float64[], margin_s = Float64[])
    for k in eachindex(sol.t)
        c = components_at(sol, p, k)
        P_cc = c.cc.P[1]
        sp = sat_props(Tsat(P_cc, tables), tables)
        lg = cc_layer_geometry(c.ccl.X[1], p.p_cc)
        h_liq, h_vap = lg.geom.h_liq, lg.geom.h_vap
        z_s    = top ? h_vap : 0.0
        z_evap = top ? h_vap + stack : -(h_liq + stack)
        dP_loop  = groove_pressure(c.grv) - P_cc
        dP_req_p = dP_loop + sp.ρ_l * g * z_evap
        dP_req_s = sp.ρ_l * g * z_s
        dP_cap_p = 2 * sp.σ / r_pore_p
        dP_cap_s = 2 * sp.σ / r_pore_s
        push!(df, (t = sol.t[k], f_v = c.ccl.X[1] / p.p_cc.V, h_liq = h_liq, h_vap = h_vap, z_s = z_s,
                   z_evap = z_evap, dP_loop = dP_loop, dP_req_p = dP_req_p, dP_cap_p = dP_cap_p,
                   margin_p = dP_cap_p - dP_req_p, dP_req_s = dP_req_s, dP_cap_s = dP_cap_s,
                   margin_s = dP_cap_s - dP_req_s))
    end
    if verbose
        println("\n毛細管圧の余裕 (orientation = ", p.p_cc.orientation, "、r_pore: primary ", r_pore_p * 1e6,
                " µm、secondary ", r_pore_s * 1e6, " µm):")
        println("  primary wick  : 必要な圧力差の最大 ", round(maximum(df.dP_req_p), digits = 1), " Pa (うち重力 ",
                round(df.dP_req_p[end] - df.dP_loop[end], digits = 1), " Pa、最終時刻)、毛細管圧 ≈ ",
                round(df.dP_cap_p[end], digits = 0), " Pa、余裕の最小 ", round(minimum(df.margin_p), digits = 1), " Pa")
        println("  secondary wick: 必要な圧力差の最大 ", round(maximum(df.dP_req_s), digits = 1), " Pa、毛細管圧 ≈ ",
                round(df.dP_cap_s[end], digits = 0), " Pa、余裕の最小 ", round(minimum(df.margin_s), digits = 1), " Pa")
        any(<(0), df.margin_p) && println("  ⚠ primary wick の余裕が負になる時刻があります (t = ",
                                          df.t[findfirst(<(0), df.margin_p)], " s〜。満液の仮定が破れている)")
        any(<(0), df.margin_s) && println("  ⚠ secondary wick の余裕が負になる時刻があります (t = ",
                                          df.t[findfirst(<(0), df.margin_s)], " s〜。満液の仮定が破れている)")
    end
    return df
end