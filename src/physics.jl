# 物性用の関数
const R = 8.31 #[m2 kg s-2 K-1 mol-1]
const M_PROPYLENE = 0.04208 #[kg mol-1]
R_PROPYLENE = R / M_PROPYLENE

# ステンレス
struct FittingParamsLog
    coeffs::SVector{9, Float64}
end

struct FittingParams{N}
    coeffs::SVector{N, Float64}
end

function evapoly(x, coeffs)
    result = coeffs[end] # 最も高い次数の係数（i.params）からスタート
    for i in (length(coeffs)-1):-1:1
        result = result * x + coeffs[i]
    end
    return result
end

function curve_fit_log(params::FittingParamsLog, T)
    return exp10(evapoly(log10(T), params.coeffs))
end

function curve_fit(params, T)
    return evapoly(T, params.coeffs)
end

# 2.2 SUS304
# 2.2.1 熱伝導率
const PARAMS_K_SUS304 = FittingParamsLog(SA[
    -1.4087,     # a0
     1.3982,     # a1
     0.2543,     # a2
    -0.6260,     # a3
     0.2334,     # a4
     0.4256,     # a5
    -0.4658,     # a6
     0.1650,     # a7
    -0.0199      # a8
])
function k_sus304(T)
    return curve_fit_log(PARAMS_K_SUS304, T)
end

# 2.2.2 比熱 (J/kg K)
const PARAMS_c_SUS304 = FittingParamsLog(SA[
      22.0061,   # a0
    -127.5528,   # a1
     303.6470,   # a2
    -381.0098,   # a3
     274.0328,   # a4
    -112.9212,   # a5
      24.7593,   # a6
      -2.239153, # a7
       0.0       # a8
])


function c_sus304(T)
    return curve_fit_log(PARAMS_c_SUS304, T)
end

# 2.2.3 密度
const PARAMS_LE_SUS304 = FittingParams(SA[
    -2.9554e2, -3.9811e-1, 9.2683e-3, -2.0261e-5, 1.7127e-8
])

function ρ_sus304(T)
    y = curve_fit(PARAMS_LE_SUS304, T)
    thermal_strain = y * 1e-5
    return 7930.0 / (1.0 + thermal_strain)^3
end

# 流体
# ==========================================
# 3. ユーザー指定の関数群（完全互換）
# ==========================================

function calc_Psat(T, FLUID)
    return FLUID_TABLES[FLUID].Psat_T(T)
end

function calc_Tsat(P, FLUID)
    return FLUID_TABLES[FLUID].Tsat_P(P)
end

function calc_Pcrit(FLUID)
    return FLUID_TABLES[FLUID].Pcrit
end

function calc_latent_heat(T, FLUID)
    return FLUID_TABLES[FLUID].latent_heat_T(T)
end

function void_correction(α)
    if α <= 0.0
        return 0.0 + 1e-10
    elseif α >= 1.0
        return 1.0 - 1e-10
    else
        return α
    end
end

function vapor_quality_correction(x)
    if x <= 0.0
        return 0.0 + 1e-10
    elseif x >= 1.0
        return 1.0 - 1e-10
    else
        return x
    end
end

function calc_ρ_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return FLUID_TABLES[FLUID].D_sat_g_T(T)
        else
            return FLUID_TABLES[FLUID].D_sat_f_T(T)
        end
    else
        return FLUID_TABLES[FLUID].D_1p_TP(T, P)
    end
end

function calc_ρ_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * FLUID_TABLES[FLUID].D_sat_g_T(T) + (1 - α_p) * FLUID_TABLES[FLUID].D_sat_f_T(T)
end 

function calc_cp_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return FLUID_TABLES[FLUID].C_sat_g_T(T)
        else
            return FLUID_TABLES[FLUID].C_sat_f_T(T)
        end
    else
        return FLUID_TABLES[FLUID].C_1p_TP(T, P)
    end
end

function calc_cp_2p(T, α, FLUID)
    α_p = void_correction(α)
    x = calc_vapor_quality(T, α_p, FLUID)
    return x * FLUID_TABLES[FLUID].C_sat_g_T(T) + (1 - x) * FLUID_TABLES[FLUID].C_sat_f_T(T)
end 

function calc_μ_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return FLUID_TABLES[FLUID].V_sat_g_T(T)
        else
            return FLUID_TABLES[FLUID].V_sat_f_T(T)
        end
    else
        return FLUID_TABLES[FLUID].V_1p_TP(T, P)
    end
end

function calc_μ_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * FLUID_TABLES[FLUID].V_sat_g_T(T) + (1 - α_p) * FLUID_TABLES[FLUID].V_sat_f_T(T)
end

function calc_k_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return FLUID_TABLES[FLUID].L_sat_g_T(T)
        else
            return FLUID_TABLES[FLUID].L_sat_f_T(T)
        end
    else
        return FLUID_TABLES[FLUID].L_1p_TP(T, P)
    end
end

function calc_k_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * FLUID_TABLES[FLUID].L_sat_g_T(T) + (1 - α_p) * FLUID_TABLES[FLUID].L_sat_f_T(T)
end

function calc_e_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return FLUID_TABLES[FLUID].U_sat_g_T(T)
        else
            return FLUID_TABLES[FLUID].U_sat_f_T(T)
        end
    else
        return FLUID_TABLES[FLUID].U_1p_TP(T, P)
    end
end

function calc_e_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * FLUID_TABLES[FLUID].U_sat_g_T(T) + (1 - α_p) * FLUID_TABLES[FLUID].U_sat_f_T(T)
end

function calc_h_1p(T, P, FLUID)
    # Psat = calc_Psat(T, FLUID)
    # if abs(P - Psat) / P < 1e-5
    #     if P < Psat
    #         return FLUID_TABLES[FLUID].H_sat_g_T(T)
    #     else
    #         return FLUID_TABLES[FLUID].H_sat_f_T(T)
    #     end
    # else
        return FLUID_TABLES[FLUID].H_1p_TP(T, P)
    # end
end

function calc_h_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * FLUID_TABLES[FLUID].H_sat_g_T(T) + (1 - α_p) * FLUID_TABLES[FLUID].H_sat_f_T(T)
end

function calc_dT_dP_sat(T, FLUID)
    return FLUID_TABLES[FLUID].dTdP_sat_T(T)
end

function calc_k_eff(T, FLUID)
    Φ = 0.56
    k_s = k_sus304(T)
    k_l = FLUID_TABLES[FLUID].L_sat_f_T(T)
    β = (1 + k_l / k_s) / (1 - k_l / k_s)
    return k_s * (β - Φ) / (β + Φ)
end

# --- 無次元数計算 ---
function calc_Re_1p(T, P, m, D_i, A_cross, FLUID)
    μ = calc_μ_1p(T, P, FLUID)
    return (m * D_i) / (μ * A_cross)
end

function calc_Re_2p(T, α, m, D_i, A_cross, FLUID)
    μ = calc_μ_2p(T, α, FLUID)
    return (m * D_i) / (μ * A_cross)
end

function calc_Pr_1p(T, P, FLUID)
    μ = calc_μ_1p(T, P, FLUID)
    cp = calc_cp_1p(T, P, FLUID)
    k = calc_k_1p(T, P, FLUID)
    return μ * cp / k
end

function calc_Pr_2p(T, α, FLUID)
    μ = calc_μ_2p(T, α, FLUID)
    cp = calc_cp_2p(T, α, FLUID)
    k = calc_k_2p(T, α, FLUID)
    return μ * cp / k
end

function calc_Re_sat(T, m_p, D_i, A_cross, phase::Symbol, FLUID)
    if phase == :liquid
        μ_p = FLUID_TABLES[FLUID].V_sat_f_T(T)
    elseif phase == :vapor
        μ_p = FLUID_TABLES[FLUID].V_sat_g_T(T)
    else
        error("invalid phase")
    end
    return (m_p * D_i) / (μ_p * A_cross)
end

# --- クオリティ・ボイド率 ---
function calc_vapor_quality(T, α, FLUID)
    α_t = void_correction(α)
    denom = (α_t * calc_ρ_2p(T, 1.0, FLUID) + (1 - α_t) * calc_ρ_2p(T, 0.0, FLUID))
    if denom < 1e-15
        return 0.0
    end
    return α_t * calc_ρ_2p(T, 1.0, FLUID) / denom
end

function calc_void_rate(T, x, FLUID)
    x_t = vapor_quality_correction(x)
    ρ_l = calc_ρ_2p(T, 0.0, FLUID)
    ρ_g = calc_ρ_2p(T, 1.0, FLUID)
    denom = x_t * ρ_l + (1 - x_t) * ρ_g
    if denom < 1e-15
        return 0.0
    end
    return x_t * ρ_l / denom
end

# --- 熱伝達率・流量 ---
function calc_h_ec_gr(T, P, m, D_i, A_cross, L, FLUID)
    Re = calc_Re_1p(T, P, m, D_i, A_cross, FLUID)
    Pr = calc_Pr_1p(T, P, FLUID)
    k = calc_k_1p(T, P, FLUID)

    if Re <= 2300
        Nu = 0.662 * (Pr^0.333) * Re^0.5
    else
        Nu = 0.0371 * Pr * Re^0.8
    end
    return Nu * k / L
end

function calc_h_wall_f_1p(T, P, m, D_i, A_cross, FLUID)
    Re = calc_Re_1p(T, P, m, D_i, A_cross, FLUID)
    Pr = calc_Pr_1p(T, P, FLUID)
    k = calc_k_1p(T, P, FLUID)

    if Re <= 2300
        return 3.66 * k / D_i
    else
        return 0.023 * (Re^0.8) * (Pr^0.3) * k / D_i
    end
end

function calc_h_condense(T, P, α, m, D_i, A_cross, FLUID) # Shah(2009)
    Re = calc_Re_2p(T, 0.0, m, D_i, A_cross, FLUID)
    Pr_l = calc_Pr_2p(T, 0.0, FLUID)
    k = calc_k_2p(T, 0.0, FLUID)
    μ_f = calc_μ_2p(T, α, FLUID)
    μ_v = calc_μ_2p(T, 1.0, FLUID)
    x = calc_vapor_quality(T, α, FLUID)
    pr = P/calc_Pcrit(FLUID)
    n = 0.0058 + 0.557 * pr
    h_LT = 0.023 * (Re^0.8) * (Pr_l^0.4) * k / D_i
    h_I = h_LT * (μ_f/(14*μ_v))^n * ((1 - x)^0.8 + (3.8 * (x^0.76) * (1 - x)^0.04) / pr^0.38)

    g = 9.8
    G = m / ((π/4)*D_i^2)
    ρ_l = calc_ρ_2p(T, 0.0, FLUID)
    ρ_v = calc_ρ_2p(T, 1.0, FLUID)
    J_g = x*G / ( g*D_i*ρ_v*(ρ_l - ρ_v) )^0.5
    Z = ((1/x - 1)^0.8) * pr^0.4
    J_bound = 0.98 * (Z + 0.263)^(-0.62)

    Re_LS = Re * (1 - x)
    h_Nu = 1.32 * Re_LS^(-1/3) * (ρ_l * (ρ_l - ρ_v) * g * k^3 / μ_f^2)^(1/3)

    h_TP = h_I + h_Nu
    if J_g > J_bound
        h_TP = h_I
    end

    # @show calc_Psat(T, FLUID)
    # @show (0.023 * (Re^0.8) * (Pr_l^0.4) * k)/D_i
    # @show (1/x-1)^0.8 * pr^0.4
    # @show J_g
    # @show h_Nu
    # return 0.023 * (Re^0.8) * (Pr_l^0.4) * k * ((1 - x)^0.8 + (3.8 * (x^0.76) * (1 - x)^0.04) / pr^0.38) / D_i
    return min(h_TP, 10000.0)
end

# println("h2p = $(calc_h_condense(300.0, 1.2e6, 0.5, 0.00015, 0.00635, π*(0.00435^2)/4, "propylene"))")

function calc_mass_flow_init(T_ev, T_ec, A, h_ev, FLUID)
    λ = calc_latent_heat(T_ev, FLUID)
    if λ < 1e-15
        return 0.0
    end
    return (h_ev * A) * (T_ec - T_ev) / λ
end

# --- 圧損 ---
function calc_τ_wall(T, P, m, D_i, A_cross, FLUID)
    Re = calc_Re_1p(T, P, m, D_i, A_cross, FLUID)
    ρ = calc_ρ_1p(T, P, FLUID)
    v = m / (ρ * A_cross)
    if Re < 2300
        f = 64 / Re
    else
        f = 0.316 / (Re^0.25)
    end
    return f * (ρ * v^2) / 8
end

function calc_ΔP_1p(T, P, m, D_i, A_cross, FLUID)
    τ = calc_τ_wall(T, P, m, D_i, A_cross, FLUID)
    return π * D_i * τ
end

function calc_τ_wall_sat(T, m_p, D_i, A_cross, phase::Symbol, FLUID)
    if phase == :liquid
        Re = calc_Re_sat(T, m_p, D_i, A_cross, :liquid, FLUID)
        ρ = FLUID_TABLES[FLUID].D_sat_f_T(T)
    elseif phase == :vapor
        Re = calc_Re_sat(T, m_p, D_i, A_cross, :vapor, FLUID)
        ρ = FLUID_TABLES[FLUID].D_sat_g_T(T)
    else
        error("invalid phase")
    end

    v = m_p / (ρ * A_cross)
    if Re < 2300
        f = 64 / Re
    else
        f = 0.316 / (Re^0.25)
    end
    return f * (ρ * v^2) / 8
end

function calc_ΔP_2p(T, α, m, D_i, A_cross, FLUID)
    x = calc_vapor_quality(T, α, FLUID)
    x = clamp(x, 0.0, 1.0)
    m_l = (1 - x) * m
    m_v = x * m
    
    Re_l = calc_Re_sat(T, m_l, D_i, A_cross, :liquid, FLUID)
    Re_v = calc_Re_sat(T, m_v, D_i, A_cross, :vapor, FLUID)

    if Re_l > 2300 && Re_v > 2300
        C = 20
    elseif Re_l > 2300 && Re_v < 2300
        C = 12
    elseif Re_l < 2300 && Re_v > 2300
        C = 10
    elseif Re_l < 2300 && Re_v < 2300
        C = 5
    else
        C = 0
        @show Re_l
        @show Re_v
        @show T
        @show m
        @show x
    end

    ΔP_l = π * D_i * calc_τ_wall_sat(T, m_l, D_i, A_cross, :liquid, FLUID)
    ΔP_v = π * D_i * calc_τ_wall_sat(T, m_v, D_i, A_cross, :vapor, FLUID)
    return ΔP_l + ΔP_v + C * sqrt(ΔP_l * ΔP_v)
end

# --- チェック & エネルギー ---
function TPcheck(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P <= 1e-1
        return :tp
    else
        if P < Psat
            return :vapor
        else
            return :liquid
        end
    end
end

function load_T_cc_only()
    # 保存時と同じファイルパス
    filename = "./res/Q_Ts/Q_Ts.csv"
    
    # ファイルが存在するか確認
    if !isfile(filename)
        error("The file not found: ", filename)
    end
    
    # CSVを読み込み、T_cc列をベクトルとして取得
    file = CSV.File(filename)
    T_cc_data = file.T_cc
    
    return T_cc_data
end

function calc_wick_energy(T, V, FLUID)
    Φ = 0.56
    E_wick = V * (1 - Φ) * ρ_sus304(T) * c_sus304(T) * T
    E_liquid = V * Φ * calc_ρ_2p(T, 0.0, FLUID) * calc_e_2p(T, 0.0, FLUID)
    return E_wick + E_liquid
end

function calc_G_pipe(k, D_o, D_i, l)
    return 2 * k * π * (D_o^2 - D_i^2) / 4l
end

function calc_Q_ev(var_old, T_wall_old, geo, condition)
    h_ev = condition.h_ev
    A_ec_wick = geo.grv.A_ec_wick

    T_ec = T_wall_old.evap
    T_ev = var_old.wick[end]

    return h_ev * A_ec_wick * (T_ec - T_ev)
end

function calc_Q_grv(m_evap, var, T_wall, geo, condition)
    FLUID = condition.fluid

    height = geo.grv.height
    A_grv_cross = geo.grv.A_grv_cross
    A_ec_grv = geo.grv.A_ec_grv
    L = geo.grv.L

    T_ec   = T_wall.evap
    T_grv  = var.grv.T
    P_grv  = var.grv.P

    h_ec_grv = calc_h_ec_gr(T_grv, P_grv, m_evap, height, A_grv_cross, L, FLUID)

    return h_ec_grv * A_ec_grv * (T_ec - T_grv)
end

function calc_Q_ec_cc(T_wall_old, geo)
    T_ec = T_wall_old.evap
    T_cc = T_wall_old.cc

    D_ec_o = geo.grv.D_o
    D_ec_i = geo.grv.D_i
    l_ec   = geo.grv.L
    k_ec   = k_sus304(T_ec)
    G_ec   = calc_G_pipe(k_ec, D_ec_o, D_ec_i, l_ec)

    D_cc_o = geo.cc.D_o
    D_cc_i = geo.cc.D_i
    l_cc   = geo.cc.L
    k_cc   = k_sus304(T_cc)
    G_cc   = calc_G_pipe(k_cc, D_cc_o, D_cc_i, l_cc)

    G_ave = 1/(1/G_ec + 1/G_cc)

    return G_ave * (T_ec - T_cc)
end

function calc_Q_ec_vl(var_old, T_wall_old, geo, condition)
    T_ec = T_wall_old.evap
    T_vl = T_wall_old.vl[1]

    D_ec_o = geo.grv.D_o
    D_ec_i = geo.grv.D_i
    l_ec   = geo.grv.L
    k_ec   = k_sus304(T_ec)
    G_ec   = calc_G_pipe(k_ec, D_ec_o, D_ec_i, l_ec)

    D_vl_o = geo.vl.D_o
    D_vl_i = geo.vl.D_i
    l_vl   = geo.vl.dl
    k_vl   = k_sus304(T_vl)
    G_vl   = calc_G_pipe(k_vl, D_vl_o, D_vl_i, l_vl)

    G_ave = 1/(1/G_ec + 1/G_vl)

    return G_ave * (T_ec - T_vl)
end

function calc_Q_wall_cond(G_adj, G, T_adj, T)
    G_ave = 1/(1/G_adj + 1/G)
    return G_ave * (T_adj - T)
end

function calc_Q_wall_amb(dl, k, D_o, D_i, T_amb, T)
    G = (2 * π * k * dl) / log(D_o / D_i)
    return G * (T_amb - T)
end

function calc_Q_conv(h, S, T_f, T_wall)
    return h * S * (T_f - T_wall)
end

function calc_Q_fluid(dl, k, D_o, D_i, T_f, T) # CCでは管内壁表面が内部流体と同じ温度になっていると仮定して計算する
    G = (2 * π * k * dl) / log(D_o / D_i)
    return G * (T_f - T)
end

function calc_Q_cc_ll(T_wall_old, geo)
    T_cc = T_wall_old.cc
    T_ll = T_wall_old.ll[end]

    D_cc_o = geo.cc.D_o
    D_cc_i = geo.cc.D_i
    l_cc   = geo.cc.L
    k_cc   = k_sus304(T_cc)
    G_cc   = calc_G_pipe(k_cc, D_cc_o, D_cc_i, l_cc)

    D_ll_o = geo.ll.D_o
    D_ll_i = geo.ll.D_i
    l_ll   = geo.ll.dl
    k_ll   = k_sus304(T_ll)
    G_ll   = calc_G_pipe(k_ll, D_ll_o, D_ll_i, l_ll)

    G_ave = 1/(1/G_cc + 1/G_ll)

    return G_ave * (T_ll - T_cc)
end

function calc_Q_wick(k_eff, l_w, D_o, D_i, T_adj, T)
    return (2 * π * k_eff * l_w) * (T_adj - T) / log(D_o/D_i)
end

function α_cropper(α)
    return clamp(α, 0.0, 1.0)
end

function h_sigmoid(bal, h_sp, h_tp; k=5.0)
    # シグモイド関数で割合 t を計算
    t = 1.0 / (1.0 + exp(-k * (bal - 1.0)))
    
    # 線形補間
    return (1.0 - t) * h_sp + t * h_tp
end