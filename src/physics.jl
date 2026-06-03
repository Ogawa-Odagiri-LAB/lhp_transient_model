# 物性用の関数

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
       0.04379   # a8
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
function calc_Psat(T, FLUID)
    return PropsSI("P", "T", T, "Q", 1.0, FLUID)
end

function calc_Tsat(P, FLUID)
    return PropsSI("T", "P", P, "Q", 1.0, FLUID)
end

function calc_latent_heat(T, FLUID)
    return PropsSI("H", "T", T, "Q", 1, FLUID) - PropsSI("H", "T", T, "Q", 0, FLUID)
end

function calc_ρ_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return PropsSI("D", "T", T, "Q", 1.0, FLUID)
        else
            return PropsSI("D", "T", T, "Q", 0.0, FLUID)
        end
    else
        return PropsSI("D", "T", T, "P", P, FLUID)
    end
end


function void_correction(α) # 仕様上，αが負や1より大きな値をとる場合がある．
    if α <= 0.0
        return 0.0 + 1e-8 # 0.0や1.0にすると，計算が破綻する
    elseif α >= 1.0
        return 1.0 - 1e-8
    else
        return α
    end
end

function vapor_quality_correction(x) # 仕様上，xが負や1より大きな値をとる場合がある．
    if x <= 0.0
        return 0.0 + 1e-8
    elseif x >= 1.0
        return 1.0 - 1e-8
    else
        return x
    end
end

function calc_ρ_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * PropsSI("D", "T", T, "Q", 1, FLUID) + (1 - α_p) * PropsSI("D", "T", T, "Q", 0, FLUID)
end 

function calc_cp_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return PropsSI("C", "T", T, "Q", 1.0, FLUID)
        else
            return PropsSI("C", "T", T, "Q", 0.0, FLUID)
        end
    else
        return PropsSI("C", "T", T, "P", P, FLUID)
    end
end

function calc_cp_2p(T, α, FLUID)
    α_p = void_correction(α)
    x = calc_vapor_quality(T, α_p, FLUID)
    return x * PropsSI("C", "T", T, "Q", 1, FLUID) + (1 - x) * PropsSI("C", "T", T, "Q", 0, FLUID)
end 

function calc_μ_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return PropsSI("V", "T", T, "Q", 1.0, FLUID)
        else
            return PropsSI("V", "T", T, "Q", 0.0, FLUID)
        end
    else
        return PropsSI("V", "T", T, "P", P, FLUID)
    end
end

function calc_μ_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * PropsSI("V", "T", T, "Q", 1, FLUID) + (1 - α_p) * PropsSI("V", "T", T, "Q", 0.0, FLUID)
end

function calc_k_1p(T, P, FLUID)
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P < 1e-5
        if P < Psat
            return PropsSI("L", "T", T, "Q", 1.0, FLUID)
        else
            return PropsSI("L", "T", T, "Q", 0.0, FLUID)
        end
    else
        return PropsSI("L", "T", T, "P", P, FLUID)
    end
end

function calc_k_2p(T, α, FLUID)
    α_p = void_correction(α)
    return α_p * PropsSI("L", "T", T, "Q", 1, FLUID) + (1 - α_p) * PropsSI("L", "T", T, "Q", 0, FLUID)
end

function calc_dT_dP_sat(T, FLUID) # クラジウス・クラペイロンの式から計算
    h_g = PropsSI("H", "T", T, "Q", 1, FLUID)
    h_f = PropsSI("H", "T", T, "Q", 0, FLUID)
    ρ_g = PropsSI("D", "T", T, "Q", 1, FLUID)
    ρ_f = PropsSI("D", "T", T, "Q", 0, FLUID)
    
    return T * (1 / ρ_g - 1 / ρ_f) / (h_g - h_f)
end

# println("dTdP=$(calc_dT_dP_sat(283.0))")
function calc_k_eff(T, FLUID)
    Φ = 0.25 # 空隙率
    k_s = k_sus304(T)
    k_l = PropsSI("L", "T", T, "Q", 0.0, FLUID)
    β = (1 + k_l / k_s) / (1 - k_l / k_s)
    return k_s * (β - Φ) / (β + Φ)
end

# 無次元数計算
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

function calc_Re_sat(T, m_p, D_i, A_cross, phase::Symbol, FLUID) # 飽和気体もしくは液体のどちらかを指定してレイノルズ数を計算する関数
    if phase == :liquid
        μ_p = PropsSI("V", "T", T, "Q", 0.0, FLUID)
    elseif phase == :vapor
        μ_p = PropsSI("V", "T", T, "Q", 1.0, FLUID)
    else
        error("invalid phase")
    end

    return (m_p * D_i) / (μ_p * A_cross)
end

# 熱伝達係数の計算
function calc_vapor_quality(T, α, FLUID) # 蒸気クオリティ計算
    α_t = void_correction(α)
    denom = (α_t * calc_ρ_2p(T, 1.0, FLUID) + (1 - α_t) * calc_ρ_2p(T, 0.0, FLUID))
    if denom < 1e-10
        return 0.0
    end
    return α_t * calc_ρ_2p(T, 1.0, FLUID) / denom
end

function calc_void_rate(T, x, FLUID) # ボイド率計算
    x_t =  vapor_quality_correction(x)
    denom = (x_t * calc_ρ_2p(T, 0.0, FLUID) + (1 - x_t) * calc_ρ_2p(T, 1.0, FLUID))
    if denom < 1e-10
        return 0.0
    end
    return x_t * calc_ρ_2p(T, 0.0, FLUID) / denom
end

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

function calc_h_wall_f_2p(T, α, m, D_i, A_cross, FLUID)
    Re = calc_Re_2p(T, α, m, D_i, A_cross, FLUID)
    Pr = calc_Pr_2p(T, α, FLUID)
    k = calc_k_2p(T, α, FLUID)
    x = calc_vapor_quality(T, α, FLUID)

    return 0.023 * (Re^0.8) * (Pr^0.4) * k * ((1 - x)^0.8 + (3.8 * (x^0.76) * (1 - x)^0.04) / Pr^0.38) / D_i
end

# 質量流量計算
function calc_mass_flow_init(T_ev, T_ec, A, h_ev, FLUID) # ウィックでの蒸発量をもとに全体の質量流量を計算 
    λ = calc_latent_heat(T_ev, FLUID) # [J/kg]
    if λ < 1e-10
        return 0.0
    end
    return (h_ev * A) * (T_ec - T_ev) / λ # [kg/s]
end

# 圧損計算
function calc_τ_wall(T, P, m, D_i, A_cross, FLUID)
    Re = calc_Re_1p(T, P, m, D_i, A_cross, FLUID)
    ρ = calc_ρ_1p(T, P, FLUID)
    v = m / (ρ * A_cross)
    if Re < 2300
        f = 64 / Re
    else
        f = 0.316 / (Re^0.25)
    end

    return f * (ρ * v^2) / 8 # 液損のほうが大きいのは，1/A ~ 1/(D^2)，1/Re ~ 1/Dの効果が大きい
end

function calc_τ_wall_sat(T, m_p, D_i, A_cross, phase::Symbol, FLUID) # 飽和気体もしくは液体のどちらかを指定してせん断応力を計算する関数
    if phase == :liquid
        Re = calc_Re_sat(T, m_p, D_i, A_cross, :liquid, FLUID)
        ρ = PropsSI("D", "T", T, "Q", 0.0, FLUID)
    elseif phase == :vapor
        Re = calc_Re_sat(T, m_p, D_i, A_cross, :vapor, FLUID)
        ρ = PropsSI("D", "T", T, "Q", 1.0, FLUID)
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

function calc_ΔP_2p(T, α, m, D_i, A_cross, FLUID) # 二相流の圧損 [Pa/m]
    x = calc_vapor_quality(T, α, FLUID)
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
    end

    ΔP_l = π * D_i * calc_τ_wall_sat(T, m_l, D_i, A_cross, :liquid, FLUID) / A_cross
    ΔP_v = π * D_i * calc_τ_wall_sat(T, m_v, D_i, A_cross, :vapor, FLUID) / A_cross
    return ΔP_l + ΔP_v + C * sqrt(ΔP_l * ΔP_v)
end

# 二相流検知関数
function TPcheck(T, P, FLUID) # 凝縮器での二相流検知
    Psat = calc_Psat(T, FLUID)
    if abs(P - Psat) / P <= 1e-1
        return :tp
    else
        if P < Psat
            return :vapor # vapor
        else
            return :liquid # liquid
        end
    end
end

# Tccをとってくる関数
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