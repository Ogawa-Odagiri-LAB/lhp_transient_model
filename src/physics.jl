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