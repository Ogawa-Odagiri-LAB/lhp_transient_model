@with_kw struct LHPCondition
    fluid::String = "propylene"
    h_ev::Float64 = 8000.0       # [W/m2 K]
    T_cond::Float64 = 15.0 + 273.15 # [K]
    Q_apply::Float64 = 40.0       # [W]
end

@with_kw struct opt_u
    T_ev::Float64
    T_ec::Float64
end

function opt_u_with_Q(Q_apply)
    T_ev = Q_apply - 10.0 + 273.15
    T_ec = Q_apply + 273.15
    return opt_u(T_ev = T_ev, T_ec = T_ec)
end

# グルーブ
@with_kw struct GrvGeometry
    N::Int = 7
    n_grv::Int = 12
    L::Float64 = 0.035
    dl::Float64 = L / N
    D_i::Float64 = 0.0015
    D_o::Float64 = 0.003
    width::Float64 = 0.001
    height::Float64 = 0.001
    depth::Float64 = 0.035
    A_ec_wick::Float64 = width * depth * n_grv
    A_ec_grv::Float64 = A_ec_wick
    A_grv_cross::Float64 = width * height
end

# 輸送管
@with_kw struct TLGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2
end

# 凝縮器
@with_kw struct CondGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2
end

@with_kw struct CCGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2
end

@with_kw struct LHPGeometry
    grv::GrvGeometry
    vl::TLGeometry
    cond::CondGeometry
    ll::TLGeometry
    bay::TLGeometry
    cc::CCGeometry
end