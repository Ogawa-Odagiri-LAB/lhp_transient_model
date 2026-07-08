@with_kw struct LHPSteadyCondition
    fluid::String = "propylene"
    h_ev::Float64 = 8000.0       # [W/m2 K]
    T_cond::Float64 = 15.0 + 273.15 # [K]
    Q_apply::Float64 = 40.0       # [W]
end

@with_kw struct LHPTransientCondition
    Q_apply::Float64
    dt::Float64 # [s]
    t_init::Float64
    t_end::Float64

    fluid::String = "ammonia"
    h_ev::Float64 = 8000.0       # [W/m2 K]
    T_cond::Float64 = 15.0 + 273.15 # [K]
    T_amb::Float64 = 25.0 + 273.15 # [K]
    Q_apply_init::Float64 = 50.0       # [W]
    ΔT_startup::Float64 = 0.0 # [K]
    T_scale::Float64 = 300
    P_scale::Float64 = 1.0e6
    α_scale::Float64 = 1.0
    m_scale::Float64 = 5.0e-4

    fluid_comp::NTuple{6, Symbol} = (:wick, :grv, :vl, :cond, :ll, :cc)
    wall_comp::NTuple{5, Symbol} = (:evap, :vl, :cond, :ll, :cc)
end

@with_kw struct opt_u
    T_ev::Float64
    T_ec::Float64
end

function opt_u_with_Q(Q_apply)
    T_ev = Q_apply - 5.0 + 273.15
    # T_ev = 20.0 + 273.15
    T_ec = Q_apply + 273.15
    # T_ec = 25.0 + 273.15
    return opt_u(T_ev = T_ev, T_ec = T_ec)
end

# ウィック
@with_kw struct WickGeometry
    N::Int = 5

    n_grv::Int = 12
    L::Float64 = 0.035

    D_i::Float64 = 0.003
    D_o::Float64 = 0.006 # グルーブ部を除いた直径
    dD::Float64 = (D_o - D_i) / N

    width::Float64 = 0.001
    height::Float64 = 0.001

end

# グルーブ
@with_kw struct GrvGeometry
    N::Int = 1
    n_grv::Int = 12
    L::Float64 = 0.035
    dl::Float64 = L / N
    D_i::Float64 = 0.008 # ケースの内径
    D_o::Float64 = 0.01 # ケースの外径
    width::Float64 = 0.001
    height::Float64 = 0.001
    A_ec_wick::Float64 = width * L * n_grv
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
    S_inner::Float64 = π * D_i * dl
    V::Float64 = A_cross * dl
end

# 凝縮器
@with_kw struct CondGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2
    S_inner::Float64 = π * D_i * dl
    V::Float64 = A_cross * dl
end

@with_kw struct CCGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2
    S_inner::Float64 = π * D_i * dl
    V::Float64 = A_cross * dl
end

@with_kw struct LHPGeometry
    wick::WickGeometry
    grv::GrvGeometry
    vl::TLGeometry
    cond::CondGeometry
    ll::TLGeometry
    # bay::TLGeometry
    cc::CCGeometry
end

@with_kw struct LHPSteadyGeometry
    wick::WickGeometry
    grv::GrvGeometry
    vl::TLGeometry
    cond::CondGeometry
    ll::TLGeometry
    bay::TLGeometry
    cc::CCGeometry
end
    