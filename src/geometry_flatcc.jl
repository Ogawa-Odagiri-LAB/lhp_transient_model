using Interpolations: interpolate, extrapolate, Gridded, Linear, Line

# -----------------------------------------------------------------------------
# 1. f_v → θ_vap の1次元テーブル
# -----------------------------------------------------------------------------

"""
    f_vapor_of_theta(θ_vap)

中心角θ_vap [rad] の気セグメントが円断面に占める面積分率 (厳密式)。
"""
f_vapor_of_theta(θ_vap) = (θ_vap - sin(θ_vap)) / (π)

"""
    build_vapor_theta_table(N=2001)

f_v → θ_vap の補間オブジェクトを構築する。θ_vap を [0,2π] で等間隔に N 点
取り、f_v(θ_vap) を計算して軸を入れ替える。f_v の格子は非等間隔なので
Gridded 補間を使う。N=2001 で f_v∈[1e-3,1-1e-3] における θ_vap の
最大誤差は約6e-6 rad (検算済み)。範囲外は線形外挿 (通常は
cc_liquid_geometry 内の soft_clamp により範囲内に収まる)。
"""
function build_vapor_theta_table(N::Int=2001)
    θ_vap_grid = collect(range(0.0, π, length=N))
    f_grid = f_vapor_of_theta.(θ_vap_grid)
    f_grid[1] = 0.0     # 丸め誤差で厳密に 0, 1 となるよう揃える
    f_grid[end] = 1.0
    return extrapolate(interpolate((f_grid,), θ_vap_grid, Gridded(Linear())), Line())
end

"""
    CC_VAPOR_THETA_TABLE

f_v → θ_vap のテーブル (ジオメトリに依存しない無次元の関係なので、
全CCで共通の1つを使い回す)。include時に1度だけ構築される。
"""
const CC_VAPOR_THETA_TABLE = build_vapor_theta_table()

"""
    theta_of_vapor_fraction(f_v)

気体積分率 f_v から中心角θ_vapを返す (補間のみ。範囲外対策はしない)。
"""
theta_of_vapor_fraction(f_v) = CC_VAPOR_THETA_TABLE(f_v)

# -----------------------------------------------------------------------------
# 2. 気体積分率
# -----------------------------------------------------------------------------

"""
    vapor_volume_fraction(ρm, x, ρ_v)

CC内の気体積分率 f_v = x ρm / ρ_v。
"""
vapor_volume_fraction(ρm, x, ρ_v) = x * ρm / ρ_v

# -----------------------------------------------------------------------------
# 3. 気液分離幾何
# -----------------------------------------------------------------------------

"""
    cc_vapor_geometry(f_v, R, dl; δ=1e-3)

気体積分率 f_v から、水平半円筒CC (内半径R, セル長さdl) の気液分離
幾何量を返す。平らな面を下 (primary wick側)、曲面を上に向けた半円筒で、
蒸気は曲面の最上部に中心角θ_vapの円冠として溜まると仮定する (θ_vap∈[0,π])。戻り値は NamedTuple:

    f_eff … クランプ後の気体積分率
    θ_vap … 蒸気セグメントの中心角 [rad]
    h_liq … 液深さ [m]
    h_vap … 蒸気層の高さ [m] (= R - h_liq)
    A_liq … 液に接する壁面積 [m2] (曲面のみ = R (π-θ_vap) dl)
    A_vap … 蒸気に接する壁面積 [m2] (= R θ_vap dl)

f_v が 0 または 1 に近づくと dθ_vap/df_v が発散する (θ_vap ∝ f_v^(1/3)) ため、
[δ, 1-δ] に滑らかにクランプしてから θ_vap を引く。δ=1e-3 のとき、

壁面積の合計は常に A_liq + A_vap = πR dl (= 半円筒の曲面積) となる。
自然対流の代表長さ L を何にするか (容器直径か、液深さ・蒸気層高さか)
は、このファイルでは決めず、h_liq, h_vap を返すだけに留める
(residuals_cc.jl / heattransfer.jl 側の判断事項)。
"""
function cc_vapor_geometry(f_v, R, dl; δ=1e-3)
    f_eff = soft_clamp(f_v, δ, 1 - δ, δ)
    θ_vap = theta_of_vapor_fraction(f_eff)
    h_liq = R * cos(θ_vap / 2)
    return (
        f_eff = f_eff,
        θ_vap = θ_vap,
        h_liq = h_liq,
        h_vap = R - h_liq,
        A_liq = R * (π - θ_vap) * dl, # 円柱側面の面積のみ（液体が触れる面積は下のprimary wick側も含まれる）
        A_vap = R * θ_vap * dl,
    )
end