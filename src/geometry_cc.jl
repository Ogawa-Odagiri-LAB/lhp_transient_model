# =============================================================================
# geometry_cc.jl (CCの気液分離幾何)
#
# 水平円筒のCC (コンペンセーションチャンバー) 内部で、重力により液が
# 下、蒸気が上に分離している (層状) と仮定したときの、液体積分率 f_l
# から液面の幾何量 (中心角θ_liq、液深さ、気液それぞれの壁面接触面積) を
# 求める。
#
# 【設計方針】
# θ_liq は DAE の未知数には含めない (types.jl 参照)。f_l → θ_liq の関係は
# ケプラー方程式型の超越方程式で解析的な逆関数を持たないため、
# オフラインでテーブル化し、実行時は補間のみを行う。
#
# 【幾何関係】 液が占める円弧の中心角を θ_liq [0,2π]、内半径を R とする。
#     液の断面積分率   f_l = (θ_liq - sin θ_liq) / (2π)
#     液深さ           h_l = R (1 - cos(θ_liq/2))
#     液接触の壁面弧長 = R θ_liq、蒸気接触の壁面弧長 = R (2π - θ_liq)
# f_l(θ_liq) は θ_liq の陽な単調増加関数なので、θ_liq を等間隔に取って f_l を
# 計算し、軸を入れ替えるだけで逆関数 θ_liq(f_l) のテーブルになる
# (二分法等の求根は不要)。θ_liq 等間隔の格子は、θ_liq ∝ f_l^(1/3) と立ち上がる
# 端点 (f_l→0,1) 付近で f_l の格子点が自然に密になる。
#
# 【依存】 smoothing.jl (soft_clamp) の後に include すること。
#         Interpolations (Gridded, Linear, extrapolate, Line) を使う。
# =============================================================================

using Interpolations: interpolate, extrapolate, Gridded, Linear, Line

# -----------------------------------------------------------------------------
# 1. f_l → θ_liq の1次元テーブル
# -----------------------------------------------------------------------------

"""
    f_liquid_of_theta(θ_liq)

中心角θ_liq [rad] の液セグメントが円断面に占める面積分率 (厳密式)。
"""
f_liquid_of_theta(θ_liq) = (θ_liq - sin(θ_liq)) / (2π)

"""
    build_theta_table(N=2001)

f_l → θ_liq の補間オブジェクトを構築する。θ_liq を [0,2π] で等間隔に N 点
取り、f_l(θ_liq) を計算して軸を入れ替える。f_l の格子は非等間隔なので
Gridded 補間を使う。N=2001 で f_l∈[1e-3,1-1e-3] における θ_liq の
最大誤差は約6e-6 rad (検算済み)。範囲外は線形外挿 (通常は
cc_liquid_geometry 内の soft_clamp により範囲内に収まる)。
"""
function build_theta_table(N::Int=2001)
    θ_liq_grid = collect(range(0.0, 2π, length=N))
    f_grid = f_liquid_of_theta.(θ_liq_grid)
    f_grid[1] = 0.0     # 丸め誤差で厳密に 0, 1 となるよう揃える
    f_grid[end] = 1.0
    return extrapolate(interpolate((f_grid,), θ_liq_grid, Gridded(Linear())), Line())
end

"""
    CC_THETA_TABLE

f_l → θ_liq のテーブル (ジオメトリに依存しない無次元の関係なので、
全CCで共通の1つを使い回す)。include時に1度だけ構築される。
"""
const CC_THETA_TABLE = build_theta_table()

"""
    theta_of_liquid_fraction(f_l)

液体積分率 f_l から中心角θ_liqを返す (補間のみ。範囲外対策はしない)。
"""
theta_of_liquid_fraction(f_l) = CC_THETA_TABLE(f_l)

# -----------------------------------------------------------------------------
# 2. 液体積分率
# -----------------------------------------------------------------------------

"""
    liquid_volume_fraction(ρm, x, ρ_l)

CC内の液体積分率 f_l = (1-x) ρm / ρ_l。

質量保存と式8 (ρm = 1/(x/ρ_g+(1-x)/ρ_l)) から、f_l = 1 - α (αは
ボイド率) と等価。CCの全質量 m=ρm V を代入すると体積Vが消えるため、
状態変数 (ρm, x) と飽和液密度 ρ_l だけで書ける。
ρ_l は sat_props(T, tables).ρ_l から渡す (このファイルは
FluidTables に依存しない)。

単相液では x=0 かつ ρm≳ρ_l (過冷液) で f_l が1をわずかに超え、
単相蒸気では x=1 で f_l=0 となる。いずれも
cc_liquid_geometry 内で範囲内にクランプされる。
"""
liquid_volume_fraction(ρm, x, ρ_l) = (1 - x) * ρm / ρ_l

# -----------------------------------------------------------------------------
# 3. 気液分離幾何
# -----------------------------------------------------------------------------

"""
    cc_liquid_geometry(f_l, R, dl; δ=1e-3)

液体積分率 f_l から、水平円筒CC (内半径R, セル長さdl) の気液分離
幾何量を返す。戻り値は NamedTuple:

    f_eff … クランプ後の液体積分率
    θ_liq … 液セグメントの中心角 [rad]
    h_liq … 液深さ [m]
    h_vap … 蒸気層の高さ [m] (= 2R - h_liq)
    A_liq … 液に接する壁面積 [m2] (= R θ_liq dl)
    A_vap … 蒸気に接する壁面積 [m2] (= R (2π-θ_liq) dl)

f_l が 0 または 1 に近づくと dθ_liq/df_l が発散する (θ_liq ∝ f_l^(1/3)) ため、
[δ, 1-δ] に滑らかにクランプしてから θ_liq を引く。δ=1e-3 のとき、
実際には液が無い極限でも液接触の壁面弧長がわずかに残る
(f_l=1e-3 で θ_liq≈0.33 rad) が、物理量としては無視できる大きさ。

壁面積の合計は常に A_liq + A_vap = 2πR dl (= 円筒側面積) となる。
自然対流の代表長さ L を何にするか (容器直径か、液深さ・蒸気層高さか)
は、このファイルでは決めず、h_liq, h_vap を返すだけに留める
(residuals_cc.jl / heattransfer.jl 側の判断事項)。
"""
function cc_liquid_geometry(f_l, R, dl; δ=1e-3)
    f_eff = soft_clamp(f_l, δ, 1 - δ, δ)
    θ_liq = theta_of_liquid_fraction(f_eff)
    h_liq = R * (1 - cos(θ_liq / 2))
    return (
        f_eff = f_eff,
        θ_liq = θ_liq,
        h_liq = h_liq,
        h_vap = 2R - h_liq,
        A_liq = R * θ_liq * dl,
        A_vap = R * (2π - θ_liq) * dl,
    )
end