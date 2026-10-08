# =============================================================================
# residuals_wall.jl (蒸気管・凝縮器・液管の管壁温度の残差)
#
# 設計ドキュメント 式7, 20, 30-32, 39-40, 51 に対応する。
#
#   C_i dTw_i/dt = Q_amb,i + Q_fluid,i + Q_cond,i       (式39)
#
#   Q_fluid = h_conv S_inner (T − Tw)        壁が流体から受け取る熱 (式31)
#   Q_amb   = h_nc A_outer (T_amb − Tw)      蒸気管・液管。周囲との自然対流 (式30)
#           = h_con (A_outer/2)(T_cond − Tw) 凝縮器。凝縮器プレートとの接触 (式51)
#   Q_cond  = 隣接セルとの軸方向熱伝導 (式32, 40)
#
# 壁はすべて微分変数で、代数変数を伴わない (types.jl 参照)。
#
# 【対象範囲】 vl, cond, ll の3本の管壁。蒸発器ケース (T_wall.evap, 式41) と
# CCケースの壁 (T_wall.cc) は、周囲の多くのコンポーネントと結合するため、
# それぞれ residuals_evapcase.jl / residuals_cc.jl で扱う。
#
# 【壁の接続 (既存コードの wall_common と同じ)】
#   vl   : 上流 = 蒸発器ケース (T_wall.evap)、下流 = cond の最初のセル
#   cond : 上流 = vl の最終セル、             下流 = ll の最初のセル
#   ll   : 上流 = cond の最終セル、           下流 = CCケースの壁 (T_wall.cc)
# 接続先の壁温度は、ゴースト値 (Tw_up, Tw_dn) として、その接続先の幾何と
# あわせて wall_ghost で渡す。
#
# 【ドキュメントとの照合結果】
#   以前の照合で見つかった次の2点は、ドキュメント更新で修正済み。
#   1. 式32/40 の係数 → k π(r_o²−r_i²)/Δz に修正済み (2π → π)。
#      本実装の kA/Δz (A=π(r_o²−r_i²)) と一致する。なお既存コードの
#      calc_G_pipe は 2kπ(D_o²−D_i²)/(4l) と半セルずつの直列合成を
#      組み合わせて同じ kA/Δz を得ており、本実装もこの結果に揃えている。
#   2. 式39/40 の添字 → Qcond,i = kA/Δz (T_{i+1} − T_i)、式39は
#      「+Qcond,i − Qcond,i−1」に修正済み。セル i の受熱は左右の隣接
#      セルからの伝導の和 kA/Δz[(T_{i−1}−T_i)+(T_{i+1}−T_i)] となり、
#      本実装と一致する。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 壁の幾何・物性
# -----------------------------------------------------------------------------

"""
    wall_params(geo) -> NamedTuple

TLGeometry / CondGeometry から、壁の残差が使う幾何量を取り出す。

    V_w     … 1セルあたりの壁の体積 [m³]
    A_in    … 内壁面積 (= S_inner) [m²]
    A_out   … 外壁面積 [m²]
    A_cross … 軸方向熱伝導に関わる環状断面積 [m²]
"""
wall_params(geo::Union{TLGeometry,CondGeometry}) = (
    N = geo.N, dl = geo.dl, D_i = geo.D_i, D_o = geo.D_o,
    V_w = π * (geo.D_o^2 - geo.D_i^2) / 4 * geo.dl,
    A_in = geo.S_inner,
    A_out = π * geo.D_o * geo.dl,
    A_cross = π * (geo.D_o^2 - geo.D_i^2) / 4,
)

"""
    wall_heat_capacity(T, wp) -> C [J/K]

セルの壁の熱容量 (SUS304の温度依存物性、physics.jl の ρ_sus304, c_sus304)。
"""
wall_heat_capacity(T, wp) = wp.V_w * ρ_sus304(T) * c_sus304(T)

"""
    wall_axial_conductance(T_a, wp_a, T_b, wp_b) -> G [W/K]

隣り合う2つの壁セル a, b の間の軸方向熱コンダクタンス。
各セルの「半セル分」の導線 (k A / (dl/2)) の直列合成で、
同一寸法・同一材料なら素直な k A / dl に一致する (検算済み)。
コンポーネント境界 (寸法や材料温度が異なる) でも同じ式を使える。
"""
function wall_axial_conductance(T_a, wp_a, T_b, wp_b)
    G_a = sus304_half_cell_conductance(T_a, wp_a)
    G_b = sus304_half_cell_conductance(T_b, wp_b)
    return 1 / (1 / G_a + 1 / G_b)
end

sus304_half_cell_conductance(T, wp) = k_sus304(T) * wp.A_cross / (wp.dl / 2)

# -----------------------------------------------------------------------------
# 2. 隣接する壁とのゴースト値
# -----------------------------------------------------------------------------

"""
    wall_ghost(T, wp) -> (T, wp)

接続先の壁セルの温度と幾何をまとめる。蒸発器ケースやCCケースのように
軸方向に分割されていない集中容量の壁も、同じ形で渡せる
(wp は、その壁の A_cross と dl に相当する値を持つ NamedTuple)。
"""
wall_ghost(T, wp) = (T = T, wp = wp)

# -----------------------------------------------------------------------------
# 3. 周囲への熱 Q_amb
# -----------------------------------------------------------------------------

"""
    Q_amb_pipe(Tw, wp, T_amb) -> Q [W]

蒸気管・液管の壁が周囲 (T_amb) と自然対流で交換する熱 (式30)。
壁が受け取る向き (T_amb − Tw) で定義している。
"""
function Q_amb_pipe(Tw, wp, T_amb)
    h = h_ambient_smooth(Tw, T_amb, wp.D_o)
    return h * wp.A_out * (T_amb - Tw)
end

"""
    Q_amb_condenser(Tw, wp, T_cond; h_con=500.0) -> Q [W]

凝縮器の壁が、凝縮器プレート (温度 T_cond) と接触して交換する熱 (式51)。
接触面積は外壁面積の半分。壁が受け取る向き (T_cond − Tw)。
"""
Q_amb_condenser(Tw, wp, T_cond; h_con=500.0) = h_con * (wp.A_out / 2) * (T_cond - Tw)

# -----------------------------------------------------------------------------
# 4. 壁の残差
# -----------------------------------------------------------------------------

"""
    residual_wall_line!(resid_Tw, du_Tw, Tw, T_fluid, h_conv, wp, up, dn, T_env, kind)

蒸気管・凝縮器・液管の壁温度 Tw (N要素) の残差を resid_Tw に書き込む。

    resid = C_i dTw_i/dt − (Q_amb + Q_fluid + Q_cond_left + Q_cond_right)

  T_fluid : その成分の流体温度 (u.vl.T など。状態変数をそのまま使う)
  h_conv  : compute_h_conv_transport! で計算した熱伝達係数。流体側の q と
            共有することで、流体が壁に失う熱と壁が流体から受け取る熱が
            厳密に一致する。
  up, dn  : wall_ghost で作った、上流・下流側の壁セルの温度と幾何
  T_env   : kind=:pipe なら周囲温度 T_amb、kind=:condenser なら T_cond
  kind    : :pipe (蒸気管・液管) または :condenser (凝縮器)
"""
function residual_wall_line!(resid_Tw, du_Tw, Tw, T_fluid, h_conv, wp, up, dn, T_env, kind::Symbol;
    h_con=500.0)

    N = wp.N
    for i in 1:N
        T = Tw[i]

        # --- 流体から受け取る熱 (式31) ---
        Q_fluid = h_conv[i] * wp.A_in * (T_fluid[i] - T)

        # --- 周囲・凝縮器プレートとの熱 (式30, 51) ---
        Q_amb = kind === :condenser ? Q_amb_condenser(T, wp, T_env; h_con=h_con) :
                                      Q_amb_pipe(T, wp, T_env)

        # --- 軸方向熱伝導 (式32, 40): 左右の隣接セルから受け取る ---
        if i == 1
            G_L = wall_axial_conductance(up.T, up.wp, T, wp)
            T_L = up.T
        else
            G_L = wall_axial_conductance(Tw[i-1], wp, T, wp)
            T_L = Tw[i-1]
        end
        if i == N
            G_R = wall_axial_conductance(dn.T, dn.wp, T, wp)
            T_R = dn.T
        else
            G_R = wall_axial_conductance(Tw[i+1], wp, T, wp)
            T_R = Tw[i+1]
        end
        Q_cond = G_L * (T_L - T) + G_R * (T_R - T)

        resid_Tw[i] = wall_heat_capacity(T, wp) * du_Tw[i] - (Q_amb + Q_fluid + Q_cond)
    end
    return nothing
end