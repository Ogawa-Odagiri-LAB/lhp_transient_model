# ==========================================
# 1. 補間オブジェクトを保持するグローバル変数
# ==========================================
struct PropyleneLookups
    # 1変数関数（飽和特性：T依存）
    Psat_T::AbstractInterpolation
    latent_heat_T::AbstractInterpolation
    D_sat_g_T::AbstractInterpolation
    D_sat_f_T::AbstractInterpolation
    C_sat_g_T::AbstractInterpolation
    C_sat_f_T::AbstractInterpolation
    V_sat_g_T::AbstractInterpolation
    V_sat_f_T::AbstractInterpolation
    L_sat_g_T::AbstractInterpolation
    L_sat_f_T::AbstractInterpolation
    U_sat_g_T::AbstractInterpolation
    U_sat_f_T::AbstractInterpolation
    H_sat_g_T::AbstractInterpolation
    H_sat_f_T::AbstractInterpolation
    dTdP_sat_T::AbstractInterpolation
    
    # 1変数関数（飽和特性：P依存）
    Tsat_P::AbstractInterpolation

    # 2変数関数（単相領域：T, P依存）
    D_1p_TP::AbstractInterpolation
    C_1p_TP::AbstractInterpolation
    V_1p_TP::AbstractInterpolation
    L_1p_TP::AbstractInterpolation
    U_1p_TP::AbstractInterpolation
    H_1p_TP::AbstractInterpolation
end

const PROPYLENE_DATA = Ref{PropyleneLookups}()

# ==========================================
# 2. テーブルの生成と初期化関数
# ==========================================
function init_propylene_lookups(;
    T_min = 230.0, T_max = 360.0, N_T = 200,  # 臨界点(364.88K)の手前まで
    P_min = 1e5,   P_max = 4.5e6, N_P = 200   # 圧力範囲（Pa）
)
    println("Start initializing props")

    FLUID = "Propylene"
    
    T_range = range(T_min, T_max, length=N_T)
    P_range = range(P_min, P_max, length=N_P)

    # 1次元データの生成 (T依存)
    Psat_vals = [PropsSI("P", "T", t, "Q", 1.0, FLUID) for t in T_range]
    latent_vals = [PropsSI("H", "T", t, "Q", 1.0, FLUID) - PropsSI("H", "T", t, "Q", 0.0, FLUID) for t in T_range]
    Dg_vals = [PropsSI("D", "T", t, "Q", 1.0, FLUID) for t in T_range]
    Df_vals = [PropsSI("D", "T", t, "Q", 0.0, FLUID) for t in T_range]
    Cg_vals = [PropsSI("C", "T", t, "Q", 1.0, FLUID) for t in T_range]
    Cf_vals = [PropsSI("C", "T", t, "Q", 0.0, FLUID) for t in T_range]
    Vg_vals = [PropsSI("V", "T", t, "Q", 1.0, FLUID) for t in T_range]
    Vf_vals = [PropsSI("V", "T", t, "Q", 0.0, FLUID) for t in T_range]
    Lg_vals = [PropsSI("L", "T", t, "Q", 1.0, FLUID) for t in T_range]
    Lf_vals = [PropsSI("L", "T", t, "Q", 0.0, FLUID) for t in T_range]
    Ug_vals = [PropsSI("U", "T", t, "Q", 1.0, FLUID) for t in T_range]
    Uf_vals = [PropsSI("U", "T", t, "Q", 0.0, FLUID) for t in T_range]
    Hg_vals = [PropsSI("H", "T", t, "Q", 1.0, FLUID) for t in T_range]
    Hf_vals = [PropsSI("H", "T", t, "Q", 0.0, FLUID) for t in T_range]
    
    dTdP_vals = map(zip(T_range, Hg_vals, Hf_vals, Dg_vals, Df_vals)) do (t, hg, hf, rg, rf)
        t * (1 / rg - 1 / rf) / (hg - hf)
    end

    # 1次元データの生成 (P依存)
    Tsat_vals = [PropsSI("T", "P", p, "Q", 1.0, FLUID) for p in P_range]

    # 2次元データの生成 (T, P依存)
    D_2d = [PropsSI("D", "T", t, "P", p, FLUID) for t in T_range, p in P_range]
    C_2d = [PropsSI("C", "T", t, "P", p, FLUID) for t in T_range, p in P_range]
    V_2d = [PropsSI("V", "T", t, "P", p, FLUID) for t in T_range, p in P_range]
    L_2d = [PropsSI("L", "T", t, "P", p, FLUID) for t in T_range, p in P_range]
    U_2d = [PropsSI("U", "T", t, "P", p, FLUID) for t in T_range, p in P_range]
    H_2d = [PropsSI("H", "T", t, "P", p, FLUID) for t in T_range, p in P_range]

    # 線形補間オブジェクトのラッパー（グリッド外アクセスは外挿 Line() で許容）
    make_interp1(grid, vals) = extrapolate(scale(interpolate(vals, BSpline(Linear())), grid), Line())
    make_interp2(g1, g2, vals) = extrapolate(scale(interpolate(vals, (BSpline(Linear()), BSpline(Linear()))), g1, g2), Line())

    PROPYLENE_DATA[] = PropyleneLookups(
        make_interp1(T_range, Psat_vals), make_interp1(T_range, latent_vals),
        make_interp1(T_range, Dg_vals), make_interp1(T_range, Df_vals),
        make_interp1(T_range, Cg_vals), make_interp1(T_range, Cf_vals),
        make_interp1(T_range, Vg_vals), make_interp1(T_range, Vf_vals),
        make_interp1(T_range, Lg_vals), make_interp1(T_range, Lf_vals),
        make_interp1(T_range, Ug_vals), make_interp1(T_range, Uf_vals),
        make_interp1(T_range, Hg_vals), make_interp1(T_range, Hf_vals),
        make_interp1(T_range, dTdP_vals), make_interp1(P_range, Tsat_vals),
        make_interp2(T_range, P_range, D_2d), make_interp2(T_range, P_range, C_2d),
        make_interp2(T_range, P_range, V_2d), make_interp2(T_range, P_range, L_2d),
        make_interp2(T_range, P_range, U_2d), make_interp2(T_range, P_range, H_2d)
    )
    println("Propylene lookup tables initialized successfully.")
end