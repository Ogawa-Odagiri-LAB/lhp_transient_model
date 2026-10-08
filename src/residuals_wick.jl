# =============================================================================
# residuals_wick.jl (ウィック)
#
# 設計ドキュメント 式46-49 に対応する。
#
# 【ノードの配置 (半径方向、内側 = CC側 → 外側 = 蒸発面)】
#   ノード 0       : CC側の境界 (ウィック内表面 D_i)。温度 T_0 = T_CC と仮定
#   ノード 1..N    : 体積を持つセル (状態変数 E, T)。中心直径 D_c,i = D_i + (i−1/2)ΔD
#   ノード N+1     : 蒸発界面 (D_o)。体積を持たない薄層。T_{N+1} = T_sat(P_grv)
#   面 k (k=1..N+1): ノード k−1 とノード k の間。質量流量 ṁ が通る。
#
# 【状態変数】 各セルに E_i (微分変数) と T_i (代数変数)
#     E_i = V_i[(1−φ)ρ_wick C_wick T_i + φ ρ_l e_l(T_i)]          (式46)
#     dE_i/dt = ṁ(⟨h⟩_i − ⟨h⟩_{i+1}) + G_i(T_{i+1}−T_i) − G_{i−1}(T_i−T_{i−1})   (式47)
# 蒸発界面は未知数を増やさず、状態から陽に計算する (wick_interface):
#     T_{N+1} = T_sat(P_grv),  Q_ec-wick = h_ev A (T_ec − T_{N+1})
#     ṁ = (Q_ec-wick − G_N (T_{N+1} − T_N)) / λ                    (式48)
# G_i = k 2πL / ln(D_c,i+1/D_c,i)                                   (式49)
#
# 【エンタルピーの風上値 ⟨h⟩】 ドナーセル法の代わりに smooth_step でブレンド
# (他のファイルと同じ方針)。ṁ はウィック内で一様なので、重みは全面で共通。
#
# 【本ファイルの判断事項 (ドキュメントに明記がなく、仮定を置いたもの)】
#   1. 式46の ρ_l は定数 (非圧縮)。ρ_l,i = ρ_l(T_i) と温度依存にすると、
#      CoolPropの基準状態の取り方で見かけの熱容量が変わってしまう
#      (アンモニア300K: ρc_v=1.67e6 に対し d(ρe)/dT=2.14e6 (基準によって
#      最大2.85e6)。ρ_l の温度変化で液の質量が変わる一方、ṁ一定の仮定で
#      質量の出入りは釣り合っているため)。ρ_l を定数 (基準温度 T_ref での
#      飽和液密度) にすれば、e_l の基準は定数項にしか効かず、d(ρe)/dT=ρc_v
#      となって物理的な熱容量に一致する。
#   2. 式48の λ は既定ではドキュメント通り h_fg。界面は体積を持たないので、
#      セルNから入る液のエンタルピー ṁ h_N と、出ていく蒸気のエンタルピー
#      ṁ h_v の差は、厳密には h_v − h_N であり、h_fg との差 ṁ c_p (T_sat−T_N)
#      だけエネルギー収支がずれる (T_N が界面より10K低いと、Q の約4%)。
#      λ_includes_sensible=true で λ = h_v − h_N (エネルギー保存型) になる。
#      どちらの場合も、グルーブへのエネルギーソース S_E は
#      「界面を通過するエネルギー」(wick_groove_source) として定義するので、
#      系全体のエネルギーは保存する。
#   3. 内側境界 (式47 i=1 の G_0, T_0 が未定義): 内表面 (D_i) からセル1中心
#      への伝導 G_0 = k 2πL / ln(D_c,1/D_i) とし、T_0 = T_CC (CC内の液温度)
#      と仮定している。この熱はCCの熱源になる (Q_leak_cc)。
#      inner_conduction=false で断熱にできる。
#   4. 液物性の圧力 P_w は、引数で与える (CC圧力を想定)。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 幾何
# -----------------------------------------------------------------------------

"""
    k_eff_maxwell(k_s, k_l, Φ) -> k_eff [W/(m K)]

ウィックの固相の実効熱伝導率 (ドキュメント 式53、マクスウェルの複合則)。
ウィックのステンレス部の熱伝導率 k_s、空隙を満たす作動流体の液の熱伝導率 k_l、
空隙率 Φ から、

    k_eff = k_s (β − Φ)/(β + Φ),     β = (1 + k_l/k_s)/(1 − k_l/k_s)

Φ→0 で k_s、Φ→1 で k_l になる (極限を数値で確認済み)。
k_l = k_s では β が発散するが、k_eff = k_s (Φに依らない) が極限なので、
その場合は k_s を返す。
"""
function k_eff_maxwell(k_s, k_l, Φ)
    r = k_l / k_s
    abs(1 - r) < 1.0e-12 && return k_s
    β = (1 + r) / (1 - r)
    return k_s * (β - Φ) / (β + Φ)
end

"""
    wick_params(geo::WickGeometry, k_s, k_l) -> NamedTuple

ウィックの幾何量を作る。k_s, k_l は、ウィックのステンレス部と作動流体の液の
熱伝導率 [W/(m K)] で、k_eff_maxwell (式53) で固相の実効熱伝導率 k_eff を
求め、熱コンダクタンス G (式52) に使う。

k_s, k_l は、基準温度 (condition.T_init) での値で、k_eff は定数として扱う
(式52, 53 は k_eff に温度の添字を付けておらず、式49 の ρ_l を定数にしたのと
同じ扱い)。k_s(T) は physics.jl の k_sus304、k_l(T) は飽和テーブルの
tables.k_l_sat を使う (residuals_main.jl の build_params)。

    Dn[k+1]  … ノード k (k=0..N+1) の直径 [m]。Dn[1]=D_i, Dn[2..N+1]=セル中心, Dn[N+2]=D_o
    V[i]     … セル i (i=1..N) の体積 (環状の殻) [m³]
    G[k+1]   … ノード k と k+1 の間の熱コンダクタンス (k=0..N) [W/K] (式52)
               G[1] は内表面→セル1中心 (G_0)、G[N+1] はセルN中心→界面 (G_N)
    k_eff    … 固相の実効熱伝導率 [W/(m K)] (式53)
"""
function wick_params(geo::WickGeometry, k_s, k_l)
    N, L, dD = geo.N, geo.L, geo.dD
    k_eff = k_eff_maxwell(k_s, k_l, geo.φ)
    Dn = vcat(geo.D_i, [geo.D_i + (i - 0.5) * dD for i in 1:N], geo.D_o)
    V = [π / 4 * ((geo.D_i + i * dD)^2 - (geo.D_i + (i - 1) * dD)^2) * L for i in 1:N]
    G = [k_eff * 2π * L / log(Dn[k+2] / Dn[k+1]) for k in 0:N]
    return (N = N, L = L, Dn = Dn, V = V, G = G, k_eff = k_eff,
            por = geo.φ, ρ_w = geo.ρ_wick, C_w = geo.Cp_wick)
end

# -----------------------------------------------------------------------------
# 2. 蒸発界面 (陽に計算。未知数は増やさない)
# -----------------------------------------------------------------------------

"""
    wick_interface(T_N, P_grv, T_ec, h_ev, A_ec_wick, G_N, P_w, tables;
                   λ_includes_sensible=false, T_cc=T_N, ΔT_startup=0.0, ε_startup=0.1)

蒸発界面 (体積を持たない薄層) の量を、状態から陽に計算する。

    T_out  = T_sat(P_grv)                      … 界面温度 (式44の T_wick,out)
    Q_ec   = h_ev A_ec-wick (T_ec − T_out)     … 蒸発器ケースから界面へ入る熱 (式44)
    Q_N    = G_N (T_out − T_N)                 … セルNが界面から受け取る熱 (式47)
    ṁ      = g_s (Q_ec − Q_N)/λ                … 蒸発量 (式48 に起動判定のゲートを掛けたもの)
    Q_extra = Q_ec − Q_N − ṁ λ = (1−g_s)(Q_ec − Q_N)

起動判定のゲート (旧コード fluid_wick の |T_ev − T_cc| > ΔT_startup に相当):
    g_s = smooth_step(|T_N − T_cc| − ΔT_startup, ε_startup)   (ΔT_startup>0 のとき)
    g_s = 1                                                      (ΔT_startup≤0: 判定なし)
旧コードは、判定を満たさないとき蒸発量を0にするだけで、界面に入る熱の行き先を
考えていなかった。本実装の界面は体積を持たない (T_out が T_sat(P_grv) に固定)
ので、そのままではエネルギーが消えてしまう。そこで、蒸発に使われなかった熱
Q_extra を、ウィックの最外セルの熱源として渡す (residual_wick! が受け取る)。
蒸発が始まる前は、最外セルが熱を受けて温まり、|T_N − T_cc| が ΔT_startup を
超えるとゲートが開く。g_s=1 では Q_extra=0 で、従来と同じ。

戻り値は NamedTuple (T_out, Q_ec, Q_N, ṁ, λ, h_N, h_out, Q_extra, g_s)。
h_N, h_out は、セルNと界面の液のエンタルピー (式47の出口側の面で使う)。
"""
function wick_interface(T_N, P_grv, T_ec, h_ev, A_ec_wick, G_N, P_w, tables::FluidTables;
    λ_includes_sensible=false, T_cc=T_N, ΔT_startup=0.0, ε_startup=0.1)

    T_out = Tsat(P_grv, tables)
    sp = sat_props(T_out, tables)
    Q_ec = h_ev * A_ec_wick * (T_ec - T_out)
    Q_N = G_N * (T_out - T_N)

    h_N = liquid_props(T_N, P_w, tables).h
    h_out = liquid_props(T_out, P_w, tables).h
    λ = λ_includes_sensible ? sp.h_v - h_N : sp.h_fg

    # 起動判定のゲート (パラメータ ΔT_startup による分岐で、状態依存の if ではない)
    g_s = ΔT_startup > 0 ? smooth_step(soft_abs(T_N - T_cc, 1.0e-3) - ΔT_startup, ε_startup) : 1.0

    ṁ = g_s * (Q_ec - Q_N) / λ
    Q_extra = Q_ec - Q_N - ṁ * λ
    return (T_out = T_out, Q_ec = Q_ec, Q_N = Q_N, ṁ = ṁ, λ = λ, h_N = h_N, h_out = h_out,
            Q_extra = Q_extra, g_s = g_s)
end

"""
    wick_groove_source(iface; ε_m=1e-8) -> (S_m_total, S_E_total)

蒸発界面からグルーブへ入る質量 [kg/s] とエネルギー [W]。
質量は ṁ。エネルギーは「界面を通過するもの」として
    S_E = ṁ ⟨h⟩_{N+1} + ṁ λ
と定義する (界面は体積を持たないので、入ったものはすべて出る)。
⟨h⟩_{N+1} は面 N+1 (セルNと界面の間) の風上エンタルピーで、
residual_wick! が使う値と同じブレンド (ṁ>0 でセルNの h_N、ṁ<0 で
界面の h_out) を使う。

起動判定のゲートがない (g_s=1) とき、ṁ λ = Q_ec − Q_N なので、
S_E = ṁ⟨h⟩ + Q_ec − Q_N と同じ。ゲートがある (g_s<1) とき、蒸発に使われなかった
熱 Q_extra は、ウィックの最外セルに渡される (residual_wick! の Q_extra) ので、
グルーブへ渡すのは蒸発に使われた分 ṁ λ だけ。これにより、λ の取り方や
ゲートの値、流れの向きによらず、蒸発器ケース・ウィック・グルーブ・CC の
合計のエネルギーが保存する (丸め誤差の範囲で閉じることをPythonで検算済み)。
グルーブの体積あたりのソースにするには、グルーブの体積で割る。
"""
function wick_groove_source(iface; ε_m=1e-8)
    wm = smooth_step(iface.ṁ, ε_m)
    hf = wm * iface.h_N + (1 - wm) * iface.h_out
    return (iface.ṁ, iface.ṁ * hf + iface.ṁ * iface.λ)
end

# -----------------------------------------------------------------------------
# 3. ウィックの残差
# -----------------------------------------------------------------------------

"""
    residual_wick!(resid, du, w, wp, iface, T_cc, P_w, tables; ρ_l_ref, ε_m=1e-8, inner_conduction=true)

ウィックの N セルの残差 (式46, 47) を resid に書き込む。
resid, du, w は (E, T) を持つ ComponentVector。

  iface  : wick_interface の戻り値 (ṁ, 界面温度 T_out, 界面の液のエンタルピーなど)
  T_cc   : CC内の液温度 (ノード0の温度 T_0、およびCCから流入する液の温度)
  P_w    : ウィック内の液物性に使う圧力 (CC圧力を想定)
  ρ_l_ref: 式46で使う液の密度 [kg/m³] (定数。判断事項1)
  ε_m    : ṁ の風上ブレンドの幅 [kg/s]。ṁ の典型値 (1e-5〜1e-4 kg/s) より十分小さく取る

戻り値は (Q_leak_cc) : ウィックからCCの液へ内表面を通って逃げる熱 [W]
(inner_conduction=true のとき G_0 (T_1 − T_CC)、CC側の熱源として使う)。

エネルギーの微分方程式とエネルギーの定義を利用した温度Tに関する代数方程式の残差を計算している
"""
function residual_wick!(resid, du, w, wp, iface, T_cc, P_w, tables::FluidTables;
    ρ_l_ref, ε_m=1e-8, inner_conduction=true)

    N = wp.N
    ṁ = iface.ṁ
    wm = smooth_step(ṁ, ε_m)          # ṁ>0 (外向き) で 1。全面で共通

    # ノード 0..N+1 の温度とエンタルピー (配列の添字は k+1)
    Tn = Vector{Float64}(undef, N + 2)
    hn = Vector{Float64}(undef, N + 2)
    Tn[1] = T_cc
    hn[1] = liquid_props(T_cc, P_w, tables).h
    for i in 1:N
        Tn[i+1] = w.T[i]
        hn[i+1] = liquid_props(w.T[i], P_w, tables).h
    end
    Tn[N+2] = iface.T_out
    hn[N+2] = iface.h_out

    # 面 k (k=1..N+1) の風上エンタルピー ⟨h⟩_k = wm h_{k-1} + (1-wm) h_k
    hf = [wm * hn[k] + (1 - wm) * hn[k+1] for k in 1:N+1]

    G0 = inner_conduction ? wp.G[1] : 0.0

    for i in 1:N
        Gm = i == 1 ? G0 : wp.G[i]          # G_{i-1} (ノード i-1 と i の間)
        Gp = wp.G[i+1]                       # G_i     (ノード i と i+1 の間)
        Q_cond = Gp * (Tn[i+2] - Tn[i+1]) - Gm * (Tn[i+1] - Tn[i])
        # 蒸発に使われなかった熱 (起動判定のゲート) は最外セルが受け取る
        Q_ex = i == N ? iface.Q_extra : 0.0
        resid.E[i] = du.E[i] - (ṁ * (hf[i] - hf[i+1]) + Q_cond + Q_ex)

        # 式46: E_i = V_i[(1-φ)ρ_w C_w T_i + φ ρ_l e_l(T_i)]
        e_l = liquid_props(w.T[i], P_w, tables).e
        resid.T[i] = w.E[i] - wp.V[i] * ((1 - wp.por) * wp.ρ_w * wp.C_w * w.T[i] +
                                         wp.por * ρ_l_ref * e_l)
    end

    return inner_conduction ? wp.G[1] * (w.T[1] - T_cc) : 0.0
end