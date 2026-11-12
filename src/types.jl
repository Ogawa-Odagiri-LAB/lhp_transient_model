# =============================================================================
# types.jl
#
# このファイルは3種類の「型」を定義する。
#
#   1. ジオメトリ構造体      … 既存設計をほぼ踏襲
#   2. 運転条件構造体        … 新しい変数セット (ρm,u,Em,P,T,x,T_wall) 用に
#                              スケーリング定数を更新
#   3. 状態ベクトルの組み立て … ComponentArray として状態を表現するための
#                              コンストラクタ関数、および微分変数/代数変数
#                              マスクのビルダー
#
# 【設計方針】
# 状態 (u, du) はカスタム struct ではなく ComponentArray として表現する。
# 理由:
#   - IDA (Sundials.jl) に渡す際、getdata(u)/getaxes(u) で平坦なVector⇔
#     階層構造を相互変換できる必要があるため (DAEProblemはVector{Float64}
#     しか受け付けない)。
#   - 既存コードもすでに ComponentArrays を利用しており (var.wick, var.grv.T
#     のようなドット記法)、この慣習をそのまま引き継げる。
#   - 「TransportState」のようなカスタム struct を別途定義すると、状態の
#     実体 (ComponentArray) と型定義が二重管理になり、同期ミスの元になる。
#
# 【状態ベクトルの階層構造 (initialize_state が組み立てるもの)】
#
#   u.vl.ρm, u.vl.u, u.vl.Em, u.vl.P, u.vl.T, u.vl.x   … 蒸気管 (N要素)
#   u.cond.*  (同上)                                     … 凝縮器
#   u.ll.*    (同上)                                     … 液管
#   u.grv.*   (同上)                                     … グルーブ
#   u.cc.ρm, u.cc.u, u.cc.Em, u.cc.P, u.cc.T, u.cc.x     … CC
#       (輸送管・グルーブと全く同じ6変数構成。気液分離幾何の角度θ_liqは
#        DAEの未知数に含めず、液体積分率 f_l から θ_liq を求める1次元
#        テーブル (以前議論したケプラー方程式型の関係) をオフラインで
#        構築し、residuals_cc.jl / heattransfer.jl 側で都度補間して
#        使う。これにより、以前懸念した「端点(空・満液)での特異性が
#        Newton反復に直接影響する」リスクを、DAE本体から切り離せる)
#   u.wick.E, u.wick.T                                   … ウィック (N要素)
#       (E: 式46の内部エネルギー, T: 温度。E=微分変数, T=代数変数)
#   u.T_wall.evap                                        … 蒸発器ケース (長さ1の配列)
#       (式41の T_ec そのもの。T_ec = u.T_wall.evap[1]。他の壁と同じ配列に
#        しているのは、ComponentArray のスカラー成分への書き込みや、
#        ./= による一括更新が、配列成分と挙動が違うため)
#   u.T_wall.vl, u.T_wall.cond, u.T_wall.ll              … 各輸送管の壁 (N要素)
#       (流体側と同じ分割・同じ添字に対応する設計)
#   u.T_wall.cc                                          … CCケースの壁 (N要素、通常N=1)
#
#   ※ グルーブ・ウィックには独立した壁構造体を持たせない
#     (グルーブは蒸発器ケース内部にあり、ウィックは式46で固体マトリクスと
#      液を一体化したエネルギーとして扱っているため)。
#
# =============================================================================

using Parameters: @with_kw
using ComponentArrays: ComponentArray

# -----------------------------------------------------------------------------
# 1. ジオメトリ構造体
# -----------------------------------------------------------------------------

# ウィック
@with_kw struct WickGeometry
    N::Int = 5                     # 半径方向の分割数

    n_grv::Int = 12                # グルーブ本数
    L::Float64 = 0.035             # 軸方向長さ [m]

    D_i::Float64 = 0.003           # 内径 [m] (CC側)
    D_o::Float64 = 0.006           # 外径 [m] (グルーブ部を除く)
    dD::Float64 = (D_o - D_i) / N  # 半径方向の分割幅 [m]

    width::Float64 = 0.001         # グルーブ溝幅 [m]
    height::Float64 = 0.001        # グルーブ溝深さ [m]

    # --- 式46 (ウィックエネルギー) で必要な材料物性 ---
    # デフォルト値は与えていない: 誤った物理定数が気づかれないまま
    # 使われることを避けるため、呼び出し側で必ず指定する。
    φ::Float64                     # 空隙率 (0〜1)
    ρ_wick::Float64                # 固体マトリクス密度 [kg/m3]
    Cp_wick::Float64               # 固体マトリクス比熱 [J/(kg・K)]
    # 固相の実効熱伝導率 k_eff は、入力にしない。ウィックのステンレス部の熱伝導率
    # (physics.jl の k_sus304) と、液の熱伝導率 (飽和テーブル) から、空隙率 φ を使って、
    # マクスウェルの複合則 (式53) で計算する (residuals_wick.jl の k_eff_maxwell)。
end

# グルーブ (蒸発器ケース内部、輸送管と同じ6変数構成で扱う)
@with_kw struct GrvGeometry
    N::Int = 1
    n_grv::Int = 12
    L::Float64 = 0.035
    dl::Float64 = L / N
    D_i::Float64 = 0.008           # ケース内径 [m]
    D_o::Float64 = 0.01            # ケース外径 [m]
    width::Float64 = 0.001
    height::Float64 = 0.001
    A_ec_wick::Float64 = width * L * n_grv   # 蒸発器ケース-ウィック界面積 [m2]
    A_ec_grv::Float64 = A_ec_wick            # 蒸発器ケース-グルーブ界面積 [m2]
    A_grv_cross::Float64 = width * height    # グルーブ1本あたりの断面積 [m2]
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

# 輸送管 (蒸気管・液管共用)
@with_kw struct TLGeometry
    N::Int
    L::Float64
    dl::Float64 = L / N
    D_i::Float64
    D_o::Float64
    A_cross::Float64 = (π / 4) * D_i^2       # 流路断面積 [m2]
    S_inner::Float64 = π * D_i * dl          # 1セルあたりの内壁表面積 [m2]
    V::Float64 = A_cross * dl                # 1セルあたりの流体体積 [m3]
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
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
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

# CC (気液分離幾何を扱うため、半径 R も明示的に保持しておく)
@with_kw struct CCGeometry
    N::Int
    L::Float64                                # 軸方向長さ [m] (geometry_cc.jl の Lc)
    dl::Float64 = L / N
    D_i::Float64                              # 内径 [m]
    D_o::Float64                              # 外径 [m]
    R::Float64 = D_i / 2                      # 内半径 [m] (水平円筒の気液分離計算用)
    A_cross::Float64 = (π / 4) * D_i^2
    S_inner::Float64 = π * D_i * dl
    V::Float64 = A_cross * dl
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

@with_kw struct LHPGeometry
    wick::WickGeometry
    grv::GrvGeometry
    vl::TLGeometry
    cond::CondGeometry
    ll::TLGeometry
    cc::CCGeometry
end

# ---- 平板型蒸発器 (金属3Dプリンター製) ----
@with_kw struct FlatWickGeometry
    N::Int = 5                     # 半径方向の分割数

    n_grv::Int = 19                # グルーブ本数
    W::Float64 = 0.030             # primary wick直径方向長さ [m]
    L::Float64 = 0.030             # primary wick長さ [m]
    H::Float64 = 0.0032            # ウィックの高さ [m]（グルーブ部を除く）
    dH::Float64 = H / N            # ウィック高さ方向の分割幅 [m]

    width::Float64 = 0.00075         # グルーブ溝幅 [m]
    height::Float64 = 0.0008        # グルーブ溝深さ [m]

    # --- 式46 (ウィックエネルギー) で必要な材料物性 ---
    # デフォルト値は与えていない: 誤った物理定数が気づかれないまま
    # 使われることを避けるため、呼び出し側で必ず指定する。
    φ::Float64                     # 空隙率 (0〜1)
    ρ_wick::Float64                # 固体マトリクス密度 [kg/m3]
    Cp_wick::Float64               # 固体マトリクス比熱 [J/(kg・K)]
end

@with_kw struct SecondaryWickGeometry
    N::Int = 1

    φ_s::Float64                     # 空隙率 (0〜1)。primary wick の φ とは別
    ρ_swick::Float64                # 固体マトリクス密度 [kg/m3]
    Cp_swick::Float64               # 固体マトリクス比熱 [J/(kg・K)]
end

@with_kw struct FlatGrvGeometry
    N::Int = 1
    n_grv::Int = 19
    W::Float64 = 0.030             # primary wick直径方向長さ [m]
    L::Float64 = 0.030             # primary wick長さ [m]
    L_grv::Float64 = 0.027         # 溝の長さ [m] (設計値 11/2)。ウィック長 L (30 mm) とは別。溝の体積・過熱面積・Re_L・出力の長さはこれを使う
    H::Float64 = 0.004             # primary wick高さ [m]
    dl::Float64 = L_grv / N        # 溝の 1 セルの長さ (溝長 L_grv を N 分割)
    W_case::Float64 = 0.032                       # 蒸発器ケース受熱面の幅 [m]
    L_case::Float64 = 0.032                       # 蒸発器ケース受熱面の長さ [m]
    H_case::Float64 = 0.0055                      # 蒸発器ケース（と考えているところまで）の高さ [m]
    T_bottom::Float64 = 0.0020               # 蒸発器ケース下面の厚さ [m] (設計値 11/2、旧 1.5 mm)。現状は計算に使われていない (C_ec_P は scripts/transient_flat.jl で直接与える)
    T_side::Float64 = 0.001                  # 蒸発器ケース側面の厚み [m]
    width::Float64 = 0.00075
    height::Float64 = 0.0008
    n_fin::Int = 18                # フィン (溝に挟まれた多孔体の壁) の本数。溝 19 本の間
    w_fin::Float64 = 0.00075       # フィンの幅 [m] (設計値 11/2。端の溝は実機 1.0 mm だがモデルは 0.75 mm のまま)
    A_ec_wick::Float64 = n_fin * w_fin * L_grv   # 蒸発器ケース-ウィック界面積 [m2] = フィン 18 本 × 0.75 mm × 27 mm = 364.5 mm2 (旧: 溝幅×30 mm×19 本 = 427.5 mm2)。h_ev A の A
    A_ec_grv::Float64 = width * L_grv * n_grv   # 蒸発器ケース-グルーブ蒸気の過熱面積 [m2] = 溝 19 本 × 0.75 mm × 27 mm = 384.75 mm2 (A_ec_wick とは別の面)
    A_grv_cross::Float64 = width * height    # グルーブ1本あたりの断面積 [m2]
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

"""
    groove_length(geo) -> L_grv [m]

溝の長さ。FlatGrvGeometry はウィック長 L とは別の L_grv、旧 GrvGeometry は L。
"""
groove_length(geo::GrvGeometry) = geo.L
groove_length(geo::FlatGrvGeometry) = geo.L_grv

# 平板型: 半円筒のCC (曲面が上、平面が下で secondary wick を介して primary wick に接する)
@with_kw struct FlatCCGeometry
    N::Int
    L::Float64 = 0.025                                              # 軸方向長さ [m] (geometry_cc.jl の Lc)
    dl::Float64 = L / N
    D_i::Float64 = 0.027                                            # 内径（流体部直径）
    D_swick::Float64 = 0.030                                        # secondary wick外径 [m]
    D_o::Float64 = 0.032                                            # 外径 [m]
    # 2 節点の CC ケース (cc_2node) の伝導断面: 図面 (外形 33 mm 幅) からの見積り。CAD で要確認。
    W_o::Float64 = 0.033                                            # ドームの外形の幅・長さ [m] (表: 蒸発器 W/L 33/33 mm)
    t_shell::Float64 = 0.0015                                       # ドームの曲面の壁厚 [m] = (33 − 30)/2 (外形と secondary wick の外径の差)
    t_end::Float64 = 0.0015                                         # ドームの端面の壁厚 [m] (実測: 肉厚はどこでも 1.5 mm。以前は (33 − 25)/2 = 4 mm としていたが、その内側 2.5 mm は secondary wick (別の経路 Q_cs)。外形は 33×33 mm)
    R::Float64 = D_i / 2                                            # 内半径 [m] (水平円筒の気液分離計算用)
    T_swick::Float64 = 0.0025                                       # secodnary wick 厚み
    A_cross::Float64 = 0.5 * (π / 4) * D_i^2                        # 流体部断面積 [m^2] (半円)
    S_inner::Float64 = π * D_i * 0.5 * dl + (π / 4) * D_i^2         # CC流体部の曲面+両端面の面積 [m^2] (平面側は含まない)
    S_swick::Float64 = π * D_swick * 0.5 * dl + (π / 4) * D_swick^2 # secondary wick表面積（蒸発器ケースとsecondary wickの接触面積）[m^2]
    A_pr_se::Float64 = D_swick^2                                    # primary wickとsecondary wickの接触面積 [m^2]
    V::Float64 = A_cross * dl
    θ::Float64 = 0.0                         # 傾斜角 [度] (水平=0)。運動量式の重力項 g sinθ に使う
    elev0::Float64 = 0.0                     # 入口 (最上流の面) の標高 [m]。Emの位置エネルギー項に使う
end

# 凝縮器プレート (A5052)。L × W の長方形 (L 方向 nx 分割、W 方向 ny 分割)、厚さ t。
# 管は短辺 (x = 0) から入り、x 方向にまっすぐ進み、半円で折り返して (間隔 pitch)、同じ辺から出る。
# コールドプレートも同じ寸法で、プレート全面で接触する (温度は一定)。
@with_kw struct CondPlateGeometry
    L::Float64 = 0.300
    W::Float64 = 0.060
    t::Float64 = 0.005
    nx::Int = 30
    ny::Int = 6
    pitch::Float64 = 0.030
end

# コールドプレート (A5052、凝縮器プレートと同じ平面寸法・同じ格子)。厚さ t、チラーの配管 (外径 D_o、
# 長さ L_tube、N 分割) は、LHP の管と反対側の短辺 (x = L) から入り、折り返して (間隔 pitch は凝縮器プレートと同じ) 戻る。
# 配管の中の冷却液 (エチレングリコール) の温度は T_cond で一定とする。
@with_kw struct ColdPlateGeometry
    t::Float64 = 0.010
    D_o::Float64 = 3.175e-3     # チラー配管 (1/8 in) の外径 [m]
    L_tube::Float64 = 0.56      # チラー配管の長さ [m] (凝縮器と同じ)
    N::Int = 50                 # 配管の分割数 (凝縮器と同じ)
end

# 平板型のLHP全体のジオメトリ。secondary wick (swick) を持つ
@with_kw struct FlatLHPGeometry
    wick::FlatWickGeometry
    swick::SecondaryWickGeometry
    grv::FlatGrvGeometry
    vl::TLGeometry
    cond::CondGeometry
    ll::TLGeometry
    cc::FlatCCGeometry
    plate::CondPlateGeometry = CondPlateGeometry()    # 凝縮器プレート
    cplate::Union{Nothing,ColdPlateGeometry} = nothing    # コールドプレート (nothing なら T_cond 一定の壁に h_cp で接する)
end

# コールドプレートの節点数 (なければ 0)
n_cplate(geo) = hasproperty(geo, :cplate) && !isnothing(geo.cplate) ? geo.plate.nx * geo.plate.ny : 0

# 円筒型・平板型の共通の上位型 (initialcondition.jl などの関数を共用するため)
abstract type AbstractLHPParams end
const AnyLHPGeometry = Union{LHPGeometry,FlatLHPGeometry}

# 平板型の CC の液相の層数と、層の厚さ (体積) の等比 (気液界面の層 1 ほど細かい)。液全体 (深さ h_liq) を N_CC_LAYER 層に、
# 層 1 の割合 g1 = (r−1)/(r^N−1)、以降 r 倍ずつ (residuals_flatcc.jl の CC_LAYER_FRAC)。f_v = 0.2 (h_liq ≈ 9 mm) で、層 1 は約 0.08 mm、層 N は約 1.8 mm。
# 旧版 (集中定数の CC + 界面の熱伝導層 δ_layer_cc) と同じ N = 16、r = 1.23。
# CC_LUMPED = true: 液を 1 層 (集中定数、液の温度は一様) にする。気相は飽和蒸気 1 セルのまま、界面の凝縮は P ≈ P_sat(T_液) に張り付くので、
#   ほぼ「気液平衡の 1 ノードの CC」(旧版) になる。secondary wick との自然対流の熱交換、secondary wick (節点 A/B)、CC 壁 (節点 A/B)、壁面凝縮は
#   そのまま使う (層数に依存しない)。false: 16 層 (界面側ほど細かい)。切り替えたら Julia を再起動して実行する (const のため)。
const CC_LUMPED = false
const N_CC_LAYER = CC_LUMPED ? 1 : 16
const CC_LAYER_RATIO = 1.23

# 平板型: 周囲との熱交換の遅い経路 (residuals_housing.jl)。断熱材 (エアロフレックス) は 3 節点 (HB 側・蒸発器ケース側・CC 側)、
# PEEK ブロックは 1 次元 N_PEEK 節点 (HB 側ほど細かい等比格子)。
const N_HI = 3
const N_PEEK = 12
const PEEK_RATIO = 1.4

"""
    initialize_state(geo::FlatLHPGeometry)

円筒型の状態に、secondary wick (u.swick.E, u.swick.T: 各1要素) を加えたもの。
"""
function initialize_state(geo::FlatLHPGeometry)
    return ComponentArray(
        vl   = _zeros_transport(geo.vl.N),
        cond = _zeros_transport(geo.cond.N),
        ll   = _zeros_transport(geo.ll.N),
        grv  = _zeros_transport(geo.grv.N),
        cc   = _zeros_transport(geo.cc.N),
        wick  = (E = zeros(geo.wick.N), T = zeros(geo.wick.N)),
        swick = (E = zeros(2), T = zeros(2)),    # secondary wick: [1] = 節点 A (平面側 + 両端面)、[2] = 節点 B (半円環)
        hi    = (T = zeros(N_HI),),      # 断熱材 (エアロフレックス) の 3 節点の温度 [K] (スケールなし)
        peek  = (T = zeros(N_PEEK),),    # HB の下の PEEK ブロックの温度 [K] (スケールなし)
        cca   = (T = zeros(1),),         # CC ケースの蒸発器側の節点 A [K] (cc_2node = true のときだけ使う。T_wall.cc は頂部側の節点 B)
        eca   = (T = zeros(1),),         # 蒸発器ケースの側面の節点 S [K] (ec_2node = true のときだけ使う。T_wall.evap は受熱面の節点 P。EC1・EC2 のセンサは S)
        cplate = (T = zeros(n_cplate(geo)),),    # コールドプレート (2次元、凝縮器プレートと同じ格子。なければ 0 要素)
        plate = (T = zeros(geo.plate.nx * geo.plate.ny),),    # 凝縮器プレート (A5052) の2次元の温度 [K] (k = (j−1) nx + i)。コールドプレート (T_cond 一定) に h_cp で接する
        # 平板型の CC (気相 1 セル + 液相 N_CC_LAYER 層、residuals_flatcc.jl):
        #   Tl: 液の各層の温度 [K] (層 1 = 気液界面側、層 N = 底 = 液管の流入・ウィックへの流出)
        #   X : 名目の空き容積 V − V_L,nom [m³] (V_L,nom = 飽和液密度で見た液の体積。気相があるときは ≈ 気相の体積。負なら液が満水で圧縮されている)
        #   mv: 気相の質量 [kg]
        # cc (ρm,u,Em,P,T,x) は、P が圧力の未知数 (代数変数)、他は ccl から決まる要約 (代数変数。診断・保存用。T は気相の飽和温度)
        ccl   = (Tl = zeros(N_CC_LAYER), X = zeros(1), mv = zeros(1)),
        nuc   = (n = zeros(geo.grv.N), nc = zeros(geo.grv.N)),     # n: 界面蒸発のゲート (r_c_wick>0 のとき。0 なら n が両方を兼ねる)、nc: グルーブの気泡核生成 (r_c_wick>0 のときだけ使う)。 グルーブの核沸騰の活性度 n (0〜1。grv_evap_gate = :nucleation のヒステリシス)
        T_wall = (
            hb   = zeros(1),
            evap = zeros(1),
            vl   = zeros(geo.vl.N),
            cond = zeros(geo.cond.N),
            ll   = zeros(geo.ll.N),
            cc   = zeros(geo.cc.N),
        ),
    )
end

"""
    initialize_differential_vars(geo::FlatLHPGeometry)

円筒型と同じ分類。secondary wick は primary wick と同じく E が微分変数、T が代数変数。
CC (cc) は、新しい CC (気相 1 + 液 N_CC_LAYER 層) では全て代数変数で、微分変数は ccl (Tl, X, mv)。
"""
function initialize_differential_vars(geo::FlatLHPGeometry)
    transport_mask(N) = (ρm=trues(N), u=trues(N), Em=trues(N),
                          P=falses(N), T=falses(N), x=falses(N))
    return ComponentArray(
        vl   = transport_mask(geo.vl.N),
        cond = transport_mask(geo.cond.N),
        ll   = transport_mask(geo.ll.N),
        grv  = transport_mask(geo.grv.N),
        cc   = (ρm = falses(geo.cc.N), u = falses(geo.cc.N), Em = falses(geo.cc.N),
                P = falses(geo.cc.N), T = falses(geo.cc.N), x = falses(geo.cc.N)),
        wick  = (E = trues(geo.wick.N), T = falses(geo.wick.N)),
        swick = (E = trues(2), T = falses(2)),
        plate = (T = trues(geo.plate.nx * geo.plate.ny),),
        cplate = (T = trues(n_cplate(geo)),),
        hi    = (T = trues(N_HI),),
        peek  = (T = trues(N_PEEK),),
        cca   = (T = trues(1),),
        eca   = (T = trues(1),),
        ccl   = (Tl = trues(N_CC_LAYER), X = trues(1), mv = trues(1)),
        nuc   = (n = trues(geo.grv.N), nc = trues(geo.grv.N)),
        T_wall = (
            hb   = trues(1),
            evap = trues(1),
            vl   = trues(geo.vl.N),
            cond = trues(geo.cond.N),
            ll   = trues(geo.ll.N),
            cc   = trues(geo.cc.N),
        ),
    )
end

# -----------------------------------------------------------------------------
# 2. 運転条件構造体
# -----------------------------------------------------------------------------

@with_kw struct LHPTransientCondition
    # --- 計算条件 ---
    Q_apply::Float64                   # ステップ後の印加熱量 [W] (最終値)
    t_init::Float64                    # 計算開始時刻 [s]
    t_end::Float64                     # 計算終了時刻 [s]
    saveat::Union{Float64,Vector{Float64}} = 1.0   # 結果を保存する時間間隔 [s]、または保存する時刻のベクトル
                                        # (IDAの内部刻み幅とは独立。出力間引き用)

    # --- 印加熱量のステップアップ (滑らかなステップ) ---
    # Q(t) = Q_apply_init + (Q_apply − Q_apply_init) s((t − t_step)/τ_step)
    #   s(x) = 0 (x≤0),  3x²−2x³ (0<x<1),  1 (x≥1)         (3次の smoothstep)
    # t ≤ t_step では厳密に Q_apply_init、t ≥ t_step+τ_step では厳密に Q_apply で、
    # その間を C¹ (1階微分が連続) につなぐ。不連続なステップは IDA の積分を難しく
    # するので、有限の立ち上がり時間 τ_step を持たせる。tanh 形と違い、t=t_init で
    # 厳密に Q_apply_init から始まり、区切り点 (t_step, t_step+τ_step) も決まるので、
    # IDA の tstops に渡せる (applied_heat_breakpoints)。
    # Q_apply_init == Q_apply なら、ステップなしの一定値になる。
    Q_apply_init::Float64 = Q_apply    # ステップ前の印加熱量 [W] (既定: ステップなし)
    t_step::Float64 = t_init           # ステップ (立ち上がり) の開始時刻 [s]
    τ_step::Float64 = 1.0              # 立ち上がりにかかる時間 [s]

    # --- 作動流体・境界条件 ---
    fluid::String = "ammonia"
    T_cond::Float64 = 15.0 + 273.15    # 凝縮器境界温度 [K]。平板型ではコールドプレート (チラー) の温度 (一定)
    # --- 凝縮器プレートのモデル (平板型。LHP の管が通る A5052 のプレート → 接触 h_cp → コールドプレート T_cond) ---
    # プレートの寸法は geo.plate (CondPlateGeometry)。管壁→プレートの熱伝達は h_con (外壁の半周)。
    f_con_area::Float64 = 0.5          # 管壁−プレートの接触面積 = f_con_area × 管の外周面積 (0.5 = 外周の半分、式51)
    h_cp::Float64 = 1000.0             # プレート−コールドプレートの接触熱伝達係数 [W/m²K]
    h_chiller::Float64 = -1.0          # チラー配管−コールドプレートの接触熱伝達係数 [W/m²K]。負なら h_con と同じ (凝縮器配管と同じ拡管)
    plate_k::Float64 = 140.0           # A5052 の熱伝導率 [W/mK]
    plate_ρ::Float64 = 2680.0          # A5052 の密度 [kg/m³]
    plate_c::Float64 = 880.0           # A5052 の比熱 [J/kgK]
    T_amb::Float64 = 25.0 + 273.15     # 周囲環境温度 [K]
    T_init::Float64 = 25.0 + 273.15    # 初期の基準温度 [K]。CC・ウィック・蒸発器ケースの温度、
                                        # ウィックの液密度の基準温度にも使う
    ΔT_startup::Float64 = 0.0          # 起動判定: |T_wick,out側 − T_CC| がこの値を超えたら蒸発を許す [K]
                                        # (0 なら判定なし。旧コードの fluid_wick の起動判定に相当)
    ε_startup::Float64 = 0.1           # 起動判定のゲート (smooth_step) の幅 [K]
    g::Float64 = 9.8                   # 重力加速度 [m/s2]

    # --- グルーブの核沸騰 (平板型。ドキュメント 2.8.1) ---
    # grv_boiling = false (既定) なら従来どおり (グルーブは蒸気だけ、h は強制対流)。
    # true にすると、グルーブが濡れている (液がある) 間は
    #   ・蒸発器ケース→グルーブの h を、自然対流 ⇔ Chen (核沸騰、サブクール・飽和とも) でブレンド
    #   ・ウィック界面の蒸発 ṁ を 0 にし、Q_ec-wick を最外層の熱源にする (流速 0)
    # 質量源は足さない (グルーブ内の液→蒸気は、均質モデルでは Em が増えて x が上がるだけ)。
    grv_boiling::Bool = false
    r_c_grv::Float64 = 2.0e-6          # 活性キャビティの開口半径 [m] (式82)。【仮の値。要設定】
    ε_onset::Float64 = 0.1             # 核沸騰開始の判定 (smooth_step) の幅 [K]
    # ウィック界面蒸発のゲートを、グルーブの気泡核生成と分ける (grv_evap_gate = :nucleation のとき):
    #   r_c_wick > 0 : 界面蒸発のゲート n は ΔT = 2σT/(ρ_v h_fg r_c_wick) (細孔半径)、幅 ε_onset_wick で開き、
    #                  グルーブの沸騰の h (Chen) は別の状態 nc (キャビティ半径 r_c_grv、幅 ε_onset) で切り替える。
    #   r_c_wick = 0 : 従来どおり、1 つの n (r_c_grv, ε_onset) が両方を動かす。
    r_c_wick::Float64 = 0.0
    # 界面蒸発のゲート n の【開始】の判定に使う温度 (r_c_wick > 0 のとき): :ec = 蒸発器ケース (受熱面 P。従来)、
    # :wick = primary wick の最外セル (x.wick.T[end]。メニスカスのある側の液の温度)。維持の判定 (n > n_lat) は、どちらも T_ec。
    onset_ref::Symbol = :ec
    ε_onset_wick::Float64 = 0.1
    x_hold_wick::Float64 = 2.0         # 界面蒸発のゲート n の x による維持 (グルーブの x がこれ以上の間は開く)。2 なら無効 (ラッチ n_lat だけで維持)
    f_wet_th::Float64 = 1.0e-3         # 濡れ判定: 液の体積分率 f_l=(1−x)ρm/ρl がこれ以下で乾き
    ε_wet::Float64 = 5.0e-4            # 濡れ判定の smooth_step の幅 (f_l 単位)
    # 沸騰の h (Chen) は、液膜が薄くなると効かなくなる。h だけ、より大きい f_l で乾きに移す
    # (液がほぼ無いのに Chen の h のまま蒸気を加熱すると、熱容量の小さい蒸気が数 ms で壁温に
    #  近づく、硬い系になる)。f_l がこの値以下で、h は蒸気側の強制対流に移る。
    f_wet_th_h::Float64 = 2.0e-2
    # グルーブの成層した出口 (液は底に溜まり、蒸気だけが出る)。液が溝を f_full 以上満たすと液も出る
    grv_stratified::Bool = false
    f_full_grv::Float64 = 0.9
    ε_full_grv::Float64 = 0.03
    ε_wet_h::Float64 = 1.0e-2
    grv_evap_gate::Symbol = :dry       # 界面蒸発を止める条件: :dry (液がある間) | :saturation (グルーブがサブクールの間だけ)
                                        #   | :onset (核沸騰開始 T_ec > T_sat + ΔT_onset(r_c_grv) まで。グルーブの沸騰と同じ判定。
                                        #            グルーブの x > x_latch_onset で開いたままにする)
                                        #   | :nucleation (:onset と同じ判定を、状態変数 nuc.n のヒステリシスつきで行う。平板型のみ)
    ΔT_sat_th::Float64 = 0.2           # :saturation のとき、サブクール度がこの値より大きいと蒸発を止める [K]
    ε_sat::Float64 = 0.05              # 同上のゲートの幅 [K] (飽和ちょうどでゲートが全開になる値)
    x_liq_gate::Float64 = 0.05         # :saturation のサブクール判定を使うのは x がこれより小さいセルだけ
    ε_x_gate::Float64 = 0.02
    # :onset のとき、界面蒸発のゲートを開いたままにする (沸騰開始のヒステリシス) グルーブの x の閾値。
    # 成層した出口ではグルーブの x は 0.005 程度までしか上がらないので、:saturation の x_liq_gate (0.05) とは別にする。
    x_latch_onset::Float64 = 1.0e-3
    ε_x_latch::Float64 = 5.0e-4
    # grv_evap_gate = :nucleation のとき (平板型。状態変数 nuc.n を使う沸騰開始のヒステリシス):
    #   dn/dt = (n_target − n)/τ_nuc,   n_target = smooth_step(T_ec − T_sat − ΔT_eff(n), ε_onset)
    #   ΔT_eff = ΔT_hold + (ΔT_onset − ΔT_hold)(1 − s_n),   s_n = smooth_step(n − n_lat_nuc, ε_n_nuc),   ΔT_hold = f_hold_onset ΔT_onset
    # n が n_lat_nuc 未満では開始の閾値 ΔT_onset、超える (沸騰中) と維持の閾値 ΔT_hold (< ΔT_onset)。
    # 界面蒸発のゲートと、グルーブの沸騰の h (Chen) の重みの両方が、この n を使う (同じ判定)。
    f_hold_onset::Float64 = 1.0
    # 維持をグルーブの蒸気量で決める: x ≥ x_hold_nuc の間は沸騰 (n → 1) を続け、x が下回ると沸騰をやめる
    # (開始は過熱度 ΔT_onset)。n_target = OR(過熱度の判定, smooth_step(x − x_hold_nuc, ε_x_nuc))。
    # 2.0 (既定) なら無効。成層した出口ではグルーブの x は 0.005 程度までなので 0.003 前後を使う。
    x_hold_nuc::Float64 = 2.0
    ε_x_nuc::Float64 = 5.0e-4
    τ_nuc::Float64 = 0.5               # n の時定数 [s]
    # ラッチ (双安定): n が n_lat_nuc (幅 ε_n_nuc) を超えると閾値が ΔT_onset から ΔT_hold = f_hold_onset ΔT_onset に切り替わる。
    # f_hold_onset < 0 なら、維持の閾値が負 (EC が飽和温度より少し低くても沸騰が続く)。
    n_lat_nuc::Float64 = 0.15
    ε_n_nuc::Float64 = 0.03
    # 【注】平板型 (FlatLHPGeometry) は、気相 1 セル + 液 16 層の新しい CC (residuals_flatcc.jl) に置き換えた。ΔP_ne は、気相の圧力 P と界面の液温 Tl[1] の
    # 飽和圧の差として結果に出る (界面の運動論的な凝縮 (97)〜(99) 式)。したがって、下の τ_ne_cc 〜 τ_layer_cc、cc_layer*、cc_meniscus、
    # ne_liq_decouple、δ_layer_cc は、平板型では使われない。cc_wall_cond は平板型でも使う (液面の凝縮と壁面の凝縮の比 = 熱浸透率 × 面積の比)。円筒型 (旧版の CC) では従来どおり。
    # CC の有限速度の凝縮 (非平衡)。CC は均質・飽和の1ノードなので、液の流入で蒸気が圧縮されると、平衡だと凝縮が即座に
    # 追いついて圧力はほとんど上がらない。実際は凝縮に有限の時間がかかり、蒸気側の圧力が飽和圧より ΔP_ne だけ高くなる:
    #   dΔP_ne/dt = K_ne_cc · q_l / V_v − ΔP_ne / τ_ne_cc ,  q_l = (ṁ_in − ṁ_out)/ρ_l (蒸気を押しのける液の体積流量)
    #   CC の圧力 P = P_sat(T) + ΔP_ne  (V_v = max(f_v, f_v_min_ne)·V_cc)
    # τ_ne_cc = 0 (既定) で無効。K_ne_cc は蒸気の体積弾性率 ≈ γP [Pa]。定常的な過剰は K q τ / V_v。
    τ_ne_cc::Float64 = 0.0
    K_ne_cc::Float64 = 1.2e6
    f_v_min_ne::Float64 = 0.05
    # 蒸気がほとんど無いとき (f_v ≪ f_v_ref_ne) は、小さな蒸気塊の凝縮は速い (界面に対して体積が小さい) ので、緩和時間を短くする:
    #   τ_eff = τ_ne_cc · f_v²/(f_v² + f_v_ref_ne²)。蒸気が消えるにつれ ΔP_ne も消える (以前は V_v の下限で押し込み続けて ΔP_ne が 20 kPa まで積み上がった)
    f_v_ref_ne::Float64 = 0.2
    # 【仮】CC の気液界面の温度層 (考察用。最終結果には使わない)。cc_layer = true のとき、ΔP_ne の式を次に置き換える
    # (τ_ne_cc > 0 のときだけ有効。同じ状態変数 ccne.dP を使う):
    #   dΔP/dt = (dP_sat/dT) · ρ_v h_fg q_l · w_v / C_layer − ΔP / τ_layer ,   w_v = f_v²/(f_v² + f_v_ref_ne²)
    # 液が蒸気を押しのけて、蒸気が界面で凝縮するときの熱 ρ_v h_fg q_l が、界面近傍の液層 (熱容量 C_layer) に入り、
    # 界面温度が上がる (ΔT = ΔP/(dP_sat/dT))。層は拡散で時定数 τ_layer で戻る。熱は T_cc には足さない (エネルギー勘定は未整合)。
    cc_layer::Bool = false
    # 【物理版】界面の液側の 1 次元熱伝導層 (cc_layer_cond = true。cc_layer より優先)。自由なフィット係数なし:
    #   凝縮熱 Q = ρ_v h_fg q_l w_v が界面 (面積 A_if = 2 sqrt(R²−h_liq²) · dl) から液側へ入り、
    #   ρ c_p ∂θ/∂t = ∂/∂z (k ∂θ/∂z) (z: 液の深さ方向、下端は断熱、厚さ δ_layer_cc)、θ(表面) = 界面の過熱度、
    #   ΔP_ne = (dP_sat/dT) θ(表面)。ρ, c_p, k は飽和液の物性。A_if は液深さ h_liq (f_v から) で変わる。
    cc_layer_cond::Bool = false
    # 【数値実験】液側 (ウィック・二次ウィック・液管出口) が見る CC の圧力を P_cc − ΔP_ne にする (蒸気側の過剰圧力 ΔP_ne を液側へ伝えない)。
    # 物理的な根拠はない (自由界面なら液の圧力は蒸気の圧力と等しい。二次ウィックの毛細管圧は 2σ/r ≈ 1.5 kPa)。ΔP_ne の影響の切り分け用
    ne_liq_decouple::Bool = false
    # 蒸気がなくなるときの ΔP_ne の消滅 (cc_layer_cond のとき): 適用する過剰圧力 = g_ex · (dP_sat/dT) θ(表面)、g_ex = f_v²/(f_v² + f_ex²)。
    # CC の蒸気がほぼ尽きると (単相の液)、蒸気だけが持つ過剰圧力はなくなる (圧力は液の状態から決まる)。f_ex は、実験の ΔP_ne の減衰
    # (+42 s: 13.9、+46 s: 10.5、+50 s: 4.5、+58 s: 0.9 kPa。ピーク 16.5 kPa) を、前線の通過時刻から逆算した蒸気体積 f_v(t) に対して
    # 見たときの値 (f_ex = 0.03〜0.05) の中央。0 で無効。
    f_ex_layer_cc::Float64 = 0.035
    cc_meniscus::Bool = false          # true: ΔP_ne の消滅の尺度を、蒸気キャップの高さ s と毛管長さ l_c から決める (g = s²/(s²+(f_men l_c)²))。f_ex_layer_cc は使われない
    f_men_layer_cc::Float64 = 1.0      # g = 1/2 になる s / l_c
    # 壁面凝縮 (cc_layer_cond のとき): true なら、凝縮熱 ρ_v h_fg q_l を液面の層と、蒸気に触れる secondary wick の壁面で、
    # 熱浸透率 √(k ρc) × 面積 の比で分け合う (面積は cc_vapor_geometry の A_if と A_vap)。液面の層に入る割合は w_v の代わりにこれを使い (f_v_ref_ne は使われない)、
    # 壁面に入る熱は CC の流体から引いて secondary wick の節点 (swick) に足す。
    cc_wall_cond::Bool = false
    cc_wall_dT::Bool = true       # 壁面凝縮に、駆動する温度差の比 g_T = (T_v − T_swick)/(T_v − T̄_液) をかける (平板型・cc_wall_cond = true のとき)。false で熱浸透率×面積の比のみ (v8c)
    wall_radial_R::Bool = false   # 管壁 (SUS304) の半径方向の熱伝導抵抗を入れる: 流体−壁の h と 壁−プレートの h_con に、管壁の中心を境に半分ずつ直列。false で壁は半径方向に等温 (従来)
    δ_layer_cc::Float64 = 9.0e-3       # 熱伝導層を解く深さ [m] (液深さ h_liq ≈ 9.2 mm にそろえる。下端は断熱で近似)
    C_layer_cc::Float64 = 17.0         # 層の熱容量 [J/K] (CC の液全体 ρ c V ≈ 16 J/K 程度)
    τ_layer_cc::Float64 = 10.0         # 層の緩和時間 [s] (実験の EC・VL1 の戻り +49 → +60 s が約 10 s)
    # 凝縮の熱伝達 (Shah) にかける係数 (1 = そのまま)。蒸気管 (VL) と凝縮器 (cond) に同じ値を使う (同じ内径 1.76 mm、同じ管)。
    # Shah の重力項 h_Nu は大径管の水平管向けで、小径管 (D = 1.76 mm、ボンド数 ≈ 0.8) では過大になる可能性がある。
    # 値は凝縮器の CON1 (管壁) の実測から決めた (0.4〜0.5 程度)。VL では実験から確認できていない (効果も小さい: 壁の熱容量で決まる)
    f_h_shah::Float64 = 1.0
    f_h_shah_ll::Float64 = -1.0        # LL だけの係数。負なら f_h_shah と同じ (切り分け用)
    ε_phase_h::Float64 = 0.05          # 流体−壁の熱伝達 h の、単相⇔二相のブレンド幅 (疑似乾き度)。既定 0.05 (従来)。小さくすると、飽和から離れた液体・蒸気が二相側の式 (Shah の凝縮) を使わなくなる。摩擦のブレンド幅は別 (0.05 のまま)
    h_xps_mode::Symbol = :P            # h のブレンド用の疑似乾き度: :P = 圧力の飽和温度基準 (過冷却液で負)、:T = 従来 (セル自身の温度基準。切り分け用)
    C_cc_extra::Float64 = 0.0          # CC ケースに加える熱容量 [J/K] (幾何の V_w に入らない分)。cc_2node = true のときは節点 A (蒸発器側) の熱容量
    # CC ケースを 2 節点にする: A (蒸発器側、ドームの下部と端面。熱容量 C_cc_extra) と B (頂部側 = T_wall.cc、ドームの曲面。幾何の V_w)。
    # 蒸発器ケース → A → B の順に熱が伝わり、B が secondary wick・液管・周囲と熱をやり取りする (RES のセンサは B)。
    # 熱コンダクタンスは k_SUS304 と図面の断面・経路 (FlatCCGeometry の W_o, t_shell, t_end) から決まり、f_ec_cc を掛ける。
    cc_2node::Bool = false
    f_ec_cc::Float64 = 1.0             # 蒸発器ケース→CC ケースの熱伝導 (Q_ec_cc) にかける係数 (漏れの大きさの調整)
    # 蒸発器ケースの足 (HB にボルト締め) を通る、HB → CC ケースの熱伝導 [W/K]。Q_foot = G_foot_cc (T_hb − T_cca) (cc_2node = true: 節点 A へ、false: CC ケースへ)。
    # HB から出て CC ケースへ入る (内部の熱の移動。エネルギーは保存)。0 なら従来どおり。値は HB から CC ケースまでの全体 (接触抵抗を含む) の見積り (実測に基づく。フィッティングではない)
    G_foot_cc::Float64 = 0.0
    # 孔内の液 (primary wick・secondary wick) の熱膨張: true なら、孔内の液の密度を ρ_l,sat(T) にして、孔内の液の質量の変化 (昇温で減る分) を CC の液へ戻す
    # (CC から wick へ出る流量 = ṁ_evap + dm_pore/dt)。false なら従来 (孔内の液の質量は一定)。wick のエネルギーの式の密度は ρ_l_ref のまま (近似)
    pore_expansion::Bool = false
    Q_cc_fluid::Float64 = 0.0          # 【診断用】CC の流体 (液) へ入る一定の熱 [W]。t = 0 から常に入る。壁ではなく流体に直接入る (圧力計の取付け口など、流体に触れる部分からの入熱を想定)。物理モデルではない
    cc_boiling::Bool = false           # CC の核沸騰 (secondary wick → CC液。プール沸騰の h = Chen)
    r_c_cc::Float64 = 5.0e-6           # CC の活性キャビティ開口半径 [m] (secondary wick の細孔半径)
    ε_onset_cc::Float64 = 0.1          # CC の沸騰開始ゲートの幅 [K]
    # --- CC の液の混合 (residuals_flatcc.jl の cc_mixing。既定 :none = 従来と完全に同じ) ---
    # :jet   = 液管からの流入の噴流 (巻き込み Q = 8πνx、浮力で落ち着く層を決める)
    # :plume = 壁・床の自然対流の境界層プルーム (熱を中立の高さの層まで運ぶ)
    # :both  = 両方。N_CC_LAYER = 1 のときは何も変わらない
    cc_mix::Symbol = :none
    cc_mix_θ::Float64 = 1.0            # プルームの出口温度 T_b = T_層 + θ (T_壁 − T_層)。θ = 1 が流量最小 = 持ち上がる高さ最大の上限 (感度パラメータ。フィッティングではない)
    cc_mix_ε_T::Float64 = 1.0e-3       # 【数値の幅。物理量ではない】層の選択の滑らかな重み (シグモイド) の温度幅 [K]
    cc_mix_D_port::Float64 = 1.76e-3   # 【実測】流入部の管の内径 [m] (1/8 in 管、外径 3.18 mm、内径 1.76 mm)
    cc_mix_z_port::Float64 = 4.5e-3    # 【実測】流入ポートの中心の高さ (secondary wick の内底面から) [m]
    cc_mix_x_jet::Float64 = 0.025      # 【仮定】噴流が壁に当たるまでの距離 [m] (ポートの反対側の端面まで = CC の長さ)
    cc_mix_ent::Float64 = 1.0          # 巻き込みにかける係数 (1 = 層流の円形噴流の式そのまま。感度試験用)
    q_evap_max::Float64 = 0.0          # 蒸発熱 ṁλ の上限 [W] (tanh で滑らかに飽和。0=無効)
    grv_init::Symbol = :vapor          # グルーブの初期状態 (:vapor または :liquid=満水)
    f_l_grv_init::Float64 = 0.0        # grv_init=:vapor のとき、グルーブに初めから入れる液の体積分率 (0 = 蒸気 x=1 のみ。飽和の二相で CC と同様に置く)

    # --- 整合初期条件 (initialcondition.jl) ---
    # 定常状態は求めず、静止 (u=0) の平衡に近い状態から始める。液・蒸気は、表の補間
    # 誤差の影響を受けやすい飽和線の上を避けるため、基準温度 T_init から少しずらす
    # (液は過冷却、蒸気は過熱)。壁は、接する流体と同じ温度から始める。
    ΔT_sub_init::Float64 = 1.0         # 液 (凝縮器・液管) の過冷却度 [K]。正の値
    ΔT_sup_init::Float64 = 1.0         # 蒸気 (グルーブ・蒸気管) の過熱度 [K]。正の値
    f_l_cc_init::Float64 = 0.5         # CCの初期の液体積分率 (0〜1)。CCは T_init の気液二相
    cond_init_phase::Symbol = :liquid  # 凝縮器の初期状態 (:liquid または :vapor)

    # --- 壁の初期温度 (実験の加熱前の壁温度を、そのまま初期値にしたいとき) ---
    # nothing (既定) なら、これまでどおり (蒸発器ケースは界面と同じ T_sat(P_grv)、CCケースは T_init、
    # 蒸気管・凝縮器・液管は接する流体と同じ温度)。値を指定した壁だけ、その温度から始める。
    # 流体の初期温度は変えない (CC は T_init の二相、液は T_init−ΔT_sub_init、蒸気は
    # T_init+ΔT_sup_init)。流体の熱容量は壁より小さいので、壁の温度を合わせれば、起動時の
    # 熱的な応答はほぼ決まる。流体と壁の温度差は、最初に流体側が壁に近づくことで解消される。
    T_wall_init_evap::Union{Nothing,Float64,Vector{Float64}} = nothing   # 蒸発器ケースの壁 [K]
    T_wall_init_cc::Union{Nothing,Float64,Vector{Float64}}   = nothing   # CCケースの壁 [K]
    T_wall_init_vl::Union{Nothing,Float64,Vector{Float64}}   = nothing   # 蒸気管の壁 [K]
    T_wall_init_cond::Union{Nothing,Float64,Vector{Float64}} = nothing   # 凝縮器の壁 [K]
    T_wall_init_ll::Union{Nothing,Float64,Vector{Float64}}   = nothing   # 液管の壁 [K]
    # 実測の初期条件 (initial_from_measurements.jl): 管のセルごとの初期の位相 (:liquid / :vapor)。
    # nothing なら、従来どおり (蒸気管 = 蒸気、凝縮器 = cond_init_phase、液管 = 液)。
    # T_wall_init_* をベクトル (セルごと) にもできる。流体の温度は壁に合わせる (飽和線から ΔT_sub/sup_init 離す)。
    init_phase_vl::Union{Nothing,Vector{Symbol}}   = nothing
    init_phase_cond::Union{Nothing,Vector{Symbol}} = nothing
    init_phase_ll::Union{Nothing,Vector{Symbol}}   = nothing    # 要素: :liquid | :vapor | :sat (飽和に近い液。下の ΔT_sub_ll_sat)
    x_vap_vl_init::Float64 = 1.0       # 蒸気管の初期蒸気セルを、飽和の気液二相 (乾き度 x) にする (1 なら使わず、表の単相蒸気 = 飽和 +0.7 K)
    f_v_ll_init::Float64 = 0.0         # 液管の CC 側の末端 n_ll_2ph セルを、飽和の気液二相 (蒸気の体積分率 f_v) にする (0 なら使わない)
    n_ll_2ph::Int = 1                  # 上の対象セル数 (CC 側から数える。1 = 最終セル = LL3 を含む)
    ΔT_sub_ll_sat::Float64 = 0.02      # init_phase が :sat のセルの過冷却度 [K] (表が単相液にならなければ、0.04, 0.06, ΔT_sub_init の順に大きくする)
    # ヒーターブロック HB (平板型): 集中熱容量の独立したノード (T_wall.hb)。ヒーターの熱 Q_apply は HB に入り、
    # 接触熱伝達 G = h_hb_ec × A_hb_ec で蒸発器ケースへ伝わる。
    #   実験の HB: A5052 53×33×10 mm ≈ 2680 kg/m³ × 17.49 cm³ × 900 J/(kg K) ≈ 42 J/K
    C_hb::Float64 = 42.0                 # HB の熱容量 [J/K]
    C_ec_extra::Float64 = 0.0            # 蒸発器ケースに加える熱容量 [J/K] (固定用の足など、幾何の V_w に入らないもの)。ec_2node = true のときは使わない
    # 蒸発器ケースを 2 節点にする: P (受熱面。T_wall.evap。HB・ウィック・グルーブ・蒸気管・CC ケースとの熱のやり取りと、
    # 核沸騰のゲートの T_ec はすべて P。熱容量 C_ec_P) と S (側面。EC1・EC2 のセンサ。熱容量 C_ec_S)。
    # P と S の間は G_ec_ps で結ぶ。S は P だけと熱をやり取りする (行き止まり)。
    ec_2node::Bool = false
    C_ec_P::Float64 = 0.0                # 節点 P の熱容量 [J/K] (ec_2node = true のとき。幾何の V_w の代わりにこの値を使う)
    C_ec_S::Float64 = 0.0                # 節点 S の熱容量 [J/K]
    G_ec_ps::Float64 = 0.0               # P–S の熱コンダクタンス [W/K]
    h_hb_ec::Float64 = 500.0             # HB–蒸発器ケースの接触熱伝達係数 [W/(m² K)]
    A_hb_ec::Float64 = 0.032 * 0.032     # 接触面積 [m²] (蒸発器ケースの受熱面 32×32 mm)
    T_wall_init_hb::Union{Nothing,Float64} = nothing   # HB の初期温度 [K] (nothing なら蒸発器ケースと同じ)

    # --- 蒸気管・液管と周囲との熱交換 ---
    #   :bare       … 裸管の自然対流 (従来)
    #   :insulated  … 断熱材 (内径 ins_D_in、外径 ins_D_out、熱伝導率 ins_k) の穴に管が通っている: 空気層 + 断熱材 + 外表面 (h_o 一定) の直列抵抗
    #   :adiabatic  … 周囲と熱交換しない
    amb_mode::Symbol = :bare
    ins_k::Float64 = 0.036             # 断熱材 (エアロフレックスなどのエラストマー系発泡) の熱伝導率 [W/(m K)]
    ins_D_in::Float64 = 0.010          # 断熱材 (エアロフレックス) の内径 [m]。管はこの穴の中に通っていて、管と穴の間は空気層
    ins_D_out::Float64 = 0.030         # 断熱材の外径 [m]
    ins_k_air::Float64 = 0.026         # 空気層の熱伝導率 [W/(m K)] (伝導のみ)
    ins_h_o::Float64 = 8.0             # 断熱材の外表面の熱伝達係数 (自然対流 + 放射) [W/(m2 K)]

    # --- 平板型: HB・蒸発器ケース・CC の上側を包む断熱材 (エアロフレックス、1 枚) と、HB の下の PEEK ブロック (residuals_housing.jl) ---
    hi_on::Bool = false                # true なら、断熱材 3 節点 (HB 側・蒸発器ケース側・CC 側) を壁につなぐ
    cc_adiabatic::Bool = false         # 【試験用】 true なら、CC ケース → 断熱材 (節点 3) の熱伝導を 0 にする (CC ケースを完全断熱にした上限の試験)。節点 3 は残す (熱容量 0 の節点は DAE が特異になる)。hi_on = true が必要
    hi_t::Float64 = 0.020              # 断熱材の厚さ [m]
    hi_rcp::Float64 = 1.1e5            # 断熱材の ρ c_p [J/(m3 K)] (ρ ≈ 60 kg/m3、c_p ≈ 1800 J/(kg K) の目安。要確認)
    hi_A_hb::Float64 = 2 * 0.032 * 0.010   # HB 側の内面積 [m2] (HB の側面 2 面 × 厚さ 10 mm。覆われるのは側面だけ)
    hi_A_ev::Float64 = 2 * 0.032 * 0.0055   # 蒸発器ケース側の内面積 [m2] (側面 2 面 × 高さ 5.5 mm)
    hi_A_cc::Float64 = π * 0.032 / 2 * 0.025   # CC 側の内面積 [m2] (半円筒の曲面 = 半周 × 軸方向長さ 25 mm。両端面は含めない)
    hi_G_hb_ev::Float64 = 5.9e-3       # 断熱材の中の横方向の熱コンダクタンス HB 側 ↔ 蒸発器ケース側 [W/K] (k t (2×32 mm)/7.75 mm の目安。要確認)
    hi_G_ev_cc::Float64 = 1.2e-3       # 同 蒸発器ケース側 ↔ CC 側 [W/K] (k t 32 mm/20 mm の目安。要確認)
    peek_on::Bool = false              # true なら、HB の下の PEEK ブロックを 1 次元の熱伝導でつなぐ (下端は T_amb)
    peek_L::Float64 = 0.060            # PEEK ブロックの高さ [m]
    peek_A::Float64 = 0.032 * 0.032    # 断面積 [m2] (HB と同じ面積)
    peek_k::Float64 = 0.25             # PEEK の熱伝導率 [W/(m K)]
    peek_rcp::Float64 = 1.74e6         # PEEK の ρ c_p [J/(m3 K)] (ρ 1300 kg/m3、c_p 1340 J/(kg K))
    # 液管の CC 側の末端 n_ll_ci セルは、CC を包む断熱材と一続きの断熱材の中にある: 周囲温度を T_amb ではなく CC 側の断熱材の節点 (hi.T[3]) の温度にし、
    # その熱は節点 3 へ入る (エネルギー保存)。0 で無効 (従来: 全長が T_amb)。hi_on = true が必要。
    n_ll_ci::Int = 0

    # --- 熱橋: 熱源 (HB または蒸発器ケース) → 液管の CC 側の末端 ---
    # 平板型は CC・蒸発器・HB が一体で、液管の末端が熱源に近い。G_br [W/K] の熱伝導で、液管の CC 側から n_br セルの壁を温める。
    G_br::Float64 = 0.0
    n_br::Int = 3                      # 熱を配分する、CC 側の末端からのセル数
    br_decay::Float64 = Inf            # 配分の減衰長 [セル]。熱は exp(−k/br_decay) で末端から減る (Inf なら一様)
    br_src::Symbol = :hb               # 熱の取り元: :hb (ヒーターブロック) | :ec (蒸発器ケース)

    # --- 熱伝達係数 (定数) ---
    h_ev::Float64  = 10000.0           # 蒸発熱伝達係数 (Kaya[14]を参考) [W/m2K]。h_ev_end が NaN でなければ「初期値」(加熱前〜初期)
    # 【試験用: h_ev の時間変化】h_ev_end が NaN (既定) なら h_ev は定数 (従来どおり)。
    # NaN でなければ h_ev(t) = h_ev から h_ev_end へ、t = h_ev_t_mid を中心に tanh で滑らかに移る (対数空間で補間 = 比で滑らかに変わる)。物理モデルではなく what-if。
    h_ev_end::Float64   = NaN          # 後半の h_ev [W/m2K]
    h_ev_t_mid::Float64 = 60.0         # 遷移の中心のモデル時刻 [s] (τ + 加熱開始 10 s)
    h_ev_width::Float64 = 15.0         # tanh の長さスケール [s] (遷移はおよそ中心 ± 1.1×この値 で 75 % 進む)
    h_con::Float64 = 500.0             # 凝縮器プレートとの接触熱伝達係数 (式51) [W/m2K]

    # --- 風上ブレンドの幅 ---
    ε_u::Float64 = 1.0e-4              # 輸送管の速度の判定幅 [m/s] (流速の典型値より十分小さく)

    # --- スケーリング定数 ---
    # IDAのNewton反復・誤差評価は無次元の残差ノルムに基づくため、
    # 変数間の桁違い (P~1e6, x~O(1), T~1e2 など) を解消しておくと
    # 収束性が大きく改善する。state_scaleのみ用い、scale_var/restore_var
    # (以前のTrustRegion実装) と同じ考え方を新変数セットに合わせて再定義する。
    ρm_scale::Float64 = 50.0           # [kg/m3]  (気液平均密度の代表値)
    u_scale::Float64  = 10.0           # [m/s]
    Em_scale::Float64 = 4.0e5          # [J/kg]   (全エネルギーの代表値)
    P_scale::Float64  = 1.0e6          # [Pa]
    T_scale::Float64  = 300.0          # [K]
    x_scale::Float64  = 1.0            # [-]      (乾き度はもともとO(1))

    # --- 反復対象コンポーネントの列挙 (ループ処理の共通化に利用) ---
    fluid_comp::NTuple{5,Symbol} = (:grv, :vl, :cond, :ll, :cc)  # TransportState系
    wall_comp::NTuple{5,Symbol}  = (:evap, :vl, :cond, :ll, :cc) # T_wall系
end

# -----------------------------------------------------------------------------
# 印加熱量 Q_apply(t)
# -----------------------------------------------------------------------------

"""
    smoothstep3(x)

3次の smoothstep: x≤0 で 0、x≥1 で 1、その間は 3x²−2x³ (C¹ 連続)。
"""
smoothstep3(x) = x <= 0 ? 0.0 : (x >= 1 ? 1.0 : x^2 * (3 - 2x))

"""
    applied_heat(t, condition) -> Q_apply(t) [W]

時刻 t の印加熱量。Q_apply_init から Q_apply へ、t_step から τ_step かけて
滑らかにステップアップする (LHPTransientCondition のコメント参照)。
引数 t は時刻 (状態変数ではないので、ここでの区間分けの分岐は、残差の
状態依存の不連続を生まない)。
"""
applied_heat(t, c::LHPTransientCondition) =
    c.Q_apply_init + (c.Q_apply - c.Q_apply_init) * smoothstep3((t - c.t_step) / c.τ_step)

"""
    applied_heat_breakpoints(condition) -> Vector{Float64}

Q_apply(t) の区切り点 (立ち上がりの開始・終了時刻)。IDA に tstops として渡すと、
ここで積分が確実に刻まれ、C¹ の接続点を刻み幅が飛び越えない。
Q_apply_init == Q_apply (ステップなし) のときは空。計算区間外の点は除く。
"""
function applied_heat_breakpoints(c::LHPTransientCondition)
    c.Q_apply_init == c.Q_apply && return Float64[]
    # 立ち上がり (τ_step) の途中にも 0.1 s 刻みで区切り点を置く。区切り点は IDA の刻みを止めるだけなので、両端だけだと、
    # 後退オイラー (次数 1) が 1 s 級の刻みで立ち上がりをまたぎ、入熱を「刻みの終わりの値 × 刻み」で数えて、
    # 台形則の積分より多く熱を入れてしまう (v3: +10 s で 3.8 J = 入熱の 13 %。0.1 s 刻みなら 0.3 J 以下)。
    n_seg = max(1, ceil(Int, c.τ_step / 0.1))
    pts = [c.t_step + c.τ_step * k / n_seg for k in 0:n_seg]
    return filter(t -> c.t_init < t <= c.t_end, pts)
end

# -----------------------------------------------------------------------------
# 3. 状態ベクトルの組み立て
# -----------------------------------------------------------------------------

"""
    _zeros_transport(N)

輸送管・グルーブ用の6変数 (ρm,u,Em,P,T,x) をすべてゼロ初期化した
NamedTuple を返す。ComponentArray の入れ子として使う。

要素型は Float64 に固定している。疎行列パターン (jac_prototype) は
自動検出ツールに頼らず手動で構築する方針としたため、Dual数対応は
不要と判断した。
"""
_zeros_transport(N) = (ρm=zeros(N), u=zeros(N), Em=zeros(N), P=zeros(N), T=zeros(N), x=zeros(N))

"""
    initialize_state(geo::LHPGeometry) -> ComponentArray

与えられたジオメトリに対応する形状の状態ベクトルを、すべてゼロで
初期化して返す。実際の値は initialcondition.jl 側の
build_consistent_ic で、環境温度・封入量などから物理的に設定する。
"""
function initialize_state(geo::LHPGeometry)
    return ComponentArray(
        vl   = _zeros_transport(geo.vl.N),
        cond = _zeros_transport(geo.cond.N),
        ll   = _zeros_transport(geo.ll.N),
        grv  = _zeros_transport(geo.grv.N),
        cc   = _zeros_transport(geo.cc.N),   # θ_liqは状態に持たず、residuals_cc.jl内で
                                              # f_l → θ_liq の1次元テーブル補間として
                                              # その場で計算する(DAEの未知数には含めない)
        wick = (E = zeros(geo.wick.N), T = zeros(geo.wick.N)),
        T_wall = (
            evap = zeros(1),
            vl   = zeros(geo.vl.N),
            cond = zeros(geo.cond.N),
            ll   = zeros(geo.ll.N),
            cc   = zeros(geo.cc.N),
        ),
    )
end

"""
    initialize_differential_vars(geo::LHPGeometry) -> ComponentArray{Bool}

initialize_state と同じ形状で、各スロットが微分変数 (true) か
代数変数 (false) かを示す Bool の ComponentArray を返す。
DAEProblem の differential_vars にそのまま渡す。

対応関係 (これまでの議論で確定した分類):
    ρm, u, Em        … 微分変数 (質量・運動量・エネルギー保存則)
                       (ただしCCのuは代数変数。CCは貯水池として扱い、運動量式を使わない)
    P, T, x          … 代数変数 (CoolPropのflash拘束)
    wick.E           … 微分変数 (式47のエネルギー保存則)
    wick.T           … 代数変数 (式46とCoolPropの整合性拘束)
    T_wall.*         … すべて微分変数 (式20/41、代数変数を伴わない純粋なODE)
"""
function initialize_differential_vars(geo::LHPGeometry)
    transport_mask(N) = (ρm=trues(N), u=trues(N), Em=trues(N),
                          P=falses(N), T=falses(N), x=falses(N))

    return ComponentArray(
        vl   = transport_mask(geo.vl.N),
        cond = transport_mask(geo.cond.N),
        ll   = transport_mask(geo.ll.N),
        grv  = transport_mask(geo.grv.N),
        # CCは貯水池 (0次元) として扱う (residuals_cc.jl)。運動量式は使わず、
        # 出口面の速度 u は代数拘束 u = ṁ_out/(ρ_l A) で決まるため、
        # ρm, Em だけが微分変数で、u は P, T, x と同じ代数変数。
        cc   = merge(transport_mask(geo.cc.N), (u = falses(geo.cc.N),)),
        wick = (E = trues(geo.wick.N), T = falses(geo.wick.N)),
        T_wall = (
            evap = trues(1),
            vl   = trues(geo.vl.N),
            cond = trues(geo.cond.N),
            ll   = trues(geo.ll.N),
            cc   = trues(geo.cc.N),
        ),
    )
end

"""
    scale_state(u, condition::LHPTransientCondition) -> ComponentArray

状態ベクトルを無次元化する。IDAに渡す前に適用する。
対応する復元関数 restore_state とペアで使うこと。
"""
function scale_state(u, condition::LHPTransientCondition)
    scaled = copy(u)
    for comp in condition.fluid_comp
        c = getproperty(scaled, comp)
        c.ρm ./= condition.ρm_scale
        c.u  ./= condition.u_scale
        c.Em ./= condition.Em_scale
        c.P  ./= condition.P_scale
        c.T  ./= condition.T_scale
        c.x  ./= condition.x_scale
    end
    scaled.wick.E ./= condition.Em_scale   # ウィックのEも同じエネルギースケールを流用
    scaled.wick.T ./= condition.T_scale
    if :swick in keys(scaled)              # 平板型: secondary wick
        scaled.swick.E ./= condition.Em_scale
        scaled.swick.T ./= condition.T_scale
    end
    for comp in condition.wall_comp
        c = getproperty(scaled.T_wall, comp)
        c ./= condition.T_scale
    end
    :plate in keys(scaled) && (scaled.plate.T ./= condition.T_scale)   # 平板型: 凝縮器プレート
    :cplate in keys(scaled) && (scaled.cplate.T ./= condition.T_scale) # 平板型: コールドプレート
    return scaled
end

"""
    restore_state(u_scaled, condition::LHPTransientCondition) -> ComponentArray

scale_state の逆変換。
"""
function restore_state(u_scaled, condition::LHPTransientCondition)
    restored = copy(u_scaled)
    for comp in condition.fluid_comp
        c = getproperty(restored, comp)
        c.ρm .*= condition.ρm_scale
        c.u  .*= condition.u_scale
        c.Em .*= condition.Em_scale
        c.P  .*= condition.P_scale
        c.T  .*= condition.T_scale
        c.x  .*= condition.x_scale
    end
    restored.wick.E .*= condition.Em_scale
    restored.wick.T .*= condition.T_scale
    if :swick in keys(restored)
        restored.swick.E .*= condition.Em_scale
        restored.swick.T .*= condition.T_scale
    end
    for comp in condition.wall_comp
        c = getproperty(restored.T_wall, comp)
        c .*= condition.T_scale
    end
    :plate in keys(restored) && (restored.plate.T .*= condition.T_scale)
    :cplate in keys(restored) && (restored.cplate.T .*= condition.T_scale)
    return restored
end