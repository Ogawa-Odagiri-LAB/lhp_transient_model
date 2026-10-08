# =============================================================================
# residuals_transport.jl (輸送管・グルーブ・CC 共通の流体残差)
#
# 設計ドキュメント 式14-19, 33-38 に対応する、IDA (完全陰的DAE) 用の
# 残差関数。蒸気管・凝縮器・液管の流体側に使う。グルーブとCCも、同じ
# 6変数 (ρm,u,Em,P,T,x) 構成で、熱の入り方と質量・エネルギーのソースだけ
# が違うので、核となる residual_transport! を共用する
# (熱流束 q とソース S_m, S_E は引数として外から与える)。
#
# 【状態変数の配置 (スタッガード格子)】
#   セル中心:  ρm[i], Em[i], P[i], T[i], x[i]   (i = 1..N)
#   セル境界:  u[i] = u_{i+1/2}  (セル i の「下流側の面」の速度)
#   セル i の上流面の速度は u[i-1]、セル1の上流面 (=入口面) の速度は
#   上流コンポーネントの最終面の速度 (ゴースト値 up.u) で与える。
#
# 【隣接コンポーネントとの接続 (ゴースト値)】
#   up   = (ρ, Em, P, u, dl)  上流コンポーネントの最終セルの値。
#                              u は入口面の速度で、この成分の断面積に換算済み
#                              (質量流量 ρuA を保つため u_up * A_up / A)。
#   down = (ρ, Em, P, u, dl)  下流コンポーネントの最初のセルの値。
#                              u はその下流側の面の速度 (u_{N+3/2}) で、同様に
#                              この成分の断面積に換算済み。
#   dl は隣のセルの長さ。セル中心間距離は (dl + dl_neighbor)/2。
#   up/down を作る補助として transport_ghost_up / transport_ghost_down を用意。
#
# 【風上差分】 ドキュメントのドナーセル法 ⟨φ⟩ は、すべて smooth_step による
# tanh ブレンドに置き換えている (smoothing.jl)。u≫ε_u で厳密に風上差分に一致
# することを検算済み。
#
# 【要確認事項 (ドキュメントとの照合で見つかったもの)】
#   1. 式34の圧力勾配 −(⟨P⟩_{i+1/2} − ⟨P⟩_{i−1/2})/Δz
#      ⟨P⟩ に風上ブレンドを適用すると、u>0 では (P_i − P_{i-1})/Δz
#      (面 i−1/2 を中心とする差分)、u<0 では (P_{i+1} − P_i)/Δz (面 i+1/2
#      を中心とする差分) となり、流れの向きで差分の中心位置が変わってしまう。
#      圧力は移流される量ではなく、スタッガード格子では面を挟む2セルの
#      差分 (P_{i+1} − P_i)/Δz が自然なため、本実装は後者を採用している。
#      ※ ドキュメントの式34は更新されておらず、⟨P⟩ のままである。
#         ⟨P⟩ を風上ブレンドすると、差分の中心位置が流れの向きで変わり
#         (u>0: 面 i−1/2、u=0: セル i、u<0: 面 i+1/2)、u=0 では
#         (P_{i+1}−P_{i−1})/2Δz という2セル飛びの差分になって、市松模様の
#         圧力を検知できなくなる。⟨P⟩ を面の単純平均と読んでも同じ中心差分に
#         なる。いずれの読み方でも節点 (面 i+1/2) からずれるため、
#         圧力勾配は面を挟む2セルの差を採用している。圧力は移流される量では
#         なく、流れの向きによらず両側の差が力として働くため、対流項と違って
#         風上化する物理的な根拠がない。ドキュメントの意図と異なる場合は要相談。
#   2. 式25 の次元 → ドキュメント更新で ΔPn = 4τw/D [Pa/m] に修正済み。
#      本実装 (pressure_drop_single_phase) もドキュメントと同じ形。
#   3. 式30 Qamb の符号 → ドキュメント更新で修正済み。
#   4. 摩擦の向き (式24 の sign(u)) → ドキュメント更新で反映済み。
#      滑らかな符号関数 u/√(u²+δ²) を使っている。
#   5. 摩擦は式34の Fw,i+1/2 に従い、面で評価している。
# =============================================================================

# -----------------------------------------------------------------------------
# 1. 幾何情報の取り出し
# -----------------------------------------------------------------------------

"""
    pipe_params(geo) -> NamedTuple

TLGeometry / CondGeometry / CCGeometry から、残差関数が使う幾何量を
取り出す。D は熱伝達・摩擦の代表直径 (内径)、S は1セルあたりの内壁面積、
V は1セルの流体体積、A は流路断面積。
"""
pipe_params(geo::Union{TLGeometry,CondGeometry,CCGeometry}) = (
    N = geo.N, dl = geo.dl, D = geo.D_i, A = geo.A_cross,
    S = geo.S_inner, V = geo.V, θ = geo.θ, elev0 = geo.elev0,
)

"""
    cell_elevation(p, i) -> 標高 [m]

セル i (1始まり) の中心の標高。直管を前提とし、入口標高 elev0 から
傾斜角 θ [度] で上がるとして z sinθ を計算する。
"""
cell_elevation(p, i) = p.elev0 + (i - 0.5) * p.dl * sind(p.θ)

# -----------------------------------------------------------------------------
# 2. 隣接コンポーネントとのゴースト値
# -----------------------------------------------------------------------------

"""
    transport_ghost_up(u_up, p_up, p) -> up

上流コンポーネント (状態 u_up, 幾何 p_up) の最終セルから、この成分
(幾何 p) の入口側ゴースト値を作る。入口面の速度は、質量流量 ρuA を
保つよう、この成分の断面積に換算する。
"""
transport_ghost_up(u_up, p_up, p) = (
    ρ = u_up.ρm[end], Em = u_up.Em[end], P = u_up.P[end],
    u = u_up.u[end] * p_up.A / p.A, dl = p_up.dl,
)

"""
    transport_ghost_down(u_dn, p_dn, p) -> down

下流コンポーネントの最初のセルから、この成分の出口側ゴースト値を作る。
u は下流成分の最初の面 (u_{N+3/2}) の速度を、この成分の断面積に換算した値。
"""
transport_ghost_down(u_dn, p_dn, p) = (
    ρ = u_dn.ρm[1], Em = u_dn.Em[1], P = u_dn.P[1],
    u = u_dn.u[1] * p_dn.A / p.A, dl = p_dn.dl,
)

# -----------------------------------------------------------------------------
# 2.5 面の風上ブレンド (コンポーネント境界でも共有する)
# -----------------------------------------------------------------------------

"""
    face_upwind(uf, qL, qR, ε_u) -> (q_face, w)

面の速度 uf に対する風上ブレンド。左セルの値 qL と右セルの値 qR から、
面の値 q_face = w qL + (1-w) qR を返す (w = smooth_step(uf, ε_u))。

residual_transport! の面の質量流束や、CCのように別の残差関数が
コンポーネント境界の面を扱うとき (residuals_cc.jl) に、同じ式・同じ
ε_u を使うことで、境界での質量の出入りが厳密に一致する。
"""
function face_upwind(uf, qL, qR, ε_u)
    w = smooth_step(uf, ε_u)
    return w * qL + (1 - w) * qR, w
end

# -----------------------------------------------------------------------------
# 3. セル中心の速度
# -----------------------------------------------------------------------------

"""
    cell_velocity(u, up, i) -> u_c

セル i の中心の速度 = 上流面と下流面の速度の平均。
"""
cell_velocity(u, up, i) = 0.5 * ((i == 1 ? up.u : u.u[i-1]) + u.u[i])

# -----------------------------------------------------------------------------
# 4. 熱伝達係数と壁からの熱流束
# -----------------------------------------------------------------------------

"""
    compute_h_conv_transport!(h, u, up, Tw, p, g, tables; kwargs...)

各セルの流体-壁間熱伝達係数 h_conv[i] [W/m²K] を計算して h に書き込む。

h_conv は、流体側の熱流束 q と、壁側の Q_fluid (residuals_wall.jl) の
両方で使われる。両者で別々に計算するとエネルギー収支が微妙にずれる
ため、全体の残差関数 (residuals_main.jl) で1回だけ計算して共有する。
"""
function compute_h_conv_transport!(h, u, up, Tw, p, g, tables::FluidTables;
    δ_u=1e-4, kwargs...)
    for i in 1:p.N
        u_c = cell_velocity(u, up, i)
        G = u.ρm[i] * soft_abs(u_c, δ_u)
        h[i] = h_transport_wall(u.ρm[i], u_c, u.Em[i], cell_elevation(p, i),
                                Tw[i], p.D, G, g, tables; kwargs...)
    end
    return h
end

"""
    transport_wall_heat_source(u, Tw, h, p) -> q [W/m³]

壁から流体への体積あたりの熱流束 (式29: q = h_conv A_inner (Twall − T) を
セル体積 V で割ったもの)。式3の q は体積あたりの量なので V で割る。
T は状態変数 u.T をそのまま使う (壁側の Q_fluid = h A (T − Twall) と
厳密に符号反転した値になる)。
"""
transport_wall_heat_source(u, Tw, h, p) =
    [h[i] * p.S / p.V * (Tw[i] - u.T[i]) for i in 1:p.N]

# -----------------------------------------------------------------------------
# 5. 流体の残差 (質量・運動量・エネルギー・状態拘束)
# -----------------------------------------------------------------------------

"""
    residual_transport!(resid, du, u, up, down, q, p, g, tables;
                        S_m=nothing, S_E=nothing, ε_u=1e-4, δ_u=1e-4, kwargs...)

1つの輸送成分 (N セル) の流体残差を resid に書き込む。
resid, du, u は、その成分の ComponentVector (ρm,u,Em,P,T,x) を想定。

  質量:     du.ρm + (F_{i+1/2} − F_{i−1/2})/Δz = S_m           (式33)
  運動量:   ⟨ρ⟩(du.u + u ∂u/∂z) + ∂P/∂z + ⟨ρ⟩ g sinθ + Fw = 0  (式34)
  エネルギー: ρ du.Em + ρ u ∂φ/∂z − (P/ρ) du.ρm = q + S_E − φ S_m (式35)
  状態拘束: P − P(ρm,em) = 0, T − T(ρm,em) = 0, x − x(ρm,em) = 0 (式36-38)

エネルギー式の S_m, S_E は、連続の式に質量ソース S_m がある場合の
非保存形 (導出は sympy で検算済み):
    ∂t(ρE) + ∂z(ρuφ) = q + S_E  を連続の式 ∂tρ + ∂z(ρu) = S_m で整理すると
    ρ ∂tE + ρu ∂zφ − (P/ρ)∂tρ = q + S_E − φ S_m
S_E は注入される質量が持ち込むエンタルピーも含めたエネルギーソース。
ソースが無ければ S_m=S_E=0 で式16そのもの (nothing を渡せばよい)。

q [W/m³] は、輸送管では transport_wall_heat_source、CCでは気液別の
自然対流、グルーブでは蒸発器ケースからの過熱、といった具合に、
呼び出し側が与える。

kwargs は friction_force_transport へ渡される (ε_Re, ε_phase)。
"""
function residual_transport!(resid, du, u, up, down, q, p, g, tables::FluidTables;
    S_m=nothing, S_E=nothing, ε_u=1e-4, δ_u=1e-4, kwargs...)

    N, dl = p.N, p.dl
    sinθ = sind(p.θ)

    # --- 面 j=0..N の風上ブレンド (配列の添字は j+1) ---
    # 面 j=0 は入口面、面 j (1..N) はセル j の下流面 (速度 u.u[j])
    F  = Vector{Float64}(undef, N + 1)   # 質量流束 ρ_face * u [kg/m²/s]
    ρf = Vector{Float64}(undef, N + 1)   # 面の風上密度 ⟨ρm⟩
    wf = Vector{Float64}(undef, N + 1)   # 風上ブレンド重み
    for j in 0:N
        uf = j == 0 ? up.u : u.u[j]
        ρL = j == 0 ? up.ρ : u.ρm[j]
        ρR = j == N ? down.ρ : u.ρm[j+1]
        ρf[j+1], wf[j+1] = face_upwind(uf, ρL, ρR, ε_u)
        F[j+1] = ρf[j+1] * uf
    end

    for i in 1:N
        # ---------------- 質量 (式33) ----------------
        Sm_i = S_m === nothing ? 0.0 : S_m[i]
        resid.ρm[i] = du.ρm[i] + (F[i+1] - F[i]) / dl - Sm_i

        # ---------------- 運動量 (式34, 面 i = 下流面) ----------------
        u_i    = u.u[i]
        u_prev = i == 1 ? up.u : u.u[i-1]
        u_next = i == N ? down.u : u.u[i+1]
        w      = wf[i+1]
        dl_fwd = i == N ? down.dl : dl
        conv   = u_i * (w * (u_i - u_prev) / dl + (1 - w) * (u_next - u_i) / dl_fwd)

        P_R  = i == N ? down.P : u.P[i+1]
        dz_c = i == N ? 0.5 * (dl + down.dl) : dl
        dPdz = (P_R - u.P[i]) / dz_c

        # 摩擦 Fw,i+1/2 (式34): 面で評価する。状態は ⟨ρ⟩ と同じ風上ブレンド
        # (ρ, Em)、速度は面の速度 u_i。friction_force_transport は式25の
        # ΔPn = 4τw/D [Pa/m] を返すので、そのまま式2の摩擦項に使える。
        # 向きは滑らかな符号関数で与える (式24: Fw = sign(u) ΔPn)。
        # Em_f はブレンドした全エネルギー。em を逆算する際の運動エネルギー・
        # 位置エネルギーのずれ (u²/2 や g Δelev: 運転条件では em の1e-5
        # 程度) は、摩擦係数の物性への影響として無視できる。
        Em_R = i == N ? down.Em : u.Em[i+1]
        Em_f = w * u.Em[i] + (1 - w) * Em_R
        elev_f = p.elev0 + i * dl * sind(p.θ)        # 面 i+1/2 の標高
        G_f = ρf[i+1] * soft_abs(u_i, δ_u)
        Fw_mag = friction_force_transport(ρf[i+1], u_i, Em_f, elev_f,
                                          p.D, G_f, g, tables; kwargs...)
        Fw = (u_i / soft_abs(u_i, δ_u)) * Fw_mag

        u_c = cell_velocity(u, up, i)   # セル中心速度 (エネルギー式・状態拘束で使う)

        resid.u[i] = ρf[i+1] * (du.u[i] + conv) + dPdz + ρf[i+1] * g * sinθ + Fw

        # ---------------- エネルギー (式35, セル i) ----------------
        φ_i = u.Em[i] + u.P[i] / u.ρm[i]
        φ_L = i == 1 ? up.Em + up.P / up.ρ : u.Em[i-1] + u.P[i-1] / u.ρm[i-1]
        φ_R = i == N ? down.Em + down.P / down.ρ : u.Em[i+1] + u.P[i+1] / u.ρm[i+1]
        dzL = i == 1 ? 0.5 * (dl + up.dl) : dl
        dzR = i == N ? 0.5 * (dl + down.dl) : dl
        wc  = smooth_step(u_c, ε_u)
        dφdz = wc * (φ_i - φ_L) / dzL + (1 - wc) * (φ_R - φ_i) / dzR

        Sm_i = S_m === nothing ? 0.0 : S_m[i]
        SE_i = S_E === nothing ? 0.0 : S_E[i]
        resid.Em[i] = u.ρm[i] * du.Em[i] + u.ρm[i] * u_c * dφdz -
                      (u.P[i] / u.ρm[i]) * du.ρm[i] - q[i] - SE_i + φ_i * Sm_i

        # ---------------- 状態拘束 (式36-38) ----------------
        em = internal_energy(u.Em[i], u_c, cell_elevation(p, i), g)
        st = flash_ρe(u.ρm[i], em, tables)
        resid.P[i] = u.P[i] - st.P
        resid.T[i] = u.T[i] - st.T
        resid.x[i] = u.x[i] - st.x
    end
    return nothing
end

# -----------------------------------------------------------------------------
# 6. 輸送管 (蒸気管・凝縮器・液管) 用のまとめ関数
# -----------------------------------------------------------------------------

"""
    residual_transport_line!(resid, du, u, up, down, Tw, h, p, g, tables; kwargs...)

蒸気管・凝縮器・液管の流体残差。壁からの熱流束 q を
transport_wall_heat_source で作り、residual_transport! に渡す。
h は compute_h_conv_transport! で事前に計算した熱伝達係数。
Tw はその成分の壁温度 (u.T_wall.vl など)。
"""
function residual_transport_line!(resid, du, u, up, down, Tw, h, p, g, tables::FluidTables; kwargs...)
    q = transport_wall_heat_source(u, Tw, h, p)
    return residual_transport!(resid, du, u, up, down, q, p, g, tables; kwargs...)
end