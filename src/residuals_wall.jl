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
function Q_amb_pipe(Tw, wp, T_amb; amb = nothing)
    if amb !== nothing
        if amb.mode === :adiabatic
            return 0.0 * Tw                      # 周囲と熱交換しない
        elseif amb.mode === :insulated
            # 断熱材 (内径 D_in、外径 D_out、熱伝導率 k) の穴に管が通っている: 管 → 空気層 (管と穴の間、伝導のみ)
            # → 断熱材 → 外表面の対流 (h_o 一定) の直列抵抗 (長さあたり)。
            #   R' = ln(r_b/r1)/(2π k_air) + ln(r_o/r_b)/(2π k) + 1/(h_o 2π r_o)    、 Q = (T_amb − Tw) L_cell / R'
            # 管が細い (r1 < k/h_o = 約 4 mm) ので、断熱材を巻いても熱交換はあまり減らない。
            r1 = wp.D_o / 2
            r_b, r_o = amb.D_in / 2, amb.D_out / 2
            R_len = log(r_b / r1) / (2π * amb.k_air) + log(r_o / r_b) / (2π * amb.k) +
                    1 / (amb.h_o * 2π * r_o)
            L_cell = wp.A_out / (π * wp.D_o)
            return (L_cell / R_len) * (T_amb - Tw)
        end
    end
    h = h_ambient_smooth(Tw, T_amb, wp.D_o)
    return h * wp.A_out * (T_amb - Tw)
end

"""
    ll_env_temps(cnd, N, T_hi3) -> Vector

液管 (N セル) の周囲温度。CC 側の末端 n_ll_ci セルは CC 側の断熱材の節点の温度 T_hi3、他は T_amb。
"""
function ll_env_temps(cnd, N::Int, T_hi3)
    T = fill(cnd.T_amb, N)
    if cnd.hi_on && cnd.n_ll_ci > 0
        for k in 0:min(cnd.n_ll_ci, N) - 1
            T[N - k] = T_hi3
        end
    end
    return T
end

"""
    thermal_bridge_ll(c, T_src, Tll) -> Vector

熱橋: 熱源側 (HB または蒸発器ケース、温度 T_src) から、液管の CC 側の末端 n_br セルの壁へ入る熱 [W] (セルごと、液管が受け取る向きが正)。
    Q_i = G_br w_i (T_src − Tll_i),   w_i ∝ exp(−k/br_decay) (k = CC 側の末端からのセル数 0,1,2,…、Σ w_i = 1)
CC 側の末端から n_br セルに配分する。br_decay = Inf なら一様。G_br = 0 (既定) なら 0。
熱源側の熱収支からは、この合計を引く。
液は CC へ向かって流れるので、センサ (LL3 は 5/6 の位置) より上流のセルに熱が入らないと、センサの温度は変わらない。
"""
function thermal_bridge_ll(c, T_src, Tll)
    Q = zeros(length(Tll))
    (c.G_br > 0 && c.n_br > 0) || return Q
    N = length(Tll)
    n = min(c.n_br, N)
    w = [isfinite(c.br_decay) ? exp(-k / c.br_decay) : 1.0 for k in 0:n-1]
    w ./= sum(w)
    for k in 0:n-1
        i = N - k
        Q[i] = c.G_br * w[k+1] * (T_src - Tll[i])
    end
    return Q
end

"""
    pipe_amb_spec(c, line::Symbol) -> NamedTuple or nothing

condition の amb_mode (:bare | :insulated | :adiabatic) から、蒸気管 (:vl)・液管 (:ll) の周囲との熱交換の指定を作る。
:bare なら nothing (従来の自然対流)。
"""
function pipe_amb_spec(c, line::Symbol)
    hasproperty(c, :amb_mode) || return nothing
    c.amb_mode === :bare && return nothing
    return (mode = c.amb_mode, k = c.ins_k, k_air = c.ins_k_air, D_in = c.ins_D_in, D_out = c.ins_D_out, h_o = c.ins_h_o)
end

"""
    Q_amb_condenser(Tw, wp, T_cond; h_con=500.0) -> Q [W]

凝縮器の壁が、凝縮器プレート (温度 T_cond) と接触して交換する熱 (式51)。
接触面積は外壁面積の f_con 倍 (既定 0.5 = 外周の半分)。壁が受け取る向き (T_cond − Tw)。
"""
Q_amb_condenser(Tw, wp, T_cond; h_con=500.0, f_con=0.5, wall_R=false) =
    wall_R ? wp.A_out / (1 / (h_con * f_con) + wall_R_out(Tw, wp)) * (T_cond - Tw) :
             h_con * (f_con * wp.A_out) * (T_cond - Tw)

"""
    wall_R_in(T, wp), wall_R_out(T, wp) -> R'' [m²K/W]

管壁の半径方向の熱伝導抵抗 (SUS304、温度依存の k_sus304)。壁の節点は管壁の中心 (半径 r_m = (r_i + r_o)/2) の温度とみなし、
  内面〜中心: R_in  = r_i ln(r_m/r_i)/k  (内面積基準。流体−壁の熱伝達 h に直列に入れる)
  中心〜外面: R_out = r_o ln(r_o/r_m)/k  (外面積基準。壁−プレートの接触 h_con に直列に入れる)
wall_radial_R = true のときだけ使う (false なら壁は半径方向に等温 = 従来)。フィッティング係数なし。
"""
wall_R_in(T, wp)  = (wp.D_i / 2) * log(((wp.D_i + wp.D_o) / 4) / (wp.D_i / 2)) / k_sus304(T)
wall_R_out(T, wp) = (wp.D_o / 2) * log((wp.D_o / 2) / ((wp.D_i + wp.D_o) / 4)) / k_sus304(T)

"""
    h_series(h, R) -> h/(1 + h R)

熱伝達 h と、直列の抵抗 R (面積基準 [m²K/W]) の合成 (h = 0 でも 0 のまま)。
"""
h_series(h, R) = h / (1 + h * R)

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
    h_con=500.0, f_con=0.5, amb=nothing, Q_extra=nothing, wall_R=false)

    N = wp.N
    for i in 1:N
        T = Tw[i]

        # --- 流体から受け取る熱 (式31) ---
        Q_fluid = h_conv[i] * wp.A_in * (T_fluid[i] - T)

        # --- 周囲・凝縮器プレートとの熱 (式30, 51) ---
        Q_amb = kind === :condenser ? Q_amb_condenser(T, wp, T_env isa AbstractVector ? T_env[i] : T_env; h_con=h_con, f_con=f_con, wall_R=wall_R) :
                                      Q_amb_pipe(T, wp, T_env isa AbstractVector ? T_env[i] : T_env; amb = amb)

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

        Q_x = Q_extra === nothing ? 0.0 : Q_extra[i]
        resid_Tw[i] = wall_heat_capacity(T, wp) * du_Tw[i] - (Q_amb + Q_fluid + Q_cond + Q_x)
    end
    return nothing
end


# -----------------------------------------------------------------------------
# 5. 凝縮器プレート (A5052、2次元)
# -----------------------------------------------------------------------------

"""
    plate_cell_params(c, pg) -> (C, Gx, Gy, Gc)

プレートの 1 節点 (dx × dy × t): 熱容量 C [J/K]、x・y 方向の隣との熱コンダクタンス Gx, Gy [W/K]、
コールドプレートとの接触コンダクタンス Gc [W/K] (h_cp dx dy)。
"""
function plate_cell_params(c, pg)
    dx = pg.L / pg.nx
    dy = pg.W / pg.ny
    return (C = c.plate_ρ * c.plate_c * dx * dy * pg.t,
            Gx = c.plate_k * pg.t * dy / dx,
            Gy = c.plate_k * pg.t * dx / dy,
            Gc = c.h_cp * dx * dy)
end

"""
    plate_tube_nodes(pg, L_tube, N) -> Vector{Int}

管の各セル (N 個、管の中心線上の弧長 s) が、プレートのどの節点 (k = (j−1) nx + i) に接するか。
経路: x = 0 の短辺から y1 = (W − pitch)/2 の線を x 方向へ進み、半径 pitch/2 の半円で折り返し、
y2 = (W + pitch)/2 の線を x = 0 まで戻る。直線部の長さは L_tube = 2 L_leg + π pitch/2 から決める。
"""
function plate_tube_nodes(pg, L_tube, N)
    R = pg.pitch / 2
    L_leg = (L_tube - π * R) / 2
    y1 = (pg.W - pg.pitch) / 2
    y2 = pg.W - y1
    dx = pg.L / pg.nx
    dy = pg.W / pg.ny
    nodes = zeros(Int, N)
    for i in 1:N
        s = (i - 0.5) * L_tube / N
        if s <= L_leg
            x, y = s, y1
        elseif s <= L_leg + π * R
            θ = (s - L_leg) / R
            x, y = L_leg + R * sin(θ), pg.W / 2 - R * cos(θ)
        else
            x, y = L_leg - (s - L_leg - π * R), y2
        end
        ii = clamp(floor(Int, x / dx + 1e-9) + 1, 1, pg.nx)
        jj = clamp(floor(Int, y / dy + 1e-9) + 1, 1, pg.ny)
        nodes[i] = (jj - 1) * pg.nx + ii
    end
    return nodes
end

"""
    residual_plate!(resid, dTp, Tp, Tw, wp, c, pg, nodes)

凝縮器プレートの温度 Tp (nx·ny 要素) の残差。
    C dTp_k/dt = Σ_隣 G (Tp_隣 − Tp_k) + Gc (T_cond − Tp_k) + Σ_{管セル i → k} h_con (A_out/2)(Tw_i − Tp_k)
外周は断熱。管壁との熱は、壁の残差 (Q_amb_condenser) と同じ式で符号だけ逆。
"""
function residual_plate!(resid, dTp, Tp, Tw, wp, c, pg, nodes; Tcp = nothing)
    pp = plate_cell_params(c, pg)
    nx, ny = pg.nx, pg.ny
    for j in 1:ny, i in 1:nx
        k = (j - 1) * nx + i
        T = Tp[k]
        Q = pp.Gc * ((Tcp === nothing ? c.T_cond : Tcp[k]) - T)    # コールドプレートが2次元ならその節点の温度、なければ T_cond 一定
        i > 1  && (Q += pp.Gx * (Tp[k-1] - T))
        i < nx && (Q += pp.Gx * (Tp[k+1] - T))
        j > 1  && (Q += pp.Gy * (Tp[k-nx] - T))
        j < ny && (Q += pp.Gy * (Tp[k+nx] - T))
        resid[k] = pp.C * dTp[k] - Q
    end
    for (ic, k) in enumerate(nodes)
        resid[k] += Q_amb_condenser(Tw[ic], wp, Tp[k]; h_con = c.h_con, f_con = c.f_con_area, wall_R = c.wall_radial_R)    # resid = C dT − Q_in で、プレートへの Q_in = −(壁が受ける熱)
    end
    return nothing
end

"""
    cplate_cell_params(c, pg, cg) -> (C, Gx, Gy, Gc)

コールドプレートの 1 節点 (dx × dy × cg.t): 熱容量 C、x・y 方向の隣との熱コンダクタンス、
凝縮器プレートとの接触コンダクタンス Gc = h_cp dx dy (グリス+ボルト締結)。材質は凝縮器プレートと同じ A5052。
"""
function cplate_cell_params(c, pg, cg)
    dx = pg.L / pg.nx
    dy = pg.W / pg.ny
    return (C = c.plate_ρ * c.plate_c * dx * dy * cg.t,
            Gx = c.plate_k * cg.t * dy / dx,
            Gy = c.plate_k * cg.t * dx / dy,
            Gc = c.h_cp * dx * dy)
end

"""
    cplate_pipe_nodes(pg, cg) -> Vector{Int}

チラー配管の各セルが接するコールドプレートの節点。LHP の管 (x = 0 の短辺から入る) と同じ経路を、
x = L の短辺から入るように左右反転したもの (LHP の管と対向流)。
"""
function cplate_pipe_nodes(pg, cg)
    nodes = plate_tube_nodes(pg, cg.L_tube, cg.N)
    return [((k - 1) ÷ pg.nx) * pg.nx + (pg.nx - (k - 1) % pg.nx) for k in nodes]
end

"チラー配管 1 セルとコールドプレートの熱コンダクタンス [W/K] = h_con f_con π D_o dl (凝縮器と同じ式・同じ値)。"
cplate_pipe_G(c, cg) = (c.h_chiller > 0 ? c.h_chiller : c.h_con) * c.f_con_area * π * cg.D_o * (cg.L_tube / cg.N)

"""
    residual_cplate!(resid, dTcp, Tcp, Tp, c, pg, cg, nodes)

コールドプレートの温度 Tcp の残差。 C dTcp_k/dt = Σ_隣 G (Tcp_隣 − Tcp_k) + Gc (Tp_k − Tcp_k) + Σ_{配管セル → k} G_pipe (T_cond − Tcp_k)。
外周は断熱。T_cond は冷却液 (チラー) の温度 (一定)。
"""
function residual_cplate!(resid, dTcp, Tcp, Tp, c, pg, cg, nodes)
    pp = cplate_cell_params(c, pg, cg)
    nx, ny = pg.nx, pg.ny
    for j in 1:ny, i in 1:nx
        k = (j - 1) * nx + i
        T = Tcp[k]
        Q = pp.Gc * (Tp[k] - T)
        i > 1  && (Q += pp.Gx * (Tcp[k-1] - T))
        i < nx && (Q += pp.Gx * (Tcp[k+1] - T))
        j > 1  && (Q += pp.Gy * (Tcp[k-nx] - T))
        j < ny && (Q += pp.Gy * (Tcp[k+nx] - T))
        resid[k] = pp.C * dTcp[k] - Q
    end
    Gp = cplate_pipe_G(c, cg)
    for k in nodes
        resid[k] -= Gp * (c.T_cond - Tcp[k])
    end
    return nothing
end