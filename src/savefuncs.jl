# =============================================================================
# savefuncs.jl (結果の保存)
#
# solve_transient_lhp! の結果 (sol, p) を、CSV に保存する。
# 旧実装は、時間ステップごとに状態を保存していたが、IDA は時間積分を内部で行うので、
# condition.saveat ごとに保存された解 (sol.t, sol.u) を、まとめて書き出す。
#
#   save_results(sol, p, tables, dir)
#     dir/summary.csv … 時系列の要約 (1 保存点 = 1 行)
#     dir/profiles.csv … 空間分布 (縦持ち: 時刻 × 成分 × セル)
#     dir/wick.csv     … ウィックの温度・エネルギーの分布 (時刻 × セル)
#     dir/geometry.csv … 各成分のセル数 N と長さ L (プロットで、位置座標を作るため)
#
# 保存した CSV は、計算をやり直さずに、plot.jl (単体で動く) で、いつでも眺められる。
#
# profiles.csv は縦持ち (long format)。成分ごとにセル数が違うので、横持ちにすると
# 列が不揃いになる。DataFrames で、成分・時刻を絞り込んで使う:
#   using CSV, DataFrames
#   prof = CSV.read("profiles.csv", DataFrame)
#   vl = prof[(prof.comp .== "vl") .& (prof.t .== 10.0), :]
# =============================================================================

using CSV, DataFrames

"""
    results_summary(sol, p, tables) -> DataFrame

時系列の要約。1 保存点 = 1 行。
  t [s], Q_apply [W], T_ec [K] (蒸発器ケース), P_cc [Pa], T_cc [K], P_grv [Pa] (グルーブの平均圧力),
  T_out [K] (蒸発界面の温度 = T_sat(P_grv)), mdot [kg/s] (CC からウィックへの流量 = 蒸発量),
  Q_evap [W] (蒸発による熱 = mdot × h_fg(T_out)。印加熱 Q_apply と比べると、熱が蒸発に使われているか分かる)
"""
function results_summary(sol, p::AbstractLHPParams, tables::FluidTables)
    df = DataFrame(t = Float64[], Q_apply = Float64[], T_ec = Float64[], P_cc = Float64[],
                   T_cc = Float64[], P_grv = Float64[], T_out = Float64[], mdot = Float64[],
                   Q_evap = Float64[])
    for k in eachindex(sol.t)
        c = components_at(sol, p, k)
        P_grv = groove_pressure(c.grv)
        ρ_l = tables.ρ_l_sat(c.cc.T[1])
        T_out = Tsat(P_grv, tables)
        mdot = c.cc.u[1] * ρ_l * p.p_cc.A
        h_fg = tables.h_v_sat(T_out) - tables.h_l_sat(T_out)
        push!(df, (t = sol.t[k], Q_apply = applied_heat(sol.t[k], p.cond),
                   T_ec = c.T_wall.evap[1], P_cc = c.cc.P[1], T_cc = c.cc.T[1],
                   P_grv = P_grv, T_out = T_out, mdot = mdot, Q_evap = mdot * h_fg))
    end
    return df
end

"""
    results_profiles(sol, p) -> DataFrame

空間分布 (縦持ち)。1 行 = (保存時刻, 成分, セル)。
  comp: "grv", "vl", "cond", "ll", "cc" (流体の成分)
  cell: 成分内のセル番号 (上流 = 1)
  rho [kg/m³], u [m/s] (セルの下流側の面の速度), Em [J/kg], P [Pa], T [K], x [-] (乾き度),
  T_wall [K] (壁温度。グルーブは壁を持たないので NaN)
"""
function results_profiles(sol, p::AbstractLHPParams)
    df = DataFrame(t = Float64[], comp = String[], cell = Int[], rho = Float64[], u = Float64[],
                   Em = Float64[], P = Float64[], T = Float64[], x = Float64[], T_wall = Float64[])
    for k in eachindex(sol.t)
        c = components_at(sol, p, k)
        for comp in (:grv, :vl, :cond, :ll, :cc)
            cm = getproperty(c, comp)
            n = length(cm.ρm)
            Tw = comp === :grv ? fill(NaN, n) : getproperty(c.T_wall, comp)
            for i in 1:n
                push!(df, (t = sol.t[k], comp = String(comp), cell = i, rho = cm.ρm[i], u = cm.u[i],
                           Em = cm.Em[i], P = cm.P[i], T = cm.T[i], x = cm.x[i], T_wall = Tw[i]))
            end
        end
    end
    return df
end

"""
    results_wick(sol, p) -> DataFrame

ウィックの分布。1 行 = (保存時刻, セル)。cell は内側 (CC側) = 1。平板型は、secondary wick を cell = 0 で含む。T [K], E [J] (式49)。
"""
function results_wick(sol, p::AbstractLHPParams)
    df = DataFrame(t = Float64[], cell = Int[], T = Float64[], E = Float64[])
    for k in eachindex(sol.t)
        c = components_at(sol, p, k)
        # 平板型: secondary wick は cell = 0 の行として保存する (primary wick の内側のノード0)
        :swick in keys(c) && push!(df, (t = sol.t[k], cell = 0, T = c.swick.T[1], E = c.swick.E[1]))
        for i in eachindex(c.wick.T)
            push!(df, (t = sol.t[k], cell = i, T = c.wick.T[i], E = c.wick.E[i]))
        end
    end
    return df
end

"""
    results_geometry(p) -> DataFrame

各成分のセル数 N と長さ L [m] (comp = "grv", "vl", "cond", "ll", "cc")。
プロットで、ループに沿った位置座標を作るのに使う (plot.jl の load_lengths が読む)。
"""
function results_geometry(p::AbstractLHPParams)
    g = p.geo
    return DataFrame(comp = ["grv", "vl", "cond", "ll", "cc"],
                     N = [g.grv.N, g.vl.N, g.cond.N, g.ll.N, g.cc.N],
                     L = [g.grv.L, g.vl.L, g.cond.L, g.ll.L, g.cc.L])
end

"""
    save_results(sol, p, tables, dir) -> dir

summary.csv, profiles.csv, wick.csv, geometry.csv を dir に書き出す (dir が無ければ作る)。
"""
function save_results(sol, p::AbstractLHPParams, tables::FluidTables, dir::AbstractString)
    mkpath(dir)
    CSV.write(joinpath(dir, "summary.csv"), results_summary(sol, p, tables))
    CSV.write(joinpath(dir, "profiles.csv"), results_profiles(sol, p))
    CSV.write(joinpath(dir, "wick.csv"), results_wick(sol, p))
    CSV.write(joinpath(dir, "geometry.csv"), results_geometry(p))
    println("結果を保存しました: ", abspath(dir), " (summary.csv, profiles.csv, wick.csv, geometry.csv)")
    return dir
end