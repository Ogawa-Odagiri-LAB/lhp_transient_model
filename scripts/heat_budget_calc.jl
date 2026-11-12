# heat_budget_calc.jl — write the loop heat budget and the condenser two-phase length of a finished run to <result_dir>/budget.csv
#
# Called from transient_flat.jl right after save_results (include + write_heat_budget). Diagnostics only; no physics is changed.
# Columns (all in W unless noted; tau = t - t_heat [s]):
#   Q_apply, Q_plate (condenser plate -> cold plate), Q_amb_vl, Q_amb_ll (vapor/liquid line -> ambient, + = lost),
#   Q_house (insulation -> ambient and PEEK -> base plate), Q_out = Q_plate + Q_amb_vl + Q_amb_ll + Q_house (the sum that leaves the system, as in heat_out)
#   dE_total, dE_HB, dE_EC, dE_CC, dE_wall, dE_fluid, dE_wick, dE_plate (storage rates, central difference), resid = Q_apply - Q_out - dE_total
#   Q_cond_amb_hyp = HYPOTHETICAL ambient gain of the condenser tube if it were wrapped like the vapor line (not in the model; + = heat gain)
#   L_2ph [m] = sum of dl of condenser cells with 0.001 < x < 0.999, L_vap_eq [m] = sum(dl * x), x_front [-] = (first cell with x < 0.5, from the inlet) position / L
#   n_2ph, n_vap = number of cells (two-phase, x > 0.999)
using CSV, DataFrames

function write_heat_budget(sol, p, tables, dir; t_heat::Real = 10.0)
    M = LHPTransientModel
    cnd = p.cond
    n = length(sol.t)
    parts = [M.stored_energy_parts(M.components_at(sol, p, k), p) for k in 1:n]
    Lc, Nc = p.geo.cond.L, p.geo.cond.N
    dlc = Lc / Nc
    amb_vl = M.pipe_amb_spec(cnd, :vl)
    cols = Dict{Symbol,Vector{Float64}}()
    for s in (:tau, :Q_apply, :Q_plate, :Q_amb_vl, :Q_amb_ll, :Q_house, :Q_out, :dE_total, :dE_HB, :dE_EC, :dE_CC, :dE_wall, :dE_fluid, :dE_wick, :dE_plate,
              :resid, :Q_cond_amb_hyp, :L_2ph, :L_vap_eq, :x_front, :n_2ph, :n_vap)
        cols[s] = fill(NaN, n)
    end
    for k in 1:n
        c = M.components_at(sol, p, k)
        T_hi3 = (cnd.hi_on && :hi in keys(c)) ? c.hi.T[3] : cnd.T_amb
        Qvl = 0.0; Qll = 0.0
        for (j, T) in enumerate(c.T_wall.vl);  Qvl -= M.Q_amb_pipe(T, p.wp_vl, cnd.T_amb; amb = amb_vl); end
        envll = M.ll_env_temps(cnd, length(c.T_wall.ll), T_hi3)
        for (j, T) in enumerate(c.T_wall.ll);  Qll -= M.Q_amb_pipe(T, p.wp_ll, envll[j]; amb = M.pipe_amb_spec(cnd, :ll)); end
        Qpl = M.n_cplate(p.geo) > 0 ? sum(M.cplate_pipe_G(cnd, p.geo.cplate) * (c.cplate.T[kk] - cnd.T_cond) for kk in p.cplate_nodes) :
              (:plate in keys(c) ? sum(M.plate_cell_params(cnd, p.geo.plate).Gc * (T - cnd.T_cond) for T in c.plate.T) : NaN)
        Qh = M.housing_heat_out(c, cnd)
        Qhyp = 0.0
        for T in c.T_wall.cond
            Qhyp += M.Q_amb_pipe(T, p.wp_cond, cnd.T_amb; amb = amb_vl)
        end
        # storage rate (central difference)
        k1, k2 = max(1, k - 1), min(n, k + 1)
        dt = sol.t[k2] - sol.t[k1]
        dE(f) = dt > 0 ? (f(parts[k2]) - f(parts[k1])) / dt : NaN
        d_hb = dE(e -> e.HB); d_ec = dE(e -> e.EC_case); d_cc = dE(e -> e.CC_case)
        d_wall = dE(e -> e.wall_vl + e.wall_cond + e.wall_ll)
        d_fl = dE(e -> e.fl_grv + e.fl_vl + e.fl_cond + e.fl_ll + e.fl_cc)
        d_wick = dE(e -> e.wick + e.swick)
        d_pl = dE(e -> e.plate + e.cplate + e.housing)
        dEt = d_hb + d_ec + d_cc + d_wall + d_fl + d_wick + d_pl
        Qa = M.applied_heat(sol.t[k], cnd)
        Qo = Qpl + Qvl + Qll + Qh
        x = c.cond.x
        two = [0.001 < xi < 0.999 for xi in x]
        ifront = findfirst(<(0.5), x)
        cols[:tau][k] = sol.t[k] - t_heat;  cols[:Q_apply][k] = Qa;  cols[:Q_plate][k] = Qpl
        cols[:Q_amb_vl][k] = Qvl;  cols[:Q_amb_ll][k] = Qll;  cols[:Q_house][k] = Qh;  cols[:Q_out][k] = Qo
        cols[:dE_total][k] = dEt;  cols[:dE_HB][k] = d_hb;  cols[:dE_EC][k] = d_ec;  cols[:dE_CC][k] = d_cc
        cols[:dE_wall][k] = d_wall;  cols[:dE_fluid][k] = d_fl;  cols[:dE_wick][k] = d_wick;  cols[:dE_plate][k] = d_pl
        cols[:resid][k] = Qa - Qo - dEt
        cols[:Q_cond_amb_hyp][k] = Qhyp
        cols[:L_2ph][k] = dlc * count(two);  cols[:L_vap_eq][k] = dlc * sum(clamp.(x, 0.0, 1.0))
        cols[:x_front][k] = ifront === nothing ? 1.0 : (ifront - 1) * dlc / Lc
        cols[:n_2ph][k] = count(two);  cols[:n_vap][k] = count(>(0.999), x)
    end
    df = DataFrame(cols)
    select!(df, [:tau, :Q_apply, :Q_plate, :Q_amb_vl, :Q_amb_ll, :Q_house, :Q_out, :dE_total, :dE_HB, :dE_EC, :dE_CC, :dE_wall, :dE_fluid, :dE_wick, :dE_plate,
                 :resid, :Q_cond_amb_hyp, :L_2ph, :L_vap_eq, :x_front, :n_2ph, :n_vap])
    mkpath(dir)
    CSV.write(joinpath(dir, "budget.csv"), df)
    println("熱収支を保存しました: ", abspath(joinpath(dir, "budget.csv")))
    return df
end