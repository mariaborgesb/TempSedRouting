def critical_strmpow_q_Q(fromArrayList, SF, kv, rho, Rrho, g, Dm_dict, di, mcoef, slope_dict, Qinterp_dict, qinterp_dict):
    #Calculation of critical stream power, critical flow and critical unitflow
    #using logarithmic flow resistance law.

    import numpy as np
    import scipy.stats

    # Outputs
    interceptregressi_dict = {}  # intercept per node 
    sloperegressi_dict = {}      # slope per node
    qcr_dict = {}                # critical unit discharge per node
    Qcr_dict = {}                # critical discharge per node 
    Strmpowcri_dict = {}         # critical streampower per node
    Strmpow_dict = {}            # streampower per node
    widthcr_dict = {}            # critical width per node

    # Small constants to avoid log/divide problems
    EPS_SLOPE = 1e-8   # minimum slope used in formulas
    EPS_Q     = 1e-12  # minimum q considered "nonzero"
    EPS_LOG   = 1e-12  # minimum positive value for log10 argument

    qcri = None

    for n in fromArrayList:
        #Critical values per fraction
        tetarm = 0.021 + (0.015 * np.exp(-20.0 * SF))
        slope_eff = max(float(slope_dict[n]), EPS_SLOPE)
        Dm = float(Dm_dict[n])

        taurm = tetarm * rho * g * Rrho * Dm
        bfunci = 0.67 / (1.0 + np.exp(1.5 - (di / Dm)))       
        tetari = taurm * ((di / Dm) ** bfunci) / (rho * Rrho * g * di)

        log_arg = (30.0 * tetari * Rrho * di) / (np.e * mcoef * slope_eff * Dm)
        log_arg = np.maximum(log_arg, EPS_LOG)
        Logi = np.log10(log_arg)

        Wcri = (2.3 / kv) * rho * (tetari * Rrho * g * di) ** 1.5 * Logi     
        qcri = Wcri / (rho * g * slope_eff)   

        n_space = len(Qinterp_dict[n])
        n_frac  = len(qcri)
        Qcr         = np.zeros((n_space, n_frac), dtype=float)
        widthcr     = np.zeros((n_space, n_frac), dtype=float)
        Strmpowcri  = np.zeros((n_space, n_frac), dtype=float)
        sloperegr   = np.zeros(n_space, dtype=float)
        interceptr  = np.zeros(n_space, dtype=float)

        #Stream power
        Qtemp = np.array(Qinterp_dict[n])            
        Strmpow = rho * g * Qtemp * slope_eff        
        Strmpow_dict[n] = Strmpow

        for l in range(n_space):
            q_series = np.asarray(qinterp_dict[n][l], dtype=float)
            Q_series = np.asarray(Qinterp_dict[n][l], dtype=float)

            mask = (q_series > EPS_Q) & np.isfinite(q_series) & np.isfinite(Q_series)

            if np.any(mask):
                widths = np.zeros_like(Q_series, dtype=float)
                np.divide(Q_series, q_series, out=widths, where=mask)

                if np.count_nonzero(mask) >= 2:
                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(
                        q_series[mask], widths[mask]
                    )
                else:
                    slope, intercept = 0.0, float(np.nanmedian(widths[mask]))
                    if not np.isfinite(intercept):
                        intercept = 0.0
            else:
                slope, intercept = 0.0, 0.0

            sloperegr[l]  = slope
            interceptr[l] = intercept

            width_l = slope * qcri + intercept
            width_l = np.maximum(width_l, 0.0)

            widthcr[l, :]    = width_l
            Qcr[l, :]        = width_l * qcri
            Strmpowcri[l, :] = Wcri * width_l

        # Collect per-node results
        interceptregressi_dict[n] = interceptr
        sloperegressi_dict[n]     = sloperegr
        Strmpowcri_dict[n]        = Strmpowcri
        Qcr_dict[n]               = Qcr
        qcr_dict[n]               = qcri
        widthcr_dict[n]           = widthcr

    return Strmpowcri_dict, Qcr_dict, qcr_dict, widthcr_dict, Strmpow_dict, qcri
