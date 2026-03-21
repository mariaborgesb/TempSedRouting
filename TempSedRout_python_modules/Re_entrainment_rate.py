def reentrainment (fromArrayList, qcri, Qinterp_dict, Fi_dict, Strmpow_dict,
                   cri_all_dict, Strmpowcri_dict, depthinterp_dict, rhos, rho, g):
    #Calculation of the rate of re-entrainment of particles of size d
    #using equation developed by Haddadchi and Rose (2022).

    import numpy as np

    ki = 0.3  # the fraction of depth through which sediment is lifted [-]
    EPS_H = 1e-12  # small depth to avoid division by zero

    ryd_dict = {}
    for n in fromArrayList:  # loop for each reachID
        rydtemp_dict = {}
        for k in range(len(qcri)):  # loop for each fraction
            ryd_temp = np.zeros_like(Qinterp_dict[n])
            for l in range(len(Qinterp_dict[n])):  # loop for each interpolated data array
                # denominator: ki * g * h
                denom = ki * g * np.maximum(depthinterp_dict[n][l, :], EPS_H)
                # excess stream power term
                excess = Strmpow_dict[n][l, :] - (cri_all_dict[n][0, k] * Strmpowcri_dict[n][l, k])
                ryd_temp[l, :] = (Fi_dict[n][0, k] * excess / denom) * (rhos / (rhos - rho))

            # change negative ryd to zero (-ryd=0 --> E<0 --> scenario 2, only deposition since critical stream power > stream power)
            ryd_temp[ryd_temp < 0] = 0.0

            rydtemp_dict[k] = ryd_temp
            del ryd_temp

        ryd_dict[n] = rydtemp_dict
        del rydtemp_dict

    return ryd_dict
