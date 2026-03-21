# %%
def disperustar(fromArrayList, depthinterp_dict, slope_dict, g, Dispersioncoeff_dict):
    # calculating dispersion term in mass conservation using the approach
    # proposed by Fischer et al. (1979)
    import numpy as np

    ustar_dict = {}
    Dispersion_dict = {}

    for n in fromArrayList:
        ustar_dict[n] = np.sqrt(g * depthinterp_dict[n] * slope_dict[n])
        Dispersion_dict[n] = Dispersioncoeff_dict[n] * depthinterp_dict[n] * ustar_dict[n]

    return (ustar_dict, Dispersion_dict)
