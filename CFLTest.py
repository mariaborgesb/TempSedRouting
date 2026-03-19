# %%
def cfl (fromArrayList, Dispersion_dict, depthinterp_dict, reach_distance_dict, xpoints_dict, dt):
    # The Courant-Friedrichs-Lewy (CFL) condition provides a stability criterion 
    # for numerical solutions of partial differential equations.
    # In solving the mass conservation equation numerically, 
    # a conservative CFL number of 0.1 was used to avoid significant numerical dispersion.
    import numpy as np

    CFL_min_dict = {}
    CFL_max_dict = {}

    for n in fromArrayList:
        dx = reach_distance_dict[n] / (xpoints_dict[n] + 2)

        if dx <= 0:
            CFL_min_dict[n] = np.nan
            CFL_max_dict[n] = np.nan
            continue

        Dispersion_min = np.min(Dispersion_dict[n])
        Dispersion_max = np.max(Dispersion_dict[n])

        CFL_min_dict[n] = (Dispersion_min * dt) / (dx ** 2)
        CFL_max_dict[n] = (Dispersion_max * dt) / (dx ** 2)

    return CFL_min_dict, CFL_max_dict
