# %%
def sedcoeff (main_dir, median_sed_file, size_class, fromArrayList, hydroDF):
    import os
    import pandas as pd
    import numpy as np

    DmDF = pd.read_csv(os.path.join(main_dir, median_sed_file), header=0, sep=',').sort_values(by=["HYDSEQ"])

    cri_all_dict = {}  # critical streampower coefficients
    y_dict = {}        # maximum deposition depth
    xpoints_dict = {}  # number of interpolation steps
    Fi_dict = {}       # F coefficients

    for n in fromArrayList:
        reach_interest = hydroDF.index[hydroDF['FROM_NODE'] == n]
        
        # Build array for critical coefficients
        cri_temp = np.zeros((1, 3))
        cri_temp[0, 0] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'critical_strmpow_coef_d1'].iloc[0]
        cri_temp[0, 1] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'critical_strmpow_coef_d2'].iloc[0]
        cri_temp[0, 2] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'critical_strmpow_coef_d3'].iloc[0]
        cri_all_dict[n] = cri_temp
        del cri_temp

        # Maximum deposition depth
        y_dict[n] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'depo-depth'].iloc[0]

    # xpoints per reach
    for n in fromArrayList:
        reach_interest = hydroDF.index[hydroDF['FROM_NODE'] == n]
        xpoints_dict[n] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'xpoints'].iloc[0]

    # F coefficients per reach
    for n in fromArrayList:
        Fi_temp = np.zeros((1, 3))
        reach_interest = hydroDF.index[hydroDF['FROM_NODE'] == n]
        Fi_temp[0, 0] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'Fd1'].iloc[0]
        Fi_temp[0, 1] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'Fd2'].iloc[0]
        Fi_temp[0, 2] = DmDF.loc[DmDF['nzsegv2'] == reach_interest[0], 'Fd3'].iloc[0]
        Fi_dict[n] = Fi_temp
        del Fi_temp

    return DmDF, cri_all_dict, Fi_dict, y_dict, xpoints_dict
