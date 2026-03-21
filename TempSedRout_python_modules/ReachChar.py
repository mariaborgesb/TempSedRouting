# %% 
def reachchar (DmDF, uniqueNodeList, fromArrayList, hydroDF, main_dir, reach_charac_file):
    # River reach characteristics including river bed gradient and distance from next reach downstream
    import os
    import pandas as pd
    import numpy as np

    slope_dict = {}
    reach_distance_dict = {}
    reachID_dict = {}
    Dm_dict = {}
    Dispersioncoeff_dict = {}

    reach_charac = pd.read_csv(os.path.join(main_dir, reach_charac_file), header=0, sep=',')

    h=hydroDF.reset_index(drop=False).copy()

    for n in uniqueNodeList:
        if fromArrayList.count(n)>0:
            reach_interest=hydroDF.index[hydroDF['FROM_NODE']==n]
            if len(reach_interest)==0:
                continue
            nzseg=reach_interest[0]

            slope=reach_charac.loc[reach_charac['nzsegv2']==nzseg, 'rch_slope_grad'].iloc[0]
            reach_distance=reach_charac.loc[reach_charac['nzsegv2']==nzseg, 'rch_length_m'].iloc[0]

            slope_dict[n]=slope
            reach_distance_dict[n]=reach_distance
            reachID_dict[n]=nzseg
            Dm_dict[n]=DmDF.loc[DmDF['nzsegv2']==nzseg, 'median_diameter'].iloc[0]
            Dispersioncoeff_dict[n]=DmDF.loc[DmDF['nzsegv2']==nzseg, 'Dispersion'].iloc[0]

        else:
            # identify sinks
            row=h[h['TO_NODE']==n]
            if not row.empty and int(row['NextDownID'].iloc[0])==0:
                print(f"This is the sink reach {int(row['HYDSEQ'].iloc[0])}")

    return slope_dict, reach_distance_dict, reachID_dict, Dm_dict, Dispersioncoeff_dict
