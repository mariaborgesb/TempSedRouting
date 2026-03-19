# %%
def hydrologydata(firstdate, lastdate, hydroDF, US_reach, fromArrayList, flow_dir_name, nsubs=None):
    #Extract unitflow, depth and discharge for each reach.

    import os
    import pandas as pd
    import numpy as np

    # number of files / reaches
    n_rows = len(hydroDF) if nsubs is None else int(nsubs)

    datearray = None
    unitflowarray = None
    deptharray = None
    Qarray = None

    first_dt = pd.to_datetime(firstdate, utc=True)
    last_dt  = pd.to_datetime(lastdate,  utc=True)

    for m in range(n_rows):
        sub_id = m + 1
        fname = f"Sub{sub_id}_3h.csv"
        fpath = os.path.join(flow_dir_name, fname)

        # read file
        df = pd.read_csv(fpath, header=0, sep=",", parse_dates=["datetime"])
        df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

        # window of selected dates
        mask = (df["datetime"] >= first_dt) & (df["datetime"] <= last_dt)
        df = df.loc[mask].copy()
        df["datetime"] = df["datetime"].dt.tz_convert(None)

        # select the columns
        t  = df["datetime"].to_numpy()
        uq = df["unit_q"].to_numpy(dtype=float)
        d  = df["depth"].to_numpy(dtype=float)
        Q  = df["flow3h"].to_numpy(dtype=float)

        if datearray is None:
            T = len(t)
            datearray = t
            unitflowarray = np.zeros((n_rows, T), dtype=float)
            deptharray = np.zeros((n_rows, T), dtype=float)
            Qarray = np.zeros((n_rows, T), dtype=float)

        # fill arrays rows
        unitflowarray[m, :] = uq
        deptharray[m, :]    = d
        Qarray[m, :]        = Q

    unitflowdict = dict(zip(fromArrayList[:n_rows], unitflowarray.tolist()))
    depthdict    = dict(zip(fromArrayList[:n_rows], deptharray.tolist()))
    Qdict        = dict(zip(fromArrayList[:n_rows], Qarray.tolist()))

    return datearray, unitflowdict, depthdict, Qdict
