def icbc(fromArrayList, Qinterp_dict, c_alpha, main_dir, sediment_conc_file):
    import os
    import numpy as np
    import pandas as pd

    frac_cols = ["SSC_clay", "SSC_silt", "SSC_coarse_silt"]

    def _fit_length(x, T):
        x = np.asarray(x, dtype=float).flatten()
        if x.size == 0:
            return np.zeros(T, dtype=float)
        if x.size >= T:
            return x[:T]
        pad = np.full(T - x.size, x[-1], dtype=float)
        return np.concatenate([x, pad])

    C_dict = {}

    for i, n in enumerate(fromArrayList, start=1):
        print(n)

        base = Qinterp_dict[n]
        nsteps, T = base.shape

        fname = f"Sub{i}_SSC.csv"
        fpath = os.path.join(main_dir, fname)

        ssc_df = pd.read_csv(fpath)

        if "node_id" in ssc_df.columns:
            ssc_df = ssc_df.loc[ssc_df["node_id"] == n].copy()

        if "datetime" in ssc_df.columns:
            ssc_df["datetime"] = pd.to_datetime(ssc_df["datetime"], utc=True, errors="coerce")
            ssc_df = ssc_df.sort_values("datetime")

        Ctemp_dict = {}
        for k in range(len(c_alpha)):
            arr = np.zeros((nsteps, T), dtype=float)
            if k < len(frac_cols) and frac_cols[k] in ssc_df.columns:
                ts = ssc_df[frac_cols[k]].to_numpy()
                bc = _fit_length(ts, T)
                arr[0, :] = bc
            Ctemp_dict[k] = arr

        C_dict[n] = Ctemp_dict

    return C_dict
