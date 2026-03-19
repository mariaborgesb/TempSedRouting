# %%
def interpolhydro(fromArrayList, toArray, unitflowdict, depthdict, Qdict, xpoints_dict):
    import numpy as np
    import scipy
    import pandas as pd

    qinterp_dict = {}      # interpolated unitflow
    depthinterp_dict = {}  # interpolated depth
    Qinterp_dict = {}      # interpolated discharge

    segment = 0
    n_to = len(toArray)  

    for n in fromArrayList:
        q_US     = np.array(unitflowdict[n])
        depth_US = np.array(depthdict[n])
        Q_US     = np.array(Qdict[n])

        ds_to_node = toArray[segment] if segment < n_to else None

        ds_from_node = None
        if ds_to_node is not None:
            try:
                idx = fromArrayList.index(ds_to_node)  # find the reach that starts at this TO_NODE
                ds_from_node = fromArrayList[idx]
            except ValueError:
                ds_from_node = None  # no downstream reach => it is a sink

        if (ds_from_node is not None) and (ds_from_node in unitflowdict) and (ds_from_node in depthdict) and (ds_from_node in Qdict):
            q_DS     = np.array(unitflowdict[ds_from_node])
            depth_DS = np.array(depthdict[ds_from_node])
            Q_DS     = np.array(Qdict[ds_from_node])
        else:
            q_DS     = q_US.copy()
            depth_DS = depth_US.copy()
            Q_DS     = Q_US.copy()

        cols = xpoints_dict[n] + 2  
        qinterp     = np.zeros((len(q_US),     cols))
        depthinterp = np.zeros((len(depth_US), cols))
        Qinterp     = np.zeros((len(Q_US),     cols))

        q_stack     = np.vstack([q_US,     q_DS]) 
        depth_stack = np.vstack([depth_US, depth_DS])
        Q_stack     = np.vstack([Q_US,     Q_DS])    

        qinterp_struct     = scipy.interpolate.interp1d([1, cols], q_stack,     kind='linear', axis=0, assume_sorted=True)
        depthinterp_struct = scipy.interpolate.interp1d([1, cols], depth_stack, kind='linear', axis=0, assume_sorted=True)
        Qinterp_struct     = scipy.interpolate.interp1d([1, cols], Q_stack,     kind='linear', axis=0, assume_sorted=True)

        for i in range(cols):
            pos = i + 1
            qinterp[:, i]     = qinterp_struct(pos)
            depthinterp[:, i] = depthinterp_struct(pos)
            Qinterp[:, i]     = Qinterp_struct(pos)

        qinterpt     = np.transpose(qinterp)
        depthinterpt = np.transpose(depthinterp)
        Qinterpt     = np.transpose(Qinterp)

        qinterp_dict[n]     = qinterpt
        depthinterp_dict[n] = depthinterpt
        Qinterp_dict[n]     = Qinterpt

        segment += 1
    
    interp_out = pd.DataFrame({'depth_interp': depthinterp_dict, 'q_interp': qinterp_dict})
    return (qinterp_dict, depthinterp_dict, Qinterp_dict, interp_out)
