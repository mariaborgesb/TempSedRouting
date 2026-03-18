#%%
def riverrouting (main_dir_routing, routing_file):
    import pandas as pd
    import numpy as np
    import os

    hydroDF = pd.read_csv(os.path.join(main_dir_routing, routing_file)).sort_values(by=["HYDSEQ"])
    US_reach = hydroDF["nzsegv2"].values
    hydroDF = hydroDF.set_index("nzsegv2")
    fromArray = hydroDF["FROM_NODE"].values
    toArray = hydroDF["TO_NODE"].values
    nodeArray = np.append(toArray, fromArray)
    uniqueNodeArray = np.unique(nodeArray)
    uniqueNodeList = uniqueNodeArray.tolist()
    entryArray = np.zeros(len(uniqueNodeList)).astype(float)
    fromArrayList = fromArray.tolist()
    uniquetoArray = np.unique(toArray)
    toArrayList = uniquetoArray.tolist()

    return hydroDF, US_reach, fromArrayList, toArray, toArrayList, uniqueNodeList
