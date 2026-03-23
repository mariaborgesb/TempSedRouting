# %%
def suspchar(main_dir, sediment_size_file):
    # fine sediment size classes that will be modelled 
    # input: directory to CSV file with sediment classifications
    # output: sediment size class, particle size [m], and concentration profile parameter
    import os
    import pandas as pd

    sediment_size = pd.read_csv(os.path.join(main_dir, sediment_size_file), header=0, sep=',')  # sediment size in micrometer
    size_class = sediment_size['class_size'].values
    di = sediment_size['size'].values / 1_000_000  # convert micrometers to meters
    c_alpha = sediment_size['alpha']  # alpha parameter for concentration profile

    return size_class, di, c_alpha
