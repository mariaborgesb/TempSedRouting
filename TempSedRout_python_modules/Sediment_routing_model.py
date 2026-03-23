# %% Sediment routing model
# one dimensional sediment routing model with the relationship describing
# conservation of sediment mass subjected to advection and dispersion

#  Read paths & filenames (edit if needed)
main_dir = r"\\file\Usersm$\mbo95\Home\My Documents\Sediment routing model\Files for the model"
main_dir_routing = main_dir
flow_dir_name = main_dir

routing_file = "routing_file.csv"
sediment_size_file = "sediment-size-classes_Paper.csv"
median_sed_file = "median_sed_file.csv"          
reach_charac_file = "reach_charac_file.csv"    
sediment_conc_file = "sediment_concentration.csv" 

firstdate = "2024-12-28 00:00:00"
lastdate  = "2025-01-01 21:00:00"

#%% import packages
import pandas as pd
import numpy as np
import os
import scipy

# Import modules
import RiverRouting            as RiverRouting
import SuspSedChar             as SuspSedChar
import SedimentCoeff           as SedimentCoeff
import ReachChar               as ReachChar
import Hydrologicalinput       as Hydrologicalinput
import InterploationHydro      as InterploationHydro   # keep original variable name
import DispersionUstar         as DispersionUstar
import CFLTest                 as CFLTest
import MatrixVars              as MatrixVars
import ICBCsediment            as ICBCsediment
import FallVelocity            as FallVelocity
import CriticalValues          as CriticalValues
import Re_entrainment_rate     as Re_entrainment_rate
import Scenario_model          as Scenario_model
import Deposition_rate         as Deposition_rate

# %% reading hydroDF file which has all nzsegv2, HYDSEQ, FROM_NODE, TO_NODE
hydroDF, US_reach, fromArrayList, toArray, toArrayList, uniqueNodeList = RiverRouting.riverrouting(
    main_dir_routing, routing_file)

_h = hydroDF.reset_index(drop=False)
SINK_REACHES = set(_h.loc[_h["NextDownID"] == 0, "HYDSEQ"])

#%% import suspended sediment sizes, and alpha coefficients for each fraction
size_class, di, c_alpha = SuspSedChar.suspchar(main_dir, sediment_size_file)

#%% import median riverbed diameter and critical stream power coefficients as dataframe
DmDF, cri_all_dict, Fi_dict, y_dict, xpoints_dict = SedimentCoeff.sedcoeff(
    main_dir, median_sed_file, size_class, fromArrayList, hydroDF)

#%% import reach characteristics data (slope and distance from do)
slope_dict, reach_distance_dict, reachID_dict, Dm_dict, Dispersioncoeff_dict = ReachChar.reachchar(
    DmDF, uniqueNodeList, fromArrayList, hydroDF, main_dir, reach_charac_file)

#%% import constants
Rrho=1.65       # R=((rohs/rho)-1),rhos=2650, rho=1000 - no unit
g=9.806         # acceleration gravity [m/s^2]
mcoef=1         # drag/roughness coefficients (1-10) [-]
kv=0.3          # von-Karman coefficient
SF=0.1          # sand fraction in proportion (0 - 1)
rho=1000        # Density of water [kg/m^3]
rhos=2650       # sediment density [kg/m^3]
Rrho=1.65       # R=((rohs/rho)-1)
g=9.805         # acceleration gravity [m/s^2]
nu=0.0000010533 # kinematic viscosity [m^2/s]
#% import time step
dt=60 * 60 * 3   # time step 3 hour [s]

#%% read flow data
datearray, unitflowdict, depthdict, Qdict = Hydrologicalinput.hydrologydata(
    firstdate, lastdate, hydroDF, US_reach, fromArrayList, flow_dir_name)

#%--------------------Calculations-----------------------
#%% Interpolation of unitflow (qinterp_dict), flow (Qinterp_dict), depth (depthinterp_dict)
qinterp_dict, depthinterp_dict, Qinterp_dict, interp_out = InterploationHydro.interpolhydro(
    fromArrayList, toArray, unitflowdict, depthdict, Qdict, xpoints_dict)
# Debug: check dtype
n0 = fromArrayList[0]
print("Qinterp dtype:", Qinterp_dict[n0].dtype)

#%% calculate dispersion coefficient [m^2/s] and ustar [m/s]
ustar_dict, Dispersion_dict = DispersionUstar.disperustar(
    fromArrayList, depthinterp_dict, slope_dict, g, Dispersioncoeff_dict)

#%% Courant-Friedrichs-Lewy (CFL) condition
CFL_min_dict, CFL_max_dict = CFLTest.cfl(
    fromArrayList, Dispersion_dict, depthinterp_dict, reach_distance_dict, xpoints_dict, dt)

#%% calculation of alpha, beta, lambda for matrix
A1_dict, A2_dict, A3_dict, M1_dict, M2_dict, M3_dict = MatrixVars.matrixvar(
    fromArrayList, dt, reach_distance_dict, depthinterp_dict, qinterp_dict, xpoints_dict, Dispersion_dict)

#%% import C
C_dict = ICBCsediment.icbc(fromArrayList, Qinterp_dict, c_alpha, main_dir, sediment_conc_file)

#%% calculate fall velocity
falveli = FallVelocity.fallvelocity(di)
falvel_alphai = falveli * c_alpha
print("falvel_alphai:", falvel_alphai)

# Check first reach, first fraction, first cell
print("C max/min:", np.max(C_dict[n0][0]), np.min(C_dict[n0][0]))

#%% Calculate critical discharge and critical streampower (regression) — SAFE version
# after:
Strmpowcri_dict, Qcr_dict, qcr_dict, widthcr_dict, Strmpow_dict, qcri = CriticalValues.critical_strmpow_q_Q(
    fromArrayList, SF, kv, rho, Rrho, g, Dm_dict, di, mcoef, slope_dict, Qinterp_dict, qinterp_dict)

#%% calculate re-entrainment rate [kg/m^2/s]
ryd_dict = Re_entrainment_rate.reentrainment(
    fromArrayList, qcri, Qinterp_dict, Fi_dict, Strmpow_dict, cri_all_dict,
    Strmpowcri_dict, depthinterp_dict, rhos, rho, g)

#%% scenario selection
scen_dict, E_dict, dyd_dict = Scenario_model.scenario_model(
    fromArrayList, qcri, Strmpow_dict, Qinterp_dict, Strmpowcri_dict)

#%% deposition rate [kg/m^2/s]
ded_dict = Deposition_rate.depositionrate(fromArrayList, qcri, Qinterp_dict, C_dict, falvel_alphai)

#%% Calculate C: Matrix solution
n_us1 = 0
n_us2 = 0
for n in fromArrayList:
    if n == n_us1 or n == n_us2:
        continue
    else:
        num_step, num_time = C_dict[n][0].shape  # (xpoints+2, time)
        tonode_temp = hydroDF.TO_NODE[hydroDF['FROM_NODE'] == n].iloc[0]
        toreach_temp = hydroDF.index[hydroDF['TO_NODE'] == tonode_temp]

        # 2-tributaries
        if len(toreach_temp) == 2:
            n_us1 = hydroDF.FROM_NODE[toreach_temp[0]]
            n_us2 = hydroDF.FROM_NODE[toreach_temp[1]]
            for k in range(len(qcri)):  # loop per fraction
                # --- n_us1 ---
                t_before = 0
                for j in range(1, num_time):
                    ded_dict[n_us1][k][:, t_before] = falvel_alphai[k] * C_dict[n_us1][k][:, t_before]
                    ded_dict[n_us1][k][:, t_before][ded_dict[n_us1][k][:, t_before] < 0] = 0
                    for l in range(len(Qinterp_dict[n_us1])):
                        E_dict[n_us1][k][l, t_before] = ryd_dict[n_us1][k][l, t_before] - ded_dict[n_us1][k][l, t_before]
                        dyd_dict[n_us1][k][l, t_before] = ((ded_dict[n_us1][k][l, t_before] - ryd_dict[n_us1][k][l, t_before]) / rhos) * dt

                    A = scipy.sparse.spdiags(
                        [np.append(A1_dict[n_us1][1:, j], [0]), A2_dict[n_us1][:, j], np.append([0], A3_dict[n_us1][:-1, j])],
                        (-1, 0, 1), xpoints_dict[n_us1] + 2, xpoints_dict[n_us1] + 2
                    ).toarray()
                    M = scipy.sparse.spdiags(
                        [np.append(M1_dict[n_us1][1:, j-1], [0]), M2_dict[n_us1][:, j-1], np.append([0], M3_dict[n_us1][:-1, j-1])],
                        (-1, 0, 1), xpoints_dict[n_us1] + 2, xpoints_dict[n_us1] + 2
                    ).toarray()
                    MC = np.matmul(M, C_dict[n_us1][k][:, t_before])
                    MCE = MC + E_dict[n_us1][k][:, t_before]
                    Cleft_array = np.linalg.solve(A, MCE)
                    Cleft_array[Cleft_array < 0] = 0
                    Cleft_array[0] = C_dict[n_us1][k][0, j]
                    C_dict[n_us1][k][:, j] = Cleft_array
                    t_before += 1
                    del Cleft_array, A, M, MC, MCE

                # --- n_us2 ---
                t_before = 0
                for j in range(1, num_time):
                    ded_dict[n_us2][k][:, t_before] = falvel_alphai[k] * C_dict[n_us2][k][:, t_before]
                    ded_dict[n_us2][k][:, t_before][ded_dict[n_us2][k][:, t_before] < 0] = 0
                    for l in range(len(Qinterp_dict[n_us2])):
                        E_dict[n_us2][k][l, t_before] = ryd_dict[n_us2][k][l, t_before] - ded_dict[n_us2][k][l, t_before]
                        dyd_dict[n_us2][k][l, t_before] = ((ded_dict[n_us2][k][l, t_before] - ryd_dict[n_us2][k][l, t_before]) / rhos) * dt

                    A = scipy.sparse.spdiags(
                        [np.append(A1_dict[n_us2][1:, j], [0]), A2_dict[n_us2][:, j], np.append([0], A3_dict[n_us2][:-1, j])],
                        (-1, 0, 1), xpoints_dict[n_us2] + 2, xpoints_dict[n_us2] + 2
                    ).toarray()
                    M = scipy.sparse.spdiags(
                        [np.append(M1_dict[n_us2][1:, j-1], [0]), M2_dict[n_us2][:, j-1], np.append([0], M3_dict[n_us2][:-1, j-1])],
                        (-1, 0, 1), xpoints_dict[n_us2] + 2, xpoints_dict[n_us2] + 2
                    ).toarray()
                    MC = np.matmul(M, C_dict[n_us2][k][:, t_before])
                    MCE = MC + E_dict[n_us2][k][:, t_before]
                    Cleft_array = np.linalg.solve(A, MCE)
                    Cleft_array[Cleft_array < 0] = 0
                    Cleft_array[0] = C_dict[n_us2][k][0, j]
                    C_dict[n_us2][k][:, j] = Cleft_array
                    t_before += 1
                    del Cleft_array, A, M, MC, MCE

                # Safe downstream BC/IC for 2-tributary case
                if tonode_temp in C_dict:
                    C_dict[tonode_temp][k][0, :] = (
                        C_dict[n_us1][k][num_step - 1, :] + C_dict[n_us2][k][num_step - 1, :]
                    )
                    C_dict[tonode_temp][k][:, 0] = (
                        C_dict[n_us1][k][0, 0] + C_dict[n_us2][k][0, 0]
                    )
                # else: tonode_temp may be 0 (true sink); skip

        # 1-tributary
        elif len(toreach_temp) == 1:
            n_us = hydroDF.FROM_NODE[toreach_temp[0]]
            for k in range(len(qcri)):
                t_before = 0
                for j in range(1, num_time):
                    ded_dict[n_us][k][:, t_before] = falvel_alphai[k] * C_dict[n_us][k][:, t_before]
                    ded_dict[n_us][k][:, t_before][ded_dict[n_us][k][:, t_before] < 0] = 0
                    for l in range(len(Qinterp_dict[n_us])):
                        E_dict[n_us][k][l, t_before] = ryd_dict[n_us][k][l, t_before] - ded_dict[n_us][k][l, t_before]
                        dyd_dict[n_us][k][l, t_before] = ((ded_dict[n_us][k][l, t_before] - ryd_dict[n_us][k][l, t_before]) / rhos) * dt

                    A = scipy.sparse.spdiags(
                        [np.append(A1_dict[n_us][1:, j], [0]), A2_dict[n_us][:, j], np.append([0], A3_dict[n_us][:-1, j])],
                        (-1, 0, 1), xpoints_dict[n_us] + 2, xpoints_dict[n_us] + 2
                    ).toarray()
                    M = scipy.sparse.spdiags(
                        [np.append(M1_dict[n_us][1:, j-1], [0]), M2_dict[n_us][:, j-1], np.append([0], M3_dict[n_us][:-1, j-1])],
                        (-1, 0, 1), xpoints_dict[n_us] + 2, xpoints_dict[n_us] + 2
                    ).toarray()
                    MC = np.matmul(M, C_dict[n_us][k][:, t_before])
                    MCE = MC + E_dict[n_us][k][:, t_before]
                    Cleft_array = np.linalg.solve(A, MCE)
                    Cleft_array[Cleft_array < 0] = 0
                    Cleft_array[0] = C_dict[n_us][k][0, j]
                    C_dict[n_us][k][:, j] = Cleft_array
                    t_before += 1
                    del Cleft_array, A, M, MC, MCE

                # Safe downstream BC/IC for 1-tributary case
                if tonode_temp in C_dict:
                    C_dict[tonode_temp][k][0, :] = C_dict[n_us][k][num_step - 1, :]
                    C_dict[tonode_temp][k][:, 0] = C_dict[n_us][k][0, 0]
                # else: sink; skip

#EXPORT: per-reach time series (upstream & downstream)
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt  # also used for plots below

# Name fractions consistently with your c_alpha setup
fraction_names = ["clay", "silt", "coarse_silt"][:len(c_alpha)]

# Output folder
out_dir = os.path.join(main_dir, "outputs")
os.makedirs(out_dir, exist_ok=True)

def reach_timeseries_df(n):
    #Build a DataFrame with upstream/downstream SSC for all fractions for reach n
    if n not in C_dict: 
        return None

    num_step, num_time = C_dict[n][0].shape
    T = num_time

    t = None
    try:
        t_try = pd.to_datetime(datearray)
        if hasattr(t_try, "__len__") and len(t_try) == T:
            t = t_try
    except Exception:
        t = None
    if t is None:
        t = pd.date_range(start=pd.to_datetime(firstdate), periods=T, freq=f"{int(dt)}s")

    data = {"datetime": t}
    for k, name in enumerate(fraction_names):
        up = C_dict[n][k][0, :]                 # upstream boundary (node 0)
        down = C_dict[n][k][num_step - 1, :]    # downstream node (final)
        data[f"up_{name}"] = up
        data[f"down_{name}"] = down

    df = pd.DataFrame(data).set_index("datetime")
    return df

# Write one CSV per reach
all_reaches = []
for n in fromArrayList:
    df = reach_timeseries_df(n)
    if df is None:
        continue
    csv_path = os.path.join(out_dir, f"reach_{n}_SSC_timeseries.csv")
    df.to_csv(csv_path, float_format="%.6g")
    all_reaches.append((n, df))
print(f"✅ Wrote {len(all_reaches)} per-reach CSVs to: {out_dir}")

# Write a single wide CSV combining all reaches (downstream only)
wide_rows = []
for n, df in all_reaches:
    row = df.filter(like="down_").copy()
    row.columns = [f"{c}_r{n}" for c in row.columns]
    wide_rows.append(row)

if wide_rows:
    wide_df = pd.concat(wide_rows, axis=1)
    wide_df.to_csv(os.path.join(out_dir, "all_reaches_SSC.csv"), float_format="%.6g")
    print("✅ Wrote combined downstream file")
