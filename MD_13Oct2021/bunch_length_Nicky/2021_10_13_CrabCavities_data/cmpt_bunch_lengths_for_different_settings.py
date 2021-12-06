import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

saveflag=False
# get my current path
my_path = os.getcwd().replace('\\', '/')


coast_list = ['coast_01', 'coast_02', 'coast_03']
setting_list = [['0', '1', '2', '3', '3b'], ['1'], ['1', '2', '3']]

# .h5 filename
filename = '2021_10_13_CrabCavities.h5'
filepath = os.path.join(my_path, filename).replace('\\', '/')

   
# Manually find the indeces that correspond to each different setting
column_names = ['time start', 'time end', 'tau [ns], initial', 'tau [ns], average'] # initial and average tau=4sigmat over each different setting
df_results = pd.DataFrame(columns = column_names, dtype=object)
index_list = ['coast1-setting0', 'coast1-setting1', 'coast1-setting2', 'coast1-setting3', 'coast1-setting3b', 'coast2-setting1', 'coast3-setting1', 'coast3-setting2', 'coast3-setting3']



### COAST 1 ####
# get data collected during the selected coast
df = pd.read_hdf(filepath, 'coast_01')
time = df['time']
tau_mean = df['tau_mean'] # 4σt in [ns]
# Found manually accoridng to the WS acquisitions, summarized in
# https://docs.google.com/spreadsheets/d/13jB-KxTiWfWefqBvGPFbR9BgdIgJZjoynBYuQXUDbqc/edit?usp=sharing
coast1_setting0_start, coast1_setting0_end = 0, 110
coast1_setting1_start, coast1_setting1_end = 113, 183
coast1_setting2_start, coast1_setting2_end = 193, 282
coast1_setting3_start, coast1_setting3_end = 289, 348
coast1_setting3b_start, coast1_setting3b_end = 349, 379

to_append_1 = [time[coast1_setting0_start], time[coast1_setting0_end], tau_mean[coast1_setting0_start], np.mean(tau_mean[coast1_setting0_start:coast1_setting0_end+1])]
to_append_2 = [time[coast1_setting1_start], time[coast1_setting1_end], tau_mean[coast1_setting1_start], np.mean(tau_mean[coast1_setting1_start:coast1_setting1_end+1])]
to_append_3 = [time[coast1_setting2_start], time[coast1_setting2_end], tau_mean[coast1_setting2_start], np.mean(tau_mean[coast1_setting2_start:coast1_setting2_end+1])]
to_append_4 = [time[coast1_setting3_start], time[coast1_setting3_end], tau_mean[coast1_setting3_start], np.mean(tau_mean[coast1_setting3_start:coast1_setting3_end+1])]
to_append_5 = [time[coast1_setting3b_start], time[coast1_setting3b_end], tau_mean[coast1_setting3b_start], np.mean(tau_mean[coast1_setting3b_start:coast1_setting3b_end+1])]


df_results.loc[0] = to_append_1
df_results.loc[1] = to_append_2
df_results.loc[2] = to_append_3
df_results.loc[3] = to_append_4
df_results.loc[4] = to_append_5


### COAST 2 ###
df = pd.read_hdf(filepath, 'coast_02')
time = df['time']
tau_mean = df['tau_mean'] # 4σt in [ns]
coast2_setting1_start, coast2_setting1_end = 1, 53 # 0 seems to be an outlier, as it is too large 2.7 ns

to_append_6 = [time[coast2_setting1_start], time[coast2_setting1_end], tau_mean[coast2_setting1_start], np.mean(tau_mean[coast2_setting1_start:coast2_setting1_end+1])]
df_results.loc[5] = to_append_6

### COAST 3 ###
df = pd.read_hdf(filepath, 'coast_03')
time = df['time']
tau_mean = df['tau_mean'] # 4σt in [ns]
coast3_setting1_start, coast3_setting1_end = 7, 78
coast3_setting2_start, coast3_setting2_end = 96, 140
coast3_setting3_start, coast3_setting3_end = 147, 175

to_append_7 = [time[coast3_setting1_start], time[coast3_setting1_end], tau_mean[coast3_setting1_start], np.mean(tau_mean[coast3_setting1_start:coast3_setting1_end+1])]
to_append_8 = [time[coast3_setting2_start], time[coast3_setting2_end], tau_mean[coast3_setting2_start], np.mean(tau_mean[coast3_setting2_start:coast3_setting2_end+1])]
to_append_9 = [time[coast3_setting3_start], time[coast3_setting3_end], tau_mean[coast3_setting3_start], np.mean(tau_mean[coast3_setting3_start:coast3_setting3_end+1])]

df_results.loc[6] = to_append_7
df_results.loc[7] = to_append_8
df_results.loc[8] = to_append_9



df_results.index = index_list

print(df_results)
print(df_results['tau [ns], average'].to_string(index=False)) # Print without the index


if saveflag:
    df_results.to_pickle('./mean_tau_summary_CC_MD_13Oct2021.pkl')
