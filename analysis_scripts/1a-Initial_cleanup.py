from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import glob as glob
# from operator import itemgetter

directory = 'Z:/Chaperone_subgroup/NickM/Miscellaneous/Amy/Example_data'  ##### This is where you put the directory with your data
numberofframes = 80  ######### can adjust the number of frames to average for your background here
output_folder = f'{directory}/' 

def csv_import(input_folder):
    # if not column_names:
    #     print('no column_names found, specify list to use')
    #     return
    filenames = glob.glob(f'{input_folder}/*.csv')
    dfs = []
    for filename in filenames:
        df = pd.read_csv(filename, header = 'infer')
        df['cell'] = filename.split('_')[-2]
        dfs.append(df)
        # dfs['cell'] = filename.split('_')[-2]
    test = pd.concat(dfs, ignore_index=True)
    # test_dfs.columns = column_names
    return pd.DataFrame(test)

dfs = csv_import(directory) 


##########
########## Clean up code by subtracting background and calculating the time in seconds
collated = []
for cell, df in dfs.groupby('cell'):
    background = df["pA"].iloc[:numberofframes].mean()
    df["Background corrected"] = df["pA"] - background
    df["peakfinding"] = df["Background corrected"]*-1
    df["Time (s)"] = df["Time (ms)"] / 1000
    collated.append(df)
collated_all = pd.concat(collated)
collated_all = pd.DataFrame(collated_all)
collated_all.to_csv(f'{output_folder}/cleaned_data.csv', index = False)

