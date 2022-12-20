import pandas as pd
from scipy.signal import find_peaks
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os as os

directory = 'Z:/Chaperone_subgroup/NickM/Miscellaneous/Amy/Example_data'  #### This is where you put your directory
filename = f'{directory}/cleaned_data.csv'
distance = 30000 #### This is where you set the distance between peaks

output_all = f'{directory}/all'
if not os.path.exists(output_all):
    os.makedirs(output_all)
output_major = f'{directory}/major'
if not os.path.exists(output_major):
    os.makedirs(output_major)

cleaned_data = pd.read_csv(filename, header="infer")


#### this code finds all peaks, not just major peak
def peak_plotter(cleaned_data, prominence, height):
    col_peaks = []
    for cell, df in cleaned_data.groupby('cell'):
        df.reset_index(inplace = True)
        cell_peak, _ = find_peaks(df["peakfinding"], prominence = prominence, height = height)
        peaks_list = pd.DataFrame([df.iloc[index][["peakfinding", 'Time (s)']] for index in cell_peak], index=cell_peak)
        peaks_list['peakfinding'] = peaks_list['peakfinding']*-1
        peaks_list['cell'] = cell
        col_peaks.append(peaks_list)
        plot1 = plt.figure(1)
        sns.set(style = "ticks", font_scale = 1)
        sns.lineplot(data = df, x = 'Time (s)', y = "Background corrected")
        sns.scatterplot(data = peaks_list, x = 'Time (s)', y = 'peakfinding', color = 'orange')
        plt.title(f'{prominence} and {height} and {cell}')
        plt.ylabel("Current (pA)")
        plot1.savefig(f'{output_all}/cell_{cell}.svg', dpi = 600)
        plt.show()
    collated_peak = pd.concat(col_peaks)
    collated_peak.to_csv(f'{output_all}/all_peaks.csv')
    return collated_peak

# peak_plotter(cleaned_data, 60, 20)

##### This code finds all the major peaks, which it does by finding the max response within a defined range of x values

def plot_function_major(cleaned_data, prominence, height, distance):
    """Finds the major peaks at every defined interval and then plots the results for each cell

    Args:
        cleaned_data (df): dataframe containing corrected data
        prominence (float): variable used to define the shape of the peak to help peak fitting
        height (float): variable to identify minimum height for peaks
        distance (variable): defined earlier, this variable is used so that peak finder scripts will look for a peak every X time units, where X is distance

    Returns:
        df: returns dataframe containing the index of each peak found by the script
    """
    col_peaks_major = []
    for cell, df in cleaned_data.groupby('cell'):
        df.reset_index(inplace = True)
        cell_peak_major, _ = find_peaks(df["peakfinding"], prominence = prominence, height = height, distance = distance)
        peaks_list_major = pd.DataFrame([df.iloc[index][["peakfinding", 'Time (s)']] for index in cell_peak_major], index=cell_peak_major)
        peaks_list_major['peakfinding'] = peaks_list_major['peakfinding']*-1
        peaks_list_major['pos_peaks'] = peaks_list_major['peakfinding']*-1
        peaks_list_major['cell'] = cell
        peaks_list_major['capacitance'] = df['Capacitance'].iloc[0]
        col_peaks_major.append(peaks_list_major)
        plot1 = plt.figure(1)
        sns.set(style = "ticks", font_scale = 1)
        sns.lineplot(data = df, x = 'Time (s)', y = "Background corrected")
        sns.scatterplot(data = peaks_list_major, x = 'Time (s)', y = 'peakfinding', color = 'orange')
        plt.title(f'{prominence} and {height} and {cell}')
        plt.ylabel("Current (pA)")
        plot1.savefig(f'{output_major}/cell_{cell}.svg', dpi = 600)
        plt.show()
    collated_peak_major = pd.concat(col_peaks_major)
    collated_peak_major['normalised_to_capacitance'] = collated_peak_major['pos_peaks'] / collated_peak_major['capacitance']
    collated_peak_major.to_csv(f'{output_major}/major_peaks.csv')
    return collated_peak_major

plot_function_major(cleaned_data, 30, 20, distance)


#### This can be used to optimize parameters for fitting and finding peaks
for prominence in range(30, 45, 5):
    for height in range(40, 60, 20):
        peak_plotter(cleaned_data, prominence, height)
# peak_data = pd.read_csv(f'{directory}/major/major_peaks.csv', header="infer")

