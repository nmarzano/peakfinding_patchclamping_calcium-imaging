from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import glob as glob
import os
# from operator import itemgetter

directory = 'Z:/Chaperone_subgroup/NickM/Miscellaneous/Amy/Calcium-imaging'
output_folder = f'{directory}/saved_data'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
#########
######### Provide directory of folder containing your datasets and combine into dataframe
#########

def csv_import(input_folder):
    """combines data from multiple csv files into a single dataframe.

    Args:
        input_folder (str): copy and paste directory of folder containing csv files

    Returns:
        df: returns a dataframe containing the collated data from all files within the provided directory
    """
    filenames = glob.glob(input_folder + "/*.csv")
    dfs = []
    for filename in filenames:
        df = pd.read_csv(filename, header = 'infer', index_col=False)
        dfs.append(df)
    test = pd.concat(dfs, ignore_index=True)
    return pd.DataFrame(test)

dfs = csv_import(directory) 
dfs2 = pd.melt(dfs, id_vars=['Drug','Axis [s]'], var_name= 'roi', value_name='ratio')

##########
########## background subtraction code so that all responses are normalised to 0 at time = 0
########## done by calcularting the average response between two x values (bottom_thresh and top_thresh)
########## and then subtracting that mean value from all others

bottom_thresh = 10  #### start point 
top_thresh = 30
bc_dfs = []
for roi, df in dfs2.groupby('roi'):
    mean_background = df[(df['Axis [s]'] > bottom_thresh) & (df['Axis [s]'] < top_thresh)]['ratio'].mean()
    df['corrected'] = df['ratio'] - mean_background
    bc_dfs.append(df)
compiled = pd.concat(bc_dfs)

######## plots all cells to ensure background subtraction has worked as expected
fig, ax = plt.subplots()
sns.set(style = 'ticks', font_scale = 1)
sns.lineplot(data = compiled, x = 'Axis [s]', y = 'corrected', hue = 'roi', palette = 'mako_r', legend = False)
plt.xlabel('Time (s)')
plt.ylabel('Corrected ratio (a.u.)')
plt.show()


def peak_plotter_all(dfs):
    """Finds the maximum peak within each drug treatment and combines into a dataframe. For each cell, the trace is then plotted and the max peak for each drug treatment is plotted overlayed

    Args:
        dfs (df): dataframe containing the background subtracted dataset for each cell

    Returns:
        df: returns a dataframe containing the value of the max peak within each drug treatment for each cell
    """
    col_peaks = []
    for roi, df in dfs.groupby('roi'):
        df.reset_index(inplace = True)
        for drug, df_2 in df.groupby('Drug'):
            cell_peak, _ = find_peaks(df_2["corrected"], distance = df_2['Axis [s]'].max())
            peaks_list = pd.DataFrame([df_2.iloc[index][["corrected", "Axis [s]"]] for index in cell_peak], index=cell_peak)
            peaks_list_index = peaks_list.index
            peaks_list['roi'] = roi
            peaks_list['drug'] = drug
            col_peaks.append(peaks_list)
        test = pd.concat(col_peaks)
        fig, ax = plt.subplots()
        sns.set(style = 'ticks', font_scale = 1)
        sns.lineplot(data = df, x = 'Axis [s]', y = "corrected", color = 'black')
        sns.scatterplot(data = test[test['roi'] == roi], x = 'Axis [s]', y = 'corrected', color = 'orange')
        plt.ylabel("Current (pA)")
        fig.savefig(f'{output_folder}/cell_{roi}.svg', dpi = 600)
        plt.show()
    collated_peak = pd.concat(col_peaks)
    collated_peak.to_csv(f'{output_folder}/all_peaks.csv')
    return collated_peak


peak_all = peak_plotter_all(compiled)
peak_drug_only = peak_all[peak_all['drug'] != 'SBS1']
peak_drug_only.to_csv(f'{output_folder}/peak_data.csv')




############## Plot the response for all cells for each condition
fig, ax = plt.subplots()
sns.set(style = 'ticks', font_scale = 1)
sns.swarmplot(data = peak_drug_only, y = 'corrected', x = 'drug', color = 'black', alpha = 0.5)
sns.violinplot(data = peak_drug_only, y = 'corrected', x = 'drug', ci = 'sd', palette = 'mako_r')
plt.xlabel('')
plt.ylabel('Response (a.u.)')
plt.show()

list_of_variable = list(pd.unique(peak_drug_only['drug']))

#############
############# Calculates max response of each treatment and then normalises the signal of each
############# treatment to that of a designated control treatment
def calculate_ratio(dfs, control, list_of_variable):
    compiled_ratio = []
    for roi, df in dfs.groupby('roi'):
        df['percent of control'] = ([df[df['drug']== treatment]['corrected'].mean() for treatment in list_of_variable]/df[df['drug']==control]['corrected'].mean())*100
        compiled_ratio.append(df)
    test = pd.concat(compiled_ratio)
    return test[test['drug'] != control]

percent_response_df_all = calculate_ratio(peak_drug_only, 'HIGHK', list_of_variable)


fig, ax = plt.subplots()
sns.set(style = 'ticks', font_scale = 1)
sns.swarmplot(data = percent_response_df_all, x = 'drug', y = 'percent of control', color = 'black', alpha = 0.5)
sns.violinplot(data = percent_response_df_all, y = 'percent of control', x = 'drug', ci = 'sd', palette = 'mako_r')
plt.xlabel('')
plt.ylabel('Response normalised to positive control (%)')
plt.show()


