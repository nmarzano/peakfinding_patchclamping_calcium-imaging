from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


output_dir = 'directory_to_save'
filename = 'directory/cleaned_data.csv'
major_peaks = pd.read_csv(filename, header="infer")


def plot_cells(dfs, xmin, xmax):
    sns.set(style = 'ticks', font_scale = 1)
    for cell, df in dfs.groupby('cell'):
        fig, ax = plt.subplots(2, 1)
        sns.lineplot(data = df, x = "Time (s)", y = "Background corrected", color = "black", ax = ax[0])
        sns.lineplot(data = df, x = "Time (s)", y = "Background corrected", color = "black", ax = ax[1])
        ax[1].set_xlim(xmin, xmax)
        fig.savefig(f'{output_dir}/raw_plot_cell{cell}.eps', dpi = 600)
        plt.show()

plot_cells(major_peaks, 70.5, 71)


