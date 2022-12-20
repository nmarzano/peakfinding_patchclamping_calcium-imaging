# A.J.Hulme PhD thesis: Scripts to find peaks within patch-clamping and calcium-imaging datasets

**Pre-requisites:** 

Patch-clamping (column 1 = Time (ms), column 2 = pA, column 3 = Capacitance)

Calcium imaging (column 1 = Drug, column 2 = Axis [s] for time, every subsequent column is the calcium response from a different cell)

**Workflow**: 

1a-initial_cleanup - Clean up patch-clamp data by subtracting background and calculating time in seconds
1b-find_peaks - Finds the peaks (either all peaks or just major peaks) and then saves data to directory + plots
1c-Plotting_patch_clamp_results - Plots a subplot of each cell (top panel is the entire poking experiment and the bottom is a zoomed view of a single peak)
2a-calcium-imaging_peaks - Finds the maximum peak for each treatment within a single cell and then plots. Also calculates the normalised response of each treatment compared to a control treatment. 

