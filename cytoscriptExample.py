# libraries
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
import os
from cytoscript import CytoScript
import plot_helper_fns

# set up some constants that are used throughout
log_left_lim = 4.7  # H_log_left_lim defines left extreme of plots
log_right_lim = 7  # R1_647_H_log_thresh defines right extreme of plots

# set up square gate for closed cow hot spot
hot_spot_log_left = 5.7  # marks left edge of hot-spot
hot_spot_log_right = 7  # marks right edge of hot-spot
hot_spot_log_low = 0.95
hot_spot_log_high = 0.975

log_FSC_H_log_SSC_H_hi_lim = (
    1.15  # marks upper limit of log_FSC_H / log_SSC_H for cell type of interest
)
_2D_hist_range = [[log_left_lim, log_right_lim], [0.8, 1]]
log_FSC_H_log_SSC_H_range = [0.8, 1.4]

# singlet_poly defines a polygon that estimates if an event is a singlet
# Points for which FSC-A vs FSC-H falls inside the polygon are presumed singlets
# points outside the polygon are presumed non-singlets, and are rejected
singlet_poly = [
    [4.4, 4.7],
    [4.54, 4.9],
    [5.3, 5.7],
    [5.6, 5.75],
    [4.7, 4.6],
    [4.4, 4.7],
]

# set up square gate for closed cow hot spot
R1_647_H_log_left=5.7            # marks left edge of hot-spot
R1_647_H_log_right=7            # marks right edge of hot-spot
R1_647_H_log_low=0.95
R1_647_H_log_high=0.975

R1_647_H_log_left_lim=4.7         # R1_647_H_log_thresh defines lef extremum of plots
R1_647_H_log_right_lim=7       # R1_647_H_log_thresh defines lef extremum of plots
log_FSC_H_log_SSC_H_hi_lim = 1.15 # marks upper limit of log_FSC_H / log_SSC_H for cell type of interest
R1_647_2D_hist_range = [[R1_647_H_log_left_lim, R1_647_H_log_right_lim], [0.8, 1]]

hist_bg_color = "#C7D8E9"
log_FSC_H_log_SSC_H_range = [0.5, 1.2]
hist_bg_color = "#C7D8E9"

# set up a path to test data.
# the target folder is expected to hold a bunch of fcs files
# .. or csv files that contain FCS data
base_path = os.path.dirname(os.path.abspath(__file__))
plot_dir = os.path.join(base_path, "plots")
test_data_dir = os.path.join(
    base_path,
    "P111 EVERGREEN #20 PLATE 1 CSV",  # "FlowCytometryTools", "tests", "data", "Plate02"
)

# create a new CytoScript object
# and set the test data directory as its working dir
cytoScriptInstance = CytoScript()
cytoScriptInstance.setWorkingDir(test_data_dir)


# iterate over the files in the test data directory
# cytoScriptInstance.workingDirFiles() limits the files to csv or fcs files ..
for filename in cytoScriptInstance.workingDirFiles():
    # use the filename to generate a "partial path" to file names that we will generate
    plotNamePartialPath = os.path.join(plot_dir, filename[0:-4])

    # load the target file. When loa
    cytoScriptInstance.load_csv_or_fcs(filename)

    # get the newly loaded fcs data as a dataframe
    df = cytoScriptInstance.df_dict["df"]

    # calculate logarithmic values for all columns
    cytoScriptInstance.log10(removeNA=False)

    # calculate singlet membership using the singlet_poly applied to the 2d set defined by the log10(BL2 PI-A) vs log10(BL2 PI-H),
    # Stash the singlets in the cytoScriptInstance df_dict using the key "df_singlet_events"
    # For convenience, get a reference to the singlet dataframe in a local variable named df_singlet_events
    cytoScriptInstance.calcPolygonGate(
        "log10(BL2 PI-H)", "log10(BL2 PI-A)", singlet_poly, "is_singlet"
    )
    df_is_singlet = cytoScriptInstance.df_dict["df_is_singlet"] = df[df["is_singlet"]]

    # find singlet events with a particular size and roughness
    # by applying the "log_FSC_H_vs_log_SSC_H_ellipse" to log10(FSC-H) vs log10(SSC-H)
    log_FSC_H_vs_log_SSC_H_ellipse = {
        "center": (6.25, 5.85),
        "width": 0.3,
        "height": 0.5,
    }
    cytoScriptInstance.calcEllipticalGate(
        "log10(FSC-H)", "log10(SSC-H)", log_FSC_H_vs_log_SSC_H_ellipse, "is_hot_spot"
    )
    df["is_hot_spot"] = np.logical_and(df["is_hot_spot"], df["is_singlet"])
    df_hot_spot = cytoScriptInstance.df_dict["df_hot_spot"] = df[df["is_hot_spot"]]

    # Draw the singlet map
    x_column = "log10(BL2 PI-H)"
    y_column = "log10(BL2 PI-A)"
    plt.scatter(df[x_column].values, df[y_column].values, color="black", s=1)
    plt.scatter(
        df_is_singlet[x_column].values,
        df_is_singlet[y_column].values,
        color="blue",
        s=1,
    )
    plt.scatter(
        df_hot_spot[x_column].values, df_hot_spot[y_column].values, color="yellow", s=1
    )
    plt.title(filename + ": Singlet Map")
    plt_xlim = (4.4, 5.7)
    plt_ylim = (4.6, 5.9)
    plt.xlim(plt_xlim)
    plt.ylim(plt_ylim)
    plt.xlabel(x_column)
    plt.ylabel(y_column)
    ax = plt.gca()
    ax.set_facecolor(hist_bg_color)
    ax.set_xticks(np.arange(plt_xlim[0], plt_xlim[1], 0.1), minor=True)
    ax.set_yticks(np.arange(plt_ylim[0], plt_ylim[1], 0.1), minor=True)
    ax.grid(which="minor", alpha=0.2)
    plot_helper_fns.draw_lines(singlet_poly, "red", ax=ax)
    # save the scatter diagram
    scatterFigName = plotNamePartialPath + "_HIPI.png"
    plt.savefig(scatterFigName)
    plt.show()

    # plot histogram of log10(R1 647-H).
    plt.hist(df_is_singlet['log10(R1 647-H)'].values, bins = 90, color = 'blue', range=R1_647_2D_hist_range[0])            
    plt.hist(df_hot_spot['log10(R1 647-H)'].values, bins = 90, color = 'yellow', range=R1_647_2D_hist_range[0])  
    #plt.hist(df_hot_spot_not_singlet['log_R1 647-H'].values, bins = 90, color = 'orange', range=R1_647_2D_hist_range[0])
    ax = plt.gca()
    ax.set_facecolor(hist_bg_color)
    plt.title(filename + ': Hist. of log(R1 647-H)')
    plt.xlabel('log(R1_647-H)')
    plt.ylabel("Count")
    # save the histogram
    histFigName = plotNamePartialPath + '_hist.png'
    plt.savefig(histFigName)
    plt.show()        

    # plot FSC vs SSC
    x_column = "log10(FSC-H)"
    y_column = "log10(FSC-A)"
    plt.scatter(df[x_column].values, df[y_column].values, color = 'blue', s=1)
    plt.scatter(df_hot_spot[x_column].values, df_hot_spot[y_column].values, color = 'yellow', s=1)
    plt.title(filename + ': log10(SSC-H) vs log10(FSC-H)')
    ax = plt.gca()
    ax.set_facecolor(hist_bg_color)
    ax.set_xticks(np.arange(4.5, 7.5, 0.1), minor=True)
    ax.set_yticks(np.arange(4.0, 7.0, 0.1), minor=True)
    ax.grid(which='minor', alpha=0.2)
    lasso=patches.Ellipse(xy=log_FSC_H_vs_log_SSC_H_ellipse['center'], width=log_FSC_H_vs_log_SSC_H_ellipse['width'], height=log_FSC_H_vs_log_SSC_H_ellipse['height'], linewidth=1, color="red", edgecolor="red", fill=False, zorder=2)
    ax.add_artist(lasso)
    plt.grid(True)
    plt.xlabel('log(FSC-H)')
    plt.ylabel('log(SSC-H)')
    plt.xlim([4.5,7.5])
    plt.ylim([4,7])          
    # save the scatter diagram
    scatterFigName = plotNamePartialPath + '_FSC_SSC_scatter.png'
    plt.savefig(scatterFigName)
    plt.show()  
    
    # plot 2d histogram of FSC vs SSC
    x_column = "log10(FSC-H)"
    y_column = "log10(FSC-A)"
    plt.hist2d(df_is_singlet[x_column].values, df_is_singlet[y_column].values, bins=(335,218), range=[[4.5,7.5], [4,7]],cmap=plt.cm.coolwarm)
    plt.title(filename + ': 2D Hist of log10(SSC-H) vs log10(FSC-H)')
    plt.xlabel(x_column)
    plt.ylabel(y_column)        
    # save the 2d histogram
    histFigName = plotNamePartialPath + '_SSC_vs_FSC_hist.png'
    plt.savefig(histFigName)
    plt.show()  

    print("")
