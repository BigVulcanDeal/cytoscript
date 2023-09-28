# import the PPK_CytoScript class from PPK_flow_cytometry_dx_api
from cytoscript import CytoScript
import numpy as np
import pandas as pd

#%% define a few constants that we'll use for our dx script

# define an elliptical gate that operates on the log10(SSC-H) vs log10(FSC-H) plane. Note that log10(FSC-H) is the x axis
FSC_H_SSC_H_elliptic_gate={"xy":(6.25, 5.85), "width":0.5, "height":0.6}

# singlet_poly_gatedefines a polygon that determines if an event is a presumed singlet
# Points for which "log10(BL2 PI-H)" (x-axis) and "log10(BL2 PI-A)" (y-axis) fall inside polygon are considered singlets
# points that fall outside the polygon are presumed non-singlets
singlet_poly_gate = [
        [4.4,4.7],
        [4.54,4.9],
        [5.3,5.7],
        [5.6,5.75],
        [4.7,4.6],
        [4.4,4.7]
        ]

#%%
# create a new PPK_CytoScript object
ppkcs = CytoScript()

# navigate to a directory that has fcs files
ppkcs.setWorkingDir('examplescripts/FCS Stuff/P111 EVERGREEN #20 PLATE 1 CSV')
                    
# create a list to hold any results
results=[]

# iterate over the files in the working dir
# Note that ppkcs.workingDirFiles() will return the list of files in alphabetical order
for file in ppkcs.workingDirFiles():
    
    # If we are able to load the current file as a csv or an fcs, then process the data
    if ppkcs.load_csv_or_fcs(file):
        
        print("Processing {}".format(file))
        
        # convert all columns to log10 values, since our signals of interest are lognormal
        # NOTE: set removeNA false, to ensure consistency with current algorithm ..
        ppkcs.log10(removeNA=False)
        
        # mark all of the singlets with a flag named "is_singlet"
        ppkcs.calcPolygonGate("log10(BL2 PI-H)", "log10(BL2 PI-A)", singlet_poly_gate, "is_singlet")
        
        # mark all of the all points that fall within the FSC_H_SSC_H_elliptic_gate ith a flag named "FSC_SSC_hot_spot"
        ppkcs.calcEllipticalGate('log10(FSC-H)','log10(SSC-H)', FSC_H_SSC_H_elliptic_gate, "FSC_SSC_hot_spot")
        
        # make subset rule to put all singlets in a set called "denominator"
        ppkcs.addSubsetRule("denominator", "[is_singlet]")
        
        # male subset rule to put all events in the "FSC_SSC_hot_spot" and with ([log10(R1 647-H)] > 5.5) in a subset named "numerator"
        ppkcs.addSubsetRule("numerator", "[is_singlet] & ([log10(R1 647-H)] > 5.5) & [FSC_SSC_hot_spot]")
        
        # add a row to the cumulative results for the folder ..
        numerator = len(ppkcs.getSubSet("numerator"))
        denominator = len(ppkcs.getSubSet("denominator"))
        rfu = 0 if numerator < 1 else int(np.mean(ppkcs.getSubSet("numerator")['log10(R1 647-H)'].values) * 100)/100
            
        result_row={
                "sample_id":ppkcs.sampleID,
                "all_events":len(ppkcs['df']),
                "denominator":denominator,
                "numerator":numerator,
                "rfu":rfu,
                "ratio_pct":int(numerator/denominator*10000)/100
                }
        results.append(result_row)

    # if we generated any results, then save them in a subpfolder named "results", in a file named "results.csv"
    if len(result_row) > 0:
        results_file_name = ppkcs.ensureSubFolderExists("ppkcs_results") + "/results.csv"
        df_results=pd.DataFrame(results).reindex(columns=list(result_row.keys()))
        df_results.to_csv(results_file_name)



