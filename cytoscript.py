# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 10:26:01 2019
Has functions that simplify the task of writing scripts for processing flow cytometry data

@author: mikeh
"""

# libraries
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
import os
from FlowCytometryTools import FCMeasurement

class CytoScript:
    # just initialize the internal dataframe dictionary
    def __init__(self):
        self.df_dict = {}
        self.selection_rules = {}
        self.sampleID = ""
        
    # functions to impersonate a dictionary
    def __setitem__(self, key, item):
        self.df_dict[key] = item

    def __getitem__(self, key):
        return self.df_dict[key]

    def __repr__(self):
        return repr(self.df_dict)

    def __len__(self):
        return len(self.df_dict)

    def __delitem__(self, key):
        del self.df_dict[key]

    def clear(self):
        return self.df_dict.clear()

    def copy(self):
        return self.df_dict.copy()

    def has_key(self, k):
        return k in self.df_dict

    def update(self, *args, **kwargs):
        return self.df_dict.update(*args, **kwargs)

    def keys(self):
        return self.df_dict.keys()

    def values(self):
        return self.df_dict.values()

    def items(self):
        return self.df_dict.items()

    def pop(self, *args):
        return self.df_dict.pop(*args)

    def __cmp__(self, dict_):
        return self.__cmp__(self.df_dict, dict_)

    def __contains__(self, item):
        return item in self.df_dict

    def __iter__(self):
        return iter(self.df_dict)

        

    # getDF_colNames returns a list containing names of columns that match colNameRegEx 
    def getDF_colNames(df, colNameRegEx):
        return list(df.filter(regex=colNameRegEx).columns)

#%% Transformation oriented functions
    # log10 calculates the log, base 10, of the named column
    # the function acts on the named dataframe, and it created a new column name, as specified by newColumnName.
    def log10(self, colNameRegEx = '\w+-(A|W|H)$', dfName='df', removeNA=True):
        """
        log10 calculates the log, base 10, for a designated column, in a designated dataframe. The call will alter that dataframe by adding a column whose name is
        "log10(original_columnName)", and that column will contain the log10 value for each corresponding value in the "original_columnName" column.
        
        :param colNameRegEx: Regular expression for column whose log10 will be calculated. If colName is ".", then we will act on all columns. Use a string like '\AFITC-H\Z' if you want an exact name match
        :param dfName: The name of the dataframe in the dataframe dictionary that we will act on. If no name is provided, we act on default dataframe "df"
        :param removeNA: Remove any rows that have NA values afer the operation
        :return: returns nothing
        """
        # get the desginated dataframe
        df = self.df_dict[dfName]
        
        # get the names of designated columns from the dataframe
        # remove the Time column if it is present .. there is no point calculating the log of that
        dfColNames = list(df.filter(regex=colNameRegEx).columns)
        
        # iterate over the columns, applying the np.log10 function
        for colName in dfColNames:
            newColName = "log10(" + colName + ")"
            df[newColName] = np.log10(df[colName].values)
        
        # if the removeNA flag is true, then prune any rows with NA values
        if removeNA:
            df.dropna(inplace=True)
        
    # setTasks loads a local list of tasks that will subsequently be executed by runTasks
    def setTasks(self, tasks):
        self.tasks = tasks
        
    # setTasks loads a local list of tasks that will subsequently be executed by runTasks
    def setScript(self, script):
        self.script = script
        
        
#%% File oriented methods  
    # ensureFolderExists looks in the folder specified by the parentFolder argument to see if the child-folder specified by the folder argument exists
    # if it doesn't exist, the specified folder will be created
    # ensureFolderExists  then returns the full path to the specified child folder
    def ensureFolderExists(parentFolder, folder):
        fullFolderPath = parentFolder + "/" + folder
        exists = os.path.isdir(fullFolderPath)
        if exists == False:
            os.mkdir(fullFolderPath)
        return fullFolderPath 
    
    
    def ensureSubFolderExists(self, subFolderName):
        return CytoScript.ensureFolderExists(self.workingDir, subFolderName)

    # filesInWorkingDir returns  list of files in the working dir
    def workingDirFiles(self, sort=True, only_csv_and_fcs = True):
        files = os.listdir(self.workingDir)

        # process in alphanumeric order, if so specified
        if sort:
            files.sort()
        
        # process only fcs or csv files, if so specified
        if only_csv_and_fcs:
            files = [fi for fi in files if fi.lower().endswith(".csv") | fi.lower().endswith(".fcs")]
            
        return files 
        
    # set the working dorectory
    def setWorkingDir(self, pathToWorkingDir):
        self.workingDir = pathToWorkingDir

    # get the stirng for the working dorectory    
    def getWorkingDir(self):
        return self.workingDir
    
    # given a file name, so a spart concatenation of the working directory to determine the actual full path to the file
    def fullFileName(self,fileName):
        return os.path.join(self.workingDir, fileName)
    
    # get the portion of an fcs or csv file that is to the left of the "." .. normally to use it as a sample  ID
    def getSampleIDFromFileName(self, fileName):
        self.sampleID = fileName.split("\\")[-1].split("/")[-1].split(".")[0]
        return self.sampleID
    
    # load the designated CSV file. The CSV is presumed to be in the current wotrking irector
    def loadCSV(self, csvFileName):
        self.df_dict = {}
        self.getSampleIDFromFileName(csvFileName)
        df=pd.read_csv(self.fullFileName(csvFileName))
        self.df_dict['df'] = df
        
    # load the designated fcs file
    def loadFCS(self, fcsFileName, apply_col_rename=True):

        # clear out the dataframe dictionary
        self.df_dict = {}

        # use the file name as a sample id .. not necessarily of any significance
        self.getSampleIDFromFileName(fcsFileName)
        
        # read in the fcs file
        sample = FCMeasurement(ID=self.sampleID, datafile=self.fullFileName(fcsFileName))
        
        # If apply_col_rename is True, then check apply renaming rules for any renamed columns
        rename_dict={}
        if apply_col_rename:
            columns = list(sample.data.columns)
            for i in range(0,len(columns)):
                key="$P{}S".format(i+1)
                if key in sample.meta:
                   rename_dict[columns[i]] = sample.meta[key]
        
        # place the dataframe in the dataframe dictionary, using the default dataframe name 'df' as the key.
        self.df_dict['df'] = sample.data.rename(index=str, columns=rename_dict)
        
    # load_csv_or_fcs will load either an fcs or a csv
    # it is intended for iterating over a folder that has multiple flow cytometry files
    # that are either in csv of fcs format
    def load_csv_or_fcs(self, fileName):
        retval = True
        if fileName.endswith(".csv"):
            self.loadCSV(fileName)
        elif fileName.endswith(".fcs"):
            self.loadFCS(fileName)
        else:
            retval = False;
        return retval

#%% Gate oriented functions    
    # applyEllipticalGate accepts an array of x,y pairs and an ellipse, and returns an array of booleans that tells which points are in the ellipse, and which are not
    # The points argument is an array of x, y pairs
    # the ellipse argument has the following form {"xy":(6.25, 5.85), "width":0.5, "height":0.6, angle=0.0} # where the angle argument is optional
    def applyEllipticalGate(points, ellipse):
        # create a matplotlib.patches.Ellipse object
        if 'angle' in ellipse.keys():
            ell = patches.Ellipse(xy=ellipse['xy'], width=ellipse['width'], height=ellipse['height'], angle=ellipse['angle'])
        else:
            ell = patches.Ellipse(xy=ellipse['center'], width=ellipse['width'], height=ellipse['height'])
            
        # create an array of booleans that tell whether each point in p is in the ellipse
        p_in = np.array([ell.contains_point(pi) for pi in points])
        return p_in
    
    # calcEllipticalGate takes the columns designated by xcol_name and ycol_name in the dataframe designated by the source_df_name argument
    # It applys the elliptical gate specified by ellipse argument, and deposits the Boolean array result in the dataframe designated by the dest_df_name argument
    # This result is placed in a column whose name is given by the result_name argument
    # If source_df_name matches dest_df_name, then the modification is made in place.
    # the function returns the array of booleans that tell whether the corrsponding row of the source array is in the ellipse or outside of the ellipse
    def calcEllipticalGate(self, xcol_name, ycol_name, ellipse, result_name, source_df_name='df', dest_df_name='df' ):
        source_df = self.df_dict[source_df_name]
        dest_df = self.df_dict[dest_df_name]
        points = source_df[[xcol_name, ycol_name]].values
        result = CytoScript.applyEllipticalGate(points, ellipse)
        dest_df[result_name] = result
        return result
    
    # applyPolygonGate accepts an array of x,y pairs and an (N, 2) array-like definition of a polygon and calculates a Boolean
    # vector that tells which x,y pairs are in the polygon, and which are outside of it
    def applyPolygonGate(points, polygon):
        # create a matplotlib.patches.Polygon object
        poly = patches.Polygon(xy=polygon)
            
        # create an array of booleans that tell whether each point in p is in the ellipse
        p_in = np.array([poly.contains_point(pi) for pi in points])
        return p_in    
 
    # calcPolygonGate takes the columns designated by xcol_name and ycol_name in the dataframe designated by the source_df_name argument
    # It applys the elliptical gate specified by ellipse argument, and deposits the Boolean array result in the dataframe designated by the dest_df_name argument
    # This result is placed in a column whose name is given by the result_name argument
    # If source_df_name matches dest_df_name, then the modification is made in place.
    # the function returns the array of booleans that tell whether the corrsponding row of the source array is in the ellipse or outside of the ellipse
    def calcPolygonGate(self, xcol_name, ycol_name, polygon, result_name, source_df_name='df', dest_df_name='df' ):
        source_df = self.df_dict[source_df_name]
        dest_df = self.df_dict[dest_df_name]
        points = source_df[[xcol_name, ycol_name]].values
        result = CytoScript.applyPolygonGate(points, polygon)
        dest_df[result_name] = result
        return result

#%% selection rules
    
    def addSubsetRule(self, ruleName, ruleText):
        self.selection_rules[ruleName] = ruleText
        
    def getSubsetRule(self, ruleName):
        return self.selection_rules[ruleName]
    
    def getSubSetRules(self):
        return self.selection_rules
    
    def subSetRuleEvalText(self, ruleName, dfName = 'df'):
        ruleText = self.selection_rules[ruleName]
        dfref = "self.df_dict['{}']".format(dfName)
        innerRule = ruleText.replace("]","']").replace("[", dfref + "['")
        evalRuleText = dfref + "[" + innerRule + "]"
        return evalRuleText
    
    def getSubSet(self, ruleName):
        return eval(self.subSetRuleEvalText(ruleName))
        
        
    

        
        
        
        
    