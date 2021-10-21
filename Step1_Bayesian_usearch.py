# -*- coding: utf-8 -*-
"""
Python: Bayesian Optimization for Usearch parameters values selection
Usearch v9

@author: mmarin@ias.csic.es
"""
# BAYESIAN OPTIMIZATION FOR USEARCH PARAMETERS

import numpy as np
from skopt import gp_minimize
from skopt.space import Real, Categorical
from skopt.utils import use_named_args
import os
import time
import argparse


# ARGUMENTS

parser = argparse.ArgumentParser()
parser.add_argument("--database", help="File fasta with database sequences. Example: /path/to/database/database.fasta.")
parser.add_argument("--file_intervals", help="File with intervals for parameters. Example in Example_intervals.txt.")
parser.add_argument("--trim_primers", help="Trim primers in reads if you use database without primers. Optios: YES | NO.")
parser.add_argument("--path_usearch_control", help="Path of usearch and control raw data separated by \",\" without white spaces. Example: /paht/to/usearch,/path/to/reads_control.")
args = parser.parse_args()

#OPEN FILE WITH INTERVALS OF VALUES FOR EACH PARAMETER

try:
    file_intervals = open(args.file_intervals, "r")
except:
    print("Error: The file " + args.file_intervals + " could not be openned.")

# READ THE FILE

try:
    file_intervals_read = file_intervals.readlines()
except:
    print("Error: The file " + args.file_intervals + " could not be read.")

# CONSTRUCTION OF INTERVALS DICTIONARY

try:
    intervals = {}

    for line in file_intervals_read:
        name = line.split("\t")[0].strip()
        value1 = float(line.split("\t")[1].split(",")[0].strip())
        value2 = float(line.split("\t")[1].split(",")[1].strip())
        intervals[name] = (value1, value2)
except:
    print("Error: The intervals dictionary could not be constructed.")

# CONSTRUCTION OF SPACE FOR BAYESIAN OPTIMIZATION

try:
    space  = [Real(intervals['diff'][0], intervals['diff'][1], name='diff'),
              Real(intervals['pct'][0], intervals['pct'][1], name='pct'),
              Real(intervals['maxee'][0], intervals['maxee'][1], name='maxee'),
              Real(intervals['amp'][0], intervals['amp'][1], name='amp'),
              Categorical([intervals['id'][0], intervals['id'][1]], name = 'id')
              ]
except:
    print("Error: The space for Bayesian Optimization could not be constructed.")

# PATH FOR USEARCH AND CONTROL RAW DATA

PATHUSEARCH = args.path_usearch_control_database.split(",")[0].strip()
PATHCONTROL = args.path_usearch_control_database.split(",")[1].strip()

# THE FUNCTION TO OPTIMIZE

class Usearch:
    def __init__(self, database):
        self.database = database
        pass

    def set_params(self, **params): #**params
        self.diff = params["diff"]
        self.pct = params["pct"]
        self.maxee = params["maxee"]
        self.amp = params["amp"]
        self.id = params["id"]

    def predict(self):
        start = time.time()
        #merge command.
        os.system(PATHUSEARCH + "/usearch9 -fastq_mergepairs " + PATHCONTROL + "/*_R1*.fastq -relabel @ -fastq_maxdiffs " + str(self.diff) + " -fastq_maxdiffpct " + str(self.pct) + " -fastqout merge.fq")
        if args.trim_primers == "YES":
            #Remove the first and last nuecleotides assuming that they are primer sequences.
            os.system(PATHUSEARCH + "/usearch9 -fastx_truncate merge.fq -stripleft 21 -stripright 20 -fastqout merge_strip.fq")
        #filter command.
        if args.trim_primers == "YES":
            file_merge = "merge_strip.fq"
        else:
            file_merge = "merge.fq"
        os.system(PATHUSEARCH + "/usearch9 -fastq_filter " + file_merge + " -fastq_maxee " + str(self.maxee) + " -fastaout filter.fa")
        #uniques command.
        os.system(PATHUSEARCH + "/usearch9 -fastx_uniques filter.fa -fastaout unicos.fa -sizeout")
        #unoise command.
        os.system(PATHUSEARCH + "/usearch9 -unoise2 unicos.fa -fastaout unoise.fa -minampsize " + str(self.amp) + " -unoise_alpha 3 -log stat_log_unoise.txt")
        #cluster command.
        os.system(PATHUSEARCH + "/usearch9 " + str(self.id) + " -db " + self.database + " -strand both -otutabout cluster_database.txt")

        #number of otus clustered in database.
        try:
            file = open("cluster_database.txt", "r")
            file_read = file.readlines()
        except:
            print("Error: The cluster database file could not be openned or read.")
        reads = 0
        for line in file_read[1:]:
            for col in range(1, len(line.split("\t"))):
                reads = reads + int(line.split("\t")[col].strip())

        #The Bayesian optimization considers the lower result as the most optimal. In our case, the best result is the highest, so we return -result.
        print("Running time: ", time.time() - start)
        return -reads


@use_named_args(space)
def objetive(**params):
    model = Usearch(args.database)
    model.set_params(**params)
    result = model.predict()
    return result

#%% OPTIMIZATION

res = gp_minimize(objetive,                # the function for the optimization
                  space,      # limits for each dimension
                  acq_func="EI",      # adquisition function
                  #acq_func="PI",
                  #acq_func="LCB",
                  n_calls=100,         # number of evaluations
                  n_random_starts=5,  # number of initial random points
                  noise=0.1**2,       # noise level
                  random_state=123)   # seed


#%% RESULTS

print("Optimal values: ", res.x)
print("Optimal function value: ", res.fun)
print("Samples or observations: ", res.x_iters) # check the length of x_iters
print("Obtained values: ", res.func_vals) # the minimum of these values ​​must match res.fun
print("Search space: ", res.space)

try:
    file_result = open("Bayesian_usearch.txt", "w+")
except:
    print("Error: The file for Bayesian results could not be created.")

file_result.write("Optimal values: " + str(res.x) + "\n")
file_result.write("Optimal function value: " + str(res.fun) + "\n")
file_result.write("Samples or observations: " + str(res.x_iters) + "\n")
file_result.write("Obtained values: " + str(res.func_vals) + "\n")
file_result.write("Search space: " + str(res.space) + "\n")


#%% CONVERGENCY

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

plt.figure(figsize=(5,5))
from skopt.plots import plot_convergence
plot_convergence(res)
plt.savefig("Bayesian.png")

#For the representation of convergence plot in R.
n_calls = len(res.x_iters)
mins = [np.min(res.func_vals[:i]) for i in range(1, n_calls + 1)]

file_convergency = open("Bayesian_data_res.txt", "w+")
file_convergency.write("Number of calls n" + "\t" + "minf(x) after n calls" + "\n")
for n in range(1, n_calls + 1):
    file_convergency.write(str(n) + "\t" + str(mins[n - 1]) + "\n")
