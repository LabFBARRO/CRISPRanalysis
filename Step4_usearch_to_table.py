# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 12:39:44 2020

@author: mmarin@ias.csic.es
"""


#For Step 4 of CRISPR analysis by amplicon technology. Create tables from TMM normalized otu_table.txt of usearch pipeline in Step 3.

import pandas as pd
import argparse
import os
import statistics

file_path = os.path.realpath(__file__)

parser = argparse.ArgumentParser()
parser.add_argument("--file_otu", help="File of TMM normalized otu_table from usearch. Remove \"#OTU\" from the first line.")
parser.add_argument("--file_group", help="Path to file of genotypes in WT and CRISPR. Example in Example_groups.txt.")
parser.add_argument("--prefix_output", help="Prefix to output name. Example: if you are working with BW208 groups: BW.")
parser.add_argument("--genotype", help="Genotype name. Example: if you are working with BW208 groups: BW208.")
args = parser.parse_args()

#%%
#Unique denoised amplicon table to frequencies.

#*********************VARIABLE**********************
#Read unique amplicons table with TMM normalization and keep it as dataframe
file_otu = pd.read_csv(args.file_otu, sep="\t", header=1)

#Guarda el numero de reads por genotipo
dict_reads_total_new = {}
#Create a dataframe to save the frequencies
new_file_otu = pd.DataFrame()

for column in file_otu:
    new_column = []
    total_reads = file_otu[column].sum()
    new_column = (file_otu[column]/total_reads)*100
    new_file_otu[column] = new_column
    dict_reads_total_new[column] = total_reads

#*********************VARIABLE**********************
#Prefix to output name and genotype name
name = args.prefix_output
genotype = args.genotype

#Frequency threshold. Default 0.3.
threshold = 0.3

new_file_otu.to_csv(file_path + "Amptable_frequency.txt", index=True, sep='\t', header=True)

count = 0
for index, row in new_file_otu.iterrows():
    for col in row:
        if col < 0.3:
            continue
        else:
            print(index)
            count += 1
            break

#%%
#Controls vs CRISPR lines.
#Amp (Otu) WT no CRTISPR: Amp (Otu) in Control (frequency > 0.3) not in CRISPR line.
OtuWTnoCRISPR = []
#Amp (Otu) CRISPR no WT: Amps (Otu) in CRISPR line (frequency > 0.3) not in WT.
OtuCRISPRnoWT = []

#*********************VARIABLE**********************
#Groups file.
#You have to create a file like Example_groups (in the Examples folder), which has information on which plants of the group you want to analyze are control and which are crispr.
group_file = open(args.file_group + "groups_" + genotype + ".txt", "r")
group_file_read = group_file.readlines()

dict_group = {}

for line in group_file_read:
    if line.split("\t")[0].strip() == "WT":
        WT_list = []
        for WT in line.split("\t")[1].split(","):
            WT_list.append(WT.strip())
        dict_group["WT"] = WT_list
    else:
        CRISPR_list = []
        for CRISPR in line.split("\t")[1].split(","):
            CRISPR_list.append(CRISPR.strip())
        dict_group["CRISPR"] = CRISPR_list

number_OTU_expressed_in_WT = 0
reads_expressed_in_WT = 0
#key: genotype, value: number of Amps (otus) expressed in genotype (frequency > 0.1) and number of reads in these Amps (otus).
dict_otus = {}

dict_otusWTnoCRISPR = {}
dict_otusWTnoCRISPRfreq = {}
dict_otusCRISPRnoWT = {}
dict_otusCRISPRnoWTfreq = {}

for index, row in new_file_otu.iterrows():
    expressed_in_WT = "NO"
    #Max value and line with max value.
    frequency_WT = (0.0, "")
    #If Amp (OTU) expressed in WT, keep max value for this Amp (OTU) among all WT.
    for WT in WT_list:
        if row[WT] > threshold:
            expressed_in_WT = "YES"
            if frequency_WT[0] < row[WT]:
                frequency_WT = (row[WT], WT)
    if expressed_in_WT == "YES":
        dict_otus["WT"] = dict_otus.get("WT", 0) + 1
        frequency_WT_list = []
        for WT in WT_list:
            frequency_WT_list.append((row[WT]))

    expressed_in_CRISPR = "NO"
    for CRISPR in CRISPR_list:
        if row[CRISPR] > threshold:
            expressed_in_CRISPR = "YES"
            frequency_CRISPR = row[CRISPR]
            dict_otus[CRISPR] = dict_otus.get(CRISPR, 0) + 1
            if expressed_in_WT == "NO":
                dict_otusCRISPRnoWT[CRISPR] = dict_otusCRISPRnoWT.get(CRISPR, 0) + 1
                dict_otusCRISPRnoWTfreq[CRISPR] = dict_otusCRISPRnoWTfreq.get(CRISPR, 0) + row[CRISPR]
        else:
            if expressed_in_WT == "YES":
                dict_otusWTnoCRISPR[CRISPR] = dict_otusWTnoCRISPR.get(CRISPR, 0) + 1
                dict_otusWTnoCRISPRfreq[CRISPR] = dict_otusWTnoCRISPRfreq.get(CRISPR, 0) + statistics.mean(frequency_WT_list)

    if expressed_in_WT == "YES" and expressed_in_CRISPR == "NO":
        OtuWTnoCRISPR.append(index)
    if expressed_in_WT == "NO" and expressed_in_CRISPR == "YES":
        OtuCRISPRnoWT.append(index)

#Comprobar que todos los CRISPR están en todos los dict.
for CRISPR in CRISPR_list:
    if CRISPR not in dict_otus:
        dict_otus[CRISPR] = 0
    if CRISPR not in dict_otusWTnoCRISPR:
        dict_otusWTnoCRISPR[CRISPR] = 0
    if CRISPR not in dict_otusCRISPRnoWT:
        dict_otusCRISPRnoWT[CRISPR] = 0
    if CRISPR not in dict_otusWTnoCRISPRfreq:
        dict_otusWTnoCRISPRfreq[CRISPR] = 0
    if CRISPR not in dict_otusCRISPRnoWTfreq:
        dict_otusCRISPRnoWTfreq[CRISPR] = 0

#%%
#For write file with results.

pd_otu = pd.DataFrame([dict_otus, dict_otusWTnoCRISPR, dict_otusWTnoCRISPRfreq, dict_otusCRISPRnoWT, dict_otusCRISPRnoWTfreq])
pd_otu = pd_otu.transpose()
pd_otu.columns = ["Amps presented", "Amps present in WT no CRISPR", "Accumulated frequency of Amps present in WT not in CRISPR", "Amps present in CRISPR not in WT", "Accumulated frequency of Amps present in CRISPR not in WT"]
pd_otu.to_csv(file_path + "Ampstable_brutes_" + name + ".txt", index=True, sep='\t', header=True)

pd_otus = pd.DataFrame([OtuWTnoCRISPR, OtuCRISPRnoWT])
pd_otus = pd_otus.transpose()
pd_otus.columns = ["Amps present in WT not in CRISPR", "Amps present in CRISPR not in WT"]
pd_otus.to_csv(file_path + "Amps_" + name + ".txt", index=None, sep='\t', header=True)
