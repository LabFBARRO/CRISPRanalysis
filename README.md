# CRISPRanalysis
InDels analysis of CRISPR lines by NGS amplicon sequencing technology for a multicopy gene family.

In this work, we present a workflow to analyze InDels from the multicopy α-gliadin gene family from wheat based on NGS data without the need to pre-viously establish a reference sequence for each genetic background. The pipeline was tested it in a multiple sample set, including three generations of edited wheat lines (T0, T1, and T2), from three different backgrounds and ploidy levels (hexaploid and tetraploid). Implementation of Bayesian optimization of Usearch parameters, inhouse Python, and bash scripts are reported.


# **Workflow:**

## ***Step1:***

Bayesian optimization was implemented to optimize Usearch v9.2.64 parameters from merge to search steps for the α-gliadin amplicons on the wild type lines.
````
Step1_Bayesian_usearch.py --database <str> --file_intervals <str> --trim_primers <str> --path_usearch_control <str>
````
\
*Help:*

--database	File fasta with database sequences. Example: /path/to/database/database.fasta.\
--file_intervals	File with intervals for parameters. Example in /Examples/Example_intervals.txt.\
--trim_primers	Trim primers in reads if you use database without primers. Optios: YES | NO.\
--path_usearch_control	Path of usearch and control raw data separated by \",\" without white spaces. Example: /paht/to/usearch,/path/to/reads_control.

\
*Outputs:*

* Bayesian_usearch.txt	File with optimal values, optimal function value, samples or observations, obatained values and search space.
* Bayesian.png	Convergence plot.
* Bayesian_data_res.txt	File with the minf(x) after n calls  in each iteration.

## ***Step 2:***

Usearch pipeline optimazed on wild type lines for studying results of optimization.
````
Step2_usearch_WT_to_DB.sh dif pct maxee amp id path_control name_dir_usearch path_database trim_primers
````
\
*Help:*

> Arguments must be disposed in the order indicated before.
* dif	Optimal value for dif Usearch parameter.
* pct	Optimal value for pct Usearch parameter.
* maxee	Optimal value for maxee Usearch parameter.
* amp	Optimal value for amp Usearch parameter.
* id	Optimal value for id Usearch parameter.
* path_control	Path of the wild type lines fastq files.
* name_dir_usearch	Path of Usearch.
* path_database	Path of alpha-gliadin amplicon database.
* trim_primers	Trim primers in reads if you use database without primers. Optios: YES | NO.

\
*Outputs:*

Usearch merge files, filter files, unique amplicons file, unique denoised amplicon (Amp/otu) file, otu table file.

## ***Step 3:***

Usearch pipeline optimazed on all lines (wild types and CRISPR lines) for studying denoised unique amplicon relative abundances.
```
Step3_usearch_ALL_LINES.sh dif pct maxee amp id path_ALL name_dir_usearch trim_primers
```
\
*Help:*

> Arguments must be disposed in the order indicated before.
* dif	Optimal value for dif Usearch parameter.
* pct	Optimal value for pct Usearch parameter.
* maxee	Optimal value for maxee Usearch parameter.
* amp	Optimal value for amp Usearch parameter.
* id	Optimal value for id Usearch parameter.
* path_ALL	Path of all lines (wild type and CRISPR lines) fastq files.
* name_dir_usearch	Path of Usearch.
* trim_primers	Trim primers in reads if you use database without primers. Optios: YES | NO.

\
*Outputs:*

Usearch merge files, filter files, unique amplicon file, unique denoised amplicon (Amp/otu) file, otu table file.

> Before Step 4, otu table file must be normalized by TMM normalization method (edgeR package in R).
Results of TMM normalized unique denoised amplicons table can be represented as heatmaps.
Unique denoised amplicons can be compared between them to detect Insertions and Deletions (InDels) in CRISPR lines.


## ***Step 4:***

Create tables with the presence or absence of unique denoised amplicons in each CRISPR line compared to the wild type lines.
```
python Step4_usearch_to_table.py --file_otu <str> --file_group <str> --prefix_output <str> --genotype <str>
```

\
*Help:*

--file_otu	File of TMM normalized otu_table from usearch. Remove \"#OTU\" from the first line.\
--file_group	Path to file of genotypes in wild type and CRISPR lines. Example in /Examples/Example_groups.txt.\
--prefix_output	Prefix to output name. Example: if you are working with BW208 groups: BW.\
--genotype	Genotype name. Example: if you are working with BW208 groups: BW208.\
Default threshold 0.3 % of frequency of each unique denoised amplicon (Amp) in each line.

\
*Outputs:*

> Substitute "name" in output names for the prefix_output string.
* Amptable_frequency.txt	Table of Amps (otus) transformed to frequencies for apply the threshold.
* Amptable_brutes_name.txt	Table with number of reads contained in the unique denoised amplicons (Amps) present in each line.
* Amps_name.txt	Table with number of unique denoised amplicons (Amps) in each line.

> Python 3.6 or later is required.
