#!/bin/bash

# Usearch pipeline for all the lines with optimized parameters from Step 1.

Dire=$(pwd) #Get the path of actual directory

#script output.txt #Para capturar toda la salida de la pantalla



#Parameters: optimized from Step 1.
dif=$1
pct=$2
maxee=$3
amp=$4
id=$5
path_ALL=$6
name_dir_usearch=$7
trim_primers=$8

name_dir="resultsStep3_${dif}_${pct}_${maxee}_${amp}" #Results directory name.
mkdir $name_dir #Create results directory.

#Usearch pipeline: merge.
echo "MERGING"
	$name_dir_usearch/usearch9 -fastq_mergepairs "$path_ALL"*_R1*.fastq -relabel @ -fastq_maxdiffs $dif -fastq_maxdiffpct $pct -fastqout "$name_dir"/readsnotrucanted-"$dif".fq
if($trim_primers == "YES"){
	$name_dir_usearch/usearch9 -fastx_truncate "$name_dir"/readsnotrucanted-"$dif".fq -stripleft 21 -stripright 20 -fastqout "$name_dir"/reads-"$dif".fq
}
#Usearch pipeline: filter.
echo "FILTERING"
if($trim_primers == "YES"){
	$name_dir_usearch/usearch9 -fastq_filter "$name_dir"/reads-"$dif".fq -fastq_maxee "$maxee" -fastaout "$name_dir"/filtered-"$dif-$maxee".fa
} else{
	$name_dir_usearch/usearch9 -fastq_filter "$name_dir"/readsnotruncated-"$dif".fq -fastq_maxee "$maxee" -fastaout "$name_dir"/filtered-"$dif-$maxee".fa
}
#Usearch pipeline: uniques.
echo "UNIQUES"
	$name_dir_usearch/usearch9 -fastx_uniques "$name_dir"/filtered-"$dif-$maxee".fa -fastaout "$name_dir"/uniques2-"$dif-$maxee".fa -sizeout
#Usearch pipeline: unoise.
echo "UNOISE"
	$name_dir_usearch/usearch9 -unoise2 "$name_dir"/uniques2-"$dif-$maxee".fa -fastaout "$name_dir"/denoised-"$dif-$maxee-$amp".fa -minampsize $amp -log "$name_dir"/log-"$dif-$maxee-$amp".txt
#Usearch pipeline: searching.
echo "SEARCHING"
	$name_dir_usearch/usearch9 -search_global "$name_dir"/reads-"$dif".fq -db "$name_dir"/denoised-"$dif-$maxee-$amp".fa -id $id -strand both -otutabout "$name_dir"/otutable_TODOS-"$dif-$maxee-$amp".txt
