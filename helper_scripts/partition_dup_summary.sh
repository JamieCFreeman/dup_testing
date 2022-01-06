#!/bin/bash

# 2021-9-10 JCF

# Project:
# Tn5 library duplicate testing

# Purpose:
# Have separately summarized PCR duplication statistics for fragment sizes <150 & >150.
# Now in samtools flagstat output format over multiple files, want to combine into a more useful table.

# Variables
# Input is a directory, containing a samtools flagstat file following the naming scheme "*_le150.stats" (fragments<=150)
#	and "*_ge151.stats" (fragments >=151)
DIR_IN="./dup_test/logs/dupl_partition"
DIR_OUT="./dup_test/logs/summary_stats" 

# Add header to file
printf "Sample\tTotal_below_150\tDuplicates_below_150\tTotal_above_150\tDuplicates_above_150\n" \
	> ${DIR_OUT}/dupl_partition_summary.txt
# Parse samtools flagstat files for duplicate counts
paste \
	<( find ./dup_test/logs/dupl_partition -iname "*_le150.stats" | sort -n -t - -k2  | 
		xargs grep "total" | sed 's/ +.*//' | sed 's/_le150.stats.*//' | sed 's ^.*/  ' ) \
	<( find ${DIR_IN} -iname "*_le150.stats" | sort -n -t - -k2  | 
		xargs grep "total" | sed 's/ +.*//' | sed 's/^.*le150.stats://' ) \
	<( find ${DIR_IN} -iname "*_le150.stats" | sort -n -t - -k2  | 
		xargs grep "duplicates" | sed 's/ +.*//' | sed 's/^.*le150.stats://' ) \
	<( find ${DIR_IN} -iname "*_ge151.stats" | sort -n -t - -k2  | 
		xargs grep "total" | sed 's/ +.*//' | sed 's/^.*ge151.stats://' ) \
	<( find ${DIR_IN} -iname "*_ge151.stats" | sort -n -t - -k2  | 
		xargs grep "duplicates" | sed 's/ +.*//' | sed 's/^.*ge151.stats://' ) >> ${DIR_OUT}/dupl_partition_summary.txt


