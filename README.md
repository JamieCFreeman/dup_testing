# 2022-1-6 JCF

# For collecting read duplicate information and summary stats from seq data.

# To run on marula, run in conda env "snakemake"

# Snakemake will use the yaml files in the "envs" folder to handle other dependencies"

# To run use:
"snakemake -c n_cores --useconda"

# Need to make sure to make input files in directory:
#       "config.yaml"
#       "units.tsv" (sample table)
#       as well as directory "envs" with includes specs for conda envs snakemake will use to execute
#  Maps w/ bwa mem, marks duplicates, and calculates some summary stats on different size partitions of the libraries
#       (b/c I have noticed that in our Tn5 library preps small fragments are the source of more than their share of duplicates


