#!/bin/bash


##
## SEACR-SPIKEIN NORMALIZIED peak calling
##


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
route_name=${script_name/%.sh/}
echo -e "\n ========== ROUTE: $route_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 3 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir treatment_sample_name control_sample_name \n" >&2
	exit 1
fi

# standard comparison route arguments
proj_dir=$(readlink -f "$1")
sample_treatment=$2
sample_control=$3

# paths
code_dir=$(dirname $(dirname "$script_path"))
sbatch_dir="${proj_dir}/logs-sbatch"

# display settings
echo
echo " * proj_dir: $proj_dir "
echo " * treatment sample: $sample_treatment "
echo " * control sample: $sample_control "
echo " * code_dir: $code_dir "
echo

# specify maximum runtime for sbatch job
# SBATCHTIME=6:00:00


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: dir $proj_dir does not exist \n" >&2
	exit 1
fi

# get deduplicated BAMs corresponding to sample names
bed_treatment="${proj_dir}/BEDS/${sample_treatment}_bowtie2.fragments.bedgraph"
bed_control="${proj_dir}/BEDS/${sample_control}_bowtie2.fragments.bedgraph"


bam_treatment=$(grep -s -m 1 "^${sample_treatment}," "${proj_dir}/samples.${segment_dedup}.csv" | cut -d ',' -f 2)

# check if the treatment BAM exists
if [ ! -s "$bed_treatment" ] ; then
	echo -e "\n $script_name ERROR: treatment BAM $bam_treatment does not exist \n" >&2
	exit 1
fi

# if the treatment and control samples are the same, ignore the control sample
if [ "$bed_treatment" == "$bed_control" ] ; then
	bam_control=""
fi


#########################

module purge
module add default-environment
module add r/4.2.2
module add bedtools/2.30.0

# the only segment so no need to change yet.


seacr="/gpfs/data/feskelab/McDermott/Shared/SEACR/SEACR_1.3.sh"
SEAR_DIR="${proj_dir}/SEACR"
mkdir -p $SEAR_DIR

bash_cmd="
 $seacr  ${bed_treatment} \
	 ${bed_control} \
     non stringent \
	 ${SEAR_DIR}/${sample_treatment}_${sample_control}_seacr_control.peaks"


echo -e "\n CMD: $bash_cmd \n"
$bash_cmd


#########################


# check that output generated

if [ ! -s "${SEAR_DIR}/${sample_treatment}_${sample_control}_seacr_control.peaks" ] ; then
	echo -e "\n $script_name ERROR: peaks 	 ${SEAR_DIR}/${sample_treatment}_${sample_control}_seacr_control.peaks.stringent.bed not generated \n" >&2
	exit 1
fi

#########################


date



# end
