#!/bin/bash


# Bismark alignment


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ $# -lt 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads FASTQ_R1 [FASTQ_R2] \n" >&2
	if [ $# -gt 0 ] ; then echo -e "\n ARGS: $* \n" >&2 ; fi
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
sample=$2
threads=$3
fastq_R1=$4
fastq_R2=$5


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$fastq_R1" ] ; then
	echo -e "\n $script_name ERROR: FASTQ $fastq_R1 DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

ref_bismark=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BISMARK);

if [ ! -d "$ref_bismark" ] || [ ! "$ref_bismark" ] ; then
	echo -e "\n $script_name ERROR: REF $ref_bismark DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

summary_dir="${proj_dir}/summary"
mkdir -p "$summary_dir"
summary_csv="${summary_dir}/${sample}.${segment_name}.csv"

samples_csv="${proj_dir}/samples.${segment_name}.csv"

bismark_bam_dir="${proj_dir}/BAM-Bismark"
mkdir -p "$bismark_bam_dir"
bismark_bam_final="${bismark_bam_dir}/${sample}.bam"

bismark_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$bismark_logs_dir"
bismark_report_final="${bismark_logs_dir}/${sample}.report.txt"
bismark_nucstats_final="${bismark_logs_dir}/${sample}.nucleotide_stats.txt"

bismark_report_dir="${proj_dir}/Bismark-report"
mkdir -p "$bismark_report_dir"

# v0.16.0: "File endings .fastq | .fq | .fastq.gz | .fq.gz are now removed from the output file"
bismark_id=$(basename "$fastq_R1")
bismark_id=${bismark_id/.fastq.gz/}

# adjust for single/paired end
if [ -n "$fastq_R2" ] ; then
	suffix_bam="_bismark_bt2_pe.bam"
	suffix_report="_bismark_bt2_PE_report.txt"
	suffix_nucstats="_bismark_bt2_pe.nucleotide_stats.txt"
else
	suffix_bam="_bismark_bt2.bam"
	suffix_report="_bismark_bt2_SE_report.txt"
	suffix_nucstats="_bismark_bt2.nucleotide_stats.txt"
fi

bismark_bam_original="${bismark_logs_dir}/${bismark_id}${suffix_bam}"
bismark_report_original="${bismark_logs_dir}/${bismark_id}${suffix_report}"
bismark_nucstats_original="${bismark_logs_dir}/${bismark_id}${suffix_nucstats}"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# exit if output exists already

if [ -s "$bismark_bam_final" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo "${sample},${bismark_bam_final}" >> "$samples_csv"
	exit 1
fi


#########################


# Bismark

# bismark can't use sorted bam at the next step, but other tools may need sorted bam

# bismark/0.22.1 loads bowtie2 and samtools (no version specified)
module add bismark/0.22.1

# "In order to work properly the current working directory must contain the sequence files to be analysed" (as of v0.14)
fastq_dir=$(dirname "$fastq_R1")
cd "$fastq_dir" || exit 1

# adjust for single/paired end
if [ -n "$fastq_R2" ] ; then
	bismark_flags="--maxins 1000"
	fastq_files="-1 $fastq_R1 -2 $fastq_R2"
else
	bismark_flags=""
	fastq_files="$fastq_R1"
fi

# number of parallel instances of Bismark to run (3-5 threads per instance depending on parameters)
# "a typical Bismark run will use several cores already (Bismark itself, Bowtie/Bowtie2, Samtools, gzip etc...)"
multicore_flag=$(( threads / 3 ))

# v0.15.0: "specifying --basename in conjuction with --multicore is currently not supported"

echo
echo " * Bismark: $(readlink -f $(which bismark)) "
echo " * Bismark version: $(bismark --version | grep -m 1 'Version' | tr -s '[:blank:]') "
echo " * bowtie2: $(readlink -f $(which bowtie2)) "
echo " * bowtie2 version: $(bowtie2 --version 2>&1 | head -1) "
echo " * samtools: $(readlink -f $(which samtools)) "
echo " * samtools version: $(samtools --version | head -1) "
echo " * Bismark ref: $ref_bismark "
echo " * FASTQ dir: $fastq_dir "
echo " * FASTQ R1: $fastq_R1 "
echo " * FASTQ R2: $fastq_R2 "
echo " * BAM original: $bismark_bam_original "
echo " * report original: $bismark_report_original "
echo " * BAM final: $bismark_bam_final "
echo " * report final: $bismark_report_final "
echo

bash_cmd="
bismark \
--bowtie2 \
--gzip \
--nucleotide_coverage \
--multicore $multicore_flag \
$bismark_flags \
--temp_dir $bismark_logs_dir \
--output_dir $bismark_logs_dir \
$ref_bismark \
$fastq_files
"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# check that output generated

if [ ! -s "$bismark_bam_original" ] ; then
	echo -e "\n $script_name ERROR: BAM $bismark_bam_original NOT GENERATED \n" >&2
	exit 1
fi

if [ ! -s "$bismark_report_original" ] ; then
	echo -e "\n $script_name ERROR: REPORT $bismark_report_original NOT GENERATED \n" >&2
	exit 1
fi


#########################


# clean up output name

bash_cmd="mv -v $bismark_bam_original $bismark_bam_final"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

bash_cmd="mv -v $bismark_report_original $bismark_report_final"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

bash_cmd="mv -v $bismark_nucstats_original $bismark_nucstats_final"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

sleep 5


#########################


# put reports in report directory for bismark2report
# BAMs for bismark2summary (may not be necessary)

# adjust for single/paired end
if [ -n "$fastq_R2" ] ; then
	bismark_report_bam="${bismark_report_dir}/${sample}${suffix_bam}"
	bismark_report_report="${bismark_report_dir}/${sample}${suffix_report}"
	bismark_report_nucstats="${bismark_report_dir}/${sample}${suffix_nucstats}"
else
	bismark_report_bam="${bismark_report_dir}/${sample}${suffix_bam}"
	bismark_report_report="${bismark_report_dir}/${sample}${suffix_report}"
	bismark_report_nucstats="${bismark_report_dir}/${sample}${suffix_nucstats}"
fi

bash_cmd="ln -sv ../$(basename "$bismark_bam_dir")/$(basename "$bismark_bam_final") $bismark_report_bam"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

bash_cmd="cp -v $bismark_report_final $bismark_report_report"
echo "CMD: $bash_cmd"
eval "$bash_cmd"

bash_cmd="cp -v $bismark_nucstats_final $bismark_report_nucstats"
echo "CMD: $bash_cmd"
eval "$bash_cmd"


#########################


# generate alignment summary

# header for summary file
echo "#SAMPLE,ALIGN INPUT READS,UNIQUELY ALIGNED,ALIGNMENT RATE" > $summary_csv

# print the relevant numbers
paste -d ',' \
<(echo "$sample") \
<(cat "$bismark_report_final" | grep -m 1 "Sequence.*analysed"                | cut -f 2) \
<(cat "$bismark_report_final" | grep -m 1 "alignments with a unique best hit" | cut -f 2) \
<(cat "$bismark_report_final" | grep -m 1 "Mapping efficiency"                | cut -f 2 | sed 's/ //g') \
>> $summary_csv

sleep 5

# combine all sample summaries
cat ${summary_dir}/*.${segment_name}.csv | LC_ALL=C sort -t ',' -k1,1 | uniq > "${proj_dir}/summary.${segment_name}.csv"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bismark_bam_final}" >> "$samples_csv"

sleep 5


#########################



# end
