#!/bin/bash


# Ecoli normalization and BIGwig creation
# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads BAM \n" >&2
	exit 0
fi

# arguments
proj_dir=$(readlink -f "$1")
sample=$2
threads=$3
bam_dd_unsorted=$4

#########################

bam_dd_dir="${proj_dir}/BAM-DD"
# bam_dd_unsorted="${bam_dd_dir}/${sample}.dd.bam"
bam_dd_sorted="${bam_dd_dir}/${sample}.sorted.dd.bam"
beds_dir="${proj_dir}/BEDS"
mkdir -p "$beds_dir"
bed_sample="${beds_dir}/${sample}_bowtie2.bed"
cleaned_bed="${beds_dir}/${sample}_bowtie2.clean.bed"
fragments_bed="${beds_dir}/${sample}_bowtie2.fragments.bed"
sample_bedgraph="${beds_dir}/${sample}_bowtie2.fragments.bedgraph"
chrom_sizes="/gpfs/data/igorlab/ref/mm10/chrom.sizes"
bigwig_dir="${proj_dir}/BIGWIG"
seacr_bw="${bigwig_dir}/${sample}.seacr.bw"

#########################
# unload all loaded modulefiles
module purge
module add default-environment
module add bedtools/2.30.0
module add samtools/1.9
module add deeptools/3.5.1 
module add ucscutils/374
module add r/4.1.1
#not needed for checking only
#########################
# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: proj dir $proj_dir does not exist \n" >&2
	exit 0
fi

if [ ! -s "$bam_dd_dir" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_dd_dir does not exist \n" >&2
	exit 0
fi

code_dir=$(dirname $(dirname "$script_path"))


chrom_sizes=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-CHROMSIZES);

if [ ! -s "$chrom_sizes" ] ; then
	echo -e "\n $script_name ERROR: chrom sizes $chrom_sizes does not exist \n" >&2
	exit 0
fi

blacklist=$(bash ${code_dir}/scripts/get-set-setting.sh "${proj_dir}/settings.txt" REF-BLACKLIST);

if [ ! -s "$blacklist" ] ; then
	echo -e "\n $script_name ERROR: blacklist $blacklist does not exist \n" >&2
	exit 0
fi

# exit if output exits already

if [ -s "$seacr_bw" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi

#########################


echo
echo " * Python: $(readlink -f $(which python)) "
echo " * Python version: $(python --version 2>&1 | head -1) "
echo " * deepTools: $(readlink -f $(which deeptools)) "
echo " * samtools: $(readlink -f $(which samtools)) "
echo " * bamCoverage: $(readlink -f $(which bamCoverage)) "
echo " * bamCoverage version: $(bamCoverage --version 2>&1 | head -1) "
echo " * BAM: $bam_dd_unsorted "
echo " * BIGWIG: $seacr_bw "
echo

#########################

samtools sort -n ${bam_dd_unsorted} --threads $threads  >  ${bam_dd_sorted}
bedtools bamtobed -i ${bam_dd_sorted} -bedpe > ${bed_sample}
awk '$1==$4 && $6-$2 < 500 {print $0}' ${bed_sample} > ${cleaned_bed}
cut -f 1,2,6  ${cleaned_bed} | sort -k1,1 -k2,2n -k3,3n  > ${fragments_bed}
segment_align="align-bowtie2-cutnrun_ecoli_spikein"
seqDepth=$(grep -s -m 1 "^${sample}," "${proj_dir}/summary.${segment_align}.csv" | cut -d ',' -f 6)
scale_factor=`echo "100000 / $seqDepth" | bc -l`
echo "Scaling factor for $sample is: $scale_factor"
bedtools genomecov -bg -scale $scale_factor -i ${fragments_bed} -g $chrom_sizes > ${sample_bedgraph}
bedGraphToBigWig ${sample_bedgraph}  $chrom_sizes  $seacr_bw

###########################################################################

# check that output generated

if [ ! -s "$${bam_dd_sorted}" ] ; then
	echo -e "\n $script_name ERROR: peaks $${bam_dd_sorted} not generated \n" >&2
	exit 0
fi

if [ ! -s "$sample_bedgraph" ] ; then
	echo -e "\n $script_name ERROR: XLS $sample_bedgraph not generated \n" >&2
	exit 0
fi

if [ ! -s "$seacr_bw" ] ; then
	echo -e "\n $script_name ERROR: XLS $seacr_bw not generated \n" >&2
	exit 0
fi

#########################


# generate summary

# num_peaks_unfiltered=$(cat "$peaks_file" | wc -l)
echo "amount of paired reads of ecoli : $seqDepth"
echo "scaling factor used for 100000: $scale_factor"


sleep 5

#########################



# end
