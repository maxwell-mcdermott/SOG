#!/bin/bash


# deepTools generate BigWig from BAM


# script filename
script_path="${BASH_SOURCE[0]}"
script_name=$(basename "$script_path")
segment_name=${script_name/%.sh/}
echo -e "\n ========== SEGMENT: $segment_name ========== \n" >&2

# check for correct number of arguments
if [ ! $# == 4 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name project_dir sample_name num_threads BAM \n" >&2
	exit 1
fi

# arguments
proj_dir=$(readlink -f "$1")
sample=$2
threads=$3
bam=$4



bam_dd_dir="${proj_dir}/BAM-DD"
# bam_dd_unsorted="${bam_dd_dir}/${sample}.dd.bam"
bam_dd_sorted="${bam_dd_dir}/${sample}.sorted.dd.bam"
bam_dd_10M_sorted="${bam_dd_dir}/${sample}.sorted.10M.bam"
bam_dd_10M_sorted="${bam_dd_dir}/${sample}.sorted.10M.bam"


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# settings and files

bigwig_dir="${proj_dir}/BIGWIG"
mkdir -p "$bigwig_dir"
bigwig="${bigwig_dir}/${sample}.10M_Downsample.bin1.rpkm.bw"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# exit if output exits already

if [ -s "$bigwig" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 0
fi

#########################
#downsample to 10M reads first

# skip if final bedgraph and bigwig exist
if [ -s "$bam_dd_10M_sorted" ] && [ -s "${bigwig}.bai" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo -e "\n $script_name ADD $sample TO $samples_csv \n" >&2
	echo "${sample},${sample_bedgraph},${bigwig}" >> "$samples_csv"
	exit 0
fi


module add bedtools/2.30.0
module add samtools/1.20
module add deeptools/3.5.1 
module add ucscutils/368
module add r/4.1.1
#not needed for checking only

# downsample the deduplicated BAM file
segment_align="bam-dedup-sambamba"
seqDepth=$(grep -s -m 1 "^${sample}," "${proj_dir}/summary.${segment_align}.csv" | cut -d ',' -f 3)
echo $seqDepth
scale_factor=`echo "10000000 / $seqDepth " | bc -l`
echo $scale_factor
samtools view -s $scale_factor -b -p ${bam} > ${bam_dd_10M_sorted} 
samtools index  ${bam_dd_10M_sorted} 


#########################
#########################
# unload all loaded modulefiles
module purge
module add default-environment

# generate bigWig using deepTools

module add deeptools/3.5.1
# deepTools bamCoverage requires bedGraphToBigWig
module add ucscutils/374

echo
echo " * Python: $(readlink -f $(which python)) "
echo " * Python version: $(python --version 2>&1 | head -1) "
echo " * deepTools: $(readlink -f $(which deeptools)) "
echo " * deepTools version: $(deepTools --version 2>&1 | head -1) "
echo " * bamCoverage: $(readlink -f $(which bamCoverage)) "
echo " * bamCoverage version: $(bamCoverage --version 2>&1 | head -1) "
echo " * BAM: $bam_dd_10M_sorted "
echo " * BIGWIG: $bigwig "
echo

bamcov_cmd="
bamCoverage \
--verbose \
--numberOfProcessors $threads \
--binSize 1 \
--normalizeUsing RPKM \
--ignoreForNormalization chrX \
--maxFragmentLength 700 \
--minFragmentLength 10 \
--outFileFormat bigwig \
--bam $bam_dd_10M_sorted \
--outFileName $bigwig
"
echo "CMD: $bamcov_cmd"
eval "$bamcov_cmd"


#########################


# check that output generated

if [ ! -s "$bigwig" ] ; then
	echo -e "\n $script_name ERROR: BIGWIG $bigwig NOT GENERATED \n" >&2
	exit 1
fi


#########################



# end