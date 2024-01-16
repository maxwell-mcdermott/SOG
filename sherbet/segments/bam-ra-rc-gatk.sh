#!/bin/bash


# GATK realignment and recalibration


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
proj_dir=$1
sample=$2
threads=$3
bam=$4


#########################


# settings and files

samples_csv="${proj_dir}/samples.${segment_name}.csv"

# should account for both dedup (.dd.bam) non-dedup (.bam) input BAMs
bam_base=$(basename "$bam")
bam_base=${bam_base/%.bam/}

bam_ra_dir="${proj_dir}/BAM-GATK-RA"
mkdir -p "$bam_ra_dir"
bam_ra="${bam_ra_dir}/${bam_base}.ra.bam"
bai_ra="${bam_ra_dir}/${bam_base}.ra.bai"

bam_ra_rc_dir="${proj_dir}/BAM-GATK-RA-RC"
mkdir -p "$bam_ra_rc_dir"
bam_ra_rc="${bam_ra_rc_dir}/${bam_base}.ra.rc.bam"
bai_ra_rc="${bam_ra_rc_dir}/${bam_base}.ra.rc.bai"

gatk_logs_dir="${proj_dir}/logs-${segment_name}"
mkdir -p "$gatk_logs_dir"
gatk_ra_intervals="${gatk_logs_dir}/${sample}.intervals"
gatk_rc_table1="${gatk_logs_dir}/${sample}.table1.txt"
gatk_rc_table2="${gatk_logs_dir}/${sample}.table2.txt"
gatk_rc_csv="${gatk_logs_dir}/${sample}.csv"
gatk_rc_pdf="${gatk_logs_dir}/${sample}.pdf"

# unload all loaded modulefiles
module purge
module add default-environment


#########################


# check for output

# skip if final BAM and BAI exist
if [ -s "$bam_ra_rc" ] && [ -s "${bam_ra_rc}.bai" ] ; then
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	echo "${sample},${bam_ra_rc}" >> "$samples_csv"
	exit 0
fi

# delete potentially incomplete files (since the corresponding BAM index was not generated)
if [ -s "$bam_ra_rc" ] ; then
	echo -e "\n $script_name WARNING: POTENTIALLY CORRUPT BAM $bam_ra_rc EXISTS \n" >&2
	# delete BAMs
	rm -fv "$bam_ra_rc"
	rm -fv "$bam_ra"
	# delete realignment files that are no longer needed
	rm -fv "$bam_ra"
	rm -fv "$bai_ra"
	rm -fv "$gatk_ra_intervals"
	# delete recalibration files that are no longer needed
	rm -fv "$gatk_rc_table1"
	rm -fv "$gatk_rc_table2"
fi

# if realigned BAM exists, assume that there is another process that is currently running and skip
if [ -s "$bam_ra" ] ; then
	echo -e "\n $script_name WARNING: $bam_ra EXISTS (SAMPLE IS POTENTIALLY BEING PROCESSED) \n" >&2
	echo -e "\n $script_name SKIP SAMPLE $sample \n" >&2
	exit 1
fi


#########################


# check that inputs exist

if [ ! -d "$proj_dir" ] ; then
	echo -e "\n $script_name ERROR: PROJ DIR $proj_dir DOES NOT EXIST \n" >&2
	exit 1
fi

if [ ! -s "$bam" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam DOES NOT EXIST \n" >&2
	exit 1
fi

code_dir=$(dirname $(dirname "$script_path"))

genome_dir=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" GENOME-DIR)

if [ ! -d "$genome_dir" ] ; then
	echo -e "\n $script_name ERROR: GENOME DIR $genome_dir DOES NOT EXIST \n" >&2
	exit 1
fi

ref_fasta=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-FASTA)

if [ ! -s "$ref_fasta" ] ; then
	echo -e "\n $script_name ERROR: FASTA $ref_fasta DOES NOT EXIST \n" >&2
	exit 1
fi

ref_dict=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" REF-DICT)

if [ ! -s "$ref_dict" ] ; then
	echo -e "\n $script_name ERROR: DICT $ref_dict DOES NOT EXIST \n" >&2
	exit 1
fi

# check for BED files in the project directory and set as project BED file if it is not already defined
found_bed=$(find "$proj_dir" -maxdepth 1 -type f -iname "*.bed" | grep -v "probes" | sort | head -1)
bed=$(bash "${code_dir}/scripts/get-set-setting.sh" "${proj_dir}/settings.txt" EXP-TARGETS-BED "$found_bed")

if [ ! -s "$bed" ] ; then
	echo -e "\n $script_name ERROR: BED $bed DOES NOT EXIST \n" >&2
	exit 1
fi


#########################


# GATK settings

module add r/3.6.1

# command
gatk_jar="/gpfs/data/igorlab/software/GenomeAnalysisTK/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar"
gatk_cmd="java -Xms8G -Xmx8G -jar ${gatk_jar}"

if [ ! -s "$gatk_jar" ] ; then
	echo -e "\n $script_name ERROR: GATK $gatk_jar DOES NOT EXIST \n" >&2
	exit 1
fi

# error log (DEBUG, INFO (default), WARN, ERROR, FATAL, OFF)
gatk_log_level_arg="--logging_level ERROR"

# known variants (may vary greatly for each genome)
genome_build=$(basename "$genome_dir")
if [[ "$genome_build" == "hg19" ]] ; then
	gatk_indel_vcf_1="${genome_dir}/gatk-bundle/1000G_phase1.indels.hg19.vcf"
	gatk_indel_vcf_2="${genome_dir}/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf"
	gatk_snp_vcf="${genome_dir}/gatk-bundle/dbsnp_138.hg19.vcf"
	gatk_ra_known_arg="-known $gatk_indel_vcf_1 -known $gatk_indel_vcf_2"
	gatk_rc_known_arg="-knownSites $gatk_indel_vcf_1 -knownSites $gatk_indel_vcf_2 -knownSites $gatk_snp_vcf"
elif [[ "$genome_build" == "hg38" ]] ; then
	gatk_indel_vcf="${genome_dir}/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
	gatk_snp_vcf="${genome_dir}/gatk-bundle/dbsnp_146.hg38.vcf.gz"
	gatk_ra_known_arg="-known $gatk_indel_vcf"
	gatk_rc_known_arg="-knownSites $gatk_indel_vcf -knownSites $gatk_snp_vcf"
elif [[ "$genome_build" == "mm10" ]] ; then
	gatk_indel_vcf="${genome_dir}/MGP/mgp.v5.indels.pass.chr.sort.vcf"
	gatk_snp_vcf="${genome_dir}/dbSNP/dbsnp.146.vcf"
	gatk_ra_known_arg="-known $gatk_indel_vcf"
	gatk_rc_known_arg="-knownSites $gatk_indel_vcf -knownSites $gatk_snp_vcf"
elif [[ "$genome_build" == "dm3" ]] ; then
	gatk_indel_vcf=""
	gatk_snp_vcf="${genome_dir}/dgrp2.min.vcf"
	gatk_ra_known_arg=""
	gatk_rc_known_arg="-knownSites $gatk_snp_vcf"
elif [[ "$genome_build" == "dm6" ]] ; then
	gatk_indel_vcf="${genome_dir}/dbSNP/dbsnp.149.indel.vcf"
	gatk_snp_vcf="${genome_dir}/dbSNP/dbsnp.149.vcf"
	gatk_ra_known_arg="-known $gatk_indel_vcf"
	gatk_rc_known_arg="-knownSites $gatk_snp_vcf"
elif [[ "$genome_build" == "canFam3" ]] ; then
	gatk_indel_vcf="${genome_dir}/dbSNP/dbsnp.151.indel.vcf"
	gatk_snp_vcf="${genome_dir}/dbSNP/dbsnp.151.vcf"
	gatk_ra_known_arg="-known $gatk_indel_vcf"
	gatk_rc_known_arg="-knownSites $gatk_snp_vcf"
else
	echo -e "\n $script_name ERROR: genome $genome_build at $genome_dir not supported \n" >&2
	exit 1
fi


#########################


# test R (GATK uses R for plotting)

echo
echo " * R: $(readlink -f $(which R)) "
echo " * R version: $(R --version | head -1) "
echo " * Rscript: $(readlink -f $(which Rscript)) "
echo " * Rscript version: $(Rscript --version 2>&1) "
echo

# basic packages
Rscript --vanilla "${code_dir}/scripts/test-package.R" getopt
Rscript --vanilla "${code_dir}/scripts/test-package.R" optparse

# required by GATK 3
Rscript --vanilla "${code_dir}/scripts/test-package.R" gsalib
Rscript --vanilla "${code_dir}/scripts/test-package.R" reshape
Rscript --vanilla "${code_dir}/scripts/test-package.R" gplots
Rscript --vanilla "${code_dir}/scripts/test-package.R" ggplot2

# required by GATK 4
Rscript --vanilla "${code_dir}/scripts/test-package.R" data.table
Rscript --vanilla "${code_dir}/scripts/test-package.R" naturalsort


#########################


# realignment

echo
echo " * GATK: $(readlink -f $gatk_jar) "
echo " * GATK version: $($gatk_cmd --version) "
echo " * BAM IN: $bam "
echo " * BAM RA: $bam_ra "
echo

gatk_ra1_cmd="
$gatk_cmd -T RealignerTargetCreator -dt NONE $gatk_log_level_arg \
-nt $threads \
--reference_sequence $ref_fasta \
$gatk_ra_known_arg \
--intervals $bed --interval_padding 10 \
--input_file $bam \
--out $gatk_ra_intervals
"
echo "CMD: $gatk_ra1_cmd"
$gatk_ra1_cmd

gatk_ra2_cmd="
$gatk_cmd -T IndelRealigner -dt NONE $gatk_log_level_arg \
--reference_sequence $ref_fasta \
--maxReadsForRealignment 50000 \
$gatk_ra_known_arg \
-targetIntervals $gatk_ra_intervals \
--input_file $bam \
--out $bam_ra
"
echo "CMD: $gatk_ra2_cmd"
$gatk_ra2_cmd

sleep 5


#########################


# check that output generated

# if BAM not generated, delete any related files since something went wrong and they might be corrupted
if [ ! -s "$bam_ra" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_ra NOT GENERATED \n" >&2
	rm -fv "$gatk_ra_intervals"
	exit 1
fi

# if BAM index not generated, delete any related files since something went wrong and they might be corrupted
if [ ! -s "$bai_ra" ] ; then
	echo -e "\n $script_name ERROR: BAM INDEX $bai_ra NOT GENERATED \n" >&2
	rm -fv "$bam_ra"
	rm -fv "$gatk_ra_intervals"
	exit 1
fi


#########################


# recalibration

# first pass recalibration table file
gatk_rc1_cmd="
$gatk_cmd -T BaseRecalibrator $gatk_log_level_arg \
-nct $threads \
-rf BadCigar \
--reference_sequence $ref_fasta \
$gatk_rc_known_arg \
--intervals $bed --interval_padding 10 \
--input_file $bam_ra \
--out $gatk_rc_table1
"
echo "CMD: $gatk_rc1_cmd"
$gatk_rc1_cmd

# second pass recalibration table file
gatk_rc2_cmd="
$gatk_cmd -T BaseRecalibrator $gatk_log_level_arg \
-nct $threads \
-rf BadCigar \
--reference_sequence $ref_fasta \
$gatk_rc_known_arg \
--intervals $bed --interval_padding 10 \
--input_file $bam_ra \
-BQSR $gatk_rc_table1 \
--out $gatk_rc_table2
"
echo "CMD: $gatk_rc2_cmd"
$gatk_rc2_cmd

# generate the plots report and also keep a copy of the csv (optional)
gatk_rc3_cmd="
$gatk_cmd -T AnalyzeCovariates $gatk_log_level_arg \
--reference_sequence $ref_fasta \
-before $gatk_rc_table1 \
-after $gatk_rc_table2 \
-csv $gatk_rc_csv \
-plots $gatk_rc_pdf
"
echo "CMD: $gatk_rc3_cmd"
$gatk_rc3_cmd

# generate recalibrated BAM
gatk_rc4_cmd="
$gatk_cmd -T PrintReads $gatk_log_level_arg \
-nct $threads \
-rf BadCigar \
--reference_sequence $ref_fasta \
-BQSR $gatk_rc_table1 \
--input_file $bam_ra \
--out $bam_ra_rc
"
echo "CMD: $gatk_rc4_cmd"
$gatk_rc4_cmd

sleep 5


#########################


# check that output generated

# if BAM not generated, delete any related files since something went wrong and they might be corrupted
if [ ! -s "$bam_ra_rc" ] ; then
	echo -e "\n $script_name ERROR: BAM $bam_ra_rc NOT GENERATED \n" >&2
	rm -fv "$gatk_rc_table1"
	rm -fv "$gatk_rc_table2"
	exit 1
fi

# if BAM index not generated, delete any related files since something went wrong and they might be corrupted
if [ ! -s "$bai_ra_rc" ] ; then
	echo -e "\n $script_name ERROR: BAM INDEX $bai_ra_rc NOT GENERATED \n" >&2
	rm -fv "$bam_ra_rc"
	rm -fv "$gatk_rc_table1"
	rm -fv "$gatk_rc_table2"
	exit 1
fi


#########################


# clean up

# move .bai (GATK convention) to .bam.bai (samtools convention) since some tools expect that
mv -v "$bai_ra_rc" "${bam_ra_rc}.bai"

# delete realignment files that are no longer needed
rm -fv "$bam_ra"
rm -fv "$bai_ra"
rm -fv "$gatk_ra_intervals"

# delete recalibration files that are no longer needed
rm -fv "$gatk_rc_table1"
rm -fv "$gatk_rc_table2"


#########################


# add sample and BAM to sample sheet
echo "${sample},${bam_ra_rc}" >> "$samples_csv"

sleep 5


#########################



# end
