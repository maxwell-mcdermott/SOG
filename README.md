# SOG: sequencing data analysis pipelines

[![DOI](https://zenodo.org/badge/66501450.svg)](https://zenodo.org/badge/latestdoi/66501450)

Automated workflows for common sequencing-based (Illumina) protocols, 

such as RNA-seq, ChIP-seq, ATAC-seq, WGBS/RRBS methylation, whole genome/exome/targeted variant detection, contaminant screening, Cut & Tag, and Cut & Run.

We build off those who come before us.

All work is based off of the original SNS created by Igordot, with changes to be used for Feske Lab projects, plus a few extra routes.  If someone somehow gets into this directory please cite him not me.  
For more information on how this was built, see the full documentation at https://igordot.github.io/sns

For those looking for speed here is the demo for SnS which here would just be changed for SOG
module add git

git clone --depth 1 https://github.com/maxwell-mcdermott/SOG

SOG/generate-settings mm10

SOG/gather-fastqs /gpfs/data/sequence/results/XXXX/XXXX-XX-XX/fastq/

SOG/run cnr_with_spikein

Or

SOG/run cnr-spikein-pairs-peaks

Or

SOG/run rna-star





