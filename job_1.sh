#!/bin/bash
#PBS -N vcf_conversion 
#PBS -S /bin/bash
#PBS -l mem=40gb
#PBS -o $HOME/group/im-lab/nas40t2/jason/projects/PrediXmod/job_out/job_1.out
#PBS -e $HOME/group/im-lab/nas40t2/jason/projects/PrediXmod/job_out/job_1.err


perl /group/im-lab/nas40t2/jason/projects/PrediXmod/1_vcf2dosage.mach_gtex_hapmapSNPs.pl /group/im-lab/nas40t2/haky/Data/dbGaP/GTEx/41400/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2014-06-13/genotypes/OMNI_arrays/GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv_allchr_genot_imput_info04_maf01_HEW1E6.vcf.gz /group/im-lab/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list /group/im-lab/nas40t2/jason/projects/PrediXmod/genotypes/   gtex 
