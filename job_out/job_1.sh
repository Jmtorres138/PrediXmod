#!/bin/bash
#PBS -N vcf_conversion 
#PBS -S /bin/bash
#PBS -l mem=40gb
#PBS -o $HOME/group/im-lab/nas40t2/jason/projects/PrediXmod/job_out/job_1_Mex.out
#PBS -e $HOME/group/im-lab/nas40t2/jason/projects/PrediXmod/job_out/job_1_Mex.err


perl /group/im-lab/nas40t2/jason/projects/PrediXmod/1_vcf2dosage.mach_gtex_MexSNPs.pl 
