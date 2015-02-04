###PrediXmod
####Pipeline to develop prediction models for GTEx data and cross-validation

-Instructions for downloading expression and genotype data from **GTEx** exchange area are 
 listed in gtex.md

-Instruction for accessing GitHub from a secure server are given in tarbell_to_github.md
 
-Descriptions of scripts that execute steps in gene prediction model development pipeline shown below: 

#### Convert GTEx genotypes to downstream file formats
- for all autosomal SNPs combined in one file: 1_vcf2dosage.mach_gtex_hapmapSNPs.pl
- for separate files per autosome: 1_vcf2dosage.mach.chr_gtex_hapmapSNPs.pl

#### Convert GTEx RNA-seq data for downstream analyses
- 2_pull_protein_coding_RPKM_rebuild.R

#### Normalize and adjust the GTEx RNA-seq data by PEER factors and genotype principal components
- 3_calc_PEER_adj_exp.r

#### Run LASSO cross-validation for GTEx tissue of interest, save R<sup>2</sup> and best betas
- 4_CV_GTEx_lasso_adjusted.r

#### Run polyscore cross-validation for GTEx tissue of interest, save R<sup>2</sup> and top betas
- 8_CV_GTEX_polyscore_adjusted.r
