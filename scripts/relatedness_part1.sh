# script to analyze relationship between geographical distance and pair-wise coefficients of relatedness

# Use the output of the vcf2smartpca.py script and awk to retain only SNPs used for PCA analysis for relatedness

# Trim first column of SNP coordinates
cut -f1 par.passing_snps_q20_SNP > par.passing_snps_q20_SNP_clip

# Format trimmed file for VCFtools
awk  '{gsub("_","\t",$0); print;}' par.passing_snps_q20_SNP_clip > par.passing_snps_q20_SNP_clip_vcfinput.txt

# Calculate kinship coefficients (hereafter relatedness) using VCFtools
vcftools --vcf D_magna_RADseq.vcf --positions par.passing_snps_q20_SNP_clip_vcfinput.txt --relatedness --out par.passing_snps_q20_SNP_rel1

####### VCFtools can be used to subset dataset, as in:
vcftools --vcf D_magna_RADseq.vcf --remove-indv MN --recode --out VQSR_noMN