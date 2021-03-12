# smartpca PCA analysis and Procrustes analysis
# generate smartpca inputs
vcf2smartpca.py D_magna_RADseq.vcf passing_snps_q20 20 popfile 
# Conversion script can be obtained from:
# De Wit P, Pespeni MH, Ladner JT, et al. (2012) The simple fool's guide to population genomics via RNA-Seq: an introduction to high-throughput sequencing data analysis. Molecular Ecology Resources 12, 1058-1067

# smartpca running commands
smartpca -p par.passing_snps_q20 > passing_snps_q20_logfile.txt

# re-run the above steps in order to reduce dataset by MN and IL individuals using
vcf2smartpca.py D_magna_RADseq.vcf passing_snps_q20 20 popfile INDIVIDUALS_TO_OMIT MN IL

# and re-run smartpca to generate appropriate PCA eigenvalue and eigenvector files to go into Procrustes analysis
