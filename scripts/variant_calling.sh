# Script to generate VCF file, adapted from:
# De Wit P, Pespeni MH, Ladner JT, et al. (2012) The simple fool's guide to population genomics via RNA-Seq: an introduction to high-throughput sequencing data analysis. Molecular Ecology Resources 12, 1058-1067

# index reference
bwa index dmagna2.4.fa

# align reads to dmagna_reference
bwa mem dmagna2.4.fa AST.fq.gz | gzip -3 > AST.sam.gz
# ... for all individuals

# convert to bam, remove low mapping reads
samtools view -q 20 -bS AST.sam.gz > AST.bam
# ... for all individauls

# sort bam
samtools sort AST.bam AST.sorted
# ... for all individuals

# create file to allow maintanence of RG info in header, with format
# @RG     ID:AST  SM:AST  PL:Illumina
# REMEMBER the tabs!!!!!!!!!!!

# merge individual bam files
samtools merge -rh rg.txt merged.bam *.bam

# index merged.bam
samtools index merged.bam

# Begin GATK pipeline, v.3.1, pulled from github on June 20th, built from source
# with java 1.8

# generate dictionaries for reference
# with samtools
samtools faidx dmagna2.4.fa

# with picard
java -jar CreateSequenceDictionary.jar R=dmagna2.4.fa O=dmagna2.4.dict

# Find candidates for re-alignment
java -Xmx2g -jar GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R ~/daphnia/ref_file/dmagna2.4.fa \
  -o merged_output.intervals \
  -I merged.bam \
  --minReadsAtLocus 3 -nt 8

#  Run realigner over intervals, no nt
java -Xmx22g -jar GenomeAnalysisTK.jar \
  -I merged.bam \
  -R ~/daphnia/ref_file/dmagna2.4.fa \
  -T IndelRealigner \
  -targetIntervals merged_output.intervals \
  -o merged_realigned.bam \
  -LOD 3.0 \
  --maxReadsInMemory 1000000 \
  --maxReadsForRealignment 100000

# Call SNPs
java -Xmx22g -jar GenomeAnalysisTK.jar \
	-R dmagna2.4.fa -T UnifiedGenotyper \
	-I merged_realigned.bam -o rawSNPS_Q30.vcf \
	-gt_mode DISCOVERY \
	-stand_call_conf 30 -stand_emit_conf 10 -nt 7

# Annotations, Standard
java -Xmx16g -jar GenomeAnalysisTK.jar \
  -T VariantAnnotator \
  -l INFO \
  -R dmagna2.4.fa \
  -I merged_realigned.bam \
  -o rawSNPS_Q30_annotated.vcf \
  -V rawSNPS_Q30.vcf --list


# Call Indels
java -Xmx16g -jar GenomeAnalysisTK.jar \
	-R dmagna2.4.fa -T UnifiedGenotyper \
	-I merged_realigned.bam -o InDels_Q30.vcf \
	-gt_mode DISCOVERY \
	-glm INDEL \
	-stand_call_conf 30 -stand_emit_conf 10 -nt 7


# Filter around InDels_Q30

java -Xmx16g -jar GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R dmagna2.4.fa \
  --mask InDels_Q30.vcf \
  -V rawSNPS_Q30_annotated.vcf \
  -o Indel_filtered_Q30.vcf

# Additional filtering

java -Xmx16g -jar GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R dmagna2.4.fa \
  -V Indel_filtered_Q30.vcf \
  -o analysis_ready_Q30.vcf \
  --filterExpression "QD < 5.0" \
  --filterName "QD" \
  --filterExpression "DP < 5 " \
  --filterName "LowCoverage"

# Save only the header and all SNPS that have passed all the filters in a new file that can be used as a truth training set for the VQSR:  

cat analysis_ready_Q30.vcf | grep 'PASS\|^#' > highQualSNPS.vcf

# Lower quality calls
java -Xmx16g -jar GenomeAnalysisTK.jar \
	-R dmagna2.4.fa -T UnifiedGenotyper \
	-I merged_realigned.bam -o rawSNPS_Q4.vcf \
	-gt_mode DISCOVERY \
	-stand_call_conf 4 -stand_emit_conf 3 -nt 7

# Annotations, not really necessary

java -Xmx16g -jar GenomeAnalysisTK.jar \
  -T VariantAnnotator \
  -l INFO \
  -R dmagna2.4.fa \
  -I merged_realigned.bam \
  -o rawSNPS_Q4_annotated.vcf \
  -V rawSNPS_Q4.vcf


#Calling InDels (needed for filtering around InDels):

java -Xmx16g -jar GenomeAnalysisTK.jar \
	-R dmagna2.4.fa -T UnifiedGenotyper \
	-I merged_realigned.bam -o InDels_Q4.vcf \
	-gt_mode DISCOVERY \
	-glm INDEL \
	-stand_call_conf 4 -stand_emit_conf 3 -nt 8

#Filtering around InDels:

java -Xmx16g -jar GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R dmagna2.4.fa \
  --mask InDels_Q4.vcf \
  -V rawSNPS_Q4_annotated.vcf \
  -o Indel_filtered_Q4.vcf

#Additional filtering:

java -Xmx16g -jar GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R dmagna2.4.fa \
  -V Indel_filtered_Q4.vcf \
  -o filtered_Q4.vcf \
  --filterExpression "QD < 5.0" \
  --filterName "QD" \
  --filterExpression "DP < 5 " \
  --filterName "LowCoverage"

# Variant recalibrator

java -Xmx16g -jar GenomeAnalysisTK.jar \
   -T VariantRecalibrator \
   -R dmagna2.4.fa \
 	-input filtered_Q4.vcf \
   -resource:concordantSet,known=true,training=true,truth=true,prior=10.0 highQualSNPS.vcf \
   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
   -recalFile VQSR.recal \
   -tranchesFile VQSR.tranches \
   -rscriptFile VQSR.plots.R \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   --ignore_filter HARD_TO_VALIDATE \
   --ignore_filter LowQual

# Variant recalibrator

java -Xmx4g -jar GenomeAnalysisTK.jar  \
   -T VariantRecalibrator \
   -R dmagna2.4.fa \
   -input filtered_Q4.vcf \
   -resource:concordantSet,known=true,training=true,truth=true,prior=10.0 highQualSNPS.vcf \
   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
   -mode SNP \
   -recalFile VQSR.recal \
   -tranchesFile VQSR.tranches \
   -rscriptFile VQSR.plots.R \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   --ignore_filter HARD_TO_VALIDATE \
   --ignore_filter LowQual

#Applying the recalibration

java -Xmx16g -jar GenomeAnalysisTK.jar \
   -T ApplyRecalibration \
   -R dmagna2.4.fa \
   -input filtered_Q4.vcf \
   --ts_filter_level 99.0 \
   --ignore_filter HARD_TO_VALIDATE \
   --ignore_filter LowQual \
   -tranchesFile VQSR.tranches \
   -recalFile VQSR.recal \
   -o recalibrated_filtered_SNPS.vcf

# save all the SNPS that have passed the VQSR filter into a new vcf file:

cat recalibrated_filtered_SNPS.vcf | grep 'VQSLOD\|^#' | grep -v TruthSensitivityTranche > D_magna_RADseq.vcf



