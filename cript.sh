#!/bin/bash

# Step 1: Define paths and variables
READ1="/home/aakash/ngs/rawData/PremMore_1.fq.gz"
READ2="/home/aakash/ngs/rawData/PremMore_2.fq.gz"
TRIM_DIR="/home/aakash/ngs/docript/trim"
FASTQC_DIR="/home/aakash/ngs/docript/output_fastqc"
SORTED_BAM="/home/aakash/ngs/docript/sorted"
MARKED_BAM="/home/aakash/ngs/docript/mrkd"
REFERENCE="/home/aakash/ngs/db/GCF_000001405.26_GRCh38_genomic.fna"
BOWTIE2_INDEX="/home/aakash/ngs/db/GCF_000001405.26_GRCh38_genomic_index"
OUTPUT_VCF="/home/aakash/ngs/docript/vcf"

# Create required directories
mkdir -p $FASTQC_DIR $FASTQC_DIR/aftertrim $SORTED_BAM $MARKED_BAM $OUTPUT_VCF $TRIM_DIR

# Step 1: Check the quality of raw reads using FastQC
echo "Step 0: Running FastQC on raw reads..."
fastqc $READ1 -o $FASTQC_DIR
if [ $? -eq 0 ]; then
    echo "FastQC on $READ1 completed successfully."
else
    echo "FastQC on $READ1 failed!" && exit 1
fi

fastqc $READ2 -o $FASTQC_DIR
if [ $? -eq 0 ]; then
    echo "FastQC on $READ2 completed successfully."
else
    echo "FastQC on $READ2 failed!" && exit 1
fi

# Step 2: Trim low-quality reads
echo "Step 1: Trimming low-quality reads..."
fastp --in1 $READ1 --in2 $READ2 --out1 $TRIM_DIR/PremMore_trim_1.fq.gz --out2 $TRIM_DIR/PremMore_trim_2.fq.gz -q 20 -u 20 -l 40 --detect_adapter_for_pe -w 4 --json $TRIM_DIR/PremMore.json --html $TRIM_DIR/PremMore.html
if [ $? -eq 0 ]; then
    echo "Trimming completed successfully."
else
    echo "Trimming failed!" && exit 1
fi

# Step 3: Re-check the quality of trimmed reads using FastQC
echo "Step 2: Running FastQC on trimmed reads..."
fastqc $TRIM_DIR/PremMore_trim_1.fq.gz -o $FASTQC_DIR/aftertrim
if [ $? -eq 0 ]; then
    echo "FastQC on trimmed $READ1 completed successfully."
else
    echo "FastQC on trimmed $READ1 failed!" && exit 1
fi

fastqc $TRIM_DIR/PremMore_trim_2.fq.gz -o $FASTQC_DIR/aftertrim
if [ $? -eq 0 ]; then
    echo "FastQC on trimmed $READ2 completed successfully."
else
    echo "FastQC on trimmed $READ2 failed!" && exit 1
fi

# Step 4: Align reads to the reference genome using Bowtie2 and sort the BAM file
echo "Step 3: Aligning reads with Bowtie2 and sorting BAM file..."
bowtie2 -x $BOWTIE2_INDEX -1 $TRIM_DIR/PremMore_trim_1.fq.gz -2 $TRIM_DIR/PremMore_trim_2.fq.gz | samtools view -b | samtools sort -o $SORTED_BAM/PremMore_sorted.bam
if [ $? -eq 0 ]; then
    echo "Alignment and sorting completed successfully."
else
    echo "Alignment and sorting failed!" && exit 1
fi

# Step 5: Add read group to the BAM file
echo "Step 4: Adding read group to BAM file..."
gatk AddOrReplaceReadGroups -I $SORTED_BAM/PremMore_sorted.bam -O $SORTED_BAM/PremMore_sorted_withRG.bam -RGID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM sample1
if [ $? -eq 0 ]; then
    echo "Read group added successfully."
else
    echo "Adding read group failed!" && exit 1
fi

# Step 6: Mark duplicates using GATK
echo "Step 5: Marking duplicates..."
gatk MarkDuplicates -I $SORTED_BAM/PremMore_sorted_withRG.bam -O $MARKED_BAM/PremMore_mrkd.bam -M /home/aakash/ngs/docript/mrkd/marked_dup_metrics.txt --CREATE_INDEX TRUE
if [ $? -eq 0 ]; then
    echo "Duplicates marked successfully."
else
    echo "Marking duplicates failed!" && exit 1
fi

# Step 7: Variant calling using GATK HaplotypeCaller
echo "Step 6: Calling variants..."
gatk HaplotypeCaller -R $REFERENCE -I $MARKED_BAM/PremMore_mrkd.bam -O $OUTPUT_VCF/PremMore_varient.vcf
if [ $? -eq 0 ]; then
    echo "Variant calling completed successfully."
else
    echo "Variant calling failed!" && exit 1
fi

echo "Pipeline completed successfully!"

