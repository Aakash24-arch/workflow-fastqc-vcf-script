#!/bin/bash

# Paths
Raw_files="/home/aakash/docript/raw_files"
Result_path="/home/aakash/docript/Result"
Reference="/home/aakash/ngs/db/GCF_000001405.26_GRCh38_genomic.fna"
Bowtie2_Index="/home/aakash/ngs/db/GCF_000001405.26_GRCh38_genomic_index"  # Update with actual Bowtie2 index prefix
# this index path was not correct in previous script = /home/aakash/ngs/db/GCF_000001405.26_GRCh38_genomic_index
# Create result directory 
mkdir -p $Result_path

echo "Starting workflow: FastQC to VCF"

# Loop through raw FASTQ files
for fq1 in $Raw_files/*_1.fq.gz
do
  # Extract sample base name
  base=$(basename $fq1 _1.fq.gz)
  fq2=$Raw_files/${base}_2.fq.gz
  
  echo "================ Processing $base ================"

  ### Step 1: Quality control for raw reads
  quality_path=$Result_path/$base/quality
  mkdir -p $quality_path
  echo "Running FastQC for raw reads..."
  fastqc -t 2 -o $quality_path $fq1 $fq2

  ### Step 2: Trimming low-quality reads
  trimmed_path=$Result_path/$base/trimmed
  mkdir -p $trimmed_path
  echo "Trimming reads with Fastp for $base..."
  fastp --in1 $fq1 --in2 $fq2 \
        --out1 $trimmed_path/${base}_trim_1.fq.gz \
        --out2 $trimmed_path/${base}_trim_2.fq.gz \
        -q 20 -u 20 -l 40 --detect_adapter_for_pe -w 4 \
        --json $trimmed_path/${base}.json \
        --html $trimmed_path/${base}.html
  
  echo "Running FastQC for trimmed reads..."
  fastqc -t 2 -o $quality_path $trimmed_path/${base}_trim_1.fq.gz $trimmed_path/${base}_trim_2.fq.gz

  ### Step 3: Alignment with Bowtie2
  alignment_path=$Result_path/$base/alignment
  mkdir -p $alignment_path
  echo "Aligning reads with Bowtie2 for $base..."
  bowtie2 -x $Bowtie2_Index -1 $trimmed_path/${base}_trim_1.fq.gz -2 $trimmed_path/${base}_trim_2.fq.gz --rg-id ${base}  --rg "SM:${base}" --rg "PL:MGI" --rg "PU:Lane1" --rg "LB:MGI"  | samtools view -b | samtools sort -o $alignment_path/${base}.align.bam --threads 6
##################-R "@RG\\tID:${base}\\tSM:${base}\\tPL:MGI\\tPU:Lane1\\tLB:MGI" it was not 
  ### Step 4: Mark duplicates with GATK
  mark_path=$Result_path/$base/duplicated
  mkdir -p $mark_path
  echo "Marking duplicates for $base..."
  gatk MarkDuplicates \
       -I $alignment_path/${base}.align.bam \
       -O $mark_path/${base}.marked.bam \
       -M $mark_path/${base}.metrics.txt
  samtools index $mark_path/${base}.marked.bam

  ### Step 5: Variant calling with GATK HaplotypeCaller
  variant_path=$Result_path/$base/variants
  mkdir -p $variant_path
  echo "Calling variants with GATK HaplotypeCaller for $base..."
  
  gatk HaplotypeCaller \
       -R $Reference \
       -I $mark_path/${base}.marked.bam \
       -O $variant_path/${base}.vcf.gz \
       --native-pair-hmm-threads 4

  # Index the VCF file
  gatk IndexFeatureFile -I $variant_path/${base}.vcf.gz

  echo "Completed processing for $base!"
done

echo "Workflow completed successfully!"

