#!bin/bash 
#
#author aakash
#
# path of directories 
#
READ1="/home/aakash/ngs/rawData/PremMore_1.fq.gz"
READ2="/home/aakash/ngs/rawData/PremMore_2.fq.gz"
FASTQC_DIR="/home/aakash/ngs/2script/Output_fastqc"
TRIM_DIR="/home/aakash/ngs/2script/trim"
BOWTIE_INDEX="/home/aakash/ngs/db/GCF_000001405.26_GRCh38_genomic_index"
SORTED="/home/aakash/ngs/2script/sorted"
MARKED="/home/aakash/ngs/2script/marked"
REFERENCE="/home/aakash/ngs/db/GCF_000001405.26_GRCh38_genomic.fna"
VARIENT_OUTPUT="/home/aakash/ngs/2script/Output_vcf"


#create directories 
#step 1 
#
mkdir -p $FASTQC $FASTQC/AFTERTRIM $TRIM $SORTED $MARKED $VARIENT_OUTPUT


##############################################################################step 1 checking the quality score via fastqc tools##################################################################################################### 
#
echo "step 1-- Read1 --  checking the quality score"
fastqc $READ1 -o $FASTQC_DIR 
if [ $? -ep 0 ];then
	echo "fastqc is completted on $Read1"
else
	echo "fastqc is failed on $Read1" && exit 1
fi

###################################################################################step 1 checking the quality score via fastqc tool#################################################################################################
#
echo "step 1 --Read2 -- checking the quality score"
fastqc $Read2 -o $FASTQC_DIR
if [ $? -ep 0 ];then
	echo "fastqc is completed on $READ2"
else
	echo "fastqc is failed on $READ2" && exit 1 
fi

#########################################################################################step 2 trimming low quality reads file via fastp tool#######################################################################################
#
echo "step 2 --trim low quality reads"
fastp --int1 $READ1 --int2 $READ2 --out1 $TRIM_DIR/PremMore_trim_1_fq.gz --out $TRIM_DIR/PremMore_trim_2_fq.gz -q 20 -u 20 -l 40 --detect_adapter_for_pe -w 4 --json $TRIM_DIR/PremMore.json --html $TRIM_DIR/PremMore.html
if [$? -ep 0 ];then
	echo "trimming of low_quality reads is completed"
else
	echo "step 2 trim_low_read  " && exit 1 
fi 


############################################################################### step 3  checking the trimmed file for quality score ###############################################################################################
#
echo "step 3 --Read 1 --- checking the quality score " 

fastqc $TRIM_DIR/PremMore_trim_1_fq.gz  -o $FASTQC_DIR/AFTERTRIM 
if [$? -ep 0 ];then
	echo "fastqc for trimmed file is compeleted"

else
	echo "quality checking of trim file 1 is failed" && exit 1 
fi


 #                  #                         #                      #          # 
 
 echo " -- read 2 -----------"

 fastqc $TRIM_DIR/PremMore_trim_1_fq.gz -o $FASTQC_DIR/AFTERTRIM
 if [$? -ep 0];then
	 echo "fastqc for trimmend file completed "
 else
	 echo"quality checking of trimmed file 2 is failed" && exit 1 
 fi



 ############################################################################################step 4 alignment of two read via bowtie2  ####################################################################

 echo "step 4 alignment WITH GENOME  AND SORT THE BAM FILE " 
 bowtie2 -x $BOWTIE_INDEX -1 $TRIM_DIR/PremMore_trim_1_fq.gz  -2 $TRIM_DIR/PremMore_trim_2_fq.gz  |samtools view -b | samtools sort -o $sorted/PremMore_sorted.bam
 if [$? -ep 0];then
	 echo "alignment and sorting  is completed"
 else
	echo "alignment and sorting  is failed" && exit 1  

 fi

 ############################################################################step 5   Add read group to the BAM file ########################################################## 
 #
 
 echo "# Step 5: Add read group to the BAM file " 
 gatk AddOrReplaceReadGroups -I $SORTED/PremMore_sorted.bam -O $SORTED/PremMore_sorted_withRG.bam -RGID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM sample1
if [$? -ep 0 ];then
	echo "addingreplacement replacementgroup is completed"
else
	echo ""addingreplacement replacementgroup is failed"
fi


