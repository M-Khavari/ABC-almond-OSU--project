
###### 
mkdir ABC_gmeBS
mkdir BF_samples
cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/1_Raw

##############
mkdir 1_Raw  2_ Pre_Trim_QC  3_Trimming  4_Post_Trim_QC  Software

#####
 mkdir Software
cd Software
# Bismark container
apptainer pull bismark.sif docker://quay.io/biocontainers/bismark:0.24.0--hdfd78af_0
# fastqc container
apptainer pull fastqc.sif docker://quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1
# multiqc container
apptainer pull multiqc.sif docker://quay.io/biocontainers/multiqc:1.13a--pyhdfd78af_1
# samtools container
apptainer pull samtools.sif docker://quay.io/biocontainers/samtools:1.16.1--h6899075_1
# BBMap binaries
wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz -O - | tar -xz
# gemBS-rs
apptainer pull gemBS.sif docker://heathsc/gembs-rs:latest
# BEDtools
apptainer pull BEDtools.sif docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6

apptainer pull vcftools.sif docker://quay.io/biocontainers/vcftools:0.1.16--he513fc3_4

# bcftools  
### https://hub.docker.com/r/staphb/bcftools/tags
apptainer pull bcftools.sif docker://quay.io/staphb/bcftools:latest


##############
# 2. Merging FASTQ files
cd ..
cd 2_Merging
mkdir R1 R2
find  /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/1_Raw/ -type f -name "*R1*gz" -exec cp {} R1/ \;
find /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/1_Raw/ -type f -name "*R2*gz" -exec cp {} R2/ \;
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/2_Merging/R1 | sed -E 's/_L.*fastq.*//' | uniq | sort > R1/list.txt
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/2_Merging/R2 | sed -E 's/_L.*fastq.*//' | uniq | sort > R2/list.txt

#####################
print list.txt
## removing extra lines
grep -n list list.txt 
sed '56d' list.txt 
####################
cd 3_Pre_Trim_QC
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/2_Merging/R1/*R1*gz | sort > list1.txt
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/2_Merging/R2/*R2*gz | sort > list2.txt
mkdir PreQC_Report

while IFS= read R1 && IFS= read R2 <&3
do
 ID1=`echo $R1 | xargs basename | sed -E 's/\.fastq\.gz//'`
 ID2=`echo $R2 | xargs basename | sed -E 's/\.fastq\.gz//'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=Pre_Trim_QC_PE'
 echo '#SBATCH --time=00:15:00'
 echo '#SBATCH --ntasks=28'
 echo '#SBATCH --exclusive'
 echo '#SBATCH --mail-type=BEGIN,END,FAIL'
 echo '#SBATCH --account=PAS0471'
echo 'cd  /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/3_Pre_Trim_QC'
echo "apptainer exec  /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/Software/fastqc.sif \
 fastqc -t 28 -o PreQC_Report \
 $R1 $R2" 
 } > fastqc_pre_$ID1.sh
sbatch fastqc_pre_$ID1.sh
done < list1.txt 3< list2.txt

# Getting compiled report of QC
apptainer exec ../Software/multiqc.sif multiqc -s -f --interactive PreQC_Report/

###########
# 4. Trimming
cd ..
cd 4_Trimming
mkdir BBMapOut
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/2_Merging/R1/*R1*gz | sort > list_R1.txt
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/2_Merging/R2/*R2*gz | sort > list_R2.txt

while IFS= read R1 && IFS= read R2 <&3
do
 ID1=`echo $R1 | xargs basename | sed -E 's/\.fastq\.gz//'`
 ID2=`echo $R2 | xargs basename | sed -E 's/\.fastq\.gz//'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=Trimming_PE'
 echo '#SBATCH --time=00:05:00'
 echo '#SBATCH --ntasks=28'
 echo '#SBATCH --exclusive'
 echo '#SBATCH --mail-type=BEGIN,END,FAIL'
 echo '#SBATCH --account=PAS1755'
 echo 'cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/4_Trimming'
 echo "../Software/bbmap/bbduk.sh in1=$R1 in2=$R2 \
 out1=BBMapOut/$ID1.fastq.gz \
 out2=BBMapOut/$ID2.fastq.gz \
 qtrim=rl trimq=20 ktrim=r k=23 mink=11 hdist=1 tpe tbo tpe t=28 \
 ref=/fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/Software/bbmap/resources/adapters.fa" 
} > trimming_$ID1.sh
 sbatch trimming_$ID1.sh
done < list_R1.txt 3< list_R2.txt
##################################################
# 5. Post Trim QC for EM-seq
cd ..
cd 5_Post_Trim_QC
mkdir PostQC_Report
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/4_Trimming/BBMapOut/*R1*gz | sort > list1.txt
ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/4_Trimming/BBMapOut/*R2*gz | sort > list2.txt

while IFS= read R1 && IFS= read R2 <&3
do
 ID1=`echo $R1 | xargs basename | sed -E 's/\.fastq\.gz//'`
 ID2=`echo $R2 | xargs basename | sed -E 's/\.fastq\.gz//'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=Post_Trim_QC'
 echo '#SBATCH --time=00:15:00'
 echo '#SBATCH --ntasks=28'
 echo '#SBATCH --exclusive'
 echo '#SBATCH --account=PAS0471'
 echo 'cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/5_Post_Trim_QC'
 echo "apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/Software/fastqc.sif \
 fastqc -t 28 -o PostQC_Report \
 $R1 $R2"
 } > fastqc_post_$ID1.sh
 sbatch fastqc_post_$ID1.sh
done < list1.txt 3< list2.txt

# Getting compiled report of QC
apptainer exec ../Software/multiqc.sif multiqc -s -f --interactive PostQC_Report/

###################
gemBS running

###########

# 9.1 Using gemBS-rs

### aaranging samples ID
grep -e '-N\|-AN\|-WN\|_N\|_WN\|_AN' list_2.txt > list_Nonp.txt
sinteractive -c 28 -A PAS1755 -t 01:00:00 J- gemBS

####################################################################
https://github.com/heathsc/gemBS-rs/blob/master/etc/config_scripts/IHEC_standard.conf
###############
cd in_gemBS
mkdir Raw_data
#### ls /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/4_Trimming/BBMapOut/*.gz | sort > renaming_list.txt
cp /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/4_Trimming/BBMapOut/*.gz Raw_data/
ls -d "$PWD"/Raw_data/* | sort > list_files.txt

# Prepare list_BH.csv
# Prepare Config.txt file

tee -a gemBS_ABC.sh <<EOF
#!/bin/bash
#SBATCH --account=PAS1755
#SBATCH --job-name=gemBS
#SBATCH --mail-type=END,FAIL
#SBATCH --time=47:59:59
#SBATCH --ntasks=48
#SBATCH --partition=hugemem

cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/6_DMC/gembs_runnning
apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/Software/gemBS.sif /gemBS-rs/gemBS prepare --config Config.txt --text-metadata ID_gembs_NBF_2.csv 
#apptainer exec gemBS.sif /gemBS-rs/gemBS index -b
#apptainer exec gemBS.sif /gemBS-rs/gemBS map --bs 
#apptainer exec gemBS.sif /gemBS-rs/gemBS call -u 
#apptainer exec gemBS.sif /gemBS-rs/gemBS extract -C -N -S 
#apptainer exec gemBS.sif /gemBS-rs/gemBS report -P 
apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/NBF_samples/Software/gemBS.sif /gemBS-rs/gemBS run
scontrol show job=$SLURM_JOB_ID

EOF
# And submit the job
sbatch gemBS_ABC.sh

##################################################################################
Configuration file 
########
# Directory definitions
#
# Note that we can use environment variables ($HOME in this case)
# and previously defined variables can be used in subsequent definitions

base = /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/6_DMC/gembs_runnning

ref_dir = ${base}/Genome
reference = ${ref_dir}/Almond.fasta

index_dir = ${base}/index
sequence_dir = ${base}/fastq/@SAMPLE    # @SAMPLE and @BARCODE are special
bam_dir = ${base}/mapping/@BARCODE      # variables that are replaced with
bcf_dir = ${base}/calls/@BARCODE        # the sample name or barcode being
extract_dir = ${base}/extract/@BARCODE  # worked on during gemBS operation
report_dir = ${base}/report

# For using CRAM
populate_cache = true
make_cram = true

# General project info
project = Almond_methyl_ABC
species = Almond

# Default parameters
threads = 48
cores = 28
jobs = 28

[mapping]

# Set names of spiked in conversion controls
underconversion_sequence = Lambda_NEB_N3011
overconversion_sequence = pUC19_Methylated

[calling]

mapq_threshold = 10
qual_threshold = 13
reference_bias = 2
left_trim = 5
right_trim = 0
keep_improper_pairs = False
keep_duplicates = False
haploid = False
conversion = 0.01,0.05
remove_individual_bcfs = True

# Contigs smaller than contig_pool_limit will be called together
contig_pool_limit = 25000000

[extract]

strand_specific = True
bigWig_strand_specific = True
phred_threshold = 10
make_cpg = True
make_non_cpg = True
make_bedmethyl = True



#############################

find /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/6_DMC/gembs_runnning/extract -type f -name '*cpg.bed.gz' -exec cp {} . \;

apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/BEDtools.sif bedtools bamtobed -i BH_C_N_BH10R_cpg.bed -g /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/6_DMC/gembs_runnning/Genome/Almond.fasta > BH_C_N_BH10R_cpg.bam

###############################

filtering 
#############

apptainer pull bcftools.sif docker://quay.io/biocontainers/bcftools:v1.9-1-deb_cv1

mkdir 7_filtering
cd 7_filtering/

find /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/6_DMC/old_gembs_runnning/calls/ -type f -name '*.bcf' -exec cp {} . \;
find /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/6_DMC/gembs_runnning/calls/ -type f -name '*.bcf.csi' -exec cp {} .\ ;
@####################

Normalization code

####
mkdir norm_dic

##########

ls *.bcf | sort > bcf_all_samples


while IFS= read R1
do
 ID1=`echo $R1 | xargs basename | sed -E 's///'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=BCF_norm'
 echo '#SBATCH --time=00:50:00'
 echo '#SBATCH --ntasks=1'
 echo '#SBATCH --mail-type=END,FAIL'
 echo '#SBATCH --account=PAS0471'

 echo '/fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/7_filtering'
 
 echo " apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools norm -m-any $R1 | apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools norm -Ov --check-ref w -f /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/6_DMC/old_gembs_runnning/Genome/Almond.fasta > norm_dic/$R1.norm.vcf"
  
} > BCF_normal_$ID1.sh
 sbatch BCF_normal_$ID1.sh
done < bcf_all_samples
##########################################

cd norm_dic
ls *.vcf | sort > list_norm.txt
mkdir run_merge

while IFS= read R1
do
 ID1=`echo $R1 | xargs basename | sed -E 's/.vcf//'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=VCF_BCF'
 echo '#SBATCH --time=00:20:00'
 echo '#SBATCH --ntasks=1'
 echo '#SBATCH --mail-type=END,FAIL'
 echo '#SBATCH --account=PAS0471'
cd 
 echo 'cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/7_filtering/norm_dic'
 echo "apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools view $R1 -o run_merge/$R1.f1.bcf"
 
} > BCF_convert_$ID1.sh
 sbatch BCF_convert_$ID1.sh
done < list_norm.txt
######################################################
################### indexing bcf files
cd run_merge
ls *.bcf.norm.vcf.f1.bcf | sort > list_bcf_index.txt

while IFS= read R1
do
 ID1=`echo $R1 | xargs basename | sed -E 's/.bcf.norm.vcf.f1.bcf//'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=BCF_index'
 echo '#SBATCH --time=00:20:00'
 echo '#SBATCH --ntasks=1'
 echo '#SBATCH --mail-type=END,FAIL'
 echo '#SBATCH --account=PAS0471'

 echo 'cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/7_filtering/norm_dic/run_merge'
 echo "apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools index $R1 -o $R1.csi"
 
} > BCF_index_$ID1.sh
 sbatch BCF_index_$ID1.sh
done < list_bcf_index.txt

###################################################

while IFS= read R1
do
 ID1=`echo $R1 | xargs basename | sed -E 's/.bcf.norm.vcf.f1.bcf//'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=BCF_index'
 echo '#SBATCH --time=00:20:00'
 echo '#SBATCH --ntasks=1'
 echo '#SBATCH --mail-type=END,FAIL'
 echo '#SBATCH --account=PAS0471'

 echo 'cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/7_filtering/norm_dic/run_merge'
 echo "apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools index $R1 -o $R1.csi"
 
} > BCF_index_$ID1.sh
 sbatch BCF_index_$ID1.sh
done < list_bcf_index.txt

###############################
# Define the input file
input_file="ABC_leaf_Failure_65_40__55555_tab_imputed.hmp.txt"

# Extract the allele columns and filter out duplicate alleles
awk -F'\t' '{for(i=12;i<=NF;i++) print $i}' $input_file | sort | uniq -d

###################################################################### separating files for leafing failure vs NBF 

# grep '^I\|^M\|^BH' list_bcf.txt | sed 's/.csi/ /g' > leafing_failure_samples.txt
 ls *.norm.vcf.f1.bcf > leafing_failure_samples.txt 
sinteractive -c 28 -A PAS1755 -t 01:00:00 J- merging_bcf

#### make a job and give it time for 

#######
#!/bin/bash
#SBATCH -J ABC_merging_samples
#SBATCH --time=12:30:00
#SBATCH --exclusive
#SBATCH --ntasks=28
#SBATCH --mem=28gb
#SBATCH --account=PAS1755
#SBATCH --mail-type=BEGIN,END,FAIL
####### step 2
cd /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/7_filtering/norm_dic/run_merge
apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools merge -o merged_all_leafing_Failure_samples.bcf --file-list leafing_failure_samples.txt --threads 28
########
apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools merge -o merged_all_leafing_Failure_samples.bcf --file-list leafing_failure_samples.txt 8 
55
##################

#ls *Stukey*bcf | sort > list_stuckey.txt
#grep '^I\|^N\|^C\|Stukey' list_bcf.txt | sed 's/.csi/ /g' > Bud_failure_samples.txt
#grep '^I\|^N\|^C' list_bcf.txt | sed 's/.csi/ /g' > Bud_failure_samples.txt
## merging Stukey and bud faliure samples
#
#sed  '/Stukey/d' Bud_failure_samples.txt > Bud_failure_not_Stukey.txt
****
#sinteractive -c 28 -A PAS1755 -t 01:00:00 J- merging_bcf
#nohup apptainer exec /fs/scratch/PAS1755/Marziye/ABC/gemBS/all_5_cultivar/9_DMR_calling/in_gemBS_all_pairing/Software/bcftools.sif bcftools merge -o merged_all_samples.bcf --file-list list_pool.txt --threads 28
2> /dev/ null < /dev/null & 

apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools merge -o merged_Leafing_fail_BH_M_I_samples.bcf --file-list leafing_failure_samples.txt --threads 28 
apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools merge -o merged_Leafing_fail_BH_M_I_samples.bcf --file-list Bud_samples.txt --threads 28 
apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools merge -o merged_Bud_failure_N_C_I_Stukey_samples.bcf --file-list Bud_failure_samples.txt --threads 28
apptainer exec /fs/scratch/PAS1755/Marziye/ABC_gmeBS/ABC_gmeBS/BF_samples/Software/bcftools.sif bcftools merge -o merged_Bud_failure_N_C_I_not_stukey_samples.bcf --file-list Bud_failure_not_Stukey.txt --threads 28


ls *.bcf | sort > list_pool.txt
######## indexing for contigs and mits

while IFS= read R1
do
 ID1=`echo $R1 | xargs basename | sed -E 's/.bcf.filtered.vcf//'`
 {
 echo '#!/bin/bash'
 echo '#SBATCH -J ondemand/sys/myjobs/basic_sequential'
 echo '#SBATCH --job-name=BCF_index'
 echo '#SBATCH --time=00:20:00'
 echo '#SBATCH --ntasks=1'
 echo '#SBATCH --mail-type=END,FAIL'
 echo '#SBATCH --account=PAS0471'

 echo 'cd /fs/scratch/PAS1755/Marziye/ABC/gemBS/all_5_cultivar/9_DMR_calling/in_gemBS_all_Raw/DMC_calling/pools_samples'
 echo "apptainer exec /fs/scratch/PAS1755/Marziye/ABC/gemBS/all_5_cultivar/9_DMR_calling/in_gemBS_all_pairing/Software/bcftools.sif bcftools index $R1 -o $R1.csi"
 
} > BCF_chr_pulling_$ID1.sh
 sbatch BCF_chr_pulling_$ID1.sh
done < list_pool.txt




