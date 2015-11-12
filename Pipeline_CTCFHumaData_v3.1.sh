echo "#################################################################################################################################################################################################"
echo "#################################################################################################################################################################################################"
echo "##############################################################        Test:CTCF Human Data for Pipeline Development      ########################################################################"
echo "#################################################################################################################################################################################################"
echo "#################################################################################################################################################################################################"
source /export/apps/rsat/RSAT_config.bashrc
ALEA=/export/storage/users/wsantana/ALEA
FileChIPseq=FileListChIPseq.txt #MISSING SPECIFICS
mkdir $ALEA
mkdir $ALEA/Data
mkdir $ALEA/Scripts
mkdir $ALEA/Software
mkdir $ALEA/variation-scan-git
mkdir $ALEA/Stats
mkdir $ALEA/Data/SRA_Files
mkdir $ALEA/Data/ChIPseq_Files
mkdir $ALEA/Data/WASP_SNPs

echo "##################################################################################################################################################################################################"
echo "#######################################################################  Part 1. Data retrieving  and Formating  #################################################################################"
echo "##################################################################################################################################################################################################"

echo "#################################################################################"
echo "#####################  Retrieve VCF files for CEU and YRI trios #################"
echo "#################################################################################"

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/trio/snps/CEU.trio.2010_03.genotypes.vcf.gz --directory-prefix=$ALEA/Data
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/trio/snps/YRI.trio.2010_03.genotypes.vcf.gz --directory-prefix=$ALEA/Data

gunzip $ALEA/Data/CEU.trio.2010_03.genotypes.vcf.gz
gunzip $ALEA/Data/YRI.trio.2010_03.genotypes.vcf.gz

echo "#################################################################################"
echo "#####################  Retrieve CTCF Files and Quality Filters  #################"
echo "#################################################################################"
sh $ALEA/Scripts/DownloadFiles_v2.0.sh
ALEA=/export/storage/users/wsantana/ALEA
while read line
do	

	NewName=$(echo $line|perl -pe 's/.+\|(.+)\|.+/$1/')
	trimmomatic SE -threads 40 -phred33 -trimlog $ALEA/Data/ChIPseq_Files/Trim20L20Q.$NewName.log $ALEA/Data/ChIPseq_Files/$NewName.fastq $ALEA/Data/ChIPseq_Files/Trim20L20Q.$NewName.fastq SLIDINGWINDOW:4:20 MINLEN:20
done < $ALEA/Data/FileListChIPseq.txt


mv $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep1_36680.fastq $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq
cat $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2_36681.fastq $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2_36792.fastq > $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq
mv $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep3_36682.fastq $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq

cat $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep1_36683.fastq $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep1_36791.fastq > $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep1.ALL.fastq
cat $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep2_36684.fastq $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep2_36793.fastq > $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep2.ALL.fastq

mv $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep1_36685.fastq $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq
mv $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep2_36686.fastq $ALEA/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq

cat $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1_36687.fastq $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1_36795.fastq > $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq
mv $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep2_36688.fastq $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq

cat $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1_36689.fastq $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1_36801.fastq > $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq
cat $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2_36690.fastq $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2_36803.fastq > $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq

mv $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Father_CTCF_Rep1_36691.fastq $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Father_CTCF_Rep1.ALL.fastq
mv $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Father_CTCF_Rep2_36692.fastq $ALEA/ChIPseq_Files/Trim20L20Q.YRI_Father_CTCF_Rep2.ALL.fastq

echo "#################################################################################"
echo "###########  Subsetting VCF data for each individual from both trios  ###########"
echo "###########                WARNING!: VCFtools needed                  ###########"
echo "#################################################################################"
vcf-subset -e -c NA12891 $ALEA/Data/CEU.trio.2010_03.genotypes.vcf > $ALEA/Data/CEU_Father.genotypes.hg18.vcf
vcf-subset -e -c NA12892 $ALEA/Data/CEU.trio.2010_03.genotypes.vcf > $ALEA/Data/CEU_Mother.genotypes.hg18.vcf
vcf-subset -e -c NA12878 $ALEA/Data/CEU.trio.2010_03.genotypes.vcf > $ALEA/Data/CEU_Daughter.genotypes.hg18.vcf
vcf-subset -e -c NA19239 $ALEA/Data/YRI.trio.2010_03.genotypes.vcf > $ALEA/Data/YRI_Father.genotypes.hg18.vcf
vcf-subset -e -c NA19238 $ALEA/Data/YRI.trio.2010_03.genotypes.vcf > $ALEA/Data/YRI_Mother.genotypes.hg18.vcf
vcf-subset -e -c NA19240 $ALEA/Data/YRI.trio.2010_03.genotypes.vcf > $ALEA/Data/YRI_Daughter.genotypes.hg18.vcf

echo "#################################################################################"
echo "##              Remapping VCFs genome coordinates from hg18 to hg19            ##"
echo "## WARNING!: CrossMap,chain hg18ToHg19.over.chain and hg19 fasta genome needed ##"
echo "#################################################################################"

CrossMap.py vcf $ALEA/Data/hg18ToHg19.over.chain.gz $ALEA/Data/CEU_Father.genotypes.hg18.vcf $ALEA/Data/hg19.fa $ALEA/Data/CEU_Father.genotypes.hg19.vcf
CrossMap.py vcf $ALEA/Data/hg18ToHg19.over.chain.gz $ALEA/Data/CEU_Mother.genotypes.hg18.vcf $ALEA/Data/hg19.fa $ALEA/Data/CEU_Mother.genotypes.hg19.vcf
CrossMap.py vcf $ALEA/Data/hg18ToHg19.over.chain.gz $ALEA/Data/CEU_Daughter.genotypes.hg18.vcf $ALEA/Data/hg19.fa $ALEA/Data/CEU_Daughter.genotypes.hg19.vcf
CrossMap.py vcf $ALEA/Data/hg18ToHg19.over.chain.gz $ALEA/Data/YRI_Daughter.genotypes.hg18.vcf $ALEA/Data/hg19.fa $ALEA/Data/YRI_Daughter.genotypes.hg19.vcf
CrossMap.py vcf $ALEA/Data/hg18ToHg19.over.chain.gz $ALEA/Data/YRI_Father.genotypes.hg18.vcf $ALEA/Data/hg19.fa $ALEA/Data/YRI_Father.genotypes.hg19.vcf
CrossMap.py vcf $ALEA/Data/hg18ToHg19.over.chain.gz $ALEA/Data/YRI_Mother.genotypes.hg18.vcf $ALEA/Data/hg19.fa $ALEA/Data/YRI_Mother.genotypes.hg19.vcf

echo "#################################################################################"
echo "######    Extracting PHASED SNPs and Filtering out the other variants      ######"
echo "#################################################################################"

awk '/#/{print;next}{if($0 !~ /\//){print}}' $ALEA/Data/CEU_Father.genotypes.hg19.vcf > $ALEA/Data/CEU_Father.genotypes.hg19.phased.vcf
awk '/#/{print;next}{if($0 !~ /\//){print}}' $ALEA/Data/CEU_Mother.genotypes.hg19.vcf > $ALEA/Data/CEU_Mother.genotypes.hg19.phased.vcf
awk '/#/{print;next}{if($0 !~ /\//){print}}' $ALEA/Data/CEU_Daughter.genotypes.hg19.vcf > $ALEA/Data/CEU_Daughter.genotypes.hg19.phased.vcf

awk '/#/{print;next}{if($0 !~ /\//){print}}' $ALEA/Data/YRI_Father.genotypes.hg19.vcf > $ALEA/Data/YRI_Father.genotypes.hg19.phased.vcf
awk '/#/{print;next}{if($0 !~ /\//){print}}' $ALEA/Data/YRI_Mother.genotypes.hg19.vcf > $ALEA/Data/YRI_Mother.genotypes.hg19.phased.vcf
awk '/#/{print;next}{if($0 !~ /\//){print}}' $ALEA/Data/YRI_Daughter.genotypes.hg19.vcf > $ALEA/Data/YRI_Daughter.genotypes.hg19.phased.vcf

awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' $ALEA/Data/CEU_Father.genotypes.hg19.phased.vcf> $ALEA/Data/CEU_Father.genotypes.hg19.phased2.vcf
awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' $ALEA/Data/CEU_Mother.genotypes.hg19.phased.vcf> $ALEA/Data/CEU_Mother.genotypes.hg19.phased2.vcf
awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' $ALEA/Data/CEU_Daughter.genotypes.hg19.phased.vcf> $ALEA/Data/CEU_Daughter.genotypes.hg19.phased2.vcf

awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' $ALEA/Data/YRI_Father.genotypes.hg19.phased.vcf> $ALEA/Data/YRI_Father.genotypes.hg19.phased2.vcf
awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' $ALEA/Data/YRI_Mother.genotypes.hg19.phased.vcf> $ALEA/Data/YRI_Mother.genotypes.hg19.phased2.vcf
awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' $ALEA/Data/YRI_Daughter.genotypes.hg19.phased.vcf> $ALEA/Data/YRI_Daughter.genotypes.hg19.phased2.vcf



echo "##################################################################################################################################################################################################"
echo "#####################################################################   Part 2. Haploid genome reconstruction  ###################################################################################"
echo "##################################################################################################################################################################################################"
ALEA=/export/storage/users/wsantana/ALEA

echo "#####################################################################################"
echo "##       Extracting Haplotypes and Reconstructing in silico Haploid Genomes        ##"
echo "#####################################################################################"
perl $ALEA/Scripts/Hg19ChromosomesFiltering_v1.2.pl -f $ALEA/Data/hg19.fa -gender female > $ALEA/Data/hg19.female.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering_v1.2.pl -f $ALEA/Data/hg19.fa -gender male > $ALEA/Data/hg19.male.fasta

perl $ALEA/variation-scan-git/vcf_to_vcfbed.pl -v $ALEA/Data/CEU_Father.genotypes.hg19.phased2.vcf > $ALEA/Data/CEU_Father.hg19.bed
perl $ALEA/variation-scan-git/vcf_to_vcfbed.pl -v $ALEA/Data/CEU_Mother.genotypes.hg19.phased2.vcf > $ALEA/Data/CEU_Mother.hg19.bed
perl $ALEA/variation-scan-git/vcf_to_vcfbed.pl -v $ALEA/Data/CEU_Daughter.genotypes.hg19.phased2.vcf > $ALEA/Data/CEU_Daughter.hg19.bed

perl $ALEA/variation-scan-git/vcf_to_vcfbed.pl -v $ALEA/Data/YRI_Father.genotypes.hg19.phased2.vcf > $ALEA/Data/YRI_Father.hg19.bed
perl $ALEA/variation-scan-git/vcf_to_vcfbed.pl -v $ALEA/Data/YRI_Mother.genotypes.hg19.phased2.vcf > $ALEA/Data/YRI_Mother.hg19.bed
perl $ALEA/variation-scan-git/vcf_to_vcfbed.pl -v $ALEA/Data/YRI_Daughter.genotypes.hg19.phased2.vcf > $ALEA/Data/YRI_Daughter.hg19.bed


perl $ALEA/variation-scan-git/vcfbed_to_fasta.pl -fasta $ALEA/Data/hg19.male.fasta -bed $ALEA/Data/CEU_Father.hg19.bed -o_ref $ALEA/Data/CEU_Father.Hap2.fasta -o $ALEA/Data/CEU_Father.Hap1.fasta 
perl $ALEA/variation-scan-git/vcfbed_to_fasta.pl -fasta $ALEA/Data/hg19.female.fasta -bed $ALEA/Data/CEU_Mother.hg19.bed -o_ref $ALEA/Data/CEU_Mother.Hap2.fasta -o $ALEA/Data/CEU_Mother.Hap1.fasta 
perl $ALEA/variation-scan-git/vcfbed_to_fasta.pl -fasta $ALEA/Data/hg19.female.fasta -bed $ALEA/Data/CEU_Daughter.hg19.bed -o_ref $ALEA/Data/CEU_Daughter.Hap2.fasta -o $ALEA/Data/CEU_Daughter.Hap1.fasta

perl $ALEA/variation-scan-git/vcfbed_to_fasta.pl -fasta $ALEA/Data/hg19.male.fasta -bed $ALEA/Data/YRI_Father.hg19.bed -o_ref $ALEA/Data/YRI_Father.Hap2.fasta -o $ALEA/Data/YRI_Father.Hap1.fasta 
perl $ALEA/variation-scan-git/vcfbed_to_fasta.pl -fasta $ALEA/Data/hg19.female.fasta -bed $ALEA/Data/YRI_Mother.hg19.bed -o_ref $ALEA/Data/YRI_Mother.Hap2.fasta -o $ALEA/Data/YRI_Mother.Hap1.fasta 
perl $ALEA/variation-scan-git/vcfbed_to_fasta.pl -fasta $ALEA/Data/hg19.female.fasta -bed $ALEA/Data/YRI_Daughter.hg19.bed -o_ref $ALEA/Data/YRI_Daughter.Hap2.fasta -o $ALEA/Data/YRI_Daughter.Hap1.fasta 

echo "#####################################################################################"
echo "##  Formatting Haploid Genomes:Extracting Fasta entries just for full chromosomes  ##"
echo "##            WARNING!:Chr1-Chr22 and ChrX only i.e.just female genomes            ##"
echo "#####################################################################################"
ALEA=/export/storage/users/wsantana/ALEA
mkdir $ALEA/Scripts
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/CEU_Daughter.Hap1.fasta > $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/CEU_Daughter.Hap2.fasta > $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/CEU_Mother.Hap1.fasta > $ALEA/Data/CEU_Mother.Chr.Hap1.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/CEU_Mother.Hap2.fasta > $ALEA/Data/CEU_Mother.Chr.Hap2.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/YRI_Daughter.Hap1.fasta > $ALEA/Data/YRI_Daughter.Chr.Hap1.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/YRI_Daughter.Hap2.fasta > $ALEA/Data/YRI_Daughter.Chr.Hap2.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/YRI_Mother.Hap1.fasta > $ALEA/Data/YRI_Mother.Chr.Hap1.fasta
perl $ALEA/Scripts/Hg19ChromosomesFiltering.pl -f $ALEA/Data/YRI_Mother.Hap2.fasta > $ALEA/Data/YRI_Mother.Chr.Hap2.fasta

echo "#####################################################################################"
echo "##                          WASP SNPs Folder construction                          ##"
echo "#####################################################################################"
mkdir $ALEA/Data/WASP_SNPs/CEU_Daughter
mkdir $ALEA/Data/WASP_SNPs/CEU_Father
mkdir $ALEA/Data/WASP_SNPs/CEU_Mother
mkdir $ALEA/Data/WASP_SNPs/YRI_Daughter
mkdir $ALEA/Data/WASP_SNPs/YRI_Father
mkdir $ALEA/Data/WASP_SNPs/YRI_Mother

awk 'BEGIN{OFS="\t"};{out="/export/storage/users/wsantana/ALEA/Data/WASP_SNPs/CEU_Daughter/"$1".snps.txt";print $3,$5,$6 > out}' $ALEA/Data/CEU_Daughter.hg19.bed ; gzip $ALEA/Data/WASP_SNPs/CEU_Daughter/*.snps.txt
awk 'BEGIN{OFS="\t"};{out="/export/storage/users/wsantana/ALEA/Data/WASP_SNPs/CEU_Father/"$1".snps.txt";print $3,$5,$6 > out}' $ALEA/Data/CEU_Father.hg19.bed; gzip $ALEA/Data/WASP_SNPs/CEU_Father/*.snps.txt
awk 'BEGIN{OFS="\t"};{out="/export/storage/users/wsantana/ALEA/Data/WASP_SNPs/CEU_Mother/"$1".snps.txt";print $3,$5,$6 > out}' $ALEA/Data/CEU_Mother.hg19.bed; gzip $ALEA/Data/WASP_SNPs/CEU_Mother/*.snps.txt
awk 'BEGIN{OFS="\t"};{out="/export/storage/users/wsantana/ALEA/Data/WASP_SNPs/YRI_Daughter/"$1".snps.txt";print $3,$5,$6 > out}' $ALEA/Data/YRI_Daughter.hg19.bed; gzip $ALEA/Data/WASP_SNPs/YRI_Daughter/*.snps.txt
awk 'BEGIN{OFS="\t"};{out="/export/storage/users/wsantana/ALEA/Data/WASP_SNPs/YRI_Father/"$1".snps.txt";print $3,$5,$6 > out}' $ALEA/Data/YRI_Father.hg19.bed; gzip $ALEA/Data/WASP_SNPs/YRI_Father/*.snps.txt
awk 'BEGIN{OFS="\t"};{out="/export/storage/users/wsantana/ALEA/Data/WASP_SNPs/YRI_Mother/"$1".snps.txt";print $3,$5,$6 > out}' $ALEA/Data/YRI_Mother.hg19.bed; gzip $ALEA/Data/WASP_SNPs/YRI_Mother/*.snps.txt



echo "##################################################################################################################################################################################################"
echo "#############################################################################   Part 3. Read Alignment   #########################################################################################"
echo "##################################################################################################################################################################################################"

echo "#####################################################################################"
echo "#########################  BWA alignment for CTCF sequences  ########################"
echo "#####################################################################################"
bwa index -a bwtsw $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta
bwa index -a bwtsw $ALEA/Data/CEU_Mother.Chr.Hap2.fasta
bwa index -a bwtsw $ALEA/Data/YRI_Daughter.Chr.Hap2.fasta
bwa index -a bwtsw $ALEA/Data/YRI_Mother.Chr.Hap2.fasta
bwa index -a bwtsw $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta
bwa index -a bwtsw $ALEA/Data/CEU_Mother.Chr.Hap1.fasta
bwa index -a bwtsw $ALEA/Data/YRI_Daughter.Chr.Hap1.fasta
bwa index -a bwtsw $ALEA/Data/YRI_Mother.Chr.Hap1.fasta

bwa aln -t 40 $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sai
bwa aln -t 40 $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sai
bwa aln -t 40 $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sai
bwa aln -t 40 $ALEA/Data/CEU_Mother.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sai
bwa aln -t 40 $ALEA/Data/CEU_Mother.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sai
bwa aln -t 40 $ALEA/Data/YRI_Daughter.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sai
bwa aln -t 40 $ALEA/Data/YRI_Daughter.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sai
bwa aln -t 40 $ALEA/Data/YRI_Mother.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sai
bwa aln -t 40 $ALEA/Data/YRI_Mother.Chr.Hap2.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sai
bwa aln -t 40 $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sai
bwa aln -t 40 $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sai
bwa aln -t 40 $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sai
bwa aln -t 40 $ALEA/Data/CEU_Mother.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sai
bwa aln -t 40 $ALEA/Data/CEU_Mother.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sai
bwa aln -t 40 $ALEA/Data/YRI_Daughter.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sai
bwa aln -t 40 $ALEA/Data/YRI_Daughter.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sai
bwa aln -t 40 $ALEA/Data/YRI_Mother.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sai
bwa aln -t 40 $ALEA/Data/YRI_Mother.Chr.Hap1.fasta $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sai

bwa samse $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sam
bwa samse $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sam
bwa samse $ALEA/Data/CEU_Daughter.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sam
bwa samse $ALEA/Data/CEU_Mother.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sam
bwa samse $ALEA/Data/CEU_Mother.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sam
bwa samse $ALEA/Data/YRI_Daughter.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sam
bwa samse $ALEA/Data/YRI_Daughter.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sam
bwa samse $ALEA/Data/YRI_Mother.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sam
bwa samse $ALEA/Data/YRI_Mother.Chr.Hap2.fasta $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sam
bwa samse $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sam
bwa samse $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sam
bwa samse $ALEA/Data/CEU_Daughter.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sam
bwa samse $ALEA/Data/CEU_Mother.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sam
bwa samse $ALEA/Data/CEU_Mother.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sam
bwa samse $ALEA/Data/YRI_Daughter.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sam
bwa samse $ALEA/Data/YRI_Daughter.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sam
bwa samse $ALEA/Data/YRI_Mother.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sam
bwa samse $ALEA/Data/YRI_Mother.Chr.Hap1.fasta $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sai $ALEA/Data/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sam

samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sorted
samtools view -@ 40 -bS $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sam | samtools sort -@ 40 - $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sorted


samtools merge -@ 40 $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sorted.bam
samtools merge -@ 40 $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sorted.bam
samtools merge -@ 40 $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sorted.bam
samtools merge -@ 40 $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sorted.bam
samtools merge -@ 40 $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sorted.bam
samtools merge -@ 40 $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sorted.bam
samtools merge -@ 40 $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sorted.bam
samtools merge -@ 40 $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sorted.bam $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sorted.bam

samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap1.sorted.sam

samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sorted.sam
samtools view $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sorted.bam > $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sorted.sam


echo "##################################################################################################################################################################################################"
echo "################################################################################  Part 4. Peak Calling  ##########################################################################################"
echo "##################################################################################################################################################################################################"

echo "###################################################################################"
echo "#########  Peak Calling with SWEMBL NO -R 0.005 for each REPLICA and SAMPLE  ######"
echo "###################################################################################"
ALEA=/export/storage/users/wsantana/ALEA
mkdir $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.Peak 
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap2.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF.Hap2.Peak 
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap1.sorted.sam -S -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF.Hap1.Peak

echo "################################################################################"
echo "#########  Peak Calling with SWEMBL -R 0.005 for each REPLICA and SAMPLE  ######"
echo "################################################################################"
ALEA=/export/storage/users/wsantana/ALEA
mkdir $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.Peak 
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap2.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF.Hap2.Peak 
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap1.sorted.sam -S -R 0.005 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF.Hap1.Peak

echo "################################################################################"
echo "#########   Peak Calling with SWEMBL -R 0.01 for each REPLICA and SAMPLE  ######"
echo "################################################################################"
ALEA=/export/storage/users/wsantana/ALEA
mkdir $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap2.Peak 
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF.Hap2.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap2.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF.Hap2.Peak 
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.CEU_Mother_CTCF.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Daughter_CTCF.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF.Hap1.Peak
$ALEA/Software/SWEMBL.3.3.1/SWEMBL -i $ALEA/Data/Trim20L20Q.YRI_Mother_CTCF.Hap1.sorted.sam -S -R 0.01 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF.Hap1.Peak

echo "##################################################################################"
echo "########  Peak Calling with MACS2 for each REPLICA and SAMPLE for n-Summits ######"
echo "##################################################################################"
ALEA=/export/storage/users/wsantana/ALEA
mkdir $ALEA/Data/MACS2_Trim20L20Q
Directory=$ALEA/Data/MACS2_Trim20L20Q
names=("`ls -lhtr $ALEA/Data/*.bam|awk '{print $9}'|perl -pe 's/.+Data\/(.+)\.sorted\.bam\n/$1 /'`")
for name in $names
do
	echo ">Peak calling of file $name"
	echo " CMD: macs2 callpeak -t $ALEA/Data/$name.sorted.bam -f BAM -g hs -n $name -B --outdir $Directory --verbose 3 --call-summits"
	macs2 callpeak -t $ALEA/Data/$name.sorted.bam -f BAM -g hs -n $name -B --outdir $Directory --verbose 3 --call-summits
done

echo "################################################################################"
echo "######  Peak Calling with MACS2 for each REPLICA and SAMPLE for 1-Summit  ######"
echo "################################################################################"
ALEA=/export/storage/users/wsantana/ALEA
mkdir $ALEA/Data/MACS2_Trim20L20Q/Summit_1
Directory=$ALEA/Data/MACS2_Trim20L20Q/Summit_1
names=("`ls -lhtr $ALEA/Data/*.bam|awk '{print $9}'|perl -pe 's/.+Data\/(.+)\.sorted\.bam\n/$1 /'`")
for name in $names
do
	echo ">Peak calling of file $name"
	echo " CMD: macs2 callpeak -t $ALEA/Data/$name.sorted.bam -f BAM -g hs -n $name.Summit_1 -B --outdir $Directory"
	macs2 callpeak -t $ALEA/Data/$name.sorted.bam -f BAM -g hs -n $name.Summit_1 -B --outdir $Directory
done


echo "##################################################################################################################################################################################################"
echo "################################################################################  Part 5. Motif discovery  #######################################################################################"
echo "##################################################################################################################################################################################################"
ALEA=/export/storage/users/wsantana/ALEA
names=("`ls -lhtr $ALEA/Data/*.bam|awk '{print $9}'|perl -pe 's/.+Data\/(.+)\.sorted\.bam\n/$1 /'`")
#                   Script name: FetchSequences_CalledPeaks_v1.0.sh          #
echo "#######################################################################"
echo "#####  Fetch sequences for SWEMBL NO -R 0.005 with FULL interval  #####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|Region" $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk 'x=$2-1{print $1"\t"x"\t"$3}'|fetch-sequences -genome hg19 -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.Full_length.sorted.fasta -v &> $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.Full_length.sorted.fas.log &
done
echo "#######################################################################"
echo "##### Fetch sequences for SWEMBL NO -R 0.005 with 300bps interval #####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|Region" $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk 'x=int($10+0.5)-151,y=int($10+0.5)+150{print $1"\t"x"\t"y}'|fetch-sequences -genome hg19 -o $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.300bp_length.sorted.fasta -v &> $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.300bp_length.sorted.fas.log &
done
echo "#######################################################################"
echo "##### Fetch sequences for SWEMBL WITH -R 0.005 with FULL interval #####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|Region" $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk 'x=$2-1{print $1"\t"x"\t"$3}'|fetch-sequences -genome hg19 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.Full_length.sorted.fasta -v &> $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.Full_length.sorted.fas.log &
done
echo "#######################################################################"
echo "#####Fetch sequences for SWEMBL WITH -R 0.005 with 300bps interval#####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|Region" $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk 'x=int($10+0.5)-151,y=int($10+0.5)+150{print $1"\t"x"\t"y}'|fetch-sequences -genome hg19 -o $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.300bp_length.sorted.fasta -v &> $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.300bp_length.sorted.fas.log &
done
echo "#######################################################################"
echo "#####  Fetch sequences for SWEMBL WITH -R 0.01 with FULL interval #####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|^$\|Region" $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk 'x=$2-1{print $1"\t"x"\t"$3}'|fetch-sequences -genome hg19 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.Full_length.sorted.fasta -v &> $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.Full_length.sorted.fas.log &
done
echo "#######################################################################"
echo "##### Fetch sequences for SWEMBL WITH -R 0.01 with 300bps interval#####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|^$\|Region" $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk 'x=int($10+0.5)-151,y=int($10+0.5)+150{print $1"\t"x"\t"y}'|fetch-sequences -genome hg19 -o $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.300bp_length.sorted.fasta -v &> $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.300bp_length.sorted.fas.log &
done
echo "#######################################################################"
echo "#####   Fetch sequences for MACS2 n-Summits with 300bps interval  #####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|start" $ALEA/Data/MACS2_Trim20L20Q/$name\_peaks.xls|sort -nr -k8|awk 'x=$5-151,y=$5+150{print $1"\t"x"\t"y}'|fetch-sequences -genome hg19 -o $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.$name.300bp_length.sorted.fasta -v &> $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.$name.300bp_length.sorted.fas.log &
done
echo "#######################################################################"
echo "#####    Fetch sequences for MACS2 1-Summit with FULL interval    #####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|start" $ALEA/Data/MACS2_Trim20L20Q/Summit_1/$name.Summit_1_peaks.xls|sort -nr -k8|awk 'x=$2-1{print $1"\t"x"\t"$3}'|fetch-sequences -genome hg19 -o $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.Full_length.sorted.fasta -v &> $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.Full_length.sorted.fas.log &
done
echo "#######################################################################"
echo "#####   Fetch sequences for MACS2 1-Summit with 300bps interval   #####"
echo "#######################################################################"
for name in $names
do
	grep -v "#\|start" $ALEA/Data/MACS2_Trim20L20Q/Summit_1/$name.Summit_1_peaks.xls|sort -nr -k8|awk 'x=$5-151,y=$5+150{print $1"\t"x"\t"y}'|fetch-sequences -genome hg19 -o $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.300bp_length.sorted.fasta -v &> $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.300bp_length.sorted.fas.log &
done
#                   Script name: PeakMotifs_CalledPeaks_v1.0.sh          #
ALEA=/export/storage/users/wsantana/ALEA
names=("`ls -lhtr $ALEA/Data/*.bam|awk '{print $9}'|perl -pe 's/.+Data\/(.+)\.sorted\.bam\n/$1 /'`")
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL NO -R 0.005 with FULL interval ALL peaks                ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.Full_length.sorted.fasta -outdir $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.Full_length.All_peaks.Sorted.$name -prefix v1.Swembl.Full_length.All_peaks.Sorted.$name -title v1.Swembl.Full_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.Full_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL NO -R 0.005 with FULL interval TOP 600 peaks            ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.Full_length.sorted.fasta -outdir $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.Full_length.Top600_peaks.Sorted.$name -prefix v1.Swembl.Full_length.Top600_peaks.Sorted.$name -title v1.Swembl.Full_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.Full_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL NO -R 0.005 with 300bps interval ALL peaks              ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.300bp_length.All_peaks.Sorted.$name -prefix v1.Swembl.300bp_length.All_peaks.Sorted.$name -title v1.Swembl.300bp_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.300bp_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL NO -R 0.005 with 300bps interval TOP 600 peaks          ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.300bp_length.Top600_peaks.Sorted.$name -prefix v1.Swembl.300bp_length.Top600_peaks.Sorted.$name -title v1.Swembl.300bp_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q/v1.Swembl.300bp_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.005 with FULL interval ALL peaks              ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.Full_length.sorted.fasta -outdir $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.Full_length.All_peaks.Sorted.$name -prefix v2.Swembl.Full_length.All_peaks.Sorted.$name -title v2.Swembl.Full_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.Full_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.005 with FULL interval TOP 600 peaks          ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.Full_length.sorted.fasta -outdir $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.Full_length.Top600_peaks.Sorted.$name -prefix v2.Swembl.Full_length.Top600_peaks.Sorted.$name -title v2.Swembl.Full_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.Full_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.005 with 300bps interval ALL peaks            ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.300bp_length.All_peaks.Sorted.$name -prefix v2.Swembl.300bp_length.All_peaks.Sorted.$name -title v2.Swembl.300bp_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.300bp_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.005 with 300bps interval TOP 600 peaks        ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.300bp_length.Top600_peaks.Sorted.$name -prefix v2.Swembl.300bp_length.Top600_peaks.Sorted.$name -title v2.Swembl.300bp_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q/v2.Swembl.300bp_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.01 with FULL interval ALL peaks              ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.Full_length.sorted.fasta -outdir $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.Full_length.All_peaks.Sorted.$name -prefix v3.Swembl.Full_length.All_peaks.Sorted.$name -title v3.Swembl.Full_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.Full_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.01 with FULL interval TOP 600 peaks          ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.Full_length.sorted.fasta -outdir $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.Full_length.Top600_peaks.Sorted.$name -prefix v3.Swembl.Full_length.Top600_peaks.Sorted.$name -title v3.Swembl.Full_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.Full_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.01 with 300bps interval ALL peaks            ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.300bp_length.All_peaks.Sorted.$name -prefix v3.Swembl.300bp_length.All_peaks.Sorted.$name -title v3.Swembl.300bp_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.300bp_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - SWEMBL WITH -R 0.01 with 300bps interval TOP 600 peaks        ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.300bp_length.Top600_peaks.Sorted.$name -prefix v3.Swembl.300bp_length.Top600_peaks.Sorted.$name -title v3.Swembl.300bp_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q/v3.Swembl.300bp_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - MACS2 1-Summit with FULL interval ALL peaks                    ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.Full_length.sorted.fasta -outdir $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.Full_length.All_peaks.Sorted.$name -prefix v1.Macs2.Full_length.All_peaks.Sorted.$name -title v1.Macs2.Full_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.Full_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - MACS2 1-Summit with FULL interval TOP 600 peaks                ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.Full_length.sorted.fasta -outdir $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.Full_length.Top600_peaks.Sorted.$name -prefix v1.Macs2.Full_length.Top600_peaks.Sorted.$name -title v1.Macs2.Full_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.Full_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - MACS2 1-Summit with 300bps interval ALL peaks                  ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.300bp_length.All_peaks.Sorted.$name -prefix v1.Macs2.300bp_length.All_peaks.Sorted.$name -title v1.Macs2.300bp_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.300bp_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - MACS2 1-Summit with 300bps interval TOP 600 peaks              ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.300bp_length.Top600_peaks.Sorted.$name -prefix v1.Macs2.300bp_length.Top600_peaks.Sorted.$name -title v1.Macs2.300bp_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/MACS2_Trim20L20Q/Summit_1/v1.Macs2.300bp_length.Top600_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - MACS2 n-Summits with 300bps interval ALL peaks                 ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.300bp_length.All_peaks.Sorted.$name -prefix v2.Macs2.300bp_length.All_peaks.Sorted.$name -title v2.Macs2.300bp_length.All_peaks.Sorted.$name -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.300bp_length.All_peaks.Sorted.$name.pm.log &
done
echo "######################################################################################################"
echo "#####              Peak motifs - MACS2 n-Summits with 300bps interval TOP 600 peaks             ######"
echo "######################################################################################################"
for name in $names
do
	peak-motifs -i $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.$name.300bp_length.sorted.fasta -outdir $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.300bp_length.Top600_peaks.Sorted.$name -prefix v2.Macs2.300bp_length.Top600_peaks.Sorted.$name -title v2.Macs2.300bp_length.Top600_peaks.Sorted.$name -top_peaks 600 -max_seq_len 2000 -nmotifs 5 -markov auto -v -motif_db cisBP_Human transfac $ALEA/Data/cisBP_Homo_sapiens_2014-10.tf -motif_db JASPAR_vertebr transfac $ALEA/Data/jaspar_core_vertebrates_2015_03.tf -disco oligos,dyads,positions,local_words &> $ALEA/Data/MACS2_Trim20L20Q/v2.Macs2.300bp_length.Top600_peaks.Sorted.$name.pm.log &
done

echo "##################################################################################################################################################################################################"
echo "################################################################################   Part 6. Variation Scan  #######################################################################################"
echo "##################################################################################################################################################################################################"
mkdir $ALEA/Data/Variation_Sequences

echo "#######################################################################"
echo "#####                Retrieve Variation Sequence                  #####"
echo "#######################################################################"

retrieve-variation-seq -v 2 -species Homo_sapiens -a_version GRCh37 -i $ALEA/Data/CEU_Daughter.hg19.bed -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Daughter.varBed -format bed
retrieve-variation-seq -v 2 -species Homo_sapiens -a_version GRCh37 -i $ALEA/Data/CEU_Father.hg19.bed -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Father.varBed -format bed
retrieve-variation-seq -v 2 -species Homo_sapiens -a_version GRCh37 -i $ALEA/Data/CEU_Mother.hg19.bed -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Mother.varBed -format bed
retrieve-variation-seq -v 2 -species Homo_sapiens -a_version GRCh37 -i $ALEA/Data/YRI_Daughter.hg19.bed -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Daughter.varBed -format bed
retrieve-variation-seq -v 2 -species Homo_sapiens -a_version GRCh37 -i $ALEA/Data/YRI_Father.hg19.bed -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Father.varBed -format bed
retrieve-variation-seq -v 2 -species Homo_sapiens -a_version GRCh37 -i $ALEA/Data/YRI_Mother.hg19.bed -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Mother.varBed -format bed

echo "##########################################################################################"
echo "#####    Background Model construction: DNase ENCODE Clusterded for ALL Cell-line    #####"
echo "##########################################################################################"

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz --output-document=ENCODE_DNaseClustered.v3.bed.gz --directory-prefix=$ALEA/Data
gunzip $ALEA/Data/ENCODE_DNaseClustered.v3.bed.gz
awk '{print $1"\t"$2"\t"$3}' $ALEA/Data/ENCODE_DNaseClustered.v3.bed|sort -k1,1 -k2,2n|bedtools merge > $ALEA/Data/ENCODE_DNase.merged.bed
awk 'x=$2-1{print $1"\t"x"\t"$3}' $ALEA/Data/ENCODE_DNase.merged.bed|fetch-sequences -genome hg19 -o $ALEA/Data/ENCODE_DNase.ALL.fasta -v
create-background-model -i $ALEA/Data/ENCODE_DNase.ALL.fasta -o $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -markov 2 -v 2

echo "#######################################################################"
echo "#####               Variation Scan: JASPAR Matrix                 #####"
echo "#######################################################################"
source /export/apps/rsat/RSAT_config.bashrc
ALEA=/export/storage/users/wsantana/ALEA
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/CEU_Daughter.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -lth score 1 -lth w_diff 1 -lth pval_ratio 10 -uth pval 1e-3 -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Daughter.JASPAR.DNaseALLmkv2.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/CEU_Father.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -lth score 1 -lth w_diff 1 -lth pval_ratio 10 -uth pval 1e-3 -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Father.DNaseALLmkv2JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/CEU_Mother.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -lth score 1 -lth w_diff 1 -lth pval_ratio 10 -uth pval 1e-3 -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Mother.DNaseALLmkv2.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/YRI_Daughter.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -lth score 1 -lth w_diff 1 -lth pval_ratio 10 -uth pval 1e-3 -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Daughter.DNaseALLmkv2.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/YRI_Father.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -lth score 1 -lth w_diff 1 -lth pval_ratio 10 -uth pval 1e-3 -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Father.DNaseALLmkv2.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/YRI_Mother.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -lth score 1 -lth w_diff 1 -lth pval_ratio 10 -uth pval 1e-3 -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Mother.DNaseALLmkv2.JASPAR.varscan
echo "#######################################################################"
echo "#####        Variation Scan NO-threshold: JASPAR Matrix           #####"
echo "#######################################################################"
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/CEU_Daughter.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Daughter.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/CEU_Father.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Father.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/CEU_Mother.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -mml 30 -o $ALEA/Data/Variation_Sequences/CEU_Mother.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/YRI_Daughter.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Daughter.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/YRI_Father.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Father.JASPAR.varscan
variation-scan -v 1 -m $ALEA/Data/JASPAR_CTCF.tf -m_format transfac -i $ALEA/Data/Variation_Sequences/YRI_Mother.varBed -bg $ALEA/Data/Variation_Sequences/DNase_ALL.mkv2.oligos.bg -mml 30 -o $ALEA/Data/Variation_Sequences/YRI_Mother.JASPAR.varscan
echo "##################################################################################################################################################################################################"
echo "#############################################################################    Part 7. PIPELINE EXECUTION   ####################################################################################"
echo "##################################################################################################################################################################################################"
#Requires py modules:
#1.- pysam
#2.- Anaconda
mkdir $ALEA/Results
echo "#######################################################################"
echo "#####                Converting peak files to bed                 #####"
echo "#######################################################################"
names=("`ls -lhtr $ALEA/Data/*.bam|awk '{print $9}'|perl -pe 's/.+Data\/(.+)\.sorted\.bam\n/$1 /'`")
old_suffix=(.Summit_1_peaks.xls _peaks.xls .Peak .Peak .Peak)
new_suffix=(.macs2.v1.bed .macs2.v2.bed .swembl.v1.bed .swembl.v2.bed .swembl.v3.bed)
PeakCalled=($ALEA/Data/MACS2_Trim20L20Q/Summit_1 $ALEA/Data/MACS2_Trim20L20Q $ALEA/Data/v1.PeakCallSwembl_Trim20L20Q $ALEA/Data/v2.PeakCallSwembl_Trim20L20Q $ALEA/Data/v3.PeakCallSwembl_Trim20L20Q)
for i in {0..4}
do
	for name in $names
	do
		grep -v "#\|^$\|start" ${PeakCalled[$i]}/$name${old_suffix[$i]}|awk '{print $1"\t"$2"\t"$3}' > ${PeakCalled[$i]}/$name${new_suffix[$i]}
	done
done
echo "#######################################################################"
echo "#####                 CEU_DAUGHTER: JASPAR MATRIX                 #####"
echo "#######################################################################"

perl $ALEA/variation-scan-git/make_all.pl -tf CTCF -b $ALEA/Data/Trim20L20Q.CEU_Daughter_CTCF.Hap2.sorted.bam -p $ALEA/Data/Variation_Sequences/Trim20L20Q.CEU_Daughter_CTCF.Hap2.Summit_1.peaks.bed -v $ALEA/Data/Variation_Sequences/CEU_Daughter.JASPAR -o $ALEA/Results/CEU_Daughter_JASPAR_test_bgAllHg19 -d 0


















































################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
#									RENAMING FILES FOR QC												       #
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################

ALEA=/export/storage/users/wsantana/ALEA/Data
mkdir $ALEA/CTCF_Trim20L20Q
cp $ALEA/ChIPseq_Files/*.fastq $ALEA/CTCF_Trim20L20Q
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM01_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM02_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM03_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep3.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Father_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM04_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap2_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Father_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM05_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap2_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM06_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM07_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM08_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM09_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Father_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM10_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap2_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Father_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM11_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap2_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM12_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM13_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep2.fastq
	
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM14_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM15_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM16_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep3.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Father_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM17_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Father_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM18_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM19_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM20_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM21_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM22_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Father_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM23_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Father_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM24_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM25_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq $ALEA/CTCF_Trim20L20Q/AM26_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep2.fastq

mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM01
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM02
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM03
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM04
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM05
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM06
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM07
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM08
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM09
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM10
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM11
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM12
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM13

mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM14
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM15
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM16
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM17
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM18
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM19
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM20
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM21
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM22
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM23
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM24
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM25
mkdir $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM26

cp $ALEA/CTCF_Trim20L20Q/AM01_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM01
cp $ALEA/CTCF_Trim20L20Q/AM02_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM02
cp $ALEA/CTCF_Trim20L20Q/AM03_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep3.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM03
cp $ALEA/CTCF_Trim20L20Q/AM04_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM04
cp $ALEA/CTCF_Trim20L20Q/AM05_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM05
cp $ALEA/CTCF_Trim20L20Q/AM06_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM06
cp $ALEA/CTCF_Trim20L20Q/AM07_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM07
cp $ALEA/CTCF_Trim20L20Q/AM08_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM08
cp $ALEA/CTCF_Trim20L20Q/AM09_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM09
cp $ALEA/CTCF_Trim20L20Q/AM10_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM10
cp $ALEA/CTCF_Trim20L20Q/AM11_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM11
cp $ALEA/CTCF_Trim20L20Q/AM12_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM12
cp $ALEA/CTCF_Trim20L20Q/AM13_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM13

mv $ALEA/CTCF_Trim20L20Q/AM01_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/AM14_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/AM02_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/AM15_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/AM03_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep3.fastq $ALEA/CTCF_Trim20L20Q/AM16_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep3.fastq
mv $ALEA/CTCF_Trim20L20Q/AM04_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/AM17_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/AM05_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/AM18_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/AM06_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/AM19_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/AM07_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/AM20_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/AM08_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/AM21_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/AM09_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/AM22_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/AM10_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/AM23_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/AM11_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/AM24_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap1_none_Rep2.fastq
mv $ALEA/CTCF_Trim20L20Q/AM12_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/AM25_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep1.fastq
mv $ALEA/CTCF_Trim20L20Q/AM13_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/AM26_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep2.fastq

mv $ALEA/CTCF_Trim20L20Q/AM14_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM14
mv $ALEA/CTCF_Trim20L20Q/AM15_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM15
mv $ALEA/CTCF_Trim20L20Q/AM16_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep3.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM16
mv $ALEA/CTCF_Trim20L20Q/AM17_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap1_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM17
mv $ALEA/CTCF_Trim20L20Q/AM18_CTCF_Lymphoblastoid_Father_CEU_hsap_Hap1_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM18
mv $ALEA/CTCF_Trim20L20Q/AM19_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM19
mv $ALEA/CTCF_Trim20L20Q/AM20_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM20
mv $ALEA/CTCF_Trim20L20Q/AM21_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM21
mv $ALEA/CTCF_Trim20L20Q/AM22_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM22
mv $ALEA/CTCF_Trim20L20Q/AM23_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap1_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM23
mv $ALEA/CTCF_Trim20L20Q/AM24_CTCF_Lymphoblastoid_Father_YRI_hsap_Hap1_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM24
mv $ALEA/CTCF_Trim20L20Q/AM25_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep1.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM25
mv $ALEA/CTCF_Trim20L20Q/AM26_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep2.fastq $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM26


sh ZipFastq.sh



mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM01/AM01_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM01/AM01_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM02/AM02_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM02/AM02_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep2.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM03/AM03_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep3.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM03/AM03_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep3.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM06/AM06_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM06/AM06_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM07/AM07_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM07/AM07_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep2.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM08/AM08_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM08/AM08_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM09/AM09_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM09/AM09_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep2.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM12/AM12_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM12/AM12_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM13/AM13_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM13/AM13_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep2.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM14/AM14_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM14/AM14_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM15/AM15_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM15/AM15_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep2.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM16/AM16_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep3.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM16/AM16_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep3.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM19/AM19_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM19/AM19_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM20/AM20_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM20/AM20_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep2.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM21/AM21_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM21/AM21_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM22/AM22_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM22/AM22_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep2.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM25/AM25_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep1.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM25/AM25_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep1.bam
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM26/AM26_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep2.bam $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM26/AM26_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep2.bam


mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM01/AM01_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM01/AM01_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM02/AM02_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM02/AM02_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep2.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM03/AM03_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap2_none_Rep3.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM03/AM03_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep3.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM06/AM06_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM06/AM06_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM07/AM07_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap2_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM07/AM07_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep2.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM08/AM08_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM08/AM08_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM09/AM09_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap2_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM09/AM09_CTCF_NA_mp07-729_hsap_Daughter_Hap2_Rep2.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM12/AM12_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM12/AM12_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM13/AM13_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap2_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM13/AM13_CTCF_NA_mp07-729_hsap_Mother_Hap2_Rep2.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM14/AM14_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM14/AM14_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM15/AM15_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM15/AM15_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep2.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM16/AM16_CTCF_Lymphoblastoid_Daughter_CEU_hsap_Hap1_none_Rep3.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM16/AM16_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep3.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM19/AM19_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM19/AM19_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM20/AM20_CTCF_Lymphoblastoid_Mother_CEU_hsap_Hap1_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM20/AM20_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep2.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM21/AM21_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM21/AM21_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM22/AM22_CTCF_Lymphoblastoid_Daughter_YRI_hsap_Hap1_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM22/AM22_CTCF_NA_mp07-729_hsap_Daughter_Hap1_Rep2.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM25/AM25_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep1.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM25/AM25_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep1.fastq.gz
mv $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM26/AM26_CTCF_Lymphoblastoid_Mother_YRI_hsap_Hap1_none_Rep2.fastq.gz $ALEA/CTCF_Trim20L20Q/Sample_Medina_AM26/AM26_CTCF_NA_mp07-729_hsap_Mother_Hap1_Rep2.fastq.gz




date
ALEA=/export/storage/users/wsantana/ALEA
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" > $ALEA/Data/hg19.genome
grep -v "_\|chrM\|chrY\|size" $ALEA/Data/hg19.genome > $ALEA/Data/hg19.female.genomesize
num=(3 2 2 2 3 2 2 2)
Id=(AM01 AM02 AM03 AM06 AM07 AM08 AM09 AM12 AM13 AM14 AM15 AM16 AM19 AM20 AM21 AM22 AM25 AM26)
Genome=(CEU_Daughter.Chr.Hap2.fasta CEU_Mother.Chr.Hap2.fasta YRI_Daughter.Chr.Hap2.fasta YRI_Mother.Chr.Hap2.fasta CEU_Daughter.Chr.Hap1.fasta CEU_Mother.Chr.Hap1.fasta YRI_Daughter.Chr.Hap1.fasta YRI_Mother.Chr.Hap1.fasta)
k=0
for i in {0..7}
do
	echo "Este es el numero ${num[$i]}"
	echo "Este es el genoma ${Genome[$i]}"
	for (( j=0;j<${num[$i]};j++ ))
	do
		echo ">Este es la muestra ${Id[$k]}"
		name=`ls -l $ALEA/Data/CTCF_Trim20L20Q/Sample_Medina_${Id[$k]}/${Id[$k]}*.bam|awk '{print $9}'|perl -pe 's/.+'${Id[$k]}'\/'${Id[$k]}'(_.+).bam\n/$1\n/'`
		echo ">Ruta $ALEA/Data/CTCF_Trim20L20Q/Sample_Medina_${Id[$k]}/${Id[$k]}$name"
		chrom_info=$ALEA/Data/${Genome[$i]}
		
		bam="$ALEA/Data/CTCF_Trim20L20Q/Sample_Medina_${Id[$k]}/${Id[$k]}$name.bam"
		bed="$ALEA/Data/CTCF_Trim20L20Q/Sample_Medina_${Id[$k]}/${Id[$k]}$name.bed"
		uniqbed="$ALEA/Data/CTCF_Trim20L20Q/Sample_Medina_${Id[$k]}/${Id[$k]}$name.uniq.bed"
		bedgraph="$ALEA/Data/CTCF_Trim20L20Q/Sample_Medina_${Id[$k]}/${Id[$k]}$name.bedGraph"
		bigwig="$ALEA/Data/CTCF_Trim20L20Q/Sample_Medina_${Id[$k]}/${Id[$k]}$name.bigWig"
		echo ">Indexing BAM file: $bam"
		samtools index $bam
		
		echo ">Generating BED from BAM: $bam"
		bamToBed -i $bam > $bed
	
		echo ">Sorting BED file: $bed"
		sort -k1,1 -k2,2n -k3,3n -k6,6 $bed > $bed".sorted"
		mv $bed".sorted" $bed

		echo ">Generating unique mapped BED file: $uniqbed"
		grep -vP '\t0\t' $bed > $uniqbed

		echo ">Generating bedGraph from BED: $bed"
		genomeCoverageBed -bg -i $bed -g $chrom_info > $bedgraph

		echo ">Generating bigWig from bedGraph: $bedGraph"
		bedGraphToBigWig $bedgraph $chrom_info $bigwig

		k=$((k+1))

	done
done
date
mkdir /export/storage/users/wsantana/ALEA/Data/CTCF_Organized/uniqbed
mkdir /export/storage/users/wsantana/ALEA/Data/CTCF_Organized/fastq
mkdir /export/storage/users/wsantana/ALEA/Data/CTCF_Organized/bam
mkdir /export/storage/users/wsantana/ALEA/Data/CTCF_Organized/bed
mkdir /export/storage/users/wsantana/ALEA/Data/CTCF_Organized/bedGraph
mkdir /export/storage/users/wsantana/ALEA/Data/CTCF_Organized/bigWig

ALEA=export/storage/users/wsantana/ALEA
sample=("01 02 03 06 07 08 09 12 13 14 15 16 19 20 21 22 25 26")
for i in $sample;do echo "AM$i"; done| xargs $ALEA/qc.1/qc/QC.py -s > qc.csv









