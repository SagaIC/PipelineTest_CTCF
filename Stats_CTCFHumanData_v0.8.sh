echo "###########################################################################"
echo "#########            Allele-specific Variant Analysis             #########"
echo "#########                  CTCF Human Data Stats                  #########"
echo "#########        Script name: Stats_CTCFHumanData_v0.8.sh         #########"
echo "###########################################################################"
date

# $PATHs#
ALEA=/export/storage/users/wsantana/ALEA
DATA=$ALEA/Data
STATS=$ALEA/Stats
SCRIPTS=$ALEA/Scripts
mkdir $STATS


#(1) Number of SNP entries in each family member of both trios#
echo "Sample Group Counts" >> $STATS/CountsGenotypes.hg18.txt
grep -v "#\|^$" $DATA/CEU_Father.genotypes.hg18.vcf|wc -l|xargs echo "Father CEU" >> $STATS/CountsGenotypes.hg18.txt
grep -v "#\|^$" $DATA/CEU_Mother.genotypes.hg18.vcf|wc -l|xargs echo "Mother CEU" >> $STATS/CountsGenotypes.hg18.txt
grep -v "#\|^$" $DATA/CEU_Daughter.genotypes.hg18.vcf|wc -l|xargs echo "Daughter CEU" >> $STATS/CountsGenotypes.hg18.txt
grep -v "#\|^$" $DATA/YRI_Father.genotypes.hg18.vcf|wc -l|xargs echo "Father YRI" >> $STATS/CountsGenotypes.hg18.txt
grep -v "#\|^$" $DATA/YRI_Mother.genotypes.hg18.vcf|wc -l|xargs echo "Mother YRI" >> $STATS/CountsGenotypes.hg18.txt
grep -v "#\|^$" $DATA/YRI_Daughter.genotypes.hg18.vcf|wc -l|xargs echo "Daughter YRI" >> $STATS/CountsGenotypes.hg18.txt
Rscript TablePlot_v1.0.R -input=$STATS/CountsGenotypes.hg18.txt,$STATS/CountsGenotypes.hg19.txt,$STATS/CountsGenotypes.hg19_phased.txt,$STATS/CountsGenotypes.hg19_filtered.txt -output=$STATS/TableSNPCounts.png -header=Group,Sample -samp_nam=Hg18,Hg19,Hg19_Phased,Hg19_Filtered -title=SNPs

#(2) Number of Raw Reads for each REPLICA
echo "Sample Group Replica Counts" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Father_CTCF_Rep1_36683.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Father_CTCF_Rep1_36791.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Father_CTCF_Rep2_36684.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Father_CTCF_Rep2_36793.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Mother_CTCF_Rep1_36685.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother CEU Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Mother_CTCF_Rep2_36686.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother CEU Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Daughter_CTCF_Rep1_36680.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Daughter_CTCF_Rep2_36681.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Daughter_CTCF_Rep2_36792.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/CEU_Daughter_CTCF_Rep3_36682.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep3" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Father_CTCF_Rep1_36691.fastq|awk 'x=$1/4{print x}'|xargs echo "Father YRI Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Father_CTCF_Rep2_36692.fastq|awk 'x=$1/4{print x}'|xargs echo "Father YRI Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Mother_CTCF_Rep1_36689.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Mother_CTCF_Rep1_36801.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Mother_CTCF_Rep2_36690.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Mother_CTCF_Rep2_36803.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep2" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Daughter_CTCF_Rep1_36687.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Daughter_CTCF_Rep1_36795.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep1" >> $STATS/CountsCTCF_RawReads.txt
wc -l $DATA/ChIPseq_Files/YRI_Daughter_CTCF_Rep2_36688.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep2" >> $STATS/CountsCTCF_RawReads.txt

#(3) Number of Raw Reads for each SAMPLE
Rscript BarPlot.R -input=$STATS/CountsCTCF_RawReads.txt -output=$STATS/RawReadsNumberSample_CTCF.jpg -title=Raw reads count

#(4) Number of Raw Reads for each REPLICA after Trimming 20L20Q
echo "Sample Group Replica Counts" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep1.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep1" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.CEU_Father_CTCF_Rep2.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep2" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep1.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother CEU Rep1" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.CEU_Mother_CTCF_Rep2.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother CEU Rep2" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep1.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep1" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep2.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep2" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.CEU_Daughter_CTCF_Rep3.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep3" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.YRI_Father_CTCF_Rep1.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Father YRI Rep1" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.YRI_Father_CTCF_Rep2.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Father YRI Rep2" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep1.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep1" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.YRI_Mother_CTCF_Rep2.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep2" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep1.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep1" >> $STATS/CountsCTCF_Trim20L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim20L20Q.YRI_Daughter_CTCF_Rep2.ALL.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep2" >> $STATS/CountsCTCF_Trim20L20Q.txt
Rscript TablePlot_v1.0.R -input=$STATS/CountsCTCF_RawReads.txt,$STATS/CountsCTCF_Trim20L20Q.txt,$STATS/CountsCTCF_Trim30L20Q.txt -output=$STATS/TableReadsCounts.png -header=Group,Sample,Replica -samp_nam=Raw,Trim20L20Q,Trim30L20Q -title=Reads

#(5) Number of Raw Reads for each SAMPLE after Trimming 20L20Q
Rscript BarPlot.R -input=$STATS/CountsCTCF_Trim20L20Q.txt -output=$STATS/Trim20L20QNumberSample_CTCF.jpg -title=Trim20L20Q reads count

#(6) Number of Raw Reads for each REPLICA after Trimming 30L20Q
echo "Sample Group Replica Counts" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Father_CTCF_Rep1_36683.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Father_CTCF_Rep1_36791.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Father_CTCF_Rep2_36684.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Father_CTCF_Rep2_36793.fastq|awk 'x=$1/4{print x}'|xargs echo "Father CEU Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Mother_CTCF_Rep1_36685.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother CEU Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Mother_CTCF_Rep2_36686.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother CEU Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Daughter_CTCF_Rep1_36680.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Daughter_CTCF_Rep2_36681.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Daughter_CTCF_Rep2_36792.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.CEU_Daughter_CTCF_Rep3_36682.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter CEU Rep3" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Father_CTCF_Rep1_36691.fastq|awk 'x=$1/4{print x}'|xargs echo "Father YRI Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Father_CTCF_Rep2_36692.fastq|awk 'x=$1/4{print x}'|xargs echo "Father YRI Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Mother_CTCF_Rep1_36689.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Mother_CTCF_Rep1_36801.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Mother_CTCF_Rep2_36690.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Mother_CTCF_Rep2_36803.fastq|awk 'x=$1/4{print x}'|xargs echo "Mother YRI Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Daughter_CTCF_Rep1_36687.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Daughter_CTCF_Rep1_36795.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep1" >> $STATS/CountsCTCF_Trim30L20Q.txt
wc -l $DATA/ChIPseq_Files/Trim30L20Q.YRI_Daughter_CTCF_Rep2_36688.fastq|awk 'x=$1/4{print x}'|xargs echo "Daughter YRI Rep2" >> $STATS/CountsCTCF_Trim30L20Q.txt
Rscript BarPlot.R -input=$STATS/CountsCTCF_Trim30L20Q.txt -output=$STATS/Trim30L20QNumberReplica_CTCF.jpg -title=Trim30L20Q reads count


#(7) Number of Raw Reads for each SAMPLE after Trimming 30L20Q
Rscript BarPlot.R -input=$STATS/CountsCTCF_Trim30L20Q.txt -output=$STATS/Trim30L20QNumberSample_CTCF.jpg -title=Trim30L20Q reads count

#(8) Number of Unmapped SNPs after CrossMap
echo "Sample Group Counts" >> $STATS/CountsGenotypes.hg19.txt
grep -v "#\|^$" $DATA/CEU_Father.genotypes.hg19.vcf|wc -l|xargs echo "Father CEU" >> $STATS/CountsGenotypes.hg19.txt
grep -v "#\|^$" $DATA/CEU_Mother.genotypes.hg19.vcf|wc -l|xargs echo "Mother CEU" >> $STATS/CountsGenotypes.hg19.txt
grep -v "#\|^$" $DATA/CEU_Daughter.genotypes.hg19.vcf|wc -l|xargs echo "Daughter CEU" >> $STATS/CountsGenotypes.hg19.txt
grep -v "#\|^$" $DATA/YRI_Father.genotypes.hg19.vcf|wc -l|xargs echo "Father YRI" >> $STATS/CountsGenotypes.hg19.txt
grep -v "#\|^$" $DATA/YRI_Mother.genotypes.hg19.vcf|wc -l|xargs echo "Mother YRI" >> $STATS/CountsGenotypes.hg19.txt
grep -v "#\|^$" $DATA/YRI_Daughter.genotypes.hg19.vcf|wc -l|xargs echo "Daughter YRI" >> $STATS/CountsGenotypes.hg19.txt
Rscript BarPlot.R -input=$STATS/CountsGenotypes.hg19.txt -output=$STATS/RemappedSNPsNumber.hg19.jpg -title=SNPs Hg19


#(9) Number of PHASED SNPs after filtering
echo "Sample Group Counts" >> $STATS/CountsGenotypes.hg19_phased.txt
grep -v "#\|^$" $DATA/CEU_Father.genotypes.hg19.phased.vcf|wc -l|xargs echo "Father CEU" >> $STATS/CountsGenotypes.hg19_phased.txt
grep -v "#\|^$" $DATA/CEU_Mother.genotypes.hg19.phased.vcf|wc -l|xargs echo "Mother CEU" >> $STATS/CountsGenotypes.hg19_phased.txt
grep -v "#\|^$" $DATA/CEU_Daughter.genotypes.hg19.phased.vcf|wc -l|xargs echo "Daughter CEU" >> $STATS/CountsGenotypes.hg19_phased.txt
grep -v "#\|^$" $DATA/YRI_Father.genotypes.hg19.phased.vcf|wc -l|xargs echo "Father YRI" >> $STATS/CountsGenotypes.hg19_phased.txt
grep -v "#\|^$" $DATA/YRI_Mother.genotypes.hg19.phased.vcf|wc -l|xargs echo "Mother YRI" >> $STATS/CountsGenotypes.hg19_phased.txt
grep -v "#\|^$" $DATA/YRI_Daughter.genotypes.hg19.phased.vcf|wc -l|xargs echo "Daughter YRI" >> $STATS/CountsGenotypes.hg19_phased.txt
Rscript BarPlot.R -input=$STATS/CountsGenotypes.hg19_phased.txt -output=$STATS/PhasedSNPsNumber.hg19.jpg -title=SNPs Hg19


#(10) Number of SNPs after Filtering out other variants
echo "Sample Group Counts" >> $STATS/CountsGenotypes.hg19_filtered.txt
grep -v "#\|^$" $DATA/CEU_Father.genotypes.hg19.phased2.vcf|wc -l|xargs echo "Father CEU" >> $STATS/CountsGenotypes.hg19_filtered.txt
grep -v "#\|^$" $DATA/CEU_Mother.genotypes.hg19.phased2.vcf|wc -l|xargs echo "Mother CEU" >> $STATS/CountsGenotypes.hg19_filtered.txt
grep -v "#\|^$" $DATA/CEU_Daughter.genotypes.hg19.phased2.vcf|wc -l|xargs echo "Daughter CEU" >> $STATS/CountsGenotypes.hg19_filtered.txt
grep -v "#\|^$" $DATA/YRI_Father.genotypes.hg19.phased2.vcf|wc -l|xargs echo "Father YRI" >> $STATS/CountsGenotypes.hg19_filtered.txt
grep -v "#\|^$" $DATA/YRI_Mother.genotypes.hg19.phased2.vcf|wc -l|xargs echo "Mother YRI" >> $STATS/CountsGenotypes.hg19_filtered.txt
grep -v "#\|^$" $DATA/YRI_Daughter.genotypes.hg19.phased2.vcf|wc -l|xargs echo "Daughter YRI" >> $STATS/CountsGenotypes.hg19_filtered.txt
Rscript BarPlot.R -input=$STATS/CountsGenotypes.hg19_filtered.txt  -output=$STATS/FilteredSNPsNumber.hg19.jpg -title=SNPs Hg19

#(11) Number of called Peaks per SAMPLE, REPLICA and METHOD
Rscript TablePlot_v1.0.R -input=$STATS/CountsCalledPeaks_v1Swembl_Trim20L20Q.txt,$STATS/CountsCalledPeaks_v2Swembl_Trim20L20Q.txt,$STATS/CountsCalledPeaks_v1Macs2_Trim20L20Q.txt,$STATS/CountsCalledPeaks_v2Macs2_Trim20L20Q.txt -output=$STATS/TableCalledPeaksCounts.png -header=Group,Sample,Replica,Haplotype -samp_nam=Swembl.v1,Swembl.v2,Macs2.v1,Macs2.v2 -title=CalledPeaks
#(11.01)v1.Swembl Full Length, All Peaks
names=("`ls -lhtr $ALEA/Data/*.bam|awk '{print $9}'|perl -pe 's/.+Data\/(.+)\.sorted\.bam\n/$1 /'`")
echo "Sample Group Replica Haplotype Counts" >> $STATS/CountsCalledPeaks_v1Swembl_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|Region" $DATA/v1.PeakCallSwembl_Trim20L20Q/$name.Peak|wc -l|xargs echo $data >> $STATS/CountsCalledPeaks_v1Swembl_Trim20L20Q.txt
done
#(11.02)v2.Swembl Full Length, All Peaks
echo "Sample Group Replica Haplotype Counts" >> $STATS/CountsCalledPeaks_v2Swembl_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|Region" $DATA/v2.PeakCallSwembl_Trim20L20Q/$name.Peak|wc -l|xargs echo $data >> $STATS/CountsCalledPeaks_v2Swembl_Trim20L20Q.txt
done
#(11.03)v1.MACS2 Full Length, All Peaks
echo "Sample Group Replica Haplotype Counts" >> $STATS/CountsCalledPeaks_v1Macs2_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|start" $DATA/MACS2_Trim20L20Q/Summit_1/$name.Summit_1_peaks.xls|wc -l|xargs echo $data >> $STATS/CountsCalledPeaks_v1Macs2_Trim20L20Q.txt
done
#(11.04)v2.MACS2 Full Length, All Peaks
echo "Sample Group Replica Haplotype Counts" >> $STATS/CountsCalledPeaks_v1Macs2_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|start" $DATA/MACS2_Trim20L20Q/$name\_peaks.xls|wc -l|xargs echo $data >> $STATS/CountsCalledPeaks_v2Macs2_Trim20L20Q.txt
done


#(12) Distribution of Peak sizes per SAMPLE, REPLICA and METHOD (violin-plot or distribution-plot)
#(12.01)v1.Swembl Full Length, All Peaks
names=("`ls -lhtr $ALEA/Data/*.bam|awk '{print $9}'|perl -pe 's/.+Data\/(.+)\.sorted\.bam\n/$1 /'`")
echo "Sample Group Replica Haplotype Length Method" >> $STATS/LengthScoreSortedCalledPeaks_ALL_Trim20L20Q.txt
echo "Sample Group Replica Haplotype Length Method" >> $STATS/LengthScoreSortedCalledPeaks_v1SwemblAllPeaks_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|Region" $DATA/v1.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv1_AllPeaks"}' >> $STATS/LengthScoreSortedCalledPeaks_v1SwemblAllPeaks_Trim20L20Q.txt
	grep -v "#\|^$\|Region" $DATA/v1.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv1_AllPeaks"}' >> $STATS/LengthScoreSortedCalledPeaks_ALL_Trim20L20Q.txt
done
#(12.02)v1.Swembl Full Length, Top 600 Peaks
echo "Sample Group Replica Haplotype Length Method" >> $STATS/LengthScoreSortedCalledPeaks_v1SwemblTop600_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|Region" $DATA/v1.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv1_Top600"}'|head -600 >> $STATS/LengthScoreSortedCalledPeaks_v1SwemblTop600_Trim20L20Q.txt
	grep -v "#\|^$\|Region" $DATA/v1.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv1_Top600"}'|head -600 >> $STATS/LengthScoreSortedCalledPeaks_ALL_Trim20L20Q.txt
done
#(12.03)v2.Swembl Full Length, All Peaks
echo "Sample Group Replica Haplotype Length Method" >> $STATS/LengthScoreSortedCalledPeaks_v2SwemblAllPeaks_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|Region" $DATA/v2.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv2_AllPeaks"}' >> $STATS/LengthScoreSortedCalledPeaks_v2SwemblAllPeaks_Trim20L20Q.txt
	grep -v "#\|^$\|Region" $DATA/v2.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv2_AllPeaks"}' >> $STATS/LengthScoreSortedCalledPeaks_ALL_Trim20L20Q.txt
done
#(12.04)v2.Swembl Full Length, Top 600 Peaks
echo "Sample Group Replica Haplotype Length Method" >> $STATS/LengthScoreSortedCalledPeaks_v2SwemblTop600_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|Region" $DATA/v2.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv2_Top600"}'|head -600 >> $STATS/LengthScoreSortedCalledPeaks_v2SwemblTop600_Trim20L20Q.txt
	grep -v "#\|^$\|Region" $DATA/v2.PeakCallSwembl_Trim20L20Q/$name.Peak|sort -nr -k7|awk -v var="$data" '{print var,$5,"Swemblv2_Top600"}'|head -600 >> $STATS/LengthScoreSortedCalledPeaks_ALL_Trim20L20Q.txt
done
#(12.05)v1.MACS2 Full Length, All Peaks
echo "Sample Group Replica Haplotype Length Method" >> $STATS/LengthScoreSortedCalledPeaks_v1Macs2AllPeaks_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|start" $DATA/MACS2_Trim20L20Q/Summit_1/$name.Summit_1_peaks.xls|sort -nr -k8|awk -v var="$data" '{print var,$4,"Macs2v1_AllPeaks"}' >> $STATS/LengthScoreSortedCalledPeaks_v1Macs2AllPeaks_Trim20L20Q.txt
	grep -v "#\|^$\|start" $DATA/MACS2_Trim20L20Q/Summit_1/$name.Summit_1_peaks.xls|sort -nr -k8|awk -v var="$data" '{print var,$4,"Macs2v1_AllPeaks"}' >> $STATS/LengthScoreSortedCalledPeaks_ALL_Trim20L20Q.txt
done
#(12.06)v1.MACS2 Full Length, Top 600 Peaks
echo "Sample Group Replica Haplotype Length Method" >> $STATS/LengthScoreSortedCalledPeaks_v1Macs2Top600_Trim20L20Q.txt
for name in $names
do
	data=("`echo $name|perl -pe 'if (/Rep/) {s/.+\.(\w+)_(\w+)_\w+_(Rep\d)\.(.+)/$2 $1 $3 $4/} else {s/.+\.(\w+)_(\w+)_\w+\.(.+)/$2 $1 NA $3/}'`")
	grep -v "#\|^$\|start" $DATA/MACS2_Trim20L20Q/Summit_1/$name.Summit_1_peaks.xls|sort -nr -k8|awk -v var="$data" '{print var,$4,"Macs2v1_Top600"}'|head -600 >> $STATS/LengthScoreSortedCalledPeaks_v1Macs2Top600_Trim20L20Q.txt
	grep -v "#\|^$\|start" $DATA/MACS2_Trim20L20Q/Summit_1/$name.Summit_1_peaks.xls|sort -nr -k8|awk -v var="$data" '{print var,$4,"Macs2v1_Top600"}'|head -600 >> $STATS/LengthScoreSortedCalledPeaks_ALL_Trim20L20Q.txt
done





echo "###########################################################################"
echo "#########                       TABLE PLOTS                       #########"
echo "###########################################################################"
Rscript $SCRIPTS/TablePlot_v1.0.R -input=$STATS/CountsGenotypes.hg18.txt,$STATS/CountsGenotypes.hg19.txt,$STATS/CountsGenotypes.hg19_phased.txt,$STATS/CountsGenotypes.hg19_filtered.txt -output=$STATS/TableSNPCounts.png -header=Group,Sample -samp_nam=Hg18,Hg19,Hg19_Phased,Hg19_Filtered -title=SNPs
Rscript $SCRIPTS/TablePlot_v1.0.R -input=$STATS/CountsCTCF_RawReads.txt,$STATS/CountsCTCF_Trim20L20Q.txt,$STATS/CountsCTCF_Trim30L20Q.txt -output=$STATS/TableReadsCounts.png -header=Group,Sample,Replica -samp_nam=Raw,Trim20L20Q,Trim30L20Q -title=Reads
Rscript $SCRIPTS/TablePlot_v1.0.R -input=$STATS/CountsCalledPeaks_v1Swembl_Trim20L20Q.txt,$STATS/CountsCalledPeaks_v2Swembl_Trim20L20Q.txt,$STATS/CountsCalledPeaks_v1Macs2_Trim20L20Q.txt,$STATS/CountsCalledPeaks_v2Macs2_Trim20L20Q.txt -output=$STATS/TableCalledPeaksCounts.png -header=Group,Sample,Replica,Haplotype -samp_nam=Swembl.v1,Swembl.v2,Macs2.v1,Macs2.v2 -title=CalledPeaks
echo "###########################################################################"
echo "#########                        BAR PLOTS                        #########"
echo "###########################################################################"




#DIDNT PLOT BECAUSE FASTA WERE 600Length (12.07)v2.MACS2 Full Length, All Peaks
#DIDNT PLOT BECAUSE FASTA WERE 600Length (12.08)v2.MACS2 Full Length, Top 600 Peaks












#####################################################################################################
#(11.01)v1.Swembl Full Length, All Peaks
#(11.02)v1.Swembl Full Length, Top 600 Peaks
#(11.03)v1.Swembl 300bp Length, All Peaks
#(11.04)v1.Swembl 300bp Length, Top 600 Peaks
#(11.05)v2.Swembl Full Length, All Peaks
#(11.06)v2.Swembl Full Length, Top 600 Peaks
#(11.07)v2.Swembl 300bp Length, All Peaks
#(11.08)v2.Swembl 300bp Length, Top 600 Peaks
#(11.09)v1.MACS2 Full Length, All Peaks
#(11.10)v1.MACS2 Full Length, Top 600 Peaks
#(11.11)v1.MACS2 300bp Length, All Peaks
#(11.12)v1.MACS2 300bp Length, Top 600 Peaks
#(11.13)v2.MACS2 300bp Length, All Peaks
#(11.14)v2.MACS2 300bp Length, Top 600 Peaks


