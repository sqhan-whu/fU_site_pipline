# name	GSE_number	samples
# GDIN1	GSM4998039	Cal62_input1
# GDIN2	GSM4998040	Cal62_input2
# GD172	GSM4998041	Cal62_output_nonmismatch1
# GD188	GSM4998042	Cal62_output_nonmismatch2
# GD272	GSM4998043	Cal62_output_mismatch1
# GD288	GSM4998044	Cal62_output_mismatch2

# 1. macs2 call peaks：
#!/bin/bash
#SBATCH -N 1 -c 16
macs2 callpeak --nomodel --slocal 200 -t GD172.sorted.bam.filter.bam -c GDIN1.filter.bam -g hs -B -f BAM -n 172_IN1 --outdir ./MACS2 
macs2 callpeak --nomodel --slocal 200 -t GD188.sorted.bam.filter.bam -c GDIN1.filter.bam -g hs -B -f BAM -n 188_IN1 --outdir ./MACS2 
macs2 callpeak --nomodel --slocal 200 -t GD272.sorted.bam.filter.bam -c GDIN2.filter.bam -g hs -B -f BAM -n 272_IN2 --outdir ./MACS2 
macs2 callpeak --nomodel --slocal 200 -t GD288.sorted.bam.filter.bam -c GDIN2.filter.bam -g hs -B -f BAM -n 288_IN2 --outdir ./MACS2 
macs2 callpeak --nomodel --slocal 200 -t GD172.sorted.bam.filter.bam -c GDIN2.filter.bam -g hs -B -f BAM -n 172_IN2 --outdir ./MACS2
macs2 callpeak --nomodel --slocal 200 -t GD188.sorted.bam.filter.bam -c GDIN2.filter.bam -g hs -B -f BAM -n 188_IN2 --outdir ./MACS2 
macs2 callpeak --nomodel --slocal 200 -t GD272.sorted.bam.filter.bam -c GDIN1.filter.bam -g hs -B -f BAM -n 272_IN1 --outdir ./MACS2 
macs2 callpeak --nomodel --slocal 200 -t GD288.sorted.bam.filter.bam -c GDIN1.filter.bam -g hs -B -f BAM -n 288_IN1 --outdir ./MACS2

# 2. merge all peaks：
cat *.xls |grep -E -v '#|start' |cut -f 1-3 |sort -k1,1 -k2,2n | awk '{print $1,$2-100,$3+100}'|tr ' ' '\t'|
grep -v '\-'|bedtools merge -i - |grep -v M > total_extend_100bp.bed

# 3. get reads in peak regions:

#for i in `ls *.bam`;do echo "bedtools intersect -a $i -b total_extend_100bp.bed > $i.region.bam";done
bedtools intersect -a GD172.sorted.bam.filter.bam -b total_extend_100bp.bed > GD172.sorted.bam.filter.bam.region.bam
bedtools intersect -a GD188.sorted.bam.filter.bam -b total_extend_100bp.bed > GD188.sorted.bam.filter.bam.region.bam
bedtools intersect -a GD272.sorted.bam.filter.bam -b total_extend_100bp.bed > GD272.sorted.bam.filter.bam.region.bam
bedtools intersect -a GD288.sorted.bam.filter.bam -b total_extend_100bp.bed > GD288.sorted.bam.filter.bam.region.bam

# 4. sam to tsv by sam2tsv
sam2tsv --reference /home/yangwei/project/00.DATABASE/hg38/bowtie2_index/hg38.fa GD172.sorted.bam.filter.bam.region.bam | awk '{if ($9 ~/M|=|X/ ) print $0}'  > GD172_both_strand_all.tsv
sam2tsv --reference /home/yangwei/project/00.DATABASE/hg38/bowtie2_index/hg38.fa GD272.sorted.bam.filter.bam.region.bam | awk '{if ($9 ~/M|=|X/ ) print $0}'  > GD272_both_strand_all.tsv
sam2tsv --reference /home/yangwei/project/00.DATABASE/hg38/bowtie2_index/hg38.fa GD188.sorted.bam.filter.bam.region.bam | awk '{if ($9 ~/M|=|X/ ) print $0}'  > GD188_both_strand_all.tsv
sam2tsv --reference /home/yangwei/project/00.DATABASE/hg38/bowtie2_index/hg38.fa GD288.sorted.bam.filter.bam.region.bam | awk '{if ($9 ~/M|=|X/ ) print $0}'  > GD288_both_strand_all.tsv

#5. Extract T->C mutation from _all.tsv considering both strands
awk '{if ($5=="C" && $8=="T") print $0; else if ( $5=="G" && $8=="A" ) print $0}'  GD172_both_strand_all.tsv > GD172_strand_all_TC.tsv
awk '{if ($5=="C" && $8=="T") print $0; else if ( $5=="G" && $8=="A" ) print $0}'  GD272_both_strand_all.tsv > GD272_strand_all_TC.tsv
awk '{if ($5=="C" && $8=="T") print $0; else if ( $5=="G" && $8=="A" ) print $0}'  GD188_both_strand_all.tsv > GD188_strand_all_TC.tsv
awk '{if ($5=="C" && $8=="T") print $0; else if ( $5=="G" && $8=="A" ) print $0}'  GD288_both_strand_all.tsv > GD288_strand_all_TC.tsv

#6. Retain T-to-C substitutions with a base Phred quality score of > 27
perl ~/script/extrac_refT_readC.pl -tsv GD172_strand_all_TC.tsv -qual 27
perl ~/script/extrac_refT_readC.pl -tsv GD188_strand_all_TC.tsv -qual 27
perl ~/script/extrac_refT_readC.pl -tsv GD272_strand_all_TC.tsv -qual 27
perl ~/script/extrac_refT_readC.pl -tsv GD288_strand_all_TC.tsv -qual 27

#7. Substract mutations from control samples

perl ~/script/background_correction.pl -bg GD172_strand_all_TC.tsv_q27.tsv -in GD188_strand_all_TC.tsv_q27.tsv

perl ~/script/background_correction.pl -bg GD272_strand_all_TC.tsv_q27.tsv -in GD288_strand_all_TC.tsv_q27.tsv



#8. get mutation number for clean expriment sample
python get_mutaion.py GD188_strand_all_TC.tsv_q27.tsv_corrected.tsv GD288_both_strand_all.tsv > GD188.mutation.raw.txt
python get_mutaion.py GD288_strand_all_TC.tsv_q27.tsv_corrected.tsv GD288_both_strand_all.tsv > GD288.mutation.raw.txt

#9. filter T_C mutaion and A_G mution
awk '{if(($7 > $4*0.15||$7 > 5) || ($8 > $4*0.15||$8 > 5) && ($7 < 0.5 * $4 && $8 < 0.5* $4))print $0}' GD188.mutation.raw.txt > tm188
awk '{if(($7 > $4*0.15||$7 > 5) || ($8 > $4*0.15||$8 > 5) && ($7 < 0.5 * $4 && $8 < 0.5* $4))print $0}' GD288.mutation.raw.txt > tm288
cat tm188 tm288 |cut -f 1-2|sort |uniq -c |awk '{if($1==2)print $2"\t"$3}' > overlap
python ~/bin_python/vlookup.py tm188 overlap > tmp188
python ~/bin_python/vlookup.py tm288 overlap > tmp288
paste tmp188 tmp288 |awk '{if($4>=10 && $12 >=10)print $0}' |sort -k1,1 -k2,2n > final.point.tsv
