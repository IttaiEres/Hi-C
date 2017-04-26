#Automated script for running a wide range of Hi-C analyses.

my=$(date +%b%g)
mkdir results_$my
cd results_$my
for line in 1 2 3 4
do
mkdir {$line,$line/hicup,$line/homer,$line/juicer,$line/bwa}
done

for dir in */; do mkdir -- "$dir"/{tmp1,foo,bar,qux}; done

#ensure it's the same as # lines in file (wc -l file)
cut -f3,8 just22.tags.tsv > 22positions

for .homer file:
grep ‘chr3.*chr3’ $i.homer > chr3pos #extract all of a single chromosome’s contacts. Has an issue! regex of chr1 and chr2 will pick up on all teen/20s chrs as well with this writing...
cut -f3,6 chr3pos > chr3positions

grep "chr22$(printf '\t').*chr22$(printf '\t')" ./21792.bed | cut -f1,2,3, > 22only_test
grep "chr22$(printf '\t').*chr22$(printf '\t')" ./21792.bed | cut -f4,5,6 >> 22only_test

#!/bin/bash

grep 'chr1\t.*chr1\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just1
grep 'chr2\t.*chr2\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just2
grep 'chr3\t.*chr3\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just3
grep 'chr4\t.*chr4\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just4
grep 'chr5\t.*chr5\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just5
grep 'chr6\t.*chr6\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just6
grep 'chr7\t.*chr7\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just7
grep 'chr8\t.*chr8\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just8
grep 'chr9\t.*chr9\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just9
grep 'chr10\t.*chr10\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just10
grep 'chr11\t.*chr11\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just11
grep 'chr12\t.*chr12\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just12
grep 'chr13\t.*chr13\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just13
grep 'chr14\t.*chr14\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just14
grep 'chr15\t.*chr15\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just15
grep 'chr16\t.*chr16\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just16
grep 'chr17\t.*chr17\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just17
grep 'chr18\t.*chr18\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just18
grep 'chr19\t.*chr19\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just19
grep 'chr20\t.*chr20\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just20
grep 'chr21\t.*chr21\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just21
grep 'chr22\t.*chr22\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/just22
grep 'chrX\t.*chrX\t' 21792.genomewide.pos | cut -f2,4 > ./chromo.positions/justX

#!/bin/bash




#To extract midpoint read positions from a .bed file:
cut -f1,2,3,4,5,6 {input}.bed | awk 'NR>0 {NF=7;$7=($2+$3)/2}1' | awk 'NR>0 {NF=8;$8=($5+$6)/2}1' | awk '{print $1 " " $7 " " $4 " " $8}' > {output}.pos

#With a full homer overlap DF, to make a bedtools-compatible file once again. First part cuts out loc1_bin column from homer overlap DF, second part grabs first part of the bin range, and 3rd part subsets down to just the digits (start of the range):
cut -f3 test.txt | cut -d ',' -f1 | grep -o '[0-9]\+'



[ittai@midway-login2]/project/gilad/ittai/roach/testing% cut -f1 test.txt >test1
[ittai@midway-login2]/project/gilad/ittai/roach/testing% cut -f3 test.txt | cut -d ',' -f1 | grep -o '[0=9]\+' >test2
[ittai@midway-login2]/project/gilad/ittai/roach/testing% paste test1 test2 | head
chr	000
chr1	0
chr1	000
chr1	000
chr1	000
chr1	0
chr1	000
chr1	000
chr1	000
chr1	000
[ittai@midway-login2]/project/gilad/ittai/roach/testing% cut -f3 test.txt | cut -d ',' -f1 | grep -o '[0-9]\+' >test2
[ittai@midway-login2]/project/gilad/ittai/roach/testing% paste test1 test2 | head                                    
chr	1
chr1	885000
chr1	1085000
chr1	1445000
chr1	1765000
chr1	2085000
chr1	2125000
chr1	2365000
chr1	2835000
chr1	2875000
[ittai@midway-login2]/project/gilad/ittai/roach/testing% head test.txt
chr	dist.bin	loc1_bin	loc2_bin	count	percentile	homer.hit	my.hits.tot	homer.hits.tot	overlap
chr1	0	[   885000,   925000)	[   925000,   965000)	18	80	0	1083	0	0
chr1	0	[  1085000,  1125000)	[  1125000,  1165000)	24	80	0	1083	0	0
chr1	0	[  1445000,  1485000)	[  1485000,  1525000)	24	80	0	1083	0	0
chr1	0	[  1765000,  1805000)	[  1805000,  1845000)	24	80	0	1083	0	0
chr1	0	[  2085000,  2125000)	[  2125000,  2165000)	26	80	0	1083	0	0
chr1	0	[  2125000,  2165000)	[  2165000,  2205000)	20	80	0	1083	0	0
chr1	0	[  2365000,  2405000)	[  2405000,  2445000)	18	80	0	1083	0	0
chr1	0	[  2835000,  2875000)	[  2875000,  2915000)	18	80	0	1083	0	0
chr1	0	[  2875000,  2915000)	[  2915000,  2955000)	18	80	0	1083	0	0
[ittai@midway-login2]/project/gilad/ittai/roach/testing% 

##For prepping a BED file for analysis with bedtools from full homer overlap DFs.
cut -f1 test.txt | tail -n +2 > chr.tmp
cut -f1 test.txt | tail -n +2 >> chr.tmp

cut -f3 test.txt | cut -d ',' -f1 | grep -o '[0-9]\+' | tail -n +2 > start.tmp
cut -f3 test.txt | cut -d ',' -f2 | grep -o '[0-9]\+' | tail -n +2 > end.tmp

cut -f4 test.txt | cut -d ',' -f1 | grep -o '[0-9]\+' | tail -n +2 >> start.tmp
cut -f4 test.txt | cut -d ',' -f2 | grep -o '[0-9]\+' | tail -n +2 >> end.tmp

paste chr.tmp start.tmp end.tmp > /project/gilad/ittai/roach/full_homer_overlap_dfs/bedtools_files/

for file in 10kb 20kb 30kb 40kb 50kb 60kb 70kb 80kb 90kb 100kb 125kb 150kb 200kb; do
	cut -f1 /project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.$file | tail -n +2 > chr.tmp
	cut -f1 /project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.$file | tail -n +2 >> chr.tmp
	
	cut -f3 /project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.$file | cut -d ',' -f1 | grep -o '[0-9]\+' | tail -n +2 > start.tmp
	cut -f3 /project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.$file | cut -d ',' -f2 | grep -o '[0-9]\+' | tail -n +2 > end.tmp
	
	cut -f4 /project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.$file | cut -d ',' -f1 | grep -o '[0-9]\+' | tail -n +2 >> start.tmp
	cut -f4 /project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.$file | cut -d ',' -f2 | grep -o '[0-9]\+' | tail -n +2 >> end.tmp
	
	paste chr.tmp start.tmp end.tmp > /project/gilad/ittai/roach/full_homer_overlap_dfs/bedtools_files/bedfile.$file
	
	rm *.tmp
done

#For extracting GC content:
for file in 10kb 30kb 40kb 50kb 60kb 80kb 200kb; do
	cut -f5 /project/gilad/ittai/roach/bedtools_files/nuc/nuc_$file.bed | tail -n +2 > /project/gilad/ittai/roach/bedtools_files/tmp.GC
	paste /project/gilad/ittai/roach/bedtools_files/info$file.bed /project/gilad/ittai/roach/bedtools_files/tmp.GC > /project/gilad/ittai/roach/bedtools_files/filled_info/infoGC_$file.bed
	rm /project/gilad/ittai/roach/bedtools_files/tmp.GC
done


#For extracting TSS content:
for file in 10kb 30kb 40kb 50kb 60kb 80kb 200kb; do
	sort -k4,4n /project/gilad/ittai/roach/bedtools_files/TSS/TSS_$file.bed > /project/gilad/ittai/roach/bedtools_files/TSS/sorted/sorted_TSS_$file.bed
	cut -f11 /project/gilad/ittai/roach/bedtools_files/TSS/sorted/sorted_TSS_$file.bed > /project/gilad/ittai/roach/bedtools_files/tmp.TSS
	paste /project/gilad/ittai/roach/bedtools_files/filled_info/infoGC_$file.bed /project/gilad/ittai/roach/bedtools_files/tmp.TSS > /project/gilad/ittai/roach/bedtools_files/filled_info/infoGC_TSS_$file.bed
	rm /project/gilad/ittai/roach/bedtools_files/tmp.TSS
done

#For extracting gene_body info content:
for file in 10kb 30kb 40kb 50kb 60kb 80kb 200kb; do
	sort -k4,4n /project/gilad/ittai/roach/bedtools_files/gene_body/GB_$file.bed > /project/gilad/ittai/roach/bedtools_files/gene_body/unsorted/unsorted_GB_$file.bed
	cut -f11 /project/gilad/ittai/roach/bedtools_files/gene_body/unsorted/unsorted_GB_$file.bed > /project/gilad/ittai/roach/bedtools_files/tmp.GB
	paste /project/gilad/ittai/roach/bedtools_files/filled_info/infoGC_TSS_$file.bed /project/gilad/ittai/roach/bedtools_files/tmp.GB > /project/gilad/ittai/roach/bedtools_files/filled_info/infoGC_TSS_GB_$file.bed
	rm /project/gilad/ittai/roach/bedtools_files/tmp.GB
done

#For re-pairing info
split -l 9571938 infoGC_TSS_GB_10kb.bed 10kb
split -l 3297298 infoGC_TSS_GB_30kb.bed 30kb
split -l 2432571 infoGC_TSS_GB_40kb.bed 40kb
split -l 1899678 infoGC_TSS_GB_50kb.bed 50kb
split -l 1551687 infoGC_TSS_GB_60kb.bed 60kb
split -l 1112807 infoGC_TSS_GB_80kb.bed 80kb
split -l 327711 infoGC_TSS_GB_200kb.bed 200kb

#For re-pairing the info and then appending it to full homer overlap DFs!
for file in 10kb 30kb 40kb 50kb 60kb 80kb 200kb; do
	paste /project/gilad/ittai/roach/bedtools_files/filled_info/one$file /project/gilad/ittai/roach/bedtools_files/filled_info/two$file | cut -f4-6,10-12 > /project/gilad/ittai/roach/bedtools_files/filled_info/$file.tmp
	cat /project/gilad/ittai/roach/bedtools_files/filled_info/header.txt /project/gilad/ittai/roach/bedtools_files/filled_info/$file.tmp > /project/gilad/ittai/roach/bedtools_files/filled_info/header.tmp
	paste /project/gilad/ittai/roach/full_homer_overlap_dfs/full_homer_overlap_df.$file /project/gilad/ittai/roach/bedtools_files/filled_info/header.tmp > /project/gilad/ittai/roach/full_homer_overlap_dfs/extracols/extra_info.$file
	rm /project/gilad/ittai/roach/bedtools_files/filled_info/$file.tmp
	rm /project/gilad/ittai/roach/bedtools_files/filled_info/header.tmp
done


#For extracting TSS content:
for file in 10kb 30kb 40kb 50kb 60kb 80kb 200kb; do
	sort -k4,4n /project/gilad/ittai/roach/bedtools_files/TSS/TSS_$file.bed > /project/gilad/ittai/roach/bedtools_files/TSS/sorted/sorted_TSS_$file.bed
	cut -f11 /project/gilad/ittai/roach/bedtools_files/TSS/sorted/sorted_TSS_$file.bed > /project/gilad/ittai/roach/bedtools_files/tmp.TSS
	paste /project/gilad/ittai/roach/bedtools_files/filled_info/infoGC_$file.bed /project/gilad/ittai/roach/bedtools_files/tmp.TSS > /project/gilad/ittai/roach/bedtools_files/filled_info/infoGC_TSS_$file.bed
	rm /project/gilad/ittai/roach/bedtools_files/tmp.TSS
done

Something VERY weird with some of the sorted info, for, say, 20kb!
head sorted_20kb.bed
10000000 10020000  10036688 
100035000 100055000  10037982 
100095000 100115000  10042090 
100155000 100175000  10037983 
10020000 10040000  10035729 
10020000 10040000  10036689 
10020000 10040000  10038814 
10020000 10040000  10040864 
10020000 10040000  10040865 
10020000 10040000  10042965 


Error: malformed BED entry at line 4309536. Start was greater than end. Exiting.
#Found in 20kb file

Error: malformed BED entry at line 253497. Start was greater than end. Exiting.
#Found in 90 kb file: head -253507 info90kb.bed | tail -20
chr4	35000	125000
chr4	0	0
chr4	125000	215000
chr4	0	0
chr4	215000	305000
chr4	0	0
chr4	305000	395000
chr4	0	0
chr4	485000	575000
chr4	0	0
chr4	575000	665000
chr4	0	0

Error: malformed BED entry at line 266886. Start was greater than end. Exiting.
#Similar at 100kb:
chr5	138960000	139060000
chr5	141060000	141160000
chr5	15000	115000
chr5	0	0
chr5	315000	415000
chr5	0	0
chr5	1015000	1115000
chr5	0	0
chr5	1315000	1415000
chr5	0	0
chr5	1415000	1515000
chr5	0	0
chr5	1615000	1715000
chr5	0	0


Error: malformed BED entry at line 322605. Start was greater than end. Exiting.
#125 kb: head -322615 info125kb.bed | tail -20
chr8	145775000	145900000
chr8	22060000	22185000
chr8	94775000	94900000
chr8	145775000	145900000
chr8	685000	810000
chr8	0	0
chr8	1435000	1560000
chr8	0	0
chr8	1685000	1810000
chr8	0	0
chr8	1935000	2060000
chr8	0	0
chr8	2185000	2310000
chr8	0	0


Error: malformed BED entry at line 256136. Start was greater than end. Exiting.
#150kb, similar attribution

#To merge .homer hicup-converted files:
cat 40300_Jan17.homer > merged.homer
awk -v OFS='\t' '{$1+=201305065}1' 40300_Mar17.homer >> merged.homer

#To get rownamed ortho files:
echo orthorownumber > rownames.txt
cat -n hum_orthos.bed | cut -f1 | head -n -1 >> rownames.txt
awk '{$1=$1};1' rownames.txt > rownames.txt
paste hum_orthos.bed rownames.txt > hum_orthos_numbered.bed

#Actually forget this, just take the gene name along with the other info! So now to obtain chr22 bed files from the ortho exo mega file:
awk '{print($3,"\t"$4"\t"$5"\t"$1)}' metaOrthoExonTrios.0.92.0.96.wExonEnsID.wSymbols.txt > humorthos.txt #for humans
awk '{print($7,"\t"$8"\t"$9"\t"$1)}' metaOrthoExonTrios.0.92.0.96.wExonEnsID.wSymbols.txt > chimporthos.txt #for chimps

#This is good, but now I realize what would actually be best is a condensed version of this file--orthologous genes, rather than exons. Can condense this by looking at unique hits for gene ID and choosing min start and max end?

#To obtain chr22 bed files from homer significance files:
grep chr22 sigs10kb.txt > just22_10kb.txt
cut -f3-5 just22_10kb.txt > pair1.tmp
cut -f9-11 just22_10kb.txt > pair2.tmp
cut -f1 just22_10kb.txt > rowids.tmp
paste pair1.tmp rowids.tmp > pair1.22.bed
paste pair2.tmp rowids.tmp > pair2.22.bed
rm *.tmp

bedtools intersect -wa -u -a pair1.22.bed -b hum_ortho_22.bed > pair1.orthos.bed
bedtools intersect -wa -u -a pair2.22.bed -b hum_ortho_22.bed > pair2.orthos.bed

join -1 4 -2 4 pair1.orthos.bed pair2.orthos.bed > join.txt

#Pull out human and chimp info (w/ strand) from Ran's file
cut -f 1,3-10 metaOrthoExonTrios.0.92.0.96.wExonEnsID.wSymbols.txt > HumChimp.txt

#Pull out gene info and human and chimp info, prep bed files for each species
cut -f1 HumChimp.txt| tail -n +2 > genes.tmp
cut -f2-5 HumChimp.txt | tail -n +2 > humbed.tmp
cut -f6-9 HumChimp.txt | tail -n +2 > chimpbed.tmp
paste humbed.tmp genes.tmp > hum_orthos.txt
paste chimpbed.tmp genes.tmp > chimp_orthos.txt
rm *.tmp

#Pull out first exon in each gene!
bedtools groupby -i hum_orthos.txt -g 5 -c 2,1,4 -o min,distinct,distinct > humstarts.txt
paste <(cut -f3 humstarts.txt) <(awk '$2-=2000' humstarts.txt | cut -d " " -f2) <(cut -f2 humstarts.txt) <(cut -f1 humstarts.txt) <(cut -f4 humstarts.txt) > hg19promoters.bed #Create the new file, BED format w/ strand and gene ID after. 2 kb upstream of first exon site to start of first exon site is region

#Pull out first exon in each gene!
bedtools groupby -i chimp_orthos.txt -g 5 -c 2,1,4 -o min,distinct,distinct > chimpstarts.txt
paste <(cut -f3 chimpstarts.txt) <(awk '$2-=2000' chimpstarts.txt | cut -d " " -f2) <(cut -f2 chimpstarts.txt) <(cut -f1 chimpstarts.txt) <(cut -f4 chimpstarts.txt) > pt3promoters.bed #Create the new file, BED format w/ strand and gene ID after. 2 kb upstream of first exon site to start of first exon site is region

#Prepping homer sig files for R analysis:
for file in *; do
	cut -f2,8,17,18 $file > $file.R
done


#Getting corrvecs:
#!/bin/bash

mkdir corrvecs

for file in *.txt; do
	grep 'chr18' $file > $file.tmp
done

join 21792_28126_50kb.corrDiff.txt.tmp 21792_28815_50kb.corrDiff.txt.tmp > join1.txt
join join1.txt 21792_28834_50kb.corrDiff.txt.tmp > join2.txt
join join2.txt 28126_28815_50kb.corrDiff.txt.tmp > join3.txt
join join3.txt 28126_28834_50kb.corrDiff.txt.tmp > join4.txt
join join4.txt 28815_28834_50kb.corrDiff.txt.tmp > joinfinal.txt

cut -d " " -f6 joinfinal.txt > corrvecs/21792_28126.corrvec
cut -d " " -f11 joinfinal.txt > corrvecs/21792_28815.corrvec
cut -d " " -f16 joinfinal.txt > corrvecs/21792_28834.corrvec
cut -d " " -f21 joinfinal.txt > corrvecs/28126_28815.corrvec
cut -d " " -f26 joinfinal.txt > corrvecs/28126_28834.corrvec
cut -d " " -f31 joinfinal.txt > corrvecs/28815_28834.corrvec

rm *.tmp
rm join*




#!/bin/bash
#SBATCH --job-name=liftover21792_mate1
#SBATCH --output=liftover21792_mate1-%j.out
#SBATCH --error=liftover21792_mate1-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=END

#!/bin/bash
#to find the orthologous region between human and chimp given a H or C file. the idea is: e.g. human to chimp:  H->C, C->H, H<->C (using reciprical chain), the second step should only have very few entries missing, caused by a large gap in H side. otherwise, the second step should already reaches a very stable state.

LIFTOVER=~/progs/miniconda2/bin/liftOver 
HCCHAIN=/project/gilad/ittai/ortho/hg19.panTro3.rbest.chain.gz
CHCHAIN=/project/gilad/ittai/ortho/panTro3.hg19.rbest.chain.gz

SPECIES=H
INPUTBED=/project/gilad/ittai/results_Mar17/committee/sig_interaxns/100kb/mate1_21792.100kb.txt.bed
BEDPLUS=3 # input bed file format


if [ "${SPECIES}" = H ]; then
	echo "\
	${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS} ${INPUTBED} ${HCCHAIN} ${INPUTBED}.panTro3 ${INPUTBED}.HtoC.unmapped ; \
	${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.panTro3 ${CHCHAIN} ${INPUTBED}.ortho.hg19 ${INPUTBED}.CtoH.unmapped ; \
	${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.ortho.hg19 ${HCCHAIN} ${INPUTBED}.ortho.panTro3 ${INPUTBED}.final.unmapped"| qsub -l h_vmem=2g -wd `pwd` -N ${SPECIES}_OrthoR

fi

if [ "${SPECIES}" = C ]; then
	echo "\
	${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS} ${INPUTBED} ${CHCHAIN} ${INPUTBED}.hg19 ${INPUTBED}.CtoH.unmapped ; \
	${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.hg19 ${HCCHAIN} ${INPUTBED}.ortho.panTro3 ${INPUTBED}.HtoC.unmapped ; \
	${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.ortho.panTro3 ${CHCHAIN} ${INPUTBED}.ortho.hg19 ${INPUTBED}.final.unmapped"| qsub -l h_vmem=2g -wd `pwd` -N ${SPECIES}_OrthoR

fi

#!/bin/bash
#SBATCH --job-name=liftover21792_mate1
#SBATCH --output=liftover21792_mate1-%j.out
#SBATCH --error=liftover21792_mate1-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=END

LIFTOVER=~/progs/miniconda2/bin/liftOver 
HCCHAIN=/project/gilad/ittai/ortho/hg19.panTro3.rbest.chain.gz
CHCHAIN=/project/gilad/ittai/ortho/panTro3.hg19.rbest.chain.gz

INPUTBED=/project/gilad/ittai/results_Mar17/committee/sig_interaxns/100kb/mate1_21792.100kb.txt.bed
BEDPLUS=3 # input bed file format

${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS} ${INPUTBED} ${HCCHAIN} ${INPUTBED}.panTro3 ${INPUTBED}.HtoC.unmapped
${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.panTro3 ${CHCHAIN} ${INPUTBED}.ortho.hg19 ${INPUTBED}.CtoH.unmapped
${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.ortho.hg19 ${HCCHAIN} ${INPUTBED}.ortho.panTro3 ${INPUTBED}.final.unmapped





#!/bin/bash

for file in *; do
	head -n1 $file > headers.tmp
	paste <(cut -f3-6 $file) <(cut -f1 $file) <(cut -f17-19 $file) | tail -n +2 > mate1_$file.bed 
	paste <(cut -f9-12 $file) <(cut -f1 $file) <(cut -f17-19 $file) | tail -n +2 > mate2_$file.bed 
	
	#Get info on whether a homer-significant hit overlaps ANY putative promoter.
	bedtools intersect -wa -wb -u -a mate1_$file.bed -b /project/gilad/ittai/results_Mar17/committee/sig_interaxns/ortho_analysis/prom.hg.18.bed > mate1overlap$file.bed
	bedtools intersect -wa -wb -u -a mate2_$file.bed -b /project/gilad/ittai/results_Mar17/committee/sig_interaxns/ortho_analysis/prom.hg.18.bed > mate2overlap$file.bed
	
	#Get info on WHICH putative promoter homer-significant hits overlap.
	bedtools intersect -wa -wb -a mate1_$file.bed -b /project/gilad/ittai/results_Mar17/committee/sig_interaxns/ortho_analysis/prom.hg.18.bed > mate1promoters$file.bed
	bedtools intersect -wa -wb -a mate2_$file.bed -b /project/gilad/ittai/results_Mar17/committee/sig_interaxns/ortho_analysis/prom.hg.18.bed > mate2promoters$file.bed
	
done



#!/bin/bash
tail -n +2 chimp25.chr18 > c25mate1.18
paste <(cut -f1,4,5 chimp25.chr18) <(cut -f2,3,6-9 chimp25.chr18) |tail -n +2 > c25mate2.18

paste c25mate1.18 <(cat -n c25mate1.18 | cut -f1) > mate1s.txt
paste c25mate2.18 <(cat -n c25mate2.18 | cut -f1) > mate2s.txt


LIFTOVER=~/progs/miniconda2/bin/liftOver 
HCCHAIN=/project/gilad/ittai/ortho/hg19.panTro3.rbest.chain.gz
CHCHAIN=/project/gilad/ittai/ortho/panTro3.hg19.rbest.chain.gz

INPUTBED=/project/gilad/ittai/results_Mar17/committee/sig_interaxns/25kb/mate1s.txt
BEDPLUS=3 # input bed file format

${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS} ${INPUTBED} ${CHCHAIN} ${INPUTBED}.hg19 ${INPUTBED}.CtoH.unmapped
${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.hg19 ${HCCHAIN} ${INPUTBED}.ortho.panTro3 ${INPUTBED}.HtoC.unmapped
${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED}.ortho.panTro3 ${CHCHAIN} ${INPUTBED}.ortho.hg19 ${INPUTBED}.final.unmapped


INPUTBED2=/project/gilad/ittai/results_Mar17/committee/sig_interaxns/25kb/mate2s.txt

${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS} ${INPUTBED2} ${CHCHAIN} ${INPUTBED2}.hg19 ${INPUTBED2}.CtoH.unmapped
${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED2}.hg19 ${HCCHAIN} ${INPUTBED2}.ortho.panTro3 ${INPUTBED}.HtoC.unmapped
${LIFTOVER} -minMatch=0.7 bedPlus=${BEDPLUS}  ${INPUTBED2}.ortho.panTro3 ${CHCHAIN} ${INPUTBED2}.ortho.hg19 ${INPUTBED}.final.unmapped


sort -k9 -n mate1s.txt.ortho.hg19 > sorted_mate1.hg19
sort -k9 -n mate2s.txt.ortho.hg19 > sorted_mate2.hg19

join -1 9 -2 9 sorted_mate1.hg19 sorted_mate2.hg19 > remaining.pairs.hg19.txt

sort -k9 -n mate1s.txt.ortho.panTro3 > sorted_mate1.PT3
sort -k9 -n mate2s.txt.ortho.panTro3 > sorted_mate2.PT3

join -1 9 -2 9 sorted_mate1.PT3 sorted_mate2.PT3 > remaining.pairs.PT3.txt

#The columns on this are human chr, start1, end1, start2, end2, #lines, Z-score; chimp chr, start1, end1, start2, end2
paste <(cut -d " " -f2-8 remaining.pairs.hg19.txt) <(cut -d " " -f2-6 remaining.pairs.PT3.txt) > H_C_remaining.txt