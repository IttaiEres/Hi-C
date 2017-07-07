##Script for dealing with orthologous exons, extracting TSS from humans and chimps (a 10kb window centered 5kb upstream of the appropriate exon--the first, for positive strand genes, and the last, for negative strand genes).

#Create new files for human and chimp, with chr, start, end, strand, ensemblID, and metaexon # as columns:
paste <(cut -f 3-6 ortho_fixed_exons_removed.txt) <(cut -f 1-2 ortho_fixed_exons_removed.txt) > human_full.tmp

paste <(cut -f 8-11 ortho_fixed_exons_removed.txt) <(cut -f 1-2 ortho_fixed_exons_removed.txt) > chimp_full.tmp

grep "\+" /project2/gilad/ittai/HiC/homersig_variance/ortho_work/human_full.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_pos.tmp
grep "\-" /project2/gilad/ittai/HiC/homersig_variance/ortho_work/human_full.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_negs.tmp


grep "\+" /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_full.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_pos.tmp
grep "\-" /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_full.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_negs.tmp

#Sort human files by ENSG ID
sort -k5,5 /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_negs.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_minus.tmp
sort -k5,5 /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_pos.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_plus.tmp

#Sort chimp files by ENSG ID 
sort -k5,5 /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_negs.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_minus.tmp
sort -k5,5 /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_pos.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_plus.tmp

#Bedtools groupby on all these files individually, grouping by ENSG ID, and taking min and max for + and - files, respectively.
bedtools groupby -i /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_plus.tmp -full -g 5 -c 2 -o min > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSpos.tmp
bedtools groupby -i /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_plus.tmp -full -g 5 -c 2 -o min > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSpos.tmp

bedtools groupby -i /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hum_minus.tmp -full -g 5 -c 3 -o max > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSneg.tmp
bedtools groupby -i /project2/gilad/ittai/HiC/homersig_variance/ortho_work/chimp_minus.tmp -full -g 5 -c 3 -o max > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSneg.tmp

#Now, take the outputs of the bedtools groupby on humans and chimps, and make two separate TSS bed files per species (+ and -) w/ ranges being a single-nucleotide interval upstream of exon start. Need to make two separate files due to wanting the TSS to fall in the last 24kb of a bin if it's on the positive strand, or the first 24kb of a bin if it's on the negative strand.
awk -v OFS="\t" '{print $1,$7-1,$7}' /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSpos.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSplus.txt
awk -v OFS="\t" '{print $1,$7,$7+1}' /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSneg.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSminus.txt

awk -v OFS="\t" '{print $1,$7-1,$7}' /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSpos.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSplus.txt
awk -v OFS="\t" '{print $1,$7,$7+1}' /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSneg.tmp > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSminus.txt

#Prepping a sigdf from my analysis as 2 bed files to run against these TSS files w/ bedtools intersect--one for negative, one for positive! Need to change some of these values to reflect this sigdf.clean file no longer looking how it once did...
tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f15 | awk -v OFS="\t" -F "-" '{print $1,$2+1000,$2+25000,NR}' > homer_sigs_H_plus.bed
tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f21 | awk -v OFS="\t" -F "-" '{print $1,$2+1000,$2+25000,NR}' >> homer_sigs_H_plus.bed

tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f15 | awk -v OFS="\t" -F "-" '{print $1,$2,$2+24000,NR}' > homer_sigs_H_minus.bed
tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f21 | awk -v OFS="\t" -F "-" '{print $1,$2,$2+24000,NR}' >> homer_sigs_H_minus.bed

tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f5 | awk -v OFS="\t" -F "-" '{print $1,$2+1000,$2+25000,NR}' > homer_sigs_C_plus.bed
tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f6 | awk -v OFS="\t" -F "-" '{print $1,$2+1000,$2+25000,NR}' >> homer_sigs_C_plus.bed

tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f5 | awk -v OFS="\t" -F "-" '{print $1,$2,$2+24000,NR}' > homer_sigs_C_minus.bed
tail -n +2 /project2/gilad/ittai/HiC/homersig_variance/sigdf.clean | cut -f6 | awk -v OFS="\t" -F "-" '{print $1,$2,$2+24000,NR}' >> homer_sigs_C_minus.bed

#Doing bedtools intersect on each of the homer_sigs files (positive and negative) to find overlap w/ TSS:
bedtools intersect -c -a /project2/gilad/ittai/HiC/homersig_variance/ortho_work/homer_sigs_H_plus.bed -b /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSplus.txt > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/H_homersigs_TSS_pos

bedtools intersect -c -a /project2/gilad/ittai/HiC/homersig_variance/ortho_work/homer_sigs_H_minus.bed -b /project2/gilad/ittai/HiC/homersig_variance/ortho_work/hTSSminus.txt > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/H_homersigs_TSS_neg

bedtools intersect -c -a /project2/gilad/ittai/HiC/homersig_variance/ortho_work/homer_sigs_C_plus.bed -b /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSplus.txt > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/C_homersigs_TSS_pos

bedtools intersect -c -a /project2/gilad/ittai/HiC/homersig_variance/ortho_work/homer_sigs_C_minus.bed -b /project2/gilad/ittai/HiC/homersig_variance/ortho_work/cTSSminus.txt > /project2/gilad/ittai/HiC/homersig_variance/ortho_work/C_homersigs_TSS_neg

#Split files using their line counts, ensuring proper re-pairing of mates.
split -l $(($(wc -l </project2/gilad/ittai/HiC/homersig_variance/H_homersigs_TSS_pos)/2)) --numeric-suffixes=1 /project2/gilad/ittai/HiC/homersig_variance/H_homersigs_TSS_pos /project2/gilad/ittai/HiC/homersig_variance/H_homersigs_TSS_pos

split -l $(($(wc -l </project2/gilad/ittai/HiC/homersig_variance/H_homersigs_TSS_neg)/2)) --numeric-suffixes=1 /project2/gilad/ittai/HiC/homersig_variance/H_homersigs_TSS_neg /project2/gilad/ittai/HiC/homersig_variance/H_homersigs_TSS_neg

split -l $(($(wc -l </project2/gilad/ittai/HiC/homersig_variance/C_homersigs_TSS_pos)/2)) --numeric-suffixes=1 /project2/gilad/ittai/HiC/homersig_variance/C_homersigs_TSS_pos /project2/gilad/ittai/HiC/homersig_variance/C_homersigs_TSS_pos

split -l $(($(wc -l </project2/gilad/ittai/HiC/homersig_variance/C_homersigs_TSS_neg)/2)) --numeric-suffixes=1 /project2/gilad/ittai/HiC/homersig_variance/C_homersigs_TSS_neg /project2/gilad/ittai/HiC/homersig_variance/C_homersigs_TSS_neg

#Now, create a file that indicates whether a TSS is found in either member of each pair (just binary, 0 for no, 1 for yes). Works by repairing and running awk on that output.
paste <(cut -f5 H_homersigs_TSS01) <(cut -f5 H_homersigs_TSS02) | awk -v OFS="\t" '{print $0,($1>0||$2>0?1:0)}' > hTSS_binary
paste <(cut -f5 C_homersigs_TSS01) <(cut -f5 C_homersigs_TSS02) | awk -v OFS="\t" '{print $0,($1>0||$2>0?1:0)}' > cTSS_binary

#Append the last column of that file (binary TSS status) onto the sigdf.clean dataframe. Now just need to incorporate all pieces of info into this dataframe, and we should be good to go! Might as well keep FULL homer info


#Now simply cut the last column of that file and paste it on to the sig_DF with all the other information from each individual line.