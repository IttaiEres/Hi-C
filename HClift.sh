#!/bin/bash

###### Script for subsetting down homer significant hits to only those that can be reciprocally lifted over.
###### Does the liftover with bed files on individual pairs, then merges on the pair interaction IDs.

homer.file=$1
out=$2
wd=$3


###Grab the chr, start, and end positions from each file; remove header line and add interaction ID as 4th column. Preps bedfiles.
cut -f3-5 $wd/${homer.file} | tail -n +2 | cat -n > $wd/${homer.file}_m1.tmp
paste <(cut -f2-4 $wd/${homer.file}_m1.tmp) <(cut -f1 $wd/${homer.file}_m1.tmp) > $wd/${homer.file}_mate1

cut -f9-11 $wd/${homer.file} | tail -n +2 | cat -n > $wd/${homer.file}_m2.tmp
paste <(cut -f2-4 $wd/${homer.file}_m2.tmp) <(cut -f1 $wd/${homer.file}_m2.tmp) > $wd/${homer.file}_mate2

rm *.tmp

#Run both sets of mate files through liftover!
bash ~/progs/HC_liftover_ortho.sh H $wd/${homer.file}_mate1
bash ~/progs/HC_liftover_ortho.sh H $wd/${homer.file}_mate2

#Sort the lifted over hg19 mate files for joining.
sort -k 4b,4 $wd/${homer.file}_mate1.ortho.hg19 > $wd/${homer.file}_mate1.hg19.sorted
sort -k 4b,4 $wd/${homer.file}_mate2.ortho.hg19 > $wd/${homer.file}_mate2.hg19.sorted

#Join hg mate files together!
join -t $'\t' -j 4 $wd/${homer.file}_mate1.hg19.sorted $wd/${homer.file}_mate2.hg19.sorted > $wd/$out.hg19.tmp1

#Sort the lifted over panTro3 mate files for joining.
sort -k 4b,4 $wd/${homer.file}_mate1.ortho.panTro3 > $wd/${homer.file}_mate1.panTro3.sorted
sort -k 4b,4 $wd/${homer.file}_mate2.ortho.panTro3 > $wd/${homer.file}_mate2.panTro3.sorted

#Join panTro3 mate files together!
join -t $'\t' -j 4 $wd/${homer.file}_mate1.panTro3.sorted $wd/${homer.file}_mate2.panTro3.sorted > $wd/$out.panTro3.tmp1

#Join the hg19 and panTro3 mate files together!
join -t $'\t' $wd/$out.hg19.tmp1 $wd/$out.panTro3.tmp1 > full.tmp1

#Get a homer.file with 'interaction' stripped from interaction IDs. Sort it and join with the joint file--sanity check, and full info! :)
tail -n +2 $wd/${homer.file} | cat -n | sort -k 1b,1 | awk -v OFS='\t' '$1=$1' > $wd/full.tmp2
join -t $'\t' $wd/full.tmp1 $wd/full.tmp2 > $wd/final.ortho.tmp

#Add the header to the final file!
echo -e 'num\tH_chr1\tstart1\tend1\tH_chr2\tstart2\tend2\tC_chr1\tstart1\tend1\tC_chr2\tstart2\tend2\tInteractionID\tPeakID(1)\tchr(1)\tstart(1)\tend(1)\tstrand(1)\tTotal_Reads(1)\tPeakID(2)\tchr(2)\tstart(2)\tend(2)\tstrand(2)\tTotal_Reads(2)\tDistance\tInteraction_Reads\tExpected_Reads\tZ-score\tLogP\tFDR\tCircos_Thickness' | cat - $wd/final.ortho.tmp > $wd/final.ortho${homer.file}

#Get rid of garbage
rm $wd/*.sorted
rm $wd/*tmp*
rm $wd/*mate*

