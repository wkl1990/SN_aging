#!/bin/bash


# run liftOver
path0=/SN_aging
path1=/ldsc/after_integra/celltype

# liftOver and merge
#wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz

# The putative enhancer regions were mapped to the human genome (hg38) using liftOver, with a strategy similar to previous reports73. Each region was required to both uniquely map to hg38, and to uniquely map back to the original region in mm10, with the requirement that >=50% of the bases in each region were mapped back to mouse after being mapped to human.

function loadavg {
    while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt 500 ]; do sleep 120; date; done;
}

cd ${path1}
mkdir raw
cd raw
ln -s /peak_calling/after_integra/final/SN_integra.final.peak.srt.bed .

finalPeak="SN_integra.final.peak.srt.bed"
i=50

liftOver ${finalPeak} /annotation/mm10ToHg38.over.chain.gz SN_integra.final.peak.srt.mm10ToHg38.bed SN_integra.final.peak.srt.mm10ToHg38.unmapped.bed -minMatch=0.${i} &
liftOver SN_integra.final.peak.srt.mm10ToHg38.bed /annotation/hg38ToMm10.over.chain.gz SN_integra.final.peak.srt.mm10ToHg38.back2mm10.bed SN_integra.final.peak.srt.mm10ToHg38.back2mm10.unmapped -minMatch=0.10 &
# check if exactly match to orginal mm10
bedtools intersect -wao -r -f 0.5 -a SN_integra.final.peak.srt.mm10ToHg38.back2mm10.bed -b SN_integra.final.peak.srt.bed | awk '$6!="."' | awk '$4==$8' | cut -f 4 | sort | uniq > SN_integra.final.peak.srt.mm10ToHg38.back2mm10.matched.peaks
join -1 1 -2 4 SN_integra.final.peak.srt.mm10ToHg38.back2mm10.matched.peaks <(sort -k4,4 SN_integra.final.peak.srt.mm10ToHg38.bed | uniq) -t$'\t' | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n | uniq > SN_integra.final.peak.srt.mm10ToHg38.map2mm10.bed

ln -s SN_integra.final.peak.srt.mm10ToHg38.map2mm10.bed SN_integra.final.peak.srt.reciprocalToHg38.bed
