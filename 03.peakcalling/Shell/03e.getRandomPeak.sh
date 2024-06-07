#!/bin/bash

#!/bin/bash

path1=/peak_calling/macs2/

peakFile="SN_integra.macs2.union.peaks.nospm.bed"
shuffleFile="SN_integra.macs2.shuffle.bed"
shuffleFilter="SN_integra.macs2.shuffle.removeovlp.bed"

# * randomly shuffled regions
bedtools shuffle -seed 10 -i ${path1}/${peakFile} -g mm10.chrom.sizes > ${path1}/${shuffleFile}

# wget https://downloads.wenglab.org/V3/mm10-cCREs.bed
# cp mm10-rDHS-Filtered.bed ${path1}
# cp mm10-cCREs.bed ${path1}/mm10-cCREs.v3.bed

bedtools intersect -v -wa -a ${path1}/${shuffleFile} -b ${path1}/${peakFile} \
    | bedtools intersect -v -wa -a - -b mm10.blacklist.bed.gz \
    | bedtools intersect -v -wa -a - -b ${path1}/mm10-rDHS-Filtered.bed \
    | bedtools intersect -v -wa -a - -b ${path1}/mm10-cCREs.v3.bed \
    | grep -v chrUn \
    | grep -v random \
    | grep -v chrM \
    | sort -k1,1 -k2,2n | uniq > ${path1}/${shuffleFilter}
sed -i -e "s/peak/random/g" ${path1}/${shuffleFilter}

#awk 'BEGIN{FS=OFS="\t"}{print $1":"$2"-"$3,$4}' ${path1}/${shuffleFilter} > mba.whole.random.coor2peak.txt

