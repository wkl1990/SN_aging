#!/bin/bash

# example: bash 03c.naiveoverlap.sh -i c1 -p /snapatac2 -o /snapatac2 

while getopts ':i:p:o:h' opt; do
  case "$opt" in
    i)
      input="$OPTARG"
      echo "Processing option 'i' with '${OPTARG}' argument"
      ;;
    p)
      path="$OPTARG"
      echo "Processing option 'p' with '${OPTARG}' argument"
      ;;
    o)
      output="$OPTARG"
      echo "Processing option 'o' with '${OPTARG}' argument"
      ;;
    h)
      echo "Usage: $(basename $0) [-i input] [-p path] [-o output]"
      exit 0
      ;;
    :)
      echo -e "option requires an argument.\nUsage: $(basename $0) [-i input] [-p path] [-o output]"
      exit 1
      ;;
    ?)
      echo -e "Invalid command option.\nUsage: $(basename $0) [-i input] [-p path] [-o output]"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"


if [ ! -d "${output}/parsePeak" ]; then
  mkdir ${output}/parsePeak
fi

#BLACKLIST=/projects/ps-renlab/yangli/genome/mm10/mm10.blacklist.bed.gz

POOLED_PEAK="${path}/peaks/${input}_peaks.narrowPeak"
POOLED_SUMMIT="${path}/peaks/${input}.summits.bed"
REP1_PEAK="${path}/peaks/${input}.rep1_peaks.narrowPeak"
REP2_PEAK="${path}/peaks/${input}.rep2_peaks.narrowPeak"
PSEUDO1_PEAK="${path}/peaks/${input}.pseudo1_peaks.narrowPeak"
PSEUDO2_PEAK="${path}/peaks/${input}.pseudo2_peaks.narrowPeak"
PooledInRep1AndRep2="${output}/parsePeak/${input}.PooledInRep1AndRep2.narrowPeak.gz"
PooledInPsRep1AndPsRep2="${output}/parsePeak/${input}.PooledInPsRep1AndPsRep2.narrowPeak.gz"
naivePeakList="${output}/parsePeak/${input}.naivePeakList.narrowPeak.gz"
naiveSummitList="${output}/parsePeak/${input}.naiveSummitList.bed"

### Naive overlap
# Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
bedtools intersect -wo -a ${POOLED_PEAK} -b ${REP1_PEAK} -nonamecheck \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | bedtools intersect -wo -a stdin -b ${REP2_PEAK} \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInRep1AndRep2}

# Find pooled peaks that overlap PooledPseudoRep1 and PooledPseudoRep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5
bedtools intersect -wo -a ${POOLED_PEAK} -b ${PSEUDO1_PEAK} -nonamecheck \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq \
| bedtools intersect -wo -a stdin -b ${PSEUDO2_PEAK} -nonamecheck \
| awk 'BEGIN{FS=OFS="\t"}($1~/chr*/){s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort -k1,1 -k2,2n | uniq | gzip -c > ${PooledInPsRep1AndPsRep2}

# Combine peak lists
zcat ${PooledInRep1AndRep2} ${PooledInPsRep1AndPsRep2}  | sort -k1,1 -k2,2n | uniq | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[0-9XY]+(?!_)' | gzip -c  > ${naivePeakList}

# Get summit 
if [ ! -f ${POOLED_SUMMIT} ] || [ ! -s ${POOLED_SUMMIT} ]; then
	cat ${POOLED_PEAK} | awk 'BEGIN{FS=OFS="\t"}{s1=$2+$10; s2=$2+$10+1}{print $1,s1,s2,$4,$9}' > ${POOLED_SUMMIT}
fi
join -1 1 -2 4 <(zcat ${naivePeakList} | cut -f 4 | sort) <(sort -k4,4 ${POOLED_SUMMIT}) -t$'\t' \
| awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1,$5}' | sort -k1,1 -k2,2n > ${naiveSummitList} 

if [ -f ${naiveSummitList} ] && [ -s ${naiveSummitList} ]; then
	echo -e "${input}\t${naiveSummitList}" >> ${output}/SN.naiveSummitList.list 
	# summary peaks in cell type
	npeak=`cat ${naiveSummitList} | wc -l`
	echo -e "${input}\t${npeak}" >> ${output}/SN.celltype.npeak4anno.txt
fi

