#!/bin/bash


# run GWAS (ldsc)
path0=/SN_aging
path1=/ldsc/after_integra/celltype

# liftOver and merge
#wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg19.over.chain.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMm10.over.chain.gz

# The putative enhancer regions were mapped to the human genome (hg19) using liftOver, with a strategy similar to previous reports73. Each region was required to both uniquely map to hg19, and to uniquely map back to the original region in mm10, with the requirement that >=50% of the bases in each region were mapped back to mouse after being mapped to human.

function loadavg {
    while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt 500 ]; do sleep 120; date; done;
}

cd ${path1}
mkdir raw
cd raw
ln -s /peak_calling/after_integra/final/SN_integra.final.peak.srt.bed .

#wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg19.over.chain.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMm10.over.chain.gz

finalPeak="SN_integra.final.peak.srt.bed"
i=50

liftOver ${finalPeak} /annotation/mm10ToHg19.over.chain.gz SN_integra.final.peak.srt.mm10ToHg19.bed SN_integra.final.peak.srt.mm10ToHg19.unmapped.bed -minMatch=0.${i} &
liftOver SN_integra.final.peak.srt.mm10ToHg19.bed /annotation/hg19ToMm10.over.chain.gz SN_integra.final.peak.srt.mm10ToHg19.back2mm10.bed SN_integra.final.peak.srt.mm10ToHg19.back2mm10.unmapped -minMatch=0.10 &
# check if exactly match to orginal mm10
bedtools intersect -wao -r -f 0.5 -a SN_integra.final.peak.srt.mm10ToHg19.back2mm10.bed -b SN_integra.final.peak.srt.bed | awk '$6!="."' | awk '$4==$8' | cut -f 4 | sort | uniq > SN_integra.final.peak.srt.mm10ToHg19.back2mm10.matched.peaks
join -1 1 -2 4 SN_integra.final.peak.srt.mm10ToHg19.back2mm10.matched.peaks <(sort -k4,4 SN_integra.final.peak.srt.mm10ToHg19.bed | uniq) -t$'\t' | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n | uniq > SN_integra.final.peak.srt.mm10ToHg19.map2mm10.bed

ln -s SN_integra.final.peak.srt.mm10ToHg19.map2mm10.bed SN_integra.final.peak.srt.reciprocalToHg19.bed


cd ${path1}
mkdir data
# overlap with liftOver peaks
ln -s /WK_aging/peak_calling/after_integra/SN.cell2anno.list SN.celltype.list
ln -s /WK_aging/peak_calling/after_integra/final/subclass.final.peak.srt .


for i in `cat SN.celltype.list`;
do echo $i;
join -1 1 -2 4 <(cut -f 4 subclass.final.peak.srt/${i}.bed | sort | uniq) <(sort -k4,4 raw/SN_integra.final.peak.srt.reciprocalToHg19.bed) | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n | uniq > data/${i}.peaks.reciprocalToHg19.bed
done

cd data
ln -s ../raw/SN_integra.final.peak.srt.mm10ToHg19.map2mm10.bed SN_integra.final.peak.srt.peaks.reciprocalToHg19.bed

# using homology peaks in default cutoff
path1=/WK_aging/ldsc/after_integra/celltype
cd $path1
mkdir qsub log pbs_log pbs_script ld_score

# run ld_score
gs=/Reference/genome/mm10/mm10.chrom.sizes
plink=/resource/gwas/1000G_phase3/1000G_EUR_Phase3_plink

function loadavg {
    while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt 500 ]; do sleep 60; date; done;
}

conda activate ldsc

for i in data/*.peaks.reciprocalToHg19.bed
do
    k=$(basename $i .peaks.reciprocalToHg19.bed)
    echo $k
    for j in {1..22}
    do
        /miniconda3/envs/ldsc/bin/python /home/kaw033/softwares/ldsc/make_annot.py \
            --bed-file $i \
            --bimfile $plink/1000G.EUR.QC.$j.bim \
            --annot-file ld_score/${k}.$j.annot.gz &
        sleep 3
        loadavg
    done
done


for i in data/*.peaks.reciprocalToHg19.bed
do
    k=$(basename $i .peaks.reciprocalToHg19.bed)
    echo $k
echo -e "
#!/bin/bash
path0=/WK_aging
path1=/WK_aging/ldsc/after_integra/celltype
cd \$path1

gs=/Reference/genome/mm10/mm10.chrom.sizes
plink=/resource/gwas/1000G_phase3/1000G_EUR_Phase3_plink

for j in {1..22}
do
    /miniconda3/envs/ldsc/bin/python /home/kaw033/softwares/ldsc/make_annot.py \
            --bed-file $i \
            --bimfile \$plink/1000G.EUR.QC.\$j.bim \
            --annot-file ld_score/${k}.\$j.annot.gz
done
" | sed '1d' > $path1/qsub/${k}.make_annot.qsub
qsub $path1/qsub/${k}.make_annot.qsub -o $path1/log/${k}.make_annot.log
done



## qsub_make.sh

DIR=/WK_aging/ldsc/after_integra/celltype
resources=/resource

for i in $(ls ld_score/*.22.annot.gz | sed 's/.22.annot.gz//g')
do
    j=$(basename $i)
    cat >pbs_script/$j.pbs <<EOF
#!/bin/bash

export PATH="/miniconda3/envs/ldsc/bin:\$PATH"
cd $DIR
source activate ldsc
for i in {1..22}
do
    /miniconda3/envs/ldsc/bin/python /home/kaw033/softwares/ldsc/ldsc.py\\
        --l2\\
        --bfile $resources/gwas/1000G_phase3/1000G_EUR_Phase3_plink/1000G.EUR.QC.\$i\\
        --ld-wind-cm 1\\
        --annot ld_score/$j.\$i.annot.gz\\
        --thin-annot\\
        --out ld_score/$j.\$i\\
        --print-snps $resources/gwas/1000G_phase3/1000G_Phase3_hapmap3/1000G_Phase3_hapmap3_print_snps.\$i.snp
done
source deactivate
EOF
done


for i in `ls pbs_script`;
do echo $i;
nohup bash pbs_script/$i 2>&1 > pbs_log/${i}.log &
sleep 3
loadavg
done;


# run ldsc_test
mkdir results

ls data/*.peaks.reciprocalToHg19.bed | grep -v "SN_integra.final.peak.srt" | sed 's/.peaks.reciprocalToHg19.bed//' | while read i
do 
    j=$(basename $i)
    printf "${j}\tld_score/${j}.,ld_score/SN_integra.final.peak.srt.\n"
done > SN_integra.final.peak.srt.ldcts

cts_name=SN_integra.final.peak.srt
DIR=/WK_aging/ldsc/after_integra/celltype
sumstats=/resource/gwas
resources=/resource

ls $sumstats/*sumstats.gz | grep sumstats | egrep "Nalls.LancetNeurol.2019.Parkinsons_disease|deLange.NatGenet.2017.Inflammatory_Bowel_Disease" | while read name
do
    i=$(basename $name .sumstats.gz)
    cat >pbs_script/$i.$cts_name.pbs <<EOF
#!/bin/bash

cd $DIR
export PATH="/miniconda3/envs/ldsc/bin:\$PATH"
source activate ldsc
/miniconda3/envs/ldsc/bin/python /home/kaw033/softwares/ldsc/ldsc.py\\
    --h2-cts $sumstats/$i.sumstats.gz \\
    --ref-ld-chr $resources/gwas/1000G_phase3/1000G_EUR_Phase3_baseline/baseline. \\
    --out results/${i}.$cts_name \\
    --ref-ld-chr-cts $cts_name.ldcts \\
    --w-ld-chr $resources/gwas/1000G_phase3/weights_hm3_no_hla/weights.
conda deactivate
EOF
done

for i in `ls pbs_script | grep .$cts_name`;
do echo $i;
qsub pbs_script/$i -o pbs_log/$i.log
done;


# for mediator run
for i in `ls pbs_script | grep .$cts_name`;
do echo $i;
nohup bash pbs_script/$i 2>&1 > pbs_log/${i}.log &
sleep 3
loadavg
done;

# summary
for set in SN_integra.final.peak.srt
do
    for i in results/*.${set}.cell_type_results.txt
    do
        j=$(basename $i .${set}.cell_type_results.txt)
        awk -v j=$j -v s=$set 'NR!=1{print $0"\t"j"\t"s}' $i
    done
done > res.summary.tsv


