#!/bin/bash

# magma
# pd gwas
awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2}' nallsEtAl2019_excluding23andMe_allVariants.tab | awk 'BEGIN{OFS="\t"}{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$9+$10}' > nallsEtAl2019_excluding23andMe_allVariants.revised.tab
cp /projects/ps-renlab/yangli/resource/GWAStraits/hapmap3.snp.bed /projects/ps-renlab2/kaw033/WK_aging/ldsc/
sed '1d' nallsEtAl2019_excluding23andMe_allVariants.revised.tab | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$3+1,$0}' | sort -k1,1 -k2,2n > nallsEtAl2019_excluding23andMe_allVariants.revised.tab.bed

bedtools intersect -a nallsEtAl2019_excluding23andMe_allVariants.revised.tab.bed \
             -b hapmap3.snp.bed \
             -wao |\
awk 'BEGIN{FS=OFS="\t"}{if($16=="."){print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}else{print $19,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}}' > nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid

sed -ie "1i SNPID\tCHR\tPOS\tA1\tA2\tfreq\tb\tse\tp\tN_cases\tN_controls\tN" nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid
awk '{sub("chr", "", $2); print}' nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid > nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid.rmchr
cp nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid.rmchr nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid.rmchr.pfile

/home/kaw033/softwares/magma/magma --annotate window=35,10 --snp-loc /projects/ps-renlab2/kaw033/WK_aging/ldsc/nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid.rmchr \
 --gene-loc ~/softwares/magma/resource/NCBI37.3/NCBI37.3.gene.loc --out nallsEtAl2019_excluding23andMe_allVariants.window

/home/kaw033/softwares/magma/magma --bfile ~/softwares/magma/resource/g1000_eur/g1000_eur \
 --pval /projects/ps-renlab2/kaw033/WK_aging/ldsc/nallsEtAl2019_excluding23andMe_allVariants.revised.tab.snpid.rmchr.pfile ncol=N \
 #use='SNP,P'
 --gene-annot nallsEtAl2019_excluding23andMe_allVariants.window.genes.annot \
 --out nallsEtAl2019_excluding23andMe_allVariants.window


/home/kaw033/softwares/magma/magma --gene-results nallsEtAl2019_excluding23andMe_allVariants.window.genes.raw \
 --gene-covar /projects/ps-renlab2/kaw033/WK_aging/rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.pseudo.subclass.avg.entrezID.txt --out after_integra/nallsEtAl2019_excluding23andMe_allVariants.window.express

# IBD gwas
cat ibd_build37_59957_20161107.txt | cut -f1-6 | awk 'BEGIN{FS="[:_\t]";OFS="\t"}{if(NR==1){print $0}else{print $1,$2,$5,$6,$7,$8,$9}}' | awk 'BEGIN{OFS="\t"}NR==1{print $1,"CHR","BP",$2,$3,$4,$5,$6,$7,"N"}NR>1{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,"59957"}' > ibd_build37_59957_20161107.revised.txt
sed '1d' ibd_build37_59957_20161107.revised.txt | awk 'BEGIN{FS=OFS="\t"}{print "chr"$2,$3,$3+1,$0}' | sort -k1,1 -k2,2n > ibd_build37_59957_20161107.revised.txt.bed

bedtools intersect -a ibd_build37_59957_20161107.revised.txt.bed \
             -b hapmap3.snp.bed \
             -wao |\
awk 'BEGIN{FS=OFS="\t"}{if($13=="."){print $4,$5,$6,$7,$8,$9,$10,$11,$12}else{print $16,$5,$6,$7,$8,$9,$10,$11,$12}}' > ibd_build37_59957_20161107.revised.txt.snpid

sed -ie "1i SNP\tCHR\tBP\tAllele1\tAllele2\tEffect\tStdErr\tP\tN" ibd_build37_59957_20161107.revised.txt.snpid

cd magma
/home/kaw033/softwares/magma/magma --annotate window=35,10 --snp-loc /projects/ps-renlab2/kaw033/WK_aging/ldsc/ibd_build37_59957_20161107.revised.txt.snpid \
 --gene-loc ~/softwares/magma/resource/NCBI37.3/NCBI37.3.gene.loc --out ibd_build37_59957_20161107.window

/home/kaw033/softwares/magma/magma --bfile ~/softwares/magma/resource/g1000_eur/g1000_eur \
 --pval /projects/ps-renlab2/kaw033/WK_aging/ldsc/ibd_build37_59957_20161107.revised.txt.snpid ncol=N \
 #use='SNP,P'
 --gene-annot ibd_build37_59957_20161107.window.genes.annot\
 --out ibd_build37_59957_20161107.window

/home/kaw033/softwares/magma/magma --gene-results ibd_build37_59957_20161107.window.genes.raw \
 --gene-covar /projects/ps-renlab2/kaw033/WK_aging/rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.pseudo.subclass.avg.entrezID.txt --out after_integra/ibd_build37_59957_20161107.window.express


