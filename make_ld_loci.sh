#!/bin/bash

# aric_afr="${page}/ARIC_AA_gwas_frz3/ARIC_AA_gwas_frz3.chr${chr}.vcf.gz"
# aric_eur="${page}/ARIC_EA_gwas_frz3/ARIC_EA_gwas_frz3.chr${chr}.vcf.gz"
# mopmap="${page}/mopmap-updatedEA/mopmap-updatedEA.chr${chr}.vcf.gz"

page="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS"
bdir="${ab_home}/bin/bcftools"
outld="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/LD/ld_final/ld_results"

chr=$1
pop=$2
[[ -z "$pop" ]] && pop=all

chr2=${chr}
[[ ${chr} -lt 10 ]] && chr2="0${chr}"
aric_afr="/proj/epi/Genetic_Data_Center/aric/aric_1kgp3/AA/chr${chr}.dose.vcf.gz"
aric_eur="/proj/epi/Genetic_Data_Center/aric/aric_1kgp3/EA/ARIC_EA_1kgp3_chr${chr}.dose.vcf.gz"
garnet="${page}/garnet-updated/garnet-updated.chr${chr}.vcf.gz"
gecco="${page}/gecco-updatedEA/gecco-updatedEA.chr${chr}.vcf.gz"
mega="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/MEGA/imputed_data/vcf/chr${chr2}.vcf.gz"
whims="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS/whims-updatedEA/whims-updatedEA.chr${chr}.vcf.gz"
mopmap="/proj/epi/CVDGeneNas/antoine/pageii_fix/mopmap/mopmap-updatedEA.chr${chr}.vcf.gz"

tfiles=""
[[ $pop == "all" ]] && studs=("aric_eur" "aric_afr" "mega" "garnet" "gecco" "whims" "mopmap")
[[ $pop == "afr" ]] && studs=("aric_afr" "mega")
[[ $pop == "eur" ]] && studs=("aric_eur" "garnet" "gecco" "whims" "mopmap")
[[ $pop == "his" ]] && studs=("mega")
[[ -z "$studs" ]] && exit

for s in ${studs[@]}
do
	echo processing ${!s}	
	tf=`mktemp`
	tfiles="${tfiles} ${tf}"
	$bdir annotate -h header_gt_gp.hdr -Ou -R infiles/${s}_ksig_pos_1e6.txt -x ID -I +'%CHROM:%POS:%REF:%ALT' ${!s} | $bdir sort -Ou | $bdir norm -m -any -Ob > $tf
	$bdir index $tf
done
# $bdir index $tfiles
gen_tmp=`mktemp`
ld_tmp=`mktemp`
ids_tmp=`mktemp`
sort infiles/pageii_ecg_${pop}_ids.txt > $ids_tmp

fn=`wc -w <<< "$tfiles"`
[[ $fn -gt 1 ]] && $bdir merge $tfiles -m none -Ou | $bdir view -S <(cut -f1 $ids_tmp) --force-samples -Ou | $bdir convert -g $gen_tmp --tag GP --chrom --vcf-ids
[[ $fn -eq 1 ]] && $bdir view -S infiles/pageii_ecg_${pop}_ids.txt --force-samples -Ou $tfiles | $bdir convert -g $gen_tmp --tag GP --chrom --vcf-ids

ids_tmp2=`mktemp`
join -j1 $ids_tmp <(cut -f1 -d' ' $gen_tmp.samples | tail -n +3) > $ids_tmp2

samp_fix=`mktemp`
head -2 $gen_tmp.samples > $samp_fix 
paste <(cut -f2 -d' ' $ids_tmp2) <(cut -f1 -d' ' $ids_tmp2) <(cut -f1-2 --complement -d' ' $gen_tmp.samples | tail -n +3) >> $samp_fix

plink --gen $gen_tmp.gen.gz --allow-extra-chr --sample $samp_fix --r2 dprime --ld-window-r2 0 --ld-window 999999 --out $ld_tmp \
--clump infiles/${pop}_ksig_vp_1e6.txt --clump-r2 0.05 --clump-kb 1000 --clump-p1 0.99 --clump-p2 1 --clump-field min

cat $ld_tmp.ld | sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | cut -f 3,6-8 > $outld/${pop}_ld_chr${chr}.ld
cat $ld_tmp.clumped | sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' > $outld/${pop}_ld_chr${chr}.clumped

rm $tfiles $gen_tmp.* $ld_tmp* $ids_tmp $ids_tmp2 $samp_fix

