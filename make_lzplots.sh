#!/bin/bash

# aric_afr="${page}/ARIC_AA_gwas_frz3/ARIC_AA_gwas_frz3.chr${chr}.vcf.gz"
# aric_eur="${page}/ARIC_EA_gwas_frz3/ARIC_EA_gwas_frz3.chr${chr}.vcf.gz"
# mopmap="${page}/mopmap-updatedEA/mopmap-updatedEA.chr${chr}.vcf.gz"

page="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS"
bdir="${ab_home}/bin/bcftools"


cd /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/LD/ld_final
snp=$1
cptid=$2
chr=$3
pos=$4
pop=$5
trait=$6
outdir=$7

[[ -z "$5" ]] && pop="all"

start=`expr $pos - 500000`
end=`expr $pos + 500000`
[[ start -lt 0 ]] && start=1

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
	$bdir annotate -h header_gt_gp.hdr -Ou -r $chr:$start-$end -x ID -I +'%CHROM:%POS:%REF:%ALT' -e 'AVG(DS) < 0.01' ${!s} | $bdir norm -m -any -Ob > $tf
	$bdir index $tf
done

# $bdir index $tfiles
gen_tmp=`mktemp`
vcf_tmp=`mktemp`
ld_tmp=`mktemp`
ids_tmp=`mktemp`
sort infiles/pageii_ecg_${pop}_ids.txt > $ids_tmp

fn=`wc -w <<< "$tfiles"`
[[ $fn -gt 1 ]] && $bdir merge $tfiles -m none -Ou | $bdir view -S <(cut -f1 $ids_tmp) --force-samples -Oz > $vcf_tmp
[[ $fn -eq 1 ]] && $bdir view -S <(cut -f1 $ids_tmp) --force-samples -Oz $tfiles > $vcf_tmp
$bdir convert -g $gen_tmp --tag GP --chrom --vcf-ids $vcf_tmp

ids_tmp2=`mktemp`
join -j1 $ids_tmp <(cut -f1 -d' ' $gen_tmp.samples | tail -n +3) > $ids_tmp2

samp_fix=`mktemp`
head -2 $gen_tmp.samples > $samp_fix 
paste <(cut -f2 -d' ' $ids_tmp2) <(cut -f1 -d' ' $ids_tmp2) <(cut -f1-2 --complement -d' ' $gen_tmp.samples | tail -n +3) >> $samp_fix

plink --gen $gen_tmp.gen.gz --allow-extra-chr --sample $samp_fix --r2 dprime --ld-window-r2 0 --ld-window 999999 --ld-snp $cptid --out $ld_tmp 

ftmp1=`mktemp`
ftmp2=`mktemp`

cat $ld_tmp.ld | sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | cut -f 3,6-8 > $ftmp1

join --nocheck-order -1 2 -2 1 <(tail -n +2 $ftmp1 | sort -k 2 ) <(cut -f1-2 infiles/${pop}_ksig_vp.txt | grep -v $'\t'$ | tail -n +2 | sort) > $ftmp2
read lines name <<< `wc -l $ftmp2`

echo -e "snp1\tsnp2\tdprime\trsquare" > topsnp_ld/${pop}_${snp}_ld.txt
paste <(cut -f5 -d' ' $ftmp2) <(yes $snp | head -${lines}) <(cut -f4 -d' ' $ftmp2) <(cut -f3 -d' ' $ftmp2) >> topsnp_ld/${pop}_${snp}_ld.txt

#cd $lz_dir
echo "Making locuszoom"
head /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/LD/ld_final/infiles/${pop}_ksig_vp.txt
head /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/LD/ld_final/topsnp_ld/${pop}_${snp}_ld.txt

[[ ! -d $outdir/$snp ]] && mkdir -pv $outdir/$snp
locuszoom --metal /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/LD/ld_final/infiles/${pop}_ksig_vp.txt --markercol rsid --pvalcol $trait --refsnp $snp \
--ld /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/LD/ld_final/topsnp_ld/${pop}_${snp}_ld.txt --flank 500kb \
--build hg19 --cache None --prefix $outdir/$snp/${trait}_${pop}_lz --plotonly --no-date

# $bdir index -t $vcf_tmp
# locuszoom --metal /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/LD/ld_final/infiles/${pop}_ksig_vp.txt --markercol rsid --pvalcol $trait --refsnp $snp \
# --ld-vcf $vcf_tmp --flank 500kb \
# --build hg19 --cache None --prefix $outdir/$snp/${trait}_${pop}_lz --plotonly --no-date

rm $tfiles $gen_tmp.* $ld_tmp.* $ftmp1 $ftmp2 $vcf_tmp

#--gwas-cat whole-cat_significant-only 
