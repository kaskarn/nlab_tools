#!/bin/bash
#Antoine Baldassari baldassa@email.unc.edu
#to call: bash thisscript study trait catalogue_trait_name. Loop for qt:
# studies=("whims" "garnet" "gecco" "mopmap" "hipfx" "share" "lls") #"aric_eur" "aric_afr" 
# for s in ${studies[@]}
# do
# sbatch --wrap="bash pageII_gwascata_check.sh $s qt"
# done

page="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS"
aric_afr="${page}/ARIC_AA_gwas_frz3"
aric_eur="${page}/ARIC_EA_gwas_frz3"
garnet="${page}/garnet-updated"
mopmap="${ab_home}/pageii_fix/mopmap"
gecco="${page}/gecco-updatedEA"
whims="${page}/whims-updatedEA"
hipfx="${page}/hipfx-updated"
share="${page}/WHI_share_Affy6.0-2015-03-05"
lls="${page}/lls-updated"

allphen="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/phenotypes/all/experimental_whi.txt"
# allphen="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/phenotypes/oldphen/WHI_1kgp3_WRONG_IDS_IN_THERE/WRONG_ID_DO_NOT_USE/ecg_whi_all_v4.txt"
sugen="/nas02/depts/epi/CVDGeneNas/sol/taor/SUGEN_github/SUGEN/SUGEN"

# idvar="analysis_id_tmpfix"
idvar="analysis_id_new2"
# idvar="analysis_id"
cov=("age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10")
# +rr_d+regnum_2+regnum_3+regnum_4+
study=$1
trait=$2
trait_name=$3
[[ -z "$trait_name" ]] && trait_name=$trait

#make tempfiles
catalog=`mktemp`; topcata=`mktemp`; topcata2=`mktemp`
chrpos=`mktemp`; tmpout=`mktemp`
tmpphen=`mktemp`; tmpvarterms=`mktemp`

#get latest gwas catalog
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -P 3306 hg19 -e "select * from gwasCatalog" > $catalog
R -e 'library(data.table); dt=fwrite(fread('"\"$catalog\""', select = c("name", "chrom", "chromEnd", "pValue", "trait"), key = "pValue")[tolower(trait) %like% '"\"$trait_name\""'][, .(rsid=unique(name),trait=.SD[which.min(pValue),trait],litp=min(pValue)), .(chrom, chromEnd)], '"\"$topcata\""', sep="\t")'
sed -e 's/chr//g' $topcata | sed -ne '/X/!p' | sed -ne '/_/!p' | sed -e 's/^om\tom/chrom\tchrom/' > $topcata2
cut -f1,2 $topcata | tail -n +2 | sed -e 's/chr//g' | sed -ne '/_/!p' | tr '\t' ':' > $chrpos

#make temporary phenotype file
sed 's/+/\n/g' <(echo $idvar+analysis_fid+$trait+$cov) > $tmpvarterms
# grep --file=$tmpvarterms -nwh <(head -1 $allphen | sed -e 's/\t/\n/g') | cut -f1 -d':' | sed 's/$/,/g' | tr -d '\n'
coln=`grep --file=$tmpvarterms -nwh <(head -1 $allphen | sed -e 's/\t/\n/g') | cut -f1 -d':' | sed ':a;N;$!ba;s/\n/,/g'`
cut -f$coln $allphen | sed -n '/\t\t/!p' | sed -n '/\t$/!p' | sed -n '/^\t/!p' > $tmpphen

formula="${trait}=$cov"
fout=${trait}_${study}_check.txt

#run checks in each chromosome and merge
bool=true
files_vcf=(${!study}/*chr*vcf.gz)
for f in ${files_vcf[@]}
do
	echo processing `basename $f`
	$sugen --pheno $tmpphen --id-col $idvar --extract-file $chrpos --family-col analysis_fid --vcf $f --formula $formula --unweighted --out-prefix $tmpout --dosage
	mv $tmpout.log check_lastlog_${study}_${trait}.log
	[[ $bool = true ]] && cat $tmpout.wald.out > $fout
	[[ $bool = false ]] && tail -n +2 $tmpout.wald.out >> $fout
	bool=false
	echo -e "done\n"
done	

fullout="full_${fout}"
R -e 'library(data.table);a=fread('"\"$fout\""', key=c("CHROM","POS"));b=fread('"\"$topcata2\""',key=c("chrom","chromEnd"));fwrite(setkey(a[b], litp),'"\"$fullout\""', sep = "\t")'
rm $catalog $topcata $chrpos $tmpout* $tmpvarterms $tmpphen $fout
