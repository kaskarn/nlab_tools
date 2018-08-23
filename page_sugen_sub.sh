#!/bin/bash

##Antoine Baldassari baldassa@email.unc.edu
##call script with no argument to summon help

sugen="/nas02/depts/epi/CVDGeneNas/sol/taor/SUGEN_github/SUGEN/SUGEN"
gwasgen="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes"
megagen="/proj/cvdgene/calico/PAGEII_genotypes/"
ab_pheno="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/phenotypes/"
## id vars

if [ $# -lt 2 ]
then
  echo -e "\n\n\tGWAS utility for PAGE studies in Genetic_DATA_Center/calico\n\t(Antoine Baldassari -- baldassa@email.unc.edu)"
  echo -e "\n\tUSAGE"
  echo -e "\n\tbash page_sugen_sub.sh -s studyname -t traitname ..."
  echo -e "\n\tResults for each study VCF file saved into ./output/study/trait/vcfname.wald.out"
  echo -e "\n\n\tREQUIRED ARGUMENTS"
  echo -e "\n\n\t-s | --study\n\n\tStudy name, or space-delimited list passed in quotes.\n\tMust be one of:\n\twhims, garnet, gecco, hipfx, mopmap, aric_afr, aric_eur, share, lls, mega_all, mega_afr, mega_his"
  echo -e "\tExamples:\t-s mega_all\n\t\t\t-s aric_eur\n\t\t\t--study \"gecco garnet whims mopmap hipfx\""
  echo -e "\n\n\t-t | --trait\n\n\tTrait name, or space-delimited list passed in quotes."
  echo -e "\tExamples:\t-t qt_interval\n\t\t\t--trait qrs\n\t\t\t-t \"pwav pr_seg tp\""
  echo -e "\n\n\tOPTIONAL ARGUMENTS"
  echo -e "\n\n\t-o | --options\n\n\tOptions passed to SUGEN, enclosed in quotes\n\n\tExample: -o \"--extract-file regcheck.txt\" will limit to SNPs in regcheck.txt"
  echo -e "\tWarning: -d option must be specified when using -o\n\tSUGEN documentation available at https://github.com/dragontaoran/SUGEN"
  echo -e "\n\n\t-d | --outdir\n\n\tDirectory where results will be saved.\n\tDefault is ./output.\n\tOption is required when arguments --formula and --options are used"
  echo -e "\n\n\t-c | --covariates\n\n\tUser-specified covariates\n\tWill overwrite study default model. specified as: cov1+cov2+cov3"
  echo -e "\tWarning: -d option must be specified when using -f"
  echo -e "\n\n\t--singlethread\n\n\tForce single-process execution of analyses.\n\n\tUseful for very short runs on small genomic sections."
  echo -e "\n\n\t--aggreg\n\n\tAgregate results into a single file for all chromosomes, using the study-appropriate python script."
  echo -e "\n\n\t--chunksize\n\n\tSubmits separate analyses for chunks of N variants, to speed up runs. \n\tRecommended in 100,000 for MEGA; 200000 for others"
  echo -e "\n\n"
  exit
fi

time="72:0:0"
while [[ $# -gt 0 ]]
do
case $1 in
	-s|--study)
	study="$2"
	shift
	;;
	-t|--trait)
	trait="$2"
	shift
	;;
	-o|--options)
	opts="$2"
	shift
	;;
	-d|--outdir)
	outdir="$2"
	shift
	;;
	-c|--covariates)
	covar="$2"
	shift
	;;
	--aggreg)
	aggreg=true
	;;
	--singlethread)
	single=true
  ;;
	--rightnow)
	single=true
	gogogo=true
	;;
  --chunksize)
  chunksize=$2
  shift
  ;;
  --norun)
  norun=true
  ;;
  --time)
  time=$2
  shift
  ;;
  *)
  echo "Unknown option $1"
  exit
	;;
esac

shift
done
[[ $single == "true" ]] || back="\""

[[ -z $trait ]] && err="ERROR: Need a list of traits. \n$err"
[[ ! -z $opts || ! -z $covar ]] && [[ -z $outdir ]] && err="ERROR: Must specify -d when using -f or -o. \n$err"
[[ -z $outdir ]] && outdir="output"

for s in $study; do
  id_col="analysis_id"
  fid_col="analysis_fid"
  case $s in
    chani)
    phen="/proj/epi/CVDGeneNas/chani/RBC_Traits/GWAS_MEGA_RBC/phenotype/MEGA_RBCtraits_phenotype.txt"
    gen="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/MEGA/imputed_data/vcf"
    d_form=" "
    ;;
    aric_afr)
    phen="$ab_pheno/all/ecg_whi_aric_sol_all.txt"
    gen="$gwasgen/GWAS/ARIC_AA_gwas_frz3"
    d_form="gender2+age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+aric_site_j+rr_d"
    ;;
    aric_eur)
    phen="$ab_pheno/all/ecg_whi_aric_sol_all.txt"
    gen="$gwasgen/GWAS/ARIC_EA_gwas_frz3"
    d_form="gender2+age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+aric_site_m+aric_site_w+rr_d"
    ;;
    garnet|whims|mopmap|gecco|hipfx|lls|share)
    # phen="$ab_pheno/all/ecg_whi_aric_sol_all.txt"
    phen="$ab_pheno/all/tmpwhiids_ecg_whi_aric_sol_all.txt"
    gen=`ls -d $gwasgen/GWAS/*${s}*`
    # d_form="age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+whi_regnum_2+whi_regnum_3+whi_regnum_4+rr_d"
    d_form="age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+whi_regnum_2+whi_regnum_3+whi_regnum_4+rr_d"
    id_col="remap_id_fix"
    id_col="remap_id_fix"
    ;;
    mega_all)
    phen="$ab_pheno/all/ecg_whi_aric_sol_all.txt"
    gen="$megagen/MEGA/imputed_data/vcf"
    d_form="age+gender2+eth_afr+rr_d+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+sol_center_s+sol_center_b+sol_center_c+sol_center_m+whi_regnum_2+whi_regnum_3+whi_regnum_4"
    modopts="--hetero-variance eth_afr"
    ;;
    mega_afr)
    phen="$ab_pheno/all/ecg_whi_aric_sol_afr.txt"
    gen="$megagen/MEGA/imputed_data/vcf"
    d_form="$age+rr_d+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+whi_regnum_2+whi_regnum_3+whi_regnum_4"
    ;;
    mega_his)
    phen="$ab_pheno/all/ecg_whi_aric_sol_his.txt"
    gen="$megagen/MEGA/imputed_data/vcf"
    d_form="age+gender2+rr_d+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+sol_center_s+sol_center_b+sol_center_c+sol_center_m+whi_regnum_2+whi_regnum_3+whi_regnum_4"
    ;;
    *)
    err="ERROR: Study name $s not recognized. \n$err"
    continue
    ;;
  esac
  [[ $s == "chani" ]] && s="mega"
  pyscript="/proj/epi/CVDGeneNas/antoine/bin/rtao/sugen_out_filter_gwas.py"
  [[ $s =~ "mega" ]] && pyscript="/proj/epi/CVDGeneNas/antoine/bin/rtao/sugen_out_filter_mega.py"
  
  #check traits exist
  traitlist=`echo $trait | tr " " "\n" | sort`
  traitfound=`grep -om 1 -wFf <(echo $traitlist | tr " " "\n") $phen | sort`
  diff=`comm -23 <(echo $traitlist | tr " " "\n") <(echo $traitfound | tr " " "\n")`
  [[ "$traitlist" == "$traitfound" ]] || err="ERROR: traits $diff  not found in $s. \n$err"


  form=$covar
  [[ -z $form ]] && form=$d_form
  for t in $trait; do
    dir_out="`pwd`/$outdir/$s/$t"
    [[ -d $dir_out ]] || mkdir -p $dir_out
    file_submit="$dir_out/submit.out"
    echo "#!/bin/bash" > $file_submit
    
    for f in $(find $gen -wholename "*vcf.gz" | grep -v multi_allele | grep -v remapped); do
      fn=`basename $f`
      if [[ -z $chunksize ]]; then
        callstr="$sugen --pheno $phen --id-col $id_col --family-col $fid_col --vcf $f --formula $t=$form --unweighted --out-prefix $dir_out/$t.$fn  $modopts $opts --dosage; read linen fnow <<< \`wc -l $dir_out/$t.$fn.wald.out\`; [[ \${linen} -eq 1 ]] && rm $dir_out/$t.$fn.wald.out" 
        if [ "$single" == "true" ]; then
          echo $callstr >> $dir_out/submit.out
        else
          echo $callstr | sed "s@^@ids=\$ids:\$(sbatch --mem=5GB -t $time -o $dir_out/$t.$fn.slurmout --wrap=\'@g" | sed "s/\$/\')/g" >> $dir_out/submit.out
        fi
      else
        chr=`echo $fn | sed s/.*chr// | sed 's/\..*//' | sed s/.*chr// | sed 's/_.*//'`
        chrnorm=`expr $chr + 0`
        inf="/proj/epi/CVDGeneNas/antoine/page_snplists/pos_$fn"
        while read start end; do 
          waldout=$dir_out/${start}_${end}_$t.$fn.wald.out
          callstr="$sugen --extract-chr $chrnorm --extract-range $start-$end --pheno $phen --id-col $id_col --family-col $fid_col --vcf $f --formula $t=$form --unweighted --out-prefix $dir_out/${start}_${end}_$t.$fn  $modopts $opts --dosage; read linen fnow <<< \`wc -l $waldout\`; [[ \${linen} -eq 1 ]] && rm $waldout" 
          if [ "$single" == "true" ]; then
            echo $callstr >> $dir_out/submit.out
          else
            echo $callstr | sed "s@^@ids=\$ids:\$(sbatch --mem=5GB -t $time -o $dir_out/${start}_${end}_$t.$fn.slurmout --wrap=\'@g" | sed "s/\$/\')/g" >> $dir_out/submit.out
          fi
        done <<< "$( paste <(sed -n "1~${chunksize}p" $inf)  <(cat <(sed -n "0~${chunksize}p" $inf)  <(tail -1 $inf)) )"
      fi
    done
    if [ "$aggreg" == "true" ]; then
      [[ "$single" == "true" ]] && echo "python $pyscript $dir_out $dir_out/${s}_${t}_allchr.txt" >> $dir_out/submit.out
      [[ -z $single ]] && echo "sbatch -t 2:0:0 --dependency=after\`echo \$ids | sed s/[[:space:]]//g | sed s/[[:alpha:]]//g\` -o $dir_out/${s}_${t}_slurmout.slurm --wrap='python $pyscript $dir_out $dir_out/${s}_${t}_allchr.txt'" >> $dir_out/submit.out  
    fi
    chmod 700 $file_submit
  done
done

[[ -z $err ]] || echo -e $err
[[ -z $err ]] || exit
[[ $norun == "true" ]] && echo "All looks good; Remove --norun option to carry out the run"
[[ -z $norun ]] || exit

for s in $study; do
  for t in $trait; do
		echo "Submitting scripts for study $s, trait $t"
    dir_out="`pwd`/$outdir/$s/$t"
    if [ "$single" == "true" ]; then
      [[ -z $gogogo ]] && ids=`sbatch --mem=6GB -o $dir_out/$s.$t.allchr.slurmout -t 24:0:0 $dir_out/submit.out`
      [[ -z $gogogo ]] || bash $dir_out/submit.out
			allids="$allids `echo $ids | sed s/[[:space:]]//g | sed s/[[:alpha:]]//g`"
    else
      allids="$allids `bash $dir_out/submit.out`"
    fi
  done
done
echo $allids
# if [[ "$s" =~ "mega" ]]; then 
        # annodir="/proj/epi/CVDGeneNas/sol/taor/mega_imputed/dosage_files/vcf_variant_list_info0.4/list"
        # anno="$annodir/$fn.list"
        # ah="echo -e position\tinfo\trs_id"
        # k=1
      # else
        # annodir="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS/short_infolists"
        # anno="$annodir/`echo $fn | sed s/.vcf.gz/.info.txt/`"
        # ah="head -1 $anno | tr \" \" \"\t\""
        # k=2
      # fi
      # outf=$dir_out/$t.$fn.wald.out
      # outa=$dir_out/$t.$fn.anno.wald.out	  
# echo "sbatch -t 2:0:0 --dependency=after\`echo \$ids | sed s/[[:space:]]//g | sed s/[[:alpha:]]//g\` -o $dir_out/${s}_${t}_slurmout.slurm --wrap=' allf=($dir_out/*.wald.out); head -1 \$allf > $dir_out/${s}_${t}_allchr.txt; tail -qn +2 \${allf[@]} >> $dir_out/${s}_${t}_allchr.txt '" >> $dir_out/submit.out
    
# if [ "$single" == "true" ]; then
  # [[ -z $aggreg ]] ||  echo "python $pyscript $dir_out $dir_out/${s}_${t}_allchr.txt" >> $dir_out/submit.out
# else
  # sed -i "1 ! s/^/ids=\$ids:\$\($cmd/g" $dir_out/submit.out
  # sed -i "1 ! s/\$/\'\)/g" $dir_out/submit.out
  # [[ -z $aggreg ]] || echo "sbatch -t 2:0:0 --dependency=after\`echo \$ids | sed s/[[:space:]]//g | sed s/[[:alpha:]]//g\` -o $dir_out/${s}_${t}_slurmout.slurm --wrap='python $pyscript $dir_out $dir_out/${s}_${t}_allchr.txt'" >> $dir_out/submit.out  
  # echo 'echo job IDs: $ids' >> $dir_out/submit.out
# fi
