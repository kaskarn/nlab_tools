#!/bin/bash
cd /proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/qt_pr/metal
xdir="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/qt_pr/gwas"
# studs=(`ls $xdir`)
studs=("whims" "mopmap" "hipfx" "gecco" "garnet" "mega_all" "aric_afr" "aric_eur")
all="whims mopmap hipfx gecco garnet mega_all aric_afr aric_eur"
afr="aric_afr mega_afr"
his="mega_his"
eur="aric_eur whims mopmap hipfx gecco garnet"
whi="whims mopmap hipfx gecco garnet"
aric="aric_eur aric_afr"

mesadir="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/mesa2/gwas"
mesagrp="CAU AFA CHN HIS"

for f in $mesadir/

# gps="all afr his eur whi aric"
# gps="afr his eur whi aric all"

# traits=(`ls $xdir/$studs/`)
# traits=("qt_interval" "pr_interval")
# traits="pwav pr_seg qrs st twav tp qt_interval pr_interval"
traits=("pwav" "pr_seg" "qrs" "st" "twav" "tp") # "qt_interval" "pr_interval")
mesatraits=("Pwave" "PRseg" "QRS" "STseg" "Twave" "TPseg")
addmesa=true
gps="all"

ids=""
for gp in $gps; do
  [[ -d "out_$gp" ]] || mkdir out_$gp
  for i in ${!traits[@]}; do
    t=${traits[$i]}
    t2=${mesatraits[$i]}
    sub="out_$gp/${gp}_metal_tmpsub_$t.metal"
    echo 'SCHEME STDERR
MAXWARNINGS 100
VERBOSE OFF
WEIGHTLABEL METAL_N

AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER ALT_AF > 0.01
ADDFILTER ALT_AF < 0.99
ADDFILTER IMPUTE2_INFO > 0.3
ADDFILTER ESS > 30

SEPARATOR TAB
MARKERLABEL CPTID
PVALUELABEL PVALUE
EFFECTLABEL BETA
STDERRLABEL SE
FREQLABEL ALT_AF
ALLELELABELS ALT REF
' > $sub
    for s in ${!gp}; do
      echo "PROCESS $xdir/$s/$t/${s}_${t}_allchr.txt" >> $sub
    done
    echo -e "\n" >> $sub
    
    if [ "$addmesa" == "true" ]; then
      echo '
REMOVEFILTERS

ADDFILTER eaf > 0.01
ADDFILTER eaf < 0.99
ADDFILTER ess > 30
ADDFILTER imputation > 0.3
SEPARATOR COMMA
MARKERLABEL cptid
PVALUELABEL pval
EFFECTLABEL beta
STDERRLABEL se
FREQLABEL eaf
ALLELELABELS other_allele effect_allele
' >> $sub
      for s in $mesagrp; do
        echo "PROCESS $mesadir/mesa_ecg_${t2}_${s}.txt" >> $sub
      done
      echo "OUTFILE out_allmesa/all_metal_${t} .txt" >> $sub
      mout="out_allmesa/all_metal_${t}1.txt"
    else
      echo "OUTFILE out_${gp}/${gp}_metal_${t} .txt" >> $sub
      mout="out_$gp/${gp}_metal_${t}1.txt"
    fi
    
    echo "ANALYZE HETEROGENEITY" >> $sub
    echo "QUIT" >> $sub
    
    [[ -f $mout ]] && mv $mout oldres
    # msg=`sbatch --mem=10GB --time=3:0:0 -o out_$gp/${gp}_${t}_metal.out --wrap="metal $sub ;sed -i 's/P-value/Pvalue/' $mout ;sed -i 's/MarkerName/Chr\tPos/' $mout ;sed -i 's/:/\t/' out_$gp/${gp}_metal_*.txt"`
    msg=`sbatch --mem=10GB --time=3:0:0 -o out_$gp/${gp}_${t}_metal.out --wrap="metal $sub ;sed -i 's/P-value/Pvalue/ ;s/MarkerName/Chr\tPos\tRef\tAlt/ ;s/:/\t/1 ;s/:/\t/1 ;s/:/\t/1' $mout"`
    idnow=`echo $msg | tr ' ' '\n' | tail -1`
    ids="${ids}:${idnow}"
  done
done

#run after previous step is complete
# for gp in $gps;do 
  # echo $gp
  # sed -i 's/P-value/Pvalue/' out_$gp/${gp}_metal_*.txt
  # sed -i 's/MarkerName/Chr\tPos/' out_$gp/${gp}_metal_*.txt
  # sed -i 's/:/\t/' out_$gp/${gp}_metal_*.txt
# done 

annoman="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/graphs/annofiles/annosig_all_min.txt"
annoqq="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/graphs/annofiles/ecg_allknown_forqq.txt"
gpath="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/qt_pr/graphs"
traits="pwav pr_seg qrs st twav tp qt_interval pr_interval"
for t in $traits; do
  f=out_all/all_metal_${t}1.txt
  outf=`echo $f | sed 's/.txt//'`
  nameout=`basename $f | sed s/_metal// | sed s/1\.txt//`
  sbatch -o "$gpath/logs/man_$nameout.out" --wrap="bash $esplot man -v Pvalue -c Chr -p Pos -n $nameout -o $gpath -f $f -a $annoman"
  sbatch -o "$gpath/logs/qq_$nameout.out" --wrap="bash $esplot qq -v Pvalue -c Chr -p Pos -n $nameout -o $gpath -f $f -r $annoqq"
done


# metout=(all_metal_*1.txt)
# annoman="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/graphs/annofiles/annosig_all_min.txt"
# annoqq="/proj/epi/CVDGeneNas/antoine/ECG_GWAS_FULL/gwas/graphs/annofiles/ecg_allknown_forqq.txt"
# for f in ${metout[@]}; do
  # echo $f
  # outf=`echo $f | sed 's/.txt//'`
  # cat <(paste <(echo -e "Chr\tPos") <(head -1 $f | sed s/P-value/Pvalue/)) <(paste <(tail -n +2 $f | cut -f1 | tr ':' '\t') <(tail -n +2 $f )) > ${outf}_fix.txt
  # nameout=`echo $f | sed s/_metal// | sed s/1\.txt//`
  # sbatch -o "man_$nameout.out" --wrap="bash $esplot man -v Pvalue -c Chr -p Pos -n plots/$nameout -f `pwd`/$f -a $annoman"
  # sbatch -o "qq_$nameout.out" --wrap="bash $esplot qq -v Pvalue -c Chr -p Pos -n plots/$nameout -f `pwd`/$f -r $annoqq"
# done

##aggregate into one
echo '
#!/bin/bash
traits=("pwav" "pr_seg" "qrs" "st" "twav" "tp" "qt_interval" "pr_interval")
# cut -f1-7 metal_pwav1_fix.txt > metal_all.txt
hcmd=`echo "paste <(head -1 out_all/all_metal_pwav1.txt | cut -f1-7)"`
agcmd=`echo "paste <(tail -n +2 out_all/all_metal_pwav1.txt | cut -f1-7)"`
for t in ${traits[@]}; do
	hcmd=`echo "$hcmd <(echo -e \"${t}_b\t${t}_p\")"`
	agcmd=`echo "$agcmd <(tail -n +2 out_all/all_metal_${t}1.txt | cut -f9,11)"`
done
echo $hcmd
echo $agcmd
eval $hcmd > out_all/metal_alltraits.txt
eval $agcmd >> out_all/metal_alltraits.txt
' > aggscript.sh
bash aggscript.sh
head metal_all.txt

 


# sbatch -o aggreg_metal.out --dependency=after$ids aggscript.sh


# for t in ${traits[@]};do 
# for s in ${studs[@]}; do
# if [[ "$s" =~ "mega" ]]; then
# sbatch --wrap="python /proj/epi/CVDGeneNas/antoine/bin/rtao/sugen_out_filter_mega.py $xdir/$s/$t $xdir/$s/$t/${s}_${t}_allchr.txt"
# else
# sbatch --wrap="python /proj/epi/CVDGeneNas/antoine/bin/rtao/sugen_out_filter_gwas.py $xdir/$s/$t $xdir/$s/$t/${s}_${t}_allchr.txt"
# fi
# done
# done