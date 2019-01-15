#!/bin/bash
#esplot=/proj/epi/CVDGeneNas/antoine/bin/easystrata/ES_mkplot.sh

type=$1
shift

tmecf=`mktemp`
tmR=`mktemp`

echo "formals(read.table)\$comment.char <- ''
assignInNamespace('read.table', read.table, 'utils')" > $tmR
[[ $type == "qq" ]]	&& echo "library(EasyQC);EasyQC(\"$tmecf\")" >> $tmR
[[ $type == "man" ]] && echo "library(EasyStrata);EasyStrata(\"$tmecf\")" >> $tmR
[[ $type == "miami" ]] && echo "library(EasyStrata);EasyStrata(\"$tmecf\")" >> $tmR

break="1"
chr="chr"
pos="pos"
pval="Pvalue"
pskip="0.1"
pline="5E-9"
pathout=`pwd`
name="myplot"
sep="TAB"
autosom="0"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
	-a|--annot)
	annot="$2"
	shift
	;;
	-b|--break)
	pbreak="$2"
	shift
	;;
	-c|--chr)
	chr="$2"
	shift
	;;
	-f|--filein)
	filein="`readlink -f $2`"
	shift
	;;
	-n|--name)
	name="$2"
	shift
	;;
	--nolambda)
	nolambda="0"
	;;
	-o|--pathout)
	pathout="$2"
	shift
	;;
	--othopts)
	othopts="$2"
	shift
	;;
	-p|--pos)
	pos="$2"
	shift
	;;
  --pline)
  pline="$2"
  shift
  ;;
	-s|--sep)
	sep="$2"
	shift
	;;
	-v|--pval)
	pval="$2"
	shift
	;;
  --pvaldown)
  pvaldown=$2
  shift
  ;;
  --pvalup)
  pvalup=$2
  shift
  ;;
	-x|--prefix)
	prefix="$2"
	shift
	;;
	-r|--remove)
	remove="$2"
	shift
	;;
	-s|--pskip)
	pskip="$2"
	shift
	;;
  --autosom)
  autosom="1"
  ;;
    *)
    ;;
esac

shift
done

echo -e "#\nDEFINE\t--strSeparator $sep" > $tmecf
[[ "$type" != "miami" ]] && echo -e "\t--acolIn ${pval};${pos};${chr}" >> $tmecf
[[ "$type" == "miami" ]] && echo -e "\t--acolIn ${pvalup};${pvaldown};${pos};${chr}" >> $tmecf
echo -e "\t--pathOut $pathout\n" >> $tmecf
echo -e "EASYIN\t--fileIn $filein" >> $tmecf
echo -e "\t--fileInShortName $type\n" >> $tmecf



if [ $type == "man" ] ; then
	echo -e "START EASYSTRATA\n" >> $tmecf
  if [ $autosom == "1" ] ; then
    echo -e "FILTER\t--rcdFilter ${chr}<=22\n\t--strFilterName tmpname_esplot\n" >> $tmecf
  fi
	echo -e "MHPLOT\t--colMHPlot $pval" >> $tmecf
	echo -e "\t--anumAddPvalLine $pline\n\t--blnPlotVerticalLines 0" >> $tmecf
	[[ ! -z "$annot" ]] && echo -e "\t--fileAnnot ${annot}\n\t--numAnnotPvalLim 5e-8" >> $tmecf
fi
if [ $type == "qq" ] ; then
	echo -e "START EASYQC\n" >> $tmecf
  if [ $autosom == "1" ] ; then
    echo -e "FILTER\t--rcdFilter ${chr}<=22\n\t--strFilterName tmpname_esplot\n" >> $tmecf
  fi
	echo -e "QQPLOT\t--acolQQPlot $pval" >> $tmecf
	[[ ! -z "$remove" ]] && echo -e "\t--fileRemove $remove\n\t--strRemovedColour blue" >> $tmecf
	[[ -z "$nolambda" ]] && echo -e "\t--blnAddLambdaGC 1" >> $tmecf
fi
if [ $type == "miami" ] ; then
	echo -e "START EASYSTRATA\n" >> $tmecf
	echo -e "MIAMIPLOT\t--colMIAMIPlotUp $pvalup\n\t--colMIAMIPlotDown $pvaldown " >> $tmecf
	echo -e "\t--anumAddPvalLine $pline\n\t--blnPlotVerticalLines 0" >> $tmecf
	[[ ! -z "$annot" ]] && echo -e "\t--fileAnnot ${annot}\n\t--numAnnotPvalLim 5e-8" >> $tmecf
fi
echo -e "\t--colInChr $chr" >> $tmecf
echo -e "\t--colInPos $pos" >> $tmecf
echo -e "\t--strPlotName $name" >> $tmecf
echo -e "\t--numPvalOffset $pskip" >> $tmecf
echo -e "\t--blnYAxisBreak $break" >> $tmecf

[[ ! -z "$othopts" ]] && echo -e "\t$othopts" >> $tmecf
[[ $type == "man" ]] && echo -e "\nSTOP EASYSTRATA" >> $tmecf
[[ $type == "miami" ]] && echo -e "\nSTOP EASYSTRATA" >> $tmecf
[[ $type == "qq" ]] && echo -e "\nSTOP EASYQC" >> $tmecf

cat $tmecf

Rscript $tmR $tmecf
tm=`basename $tmecf`

rm $tmR $tmecf $pathout/$tm*

