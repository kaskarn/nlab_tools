#!/bin/bash

homedir=`pwd`
script=$(readlink -f "$0")
rundir=$(dirname "$script")

cmd=$1
shift

case $cmd in
  setup)
  wrap="bash $rundir/setup.sh $@"
  ;;
  env_tools)
  wrap="bash $rundir/envtools.sh $@"
  ;;
  aspu)
  wrap="bash $rundir/aspu_julia.sh $@"
  ;;
  easyplot)
  bash $rundir/ES_mkplot.sh $@
  ;;
  easyld)
  wrap="bash $rundir/make_ld_loci.sh $@"
  ;;
  easylz)
  wrap="bash $rundir/make_lzplots.sh $@"
  ;;
  ncdf_whi_gwas)
  wrap="Rscript $rundir/ncgwas_script.R $@"
  ;;
  ftranspose)
  wrap="$rundir/ftranspose $@"
  ;;
  easy_sugen)
  wrap="bash $rundir/page_sugen_sub.sh $@"
  ;;
  easy_vcf_lm)
  wrap="bash $rundir/vcf_gwas.sh $@"
  ;;
  easy_sugen_filter)
  wrap="bash $rundir/sugen_filter.sh $@"
  ;;
  easy_metal)
  wrap="bash $rundir/metal_gwas_ecgsubmit.sh $@"
  ;;
  selectcols)
  bash $rundir/selectcols.sh $@
  ;;
  compute_covar)
  julia $rundir/compute_covar.jl $@
  ;;
esac

exit
