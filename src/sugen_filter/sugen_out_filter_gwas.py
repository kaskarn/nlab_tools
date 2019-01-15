#!/usr/bin/python

import sys
import glob
import re
import os

dir_res = sys.argv[1]; print "dir_res:", dir_res
osin = os.path.splitext(sys.argv[2])
if len(osin) > 1:
  ext = osin[1]
else:
  ext = ".txt"
file_out = osin[0]+ext
file_filt = osin[0]+"_excluded"+ext
print "file_out:", file_out
print "file_filt:", file_filt
dir_keep_list = "/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/GWAS/short_infolists"

files_input = sorted(glob.glob(dir_res+"/*.wald.out"))
flag_wald = True
if not files_input:
  files_input = sorted(glob.glob(dir_res+"/*.score.snp.out"))
  flag_wald = False
if not files_input:
  print dir_res, "does not contain any SUGEN output files!"
  exit(0)

trait = re.search("/[\w.]+$", dir_res).group(0).lstrip('/')
print "trait:", trait
  
output_file = open(file_out, 'w')
output_filt = open(file_filt, 'w')

actg=["a", "c", "g", "t", "A", "C", "T", "G"]
nfile = 0
pos_col = -1
chr_col = -1
n_col = -1
af_col = -1
alt_col = -1
ref_col = -1
for file_input in files_input:
  lastpos = 0
  input_results = open(file_input, 'r')
  nfile += 1;
  
  file_keep_list = re.sub(dir_res+"/", "", file_input)
  # file_keep_list = re.sub(trait+".", "", file_keep_list)
  file_keep_list = re.sub("^[/\w_]*.", "", file_keep_list)
  if flag_wald:
    file_keep_list = file_keep_list.rstrip("vcf.gz.wald.out")   
  else:
    file_keep_list = file_keep_list.rstrip(".score.snp.out")
  file_keep_list = file_keep_list+".info.txt"
  print "file_keep_list:", file_keep_list
  input_keep_list = open(dir_keep_list+"/"+file_keep_list, 'r') 
  
  snp = input_keep_list.readline()
  snp = input_keep_list.readline()
  snp = snp.rstrip()
  snp = snp.split(' ')
  pos = int(snp[1])
  info = float(snp[2])
  rsid = snp[0]
  # print snp
  # print pos
  nrow = 0
  for line in input_results:
    cantfind = 0
    # print line
    nrow += 1
    if nrow == 1:
      if nfile == 1:
        line = line.rstrip()
        arr = line.split('\t')
        for i in range(0, len(arr)):
          if arr[i] == "POS":
            pos_col = i
          elif arr[i] == "CHROM":
            chr_col = i
          elif arr[i] == "ALT_AF":
            af_col = i
          elif arr[i] == "N_INFORMATIVE":
            n_col = i
          elif arr[i] == "REF":
            ref_col = i
          elif arr[i] == "ALT":
            alt_col = i
        output_file.write("CHRPOS\tCPTID\t"+line+"\tIMPUTE2_INFO\tRS_ID\tESS\tis_snp\n")
        output_filt.write("CHRPOS\tCPTID\t"+line+"\tIMPUTE2_INFO\tRS_ID\tESS\tis_snp\n")
    else:
      line = line.rstrip()
      arr = line.split('\t')
      if int(arr[pos_col]) < pos:
        if cantfind == 0 and lastpos > int(arr[pos_col]):
          pos = 0
          # print("rewinding")
          input_keep_list.seek(0,0)
          cantfind = 1
      
      while int(arr[pos_col]) > pos:
        snp = input_keep_list.readline()
        if snp == "":
          break
        elif snp != "":
          snp = snp.rstrip()
          snp = snp.split(' ')
          if snp[1] == "position":
            continue
          pos = int(snp[1])
          info = float(snp[2])
          rsid = snp[0]
      
      cptid = arr[chr_col]+":"+arr[pos_col]+":"+arr[ref_col]+":"+arr[alt_col]
      chrpos = arr[chr_col]+":"+arr[pos_col]
      maf = float(arr[af_col])
      if maf > 0.5:
        maf = 1.-maf
      n = float(arr[n_col])
      mac = 2.*n*maf
      ess = mac*(1-maf)*info
      
      issnp = (arr[ref_col] in actg and arr[alt_col] in actg) + 0
      
      if int(arr[pos_col]) == pos:
        lastpos = int(arr[pos_col])
        outstring = chrpos+'\t'+cptid+"\t"+line+'\t'+str(info)+'\t'+rsid+'\t'+str(ess)+'\t'+str(issnp)+'\n'
        if ess > 30 and info > 0.3:
          output_file.write(outstring)
        else:
          output_file.write(outstring)
      elif snp == "":
        break
      else:
        # sys.exit("Error: Variant "+arr[chr_col]+":"+arr[pos_col]+" is not in the metric file!")
        # print("Variant "+arr[chr_col]+":"+arr[pos_col]+" is not in the metric file!")
        # input_keep_list.seek(0,0)
        if mac > 30:
          output_filt.write(chrpos+'\t'+cptid+"\t"+line+'\t\t\t\t'+str(issnp)+'\n')
        else:
          output_filt.write(chrpos+'\t'+cptid+"\t"+line+'\t\t\t\t'+str(issnp)+'\n')
        
  input_keep_list.close() 
  input_results.close()
  
output_file.close()
