for file in /home/mama5390/genome_analysis/data/raw_data/rna_seq/trimmed/trimmed/;
 do
  file_prefix = ${file%%.trim*}
  echo $file_prefix
 done

