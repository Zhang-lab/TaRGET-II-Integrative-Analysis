count_cutoff=10

for DSS_file in DSS_files/DSS_*
do
  tag=${DSS_file##*/}

  tag=${tag//"_sorted.txt_CpG"/}
  tag=${tag//"DSS_"/}


  cmd="cat "$DSS_file" | awk 'NR>1 {if(\$3 >= ${count_cutoff}) print(\$1\"\t\"\$2\"\t\"\$2+1\"\t\"\$4/\$3)}' | bedtools sort -i - > test_"$tag" "
  echo $cmd
  eval $cmd

done
