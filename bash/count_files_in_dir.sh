# OPTION 1: USE THIS FOR CT_TMP
# -bash-4.1$ pwd
# /cresis/snfs1/dataproducts/ct_data/ct_tmp
#find . -maxdepth 3 -mindepth 3 -type d -print0 | sort | while read -d '' -r dir; do

# OPTION 2: USE THIS FOR OUTPUT
# -bash-4.1$ pwd
# /cresis/snfs1/dataproducts/ct_data/accum
# /cresis/snfs1/dataproducts/ct_data/kuband
# /cresis/snfs1/dataproducts/ct_data/rds
# /cresis/snfs1/dataproducts/ct_data/snow
find . -maxdepth 2 -mindepth 2 -iname "*tmp" -type d -print0 | sort | while read -d '' -r dir; do
  #files=("$dir"/*)
  num=$(find "$dir" -ls | wc -l);
  printf "%5d files in directory %s\n" "$num" "$dir"
done
