for FILE in *_withNormalizedCounts.tsv
do
  echo "$FILE"
  grep "NODE" "$FILE"  | awk '{if ($6 < 0.05 && sqrt($3^2)){print}}' | wc -l
done
