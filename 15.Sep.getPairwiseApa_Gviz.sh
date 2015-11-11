awk 'BEGIN{gentline;} {a=$0;getline; b=$0;print a"\t"b;}' PairwiseApa.bed > tmp
awk '{l=($2<$8?$2:$8);r=($3>$9?$3:$9);nm=split($4,a,"_pro");printf($1"\t"l"\t"r"\t"a[1]"\t"$5"\t"$6"\n")}' tmp > PairwiseApa_GViz.bed
rm tmp
