awk '{p=$3-$2;getline;d=$3-$2;printf($4"\t"p"\t"d"\n")}' PairwiseApa.bed|sed 's/_dis//g'>PairwiseApa_Length.bed
