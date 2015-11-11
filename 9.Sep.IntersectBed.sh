# 2015-9-5
# Cuihua Gu
# Assign APA peaks to refSeq transcripts.
# INPUT
	# 1. bam file
	# 2. bed file
	# 3. output prefix
# OUTPUT

# USAGE
	# sh script.sh


bam=$1
bed=$2
outp=$3
dir=`dirname $outp`
pref=`basename $outp`
mkdir -p $dir/tmp

# run intersectBed
intersectBed -abam $bam -b $bed -bed -wo -s > $dir/tmp/${pref}.itst.bed

cut -f16 $dir/tmp/${pref}.itst.bed|sort -n|uniq -c|sed 's/^ *//g'|awk -F ' ' '{printf($2"\t"$1"\n")}' > $dir/${pref}.itstBed.wo.freq.tmp
awk -f merge_two_files.awk $dir/${pref}.itstBed.wo.freq.tmp /data/guc1/UTR/Output/2_ApaPair/PairwieApaList.txt > $dir/${pref}.itstBed.wo.freq.tmp1
awk 'NR%2==1' $dir/${pref}.itstBed.wo.freq.tmp1 > $dir/${pref}.tmp1
awk 'NR%2!=1' $dir/${pref}.itstBed.wo.freq.tmp1 > $dir/${pref}.tmp2
paste $dir/${pref}.tmp1 $dir/${pref}.tmp2 | sed 's/ /\t/g' > $dir/${pref}.itstBed.freq
rm $dir/${pref}.tmp1 $dir/${pref}.tmp2 $dir/tmp/${pref}.itst.bed $dir/${pref}.itstBed.wo.freq.tmp*


