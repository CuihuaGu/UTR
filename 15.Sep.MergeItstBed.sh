# 2015-9-15
# Cuihua Gu
# Merge intersectBed output into one file.
# INPUT
	# 1. input data pref(./KIRP/KIRP)
	# 2. numStart, in which ID the file name start with
	# 3. numEnd
# OUTPUT
	# 1. A table contains apapair name and read counts

pref=$1
IDstart=$2
IDend=$3

for i in `seq $IDstart $IDend`
do
	cut -f2,4 ${pref}.T${i}.itstBed.freq > ${pref}.T${i}.itstBed.freq2
	cut -f2,4 ${pref}.N${i}.itstBed.freq > ${pref}.N${i}.itstBed.freq2
	paste ${pref}.T${i}.itstBed.freq2 ${pref}.N${i}.itstBed.freq2 | sed 's/\t/,/g'> ${pref}.${i}.DP2
	sed "1i\\p$i" ${pref}.${i}.DP2 > ${pref}.${i}.DPtmp
done
paste ${pref}.*.DPtmp > ${pref}.all.DPtmp
paste /data/guc1/UTR/Output/2_ApaPair/ApaNameList.txt ${pref}.all.DPtmp > ${pref}.TdTpNdNp.txt
tmp=`dirname ${pref}`
rm ${tmp}/*.DP2 ${tmp}/*.DPtmp
rm ${tmp}/*.freq2
