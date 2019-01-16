declare -a arr=("COI" "16S" "18S" "28S")
for i in "${arr[@]}"
do
   echo "$i"
   rm -rf $i
   mkdir $i
   cd $i
   dbcAmplicons abundance -S ../workshopSamplesheet.txt -O abundance_output -F ../6.create_fixrank/$i.fixrank  --biom > $i.abundance.log
   cd ..
done
