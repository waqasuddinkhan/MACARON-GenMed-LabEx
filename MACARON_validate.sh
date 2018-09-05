#!/bin/bash

#sub1.chr22_21349676-21349677.sample02.bam

bam=$1;
name=$(echo $bam | cut -d'.' -f1)
region=$(echo $bam | cut -d'.' -f2)
sample=$(echo $bam | cut -d'.' -f3)

chr=$(echo $region | cut -d'_' -f1)
pos=$(echo $region | cut -d'_' -f2)

beg=$(echo $pos | cut -d'-' -f1)
end=$(echo $pos | cut -d'-' -f2)

#beg=$(($beg+50))
#end=$((end-50))

#echo "name $name"
#echo "region $region"
#echo "sample $sample"

echo "$name $chr:$beg-$end $sample"
echo  "$name $chr:$beg-$end $sample" >> MACARON_validate.txt

samtools view $bam | awk -v start=$beg -v stop=$end '
{
  pos=$4;
  if(pos < start){
    cstart=1+start-pos;
    cstop=1+stop-pos;
    print cstart"."cstop"."$10;
  }
}' > tmp

while read line;
do
  cb=$(echo $line | cut -d'.' -f1);
  ce=$(echo $line | cut -d'.' -f2);
  seq=$(echo $line | cut -d'.' -f3);
  echo $seq | cut -c$cb-$ce;
done <tmp | sort | uniq -c | sort -n >> MACARON_validate.txt
rm tmp




exit 0

