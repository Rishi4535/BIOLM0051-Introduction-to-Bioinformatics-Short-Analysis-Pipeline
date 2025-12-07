#Looking into all the files to find anomalies and edit them out accordingly
for file in *.FASTQ; do
  echo "=== $file ==="
  cat "$file"
done

#For cleaning the FASTQ files before conversion

#Making all the headers homogenous in Sample A files
for file in sampleA_part1.fastq sampleA_part2.fastq sampleA_part3.fastq; do
  awk 'NR==1{count=0; for(i=1;i<=length($0);i++){c=substr($0,i,1); if(c=="_"){count++; if(count==2)continue} printf "%s",c} print ""} NR>1' "$file" > temp && mv temp "$file"
done


#Deleting the duplicate headers in Sample B files
sed -i '1d' sampleB_part1.fastq
