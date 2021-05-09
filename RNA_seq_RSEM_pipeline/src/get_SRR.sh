#!/bin/bash

listGSM=$(sed 1d $1 | awk '{print $1}');
tab=$'\t'

for i in ${listGSM}; 
do
  echo 'Retrieving SRR numbers from server...'
  echo $i
  GSM=${i};
  SRR=$(esearch -db sra -query ${i} | efetch -format runinfo | cut -d "," -f1 | awk 'NR==2');
  echo $GSM${tab}$SRR >> $2
done
wait
wait
