#!/bin/bash


input="tench.txt"


while IFS= read -r line
do
  printf "\n\nDownloading " ; printf $line; echo " . . . \n"
  fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $line
done < "$input"