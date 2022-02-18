#!/bin/sh



for file in $(seq 869 904 | sed 's/[^ ]* */SRR11487&/g'); do
	printf "\n\nDownloading " ; printf $file; echo " . . . \n\n"
	fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $file
done

