#!/bin/sh



for file in $(seq 705 878 | sed 's/[^ ]* */ERR2129&/g'); do
	printf "\nDownloading " ; printf $file; echo " . . . "
	fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $file
done

#ERR2129705 ERR2129878





