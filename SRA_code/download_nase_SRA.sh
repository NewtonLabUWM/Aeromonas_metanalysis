#!/bin/sh



for file in $(seq 4 7 | sed 's/[^ ]* */ERR345855&/g'); do
	printf "\nDownloading " ; printf $file; echo " . . . "
	fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $file
done

#ERR3458554 ERR3458557




for file in $(seq 90 93 | sed 's/[^ ]* */ERR34585&/g'); do
	printf "\nDownloading " ; printf $file; echo " . . . "
	fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $file
done

#ERR3458590 ERR3458593



for file in $(seq 2 5 | sed 's/[^ ]* */ERR345864&/g'); do
	printf "\nDownloading " ; printf $file; echo " . . . "
	fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $file
done

#ERR3458642 ERR3458645



for file in $(seq 4 7 | sed 's/[^ ]* */ERR345867&/g'); do
	printf "\nDownloading " ; printf $file; echo " . . . "
	fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $file
done

#ERR3458674 ERR3458677




for file in $(seq 6 9 | sed 's/[^ ]* */ERR345874&/g'); do
	printf "\nDownloading " ; printf $file; echo " . . . "
	fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip $file
done

#ERR3458746 ERR3458749