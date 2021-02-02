#!/bin/bash

usage() { echo "Usage: $0 [-p <ddseq|10x>] (-n <number_of_cores>)" 1>&2; exit 1; }


cores=8
while getopts ":p:n:" o; do
	case "${o}" in
		p)
			p=${OPTARG}
			if [[ $p != "ddseq" && $p != "10x" ]]; then
				usage
			fi
			;;
		n)	cores=${OPTARG}
			;;
		*)
			usage
			;;
	esac
done

if [ -z "${p}" ];then
	usage
fi

if [ -d Samples/transformed ]; then
	rm -rf Samples/transformed
fi
mkdir Samples/transformed

echo $file

for file in ./Samples/*_R2*.fastq; do
	# Extract sample name, find corresponding R1.fastq (for 10X barcodes) and Extract barcode length
	SAMPLE=$(basename $file .fastq)
	READ1=$(echo Samples/$SAMPLE | sed s/R2/R1/g).fastq
	echo $READ1

	if [[ $p == "10x" ]]; then 
		echo "10x"
		Scripts/read_transform_10x.py $READ1 $file Samples/transformed/${SAMPLE}.fastq
	elif [[ $p == "ddseq"  ]]; then
		echo "ddSeq"
		cp $file Samples/transformed/${SAMPLE}.fastq
	fi
done

snakemake -kpj ${cores}
