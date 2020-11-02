#!/usr/bin/env python3


# Importing modules
import sys
import re
import os
import json


# Reading-in BAM file
sample = sys.argv[1]


# Getting name of the sample
sampleName = os.path.basename(sample)
sampleName = re.sub('_Log.final.out', '', sampleName)


# Dictionary for stats
data = {}
data['Sample'] = []
data['Number_of_input_reads'] = []
data['Average_input_read_length'] = []
data['Number_of_uniquely_mapped_reads'] = []
data['Uniquely_mapped_reads_%'] = []
data['Average_mapped_length'] = []
data['Mismatch_rate_per_base_%'] = []
data['Number_of_reads_mapped_to_multiple_loci'] = []
data['%_of_reads_mapped_to_multiple_loci'] = []
data['%_of_reads_unmapped'] = []
data['Number_of_reads_unmapped'] = []


# Opening file and reading-in lines
if sample.endswith("_Log.final.out"):
    with open (sample, 'r') as wF:
        for count, line in enumerate(wF):

            # Number of input reads
            if count == 5:
                inputReads = line.strip().split("\t")

            # Average input read length
            elif count == 6:
                average = line.strip().split("\t")

            # Uniquely mapped reads number
            elif count == 8:
                uniqueNum = line.strip().split("\t")

            # Uniquely mapped reads percent
            elif count == 9:
                uniquePer = line.strip().split("\t")

            # Average mapped length
            elif count == 10:
                averLen = line.strip().split("\t")

            # Mismatch rate per base
            elif count == 17:
                mismatch = line.strip().split("\t")

            # Number of read mapped to multiple loci
            elif count == 23:
                multiLoci = line.strip().split("\t")

            # % of reads mapped tp multiple loci
            elif count == 24:
                multiPer = line.strip().split("\t")
                multiPer = multiPer[1].split("%")

            # Number of reads mapped to too many loci
            elif count == 25:
                manyLoci = line.strip().split("\t")
                multiNum = int(multiLoci[1]) + int(manyLoci[1])

            # % of reads mapped to too many loci
            elif count == 26:
                manyPer = line.strip().split("\t")
                manyPer = manyPer[1].split("%")

                # Multi-mapping reads
                mmPer = float(multiPer[0]) + float(manyPer[0])

            # % of reads unmapped: too many mismatches
            elif count == 29:
                unmapMis = line.strip().split("\t")
                unmapMis = unmapMis[1].split("%")

            # % of reads unmapped: too short
            elif count == 31:
                unmapShort = line.strip().split("\t")
                unmapShort = unmapShort[1].split("%")

            # % of reads unmapped: other
            elif count == 33:
                unmapOth = line.strip().split("\t")
                unmapOth = unmapOth[1].split("%")

                # % of reads unmapped
                unmapped = float(unmapMis[0]) + float(unmapShort[0]) + float(unmapOth[0])

            elif count > 29:
                unmapNum = int(inputReads[1]) - (int(uniqueNum[1]) + int(multiNum))

        # Writing stats in dictionary
        data['Sample'] = sampleName
        data['Number_of_input_reads'] = int(inputReads[1])
        data['Average_input_read_length'] = int(average[1])
        data['Number_of_uniquely_mapped_reads'] = int(uniqueNum[1])
        data['Uniquely_mapped_reads_%'] = uniquePer[1]
        data['Average_mapped_length'] = float(averLen[1])
        data['Mismatch_rate_per_base_%'] = mismatch[1]
        data['Number_of_reads_mapped_to_multiple_loci'] = int(multiNum)
        data['%_of_reads_mapped_to_multiple_loci'] = str(mmPer) + "%"
        data['%_of_reads_unmapped'] = str(unmapped) + "%"
        data['Number_of_reads_unmapped'] = int(unmapNum)

# Path for working directory and file name
workingDire = os.path.join("Results/" + sampleName + "/JSON_Files/" + sampleName + "_Mapping_Rates_Sample.json")

# Writing dictionary in JSON file
with open(workingDire, 'w') as working:
    json.dump(data, working, indent=2)
