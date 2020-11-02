#!/usr/bin/env python3

import sys
import json
import os
import pysam


# Reading files
cellcounter = sys.argv[1]
mapping = sys.argv[2]
validBc = pysam.AlignmentFile(sys.argv[3], 'rb')
sampleName = sys.argv[4]


# Opening file with number of counted cell barcodes
with open(cellcounter, 'r') as cFile:

    # Getting the number of counted cell barcodes
    cellbarcodes = 0
    cellbarcodes = cFile.readline().strip().split("\t")
    #print(cellbarcodes[1])


# Opening STAR log file
with open(mapping, 'r') as mFile:
    rawReads=0

    # Reading number of raw reads
    for count, line in enumerate(mFile):

        if count == 5:

            rawReads = line.strip().split("\t")
            #print(rawReads)


# Opening file with all reads with valid barcodes
# Reading every line
bcReads = 0
for line in validBc.fetch(until_eof=True):

    read = str(line).split("\t")
    #print(read[1])

    # Check if read is not a multimapping
    if int(read[1]) != 256 and int(read[1]) != 272:

        bcReads = bcReads + 1

#print(bcReads)


# Saving all stats in a dictionary
sampleCells = {"Sample":sampleName, "Cell_barcodes":int(cellbarcodes[1]), "Raw_reads":int(rawReads[1]), "Reads_with_valid_barcode":int(bcReads)}
#print(sampleCells)


# Creating a name for JSON file
jsonFile = os.path.join("Results/" + sampleName + "/JSON_Files/" + sampleName + "_Cell_Numbers.json")


# Writing the dictionary in a JSON file
with open(jsonFile, 'w') as jFile:
    json.dump(sampleCells, jFile)
