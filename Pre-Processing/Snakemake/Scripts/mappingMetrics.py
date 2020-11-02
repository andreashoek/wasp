#!/usr/bin/env python3

import sys
import json
import os
import pysam

# Reading file in
workingFile = pysam.AlignmentFile(sys.argv[1], "rb")
workingSample = sys.argv[2]

# Dictionaries
barcodeFlag = {}
barcodeMap = {}
data = {}

# Cellbarcode stats
data['Cell_barcode'] = []
data['Number_of_uniquely_mapped_reads'] = []
data['Number_of_reads_unmapped'] = []
data['Number_of_reads_mapped_to_multiple_loci'] = []


# Reading every line
for line in workingFile.fetch(until_eof=True):
    #print(line)

    # Taking only the barcode
    cellBc = line.qname.split("_")
    propFlag = line.flag
    mapFlag = line.mapping_quality
    #print(mapFlag)

    # Write each barcode in a dictionary as key if it is not a multimapping
    if int(propFlag) != 256 and int(propFlag) != 272:
        #print(cellBc)

        if cellBc[1] in barcodeFlag.keys():

            barcodeFlag[cellBc[1]].append(propFlag)
            barcodeMap[cellBc[1]].append(mapFlag)

        else:

            barcodeFlag[cellBc[1]] = [propFlag]
            barcodeMap[cellBc[1]] = [mapFlag]
            #print(str(cellBc[0])+str(cellBc[6]))

#print(len(barcodeFlag.keys()))
#print(len(barcodeMap.keys()))

# Counting mapping metrics
for cellbarcode in barcodeFlag.keys():
    countedReads = 0
    unique = 0
    unmapped = 0
    multimapped = 0

    # Number of all reads in a cellbarcode
    countedReads = len(barcodeFlag[cellbarcode])
    #print(countedReads)

    # Counting mapping quality of reads in cellbarcode
    for mapq in barcodeMap[cellbarcode]:

        if int(mapq) == 255:

            unique = unique + 1

    # Counting unmapped reads in cellbarcode
    for flag in barcodeFlag[cellbarcode]:

        if int(flag) == 4:

            unmapped = unmapped + 1

        multimapped = int(countedReads) - (int(unmapped) + int(unique))


    #print(str(cellbarcode) + ":" + "\tCountedReads " + str(countedReads) + "\tUnique " + str(unique) + "\tUnmapped " + str(unmapped) + "\tMultimapped " + str(multimapped))

    if not int(unique) == 0:

        # Saving all metrics for each cellbarcode in directory
        data['Cell_barcode'].append(workingSample + "_Zelle_" + str(cellbarcode))
        data['Number_of_uniquely_mapped_reads'].append(int(unique))
        data['Number_of_reads_unmapped'].append(int(unmapped))
        data['Number_of_reads_mapped_to_multiple_loci'].append(int(multimapped))


# Writing directory in JSON file
workingDire = os.path.join("Results/" + workingSample + "/JSON_Files/" + workingSample + "_Mapping_Rates.json")

with open(workingDire, 'w') as working:
    json.dump(data, working, indent=2)
