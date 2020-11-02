#!/usr/bin/env python3

import sys
import os
import csv

# Reading file in
umiFile = sys.argv[1]
sampleName = sys.argv[2]

# Dictionary
bcGeneDict = {}

# Dictionary for UMI stats
data = {}
data['Cell_barcode'] = []
data['UMI_counts'] = []
data['Gene_counts'] = []


# Opening file
with open(umiFile) as file:

    # File to save the number of counted cell barcodes
    demultiJson = os.path.join("Temporary/" + sampleName + "/Demultiplexed_" + sampleName + "/CellCounter.txt")

    # Reading each line
    for num, line in enumerate(file):

        if num > 0:
            barcode = line.strip().split("\t")

            # Create dictionary structure, if it is not already existing
            bcPath = os.path.join("Temporary/" + sampleName + "/Demultiplexed_" + sampleName + "/" + barcode[1][0:2] + "/" + barcode[1][
                                                                                                         2:4] + "/")
            #print (bcPath)

            if not os.path.exists(bcPath):
                os.makedirs(bcPath)

            # Create a TSV file for each cell barcode
            cellPath = os.path.join(bcPath + sampleName + "_Zelle_" + barcode[1] + ".tsv")

            # Save cellPath in dictionary, if it is not already saved
            if cellPath not in bcGeneDict:
                bcGeneDict[cellPath] = {}
                bcGeneDict[cellPath][barcode[0]] = barcode[2]
            else:
                bcGeneDict[cellPath][barcode[0]] = (barcode[2])


# For each cell barcode: Write a TSV file with all genes and UMI counts
for cBc in bcGeneDict:
    umiCount = 0
    geneCount = 0

    for gene, umi in bcGeneDict[cBc].items():
        umiCount = umiCount + int(umi)
        geneCount = geneCount + 1


    if umiCount >= 5:

        with open(cBc, "w") as tsvFile:
            writer = csv.writer(tsvFile, delimiter="\t")
            writer.writerow(["Gene", "Count"])
            for gene, umi in bcGeneDict[cBc].items():
                writer.writerow([str(gene), int(umi)])


# Writing the number of cell barcodes in a TSV file
with open(demultiJson, "w") as jFile:
    writer = csv.writer(jFile, delimiter="\t")
    writer.writerow(["Cellbarcodes", len(bcGeneDict.keys())])
