#!/usr/bin/env python3

import sys
import os
import json

# Reading files in
umiPath = os.path.dirname(sys.argv[1])
sampleName = sys.argv[2]

# Dictionary for UMI stats
data = {}
data['Cell_barcode'] = []
data['UMI_counts'] = []
data['Gene_counts'] = []


# Taking all files in directories
for dir, subdir, files in os.walk(umiPath):

    # Get every file
    for file in files:

        # Control if is an TSV file
        if file.endswith(".tsv"):
            #print(file)

            cellbarcode = os.path.splitext(file)
            #print(cellbarcode)
            filePath = os.path.join(dir, file)
            #print(filePath)

            # Opening files
            with open(filePath, 'r') as uFile:
                umiCounter = 0
                geneCounter = 0

                # Count UMIs and genes
                for count, line in enumerate(uFile):

                    if int(count) > 0:
                        #print(line)

                        # Number of genes in cell barcode
                        geneCounter = geneCounter + 1

                        #Number of unique UMIs in cell barcode
                        spLine = line.strip().split("\t")
                        umiCounter = umiCounter + int(spLine[1])

                #print(file + " " + str(geneCounter) + " " + str(umiCounter))
                data['Cell_barcode'].append(cellbarcode[0])
                data['UMI_counts'].append(umiCounter)
                data['Gene_counts'].append(geneCounter)


# Creating a path for JSON file
workingDir = os.path.join("Results/" + sampleName + "/JSON_Files/" + sampleName + "_UMI_Counts.json")

# Writing the directory in a JSON file
with open(workingDir, 'w') as w:
    json.dump(data, w, indent=2)
