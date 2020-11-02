#!/usr/bin/env python3

import sys
import json
import os
import pysam

# Reading file in
workingFile = pysam.AlignmentFile(sys.argv[1], "rb")
workingSample = sys.argv[2]

# Dictionary
barcodeDict = {}
data = {}

# Cellbarcode stats
data["Cell_barcode"] = []
data["Assigned"] = []
data["Unassigned_unmapped"] = []
data["Unassigned_no_features"] = []
data["Unassigned_ambiguity"] = []
data["Unassigned_multimapping"] = []


# Opening file
# Read every line
for line in workingFile.fetch(until_eof=True):

    # Taking the barcode and splitting the read
    cellBc = line.query_name.split("_")
    read = str(line).strip().split("\t")

    # Write each barcode in a dictionary as key if it is not a multimapping with the assignment
    if int(read[1]) != 256 and int(read[1]) != 272:

        assigTag = line.get_tag("XS")

        if cellBc[1] in barcodeDict.keys():

            barcodeDict[cellBc[1]].append(assigTag)

        else:

            barcodeDict[cellBc[1]] = [assigTag]

#print(barcodeDict)

# Counting metrics for each cell barcode
for cellBar in barcodeDict.keys():
    assigned = 0
    unmapped = 0
    noFeatures = 0
    ambiguity = 0
    multimapping = 0

    for element in barcodeDict[cellBar]:

        # Counting assigned reads
        if "Assigned" in element:

            assigned = assigned + 1
            #print(element)

        # Counting unmapped reads
        elif "Unmapped" in element:

            unmapped = unmapped + 1
            #print(element)

        # Counting reads with no features
        elif "NoFeatures" in element:

            noFeatures = noFeatures + 1
            #print(element)

        # Counting ambiguous reads
        elif "Ambiguity" in element:

            ambiguity = ambiguity + 1
            #print(element)

        # Counting multimapping
        else:

            multimapping = multimapping + 1
            #print(element)

    # Saving all metrics for each cellbarcode in directory
    if not int(assigned) == 0:

        data["Cell_barcode"].append(workingSample + "_Zelle_" + str(cellBar))
        data["Assigned"].append(int(assigned))
        data["Unassigned_unmapped"].append(int(unmapped))
        data["Unassigned_no_features"].append(int(noFeatures))
        data["Unassigned_ambiguity"].append(int(ambiguity))
        data["Unassigned_multimapping"].append(int(multimapping))

# Writing directory in JSON file
workingDire = os.path.join("Results/" + workingSample + "/JSON_Files/" + workingSample + "_Gene_Counts.json")

with open(workingDire, "w") as wD:
    json.dump(data, wD, indent=2)
