#!/usr/bin/env python3

# Importing modules
import sys
import os
import json


# Reading-in BAM file
sample = sys.argv[1]
sampleName = sys.argv[2]


# Dictionary for stats
data = {}
data['Sample'] = []
data['Assigned'] = []
data['Unassigned_unmapped'] = []
data['Unassigned_multimapping'] = []
data['Unassigned_no_features'] = []
data['Unassigned_ambiguity'] = []


# Opening file and reading-in lines
if sample.endswith(".summary"):
    with open (sample, 'r') as wF:
        for line in wF:

            # Number of assigned reads
            if line.startswith("Assigned"):
                assigned = line.strip().rsplit("\t")

            # Number of unmapped reads
            elif "Unassigned_Unmapped" in line:
                unmapped = line.strip().rsplit("\t")

            # Number of unassigned reads: MultiMapping
            elif "Unassigned_MultiMapping" in line:
                unassiMulti = line.strip().rsplit("\t")

            # Number of unassigned reads: NoFeatures
            elif "Unassigned_NoFeatures" in line:
                unassiNoFeat = line.strip().rsplit("\t")

            # Number of unassigned reads: Ambiguity
            elif "Unassigned_Ambiguity" in line:
                unassiAmbig = line.strip().rsplit("\t")

        # Writing stats in dictionary
        data['Sample'] = sampleName
        data['Assigned'] = int(assigned[1])
        data['Unassigned_unmapped'] = int(unmapped[1])
        data['Unassigned_multimapping'] = int(unassiMulti[1])
        data['Unassigned_no_features'] = int(unassiNoFeat[1])
        data['Unassigned_ambiguity'] = int(unassiAmbig[1])


# Path for working directory and file name
workingDire = os.path.join("Results/" + sampleName + "/JSON_Files/" + sampleName + "_Gene_Counts_Sample.json")

# Writing dictionary in JSON file
with open (workingDire, 'w') as working:
    json.dump(data, working, indent=2)
