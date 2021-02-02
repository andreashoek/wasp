#!/usr/bin/env python3

import sys

read1File = sys.argv[1]
fileForDetection = sys.argv[1]
read2File = sys.argv[2]
newRead2File = sys.argv[3]

bclist10xv2 = set(line.strip() for line in open('Scripts/whitelist_10X_737K-august-2016.txt'))
bclist10xv3 = set(line.strip() for line in open('Scripts/whitelist_10X_3M-february-2018.txt'))

protocol = ""

bcUmi10x = []

lineCnt = 0
readCnt = 0

with open(fileForDetection, "r") as file:
    for line in file:
        bc = line.strip()[0:16]
        if bc in bclist10xv2:
            if bc in bclist10xv3:
                continue
            else:
                protocol = "10xv2"
                print(bc)
                break
        else:
            if bc in bclist10xv3:
                protocol = "10xv3"
                print(bc)
                break
print("Protocol:", protocol)

with open(read1File, "r") as file:
    for line in file:
        lineCnt += 1
        if lineCnt % 4 == 2:
            if protocol == "10xv2":
                bcUmi10x.append(line.strip()[0:26])
            elif protocol == "10xv3":
                bcUmi10x.append(line.strip()[0:28])
file.close()


with open(read2File, "r") as file, open(newRead2File, "a") as newFile:
    for line in file:
        lineCnt += 1
        if lineCnt % 4 == 1:
            headerSplit = line.split()
            headerParts = headerSplit[0].split(":")
            if len(headerParts) <= 7:
                for x in range(len(headerParts), 7):
                    headerParts.append("X")
            newHeader = ':'.join(headerParts) + ":" + bcUmi10x[readCnt] + " " + headerSplit[1]+"\n"
            newFile.write(newHeader)
            readCnt += 1
        else:
            newFile.write(line)

file.close()
newFile.close()
