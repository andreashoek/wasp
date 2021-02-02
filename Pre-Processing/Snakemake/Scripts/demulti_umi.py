#!/usr/bin/env python3

# Import of modules
import pysam
import sys
import Levenshtein
import py


# Reading-in BAM file of featureCounts analysis
bamFiles = pysam.AlignmentFile(sys.argv[1], "rb")
fileForDetection = pysam.AlignmentFile(sys.argv[1], "rb")
headerFile = pysam.AlignmentFile(sys.argv[2], "wb", template=bamFiles)
bclist10xv2 = set(line.strip() for line in open('Scripts/whitelist_10X_737K-august-2016.txt'))
bclist10xv3 = set(line.strip() for line in open('Scripts/whitelist_10X_3M-february-2018.txt'))

# Barcode whitelist
barcodeList = ["AACAGC", "ATACTT", "CCTCTA", "AACGTG", "ATAGCG", "CGAAAG", "GAGGCC", "GTACAG", "TCTAGC", "AAGCCA", "ATATAC", "CGAGCA", "GAGTGA", "AAAGAA", "AGTCTG", "GAGCTT", "AAGTAT", "ATCCGG", "CGCATA", "GATCAA", "AATTGG", "ATGAAG", "ACAAGG", "ATTAGT", "CCGTAA", "GACTCG", "GGTAGG", "TCGCCT", "GGTGCT", "TCGGGA", "GTCCTA", "TGAATT","GTCGGC", "TGAGAC", "CGGCGT", "GCCAGA", "GTGGTG", "TGCGGT", "CGGTCC", "GCCGTT", "GTTAAC", "TGCTAA", "GTTTCA", "TGGCAG", "ACCCAA", "CAACCG", "CGTTAT", "GCGAAT", "ACCTTC","GCGCGG", "TAAGCT", "TGTGTA", "ACGGAC", "CACCAC", "CTATTA", "GCTCCC", "TAATAG", "TGTTCG", "ACTGCA", "CACTGT","CTCAAT", "GCTGAG", "TACCGA", "TTAAGA", "AGACCC","CAGACT", "CTGTGG", "GCTTGT", "AGATGT", "CAAGTC", "CTAGGT", "CAGGAG", "CTTACG", "AGCACG", "CATAGA", "CTTGAA", "TAGAGG", "TTCGCA", "GGACGA", "TATTTC", "TTCTTG", "GGATTG", "TCAGTG", "TTGCTC","GGCCAT", "TCATCA", "TTGGAT", "AGGTTA", "CCACGC", "GAAATA", "AGTAAA", "CCGATG", "GAAGGG", "GGGATC", "TCCAAG", "TTTGGG"]


# Function to compare two sequences, allowing 1 ED
linker1 = 'TAGCCATCGCATTGC'
linker2 = 'TACCTCTGAGCTGAA'

def linkerPos(read, linker):
    min = 1
    pos = 0
    hit = False
    for i in range(len(read) - len(linker)+1):
        substr = read[i:i+len(linker)]
        hamDist = Levenshtein.hamming(substr, linker)
        if hamDist <= min:
            min = hamDist
            pos = i+1
            hit = True
    return (hit,pos)


# Function to compare barcodes of read with barcodes of whitelist
def validBarcodes(rbc, bcList):
    newBarcode = ""
    for bcL in bcList:
        if Levenshtein.hamming(rbc, bcL) <= 1:
            newBarcode = bcL
    return newBarcode


def checkBarcode10xv2():
    for read in bamFiles.fetch(until_eof=True):
        # Splitting barcode of read
        readSeq = read.qname.split(":")
        rawBarcode = readSeq[7]
        rawBarcode = rawBarcode.strip()
        readBC = rawBarcode[0:16]

        #if not rawBarcode.startswith("*"):
        if readBC in bclist10xv2:
            validBc = readBC
            readUmi = rawBarcode[16:26]
            if "Assigned" in str(read):
                bcRead = str(read)
                spRead = bcRead.split("\t")
                a = pysam.AlignedSegment()
                a.query_name = str(
                    readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[4] + ":" +
                    readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
                a.mapping_quality = int(spRead[4])
                a.query_sequence = spRead[9]
                a.flag = int(spRead[1])
                a.reference_start = int(spRead[3])
                a.reference_id = int(spRead[2])
                a.tags = (("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")), ("XT", read.get_tag("XT")))

                headerFile.write(a)
            else:
                bcRead = str(read)
                spRead = bcRead.split("\t")
                a = pysam.AlignedSegment()
                a.query_name = str(
                    readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[4] + ":" +
                    readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
                a.mapping_quality = int(spRead[4])
                a.query_sequence = spRead[9]
                a.flag = int(spRead[1])
                a.reference_start = int(spRead[3])
                a.reference_id = int(spRead[2])
                a.tags = (("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")))

                headerFile.write(a)


def checkBarcode10xv3():
    for read in bamFiles.fetch(until_eof=True):
        # Splitting barcode of read
        readSeq = read.qname.split(":")
        rawBarcode = readSeq[7]
        rawBarcode = rawBarcode.strip()
        readBC = rawBarcode[0:16]
        #if not rawBarcode.startswith("*"):
        if readBC in bclist10xv3:
            validBc = readBC
            readUmi = rawBarcode[16:29]
            if "Assigned" in str(read):
                bcRead = str(read)
                spRead = bcRead.split("\t")
                a = pysam.AlignedSegment()
                a.query_name = str(readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[4] + ":" + readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
                a.mapping_quality = int(spRead[4])
                a.query_sequence = spRead[9]
                a.flag = int(spRead[1])
                a.reference_start = int(spRead[3])
                a.reference_id = int(spRead[2])
                a.tags = (("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")), ("XT", read.get_tag("XT")))

                headerFile.write(a)
            else:
                bcRead = str(read)
                spRead = bcRead.split("\t")
                a = pysam.AlignedSegment()
                a.query_name = str(readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[4] + ":" + readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
                a.mapping_quality = int(spRead[4])
                a.query_sequence = spRead[9]
                a.flag = int(spRead[1])
                a.reference_start = int(spRead[3])
                a.reference_id = int(spRead[2])
                a.tags = (("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")))

                headerFile.write(a)

            # print("Barcode: ", rawBarcode[0:17])
            # print("UMI: ", rawBarcode[17:29])

def checkddSeq():
    # Illumina algorithm follows:
    # Reading-in each read from BAM file
    for read in bamFiles.fetch(until_eof=True):

        # Splitting barcode of read
        readSeq = read.qname.split(":")
        rawBarcode = readSeq[7]
        # print(readSeq)

        # Determine position of linker 1
        posLinker1 = linkerPos(rawBarcode, linker1)
        # print posLinker1

        # Position of linker 1 has to be bigger than 6, for BC1 to be complete
        if posLinker1[1] > 5:  # and int(posLinker1[1])<14:
            # print(posLinker1)

            # Determine position of linker 2
            posLinker2 = linkerPos(rawBarcode, linker2)
            # print(posLinker2)
            # if posLinker2[0]:
            # print(posLinker2)

            # Position of linker 2 has to be exactly position of linker 1 + 21 bases (no insertion or deletion allowed)
            if int(posLinker1[1]) + 21 == int(posLinker2[1]):
                # print(posLinker2)

                # Cutting of phase block at the beginning of the read
                matchLinker = rawBarcode[int(posLinker1[1]) - 7:]
                # print(matchLinker)

                # Cutting of access bases at the end of the read
                matchLinker = matchLinker[:int(posLinker2[1] + 34)]
                # print(matchLinker)

                # Read has to be at least 62 bases to be complete
                if len(matchLinker) >= 62:
                    # print(matchLinker)

                    # Extracting ACG and GAC sequences
                    acg = matchLinker[48:51]
                    # print(acg)
                    gac = matchLinker[59:62]
                    # print(gac)

                    # Controlling sequences while allowing 1 ED
                    if Levenshtein.hamming(acg, "ACG") <= 1 and Levenshtein.hamming(gac, "GAC") <= 1:
                        # print(matchLinker)

                        # Comparing barcodes with witelist and extracting UMI sequence
                        bc1 = validBarcodes(matchLinker[0:6], barcodeList)
                        # print bc1
                        bc2 = validBarcodes(matchLinker[21:27], barcodeList)
                        # print bc2
                        bc3 = validBarcodes(matchLinker[42:48], barcodeList)
                        # print bc3
                        validBc = bc1 + bc2 + bc3
                        # print validBc
                        readUmi = matchLinker[51:59]
                        # print(validBc + "\t" + readUmi)

                        # Composed barcodes have to have a length of exactly 18 bases
                        if len(validBc) == 18 and "Assigned" in str(read):

                            bcRead = str(read)
                            spRead = bcRead.split("\t")

                            a = pysam.AlignedSegment()
                            a.query_name = str(
                                readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[
                                    4] + ":" + readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
                            a.mapping_quality = int(spRead[4])
                            a.query_sequence = spRead[9]
                            a.flag = int(spRead[1])
                            a.reference_start = int(spRead[3])
                            a.reference_id = int(spRead[2])
                            a.tags = (
                            ("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")), ("XT", read.get_tag("XT")))

                            headerFile.write(a)

                        elif len(validBc) == 18:

                            bcRead = str(read)
                            spRead = bcRead.split("\t")

                            a = pysam.AlignedSegment()
                            a.query_name = str(
                                readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[
                                    4] + ":" + readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
                            a.mapping_quality = int(spRead[4])
                            a.query_sequence = spRead[9]
                            a.flag = int(spRead[1])
                            a.reference_start = int(spRead[3])
                            a.reference_id = int(spRead[2])
                            a.tags = (("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")))

                            headerFile.write(a)






for read in fileForDetection.fetch(until_eof=True):
    readSeq = read.qname.split(":")
    rawBarcode = readSeq[7]

    print(rawBarcode)

    if rawBarcode.startswith("*"):
        continue
    else:
        if len(rawBarcode) == 26:
            print("10X v2")
            checkBarcode10xv2()
            break
        elif len(rawBarcode) == 28:
            print("10X v3")
            checkBarcode10xv3()
            break
        else:
            print("ddSeq")
            checkddSeq()
            break
fileForDetection.close()
bamFiles.close()



# Reading-in each read from BAM file
# for read in bamFiles.fetch(until_eof=True):
    # Splitting barcode of read
    # readSeq = read.qname.split(":")
    # rawBarcode = readSeq[7]
    #print(readSeq)

    #print(rawBarcode)
    #break


    #if rawBarcode.startswith("*"):



    # # Determine position of linker 1
    # posLinker1 = linkerPos(rawBarcode, linker1)
    # #print posLinker1
    #
    # # Position of linker 1 has to be bigger than 6, for BC1 to be complete
    # if posLinker1[1] > 5: #and int(posLinker1[1])<14:
    #     #print(posLinker1)
    #
    #     # Determine position of linker 2
    #     posLinker2 = linkerPos(rawBarcode, linker2)
    #     #print(posLinker2)
    #     #if posLinker2[0]:
    #         #print(posLinker2)
    #
    #     # Position of linker 2 has to be exactly position of linker 1 + 21 bases (no insertion or deletion allowed)
    #     if int(posLinker1[1]) + 21 == int(posLinker2[1]):
    #         #print(posLinker2)
    #
    #         # Cutting of phase block at the beginning of the read
    #         matchLinker = rawBarcode[int(posLinker1[1])-7:]
    #         #print(matchLinker)
    #
    #         # Cutting of access bases at the end of the read
    #         matchLinker = matchLinker[:int(posLinker2[1]+34)]
    #         #print(matchLinker)
    #
    #         # Read has to be at least 62 bases to be complete
    #         if len(matchLinker) >= 62:
    #             #print(matchLinker)
    #
    #             # Extracting ACG and GAC sequences
    #             acg = matchLinker[48:51]
    #             #print(acg)
    #             gac = matchLinker[59:62]
    #             #print(gac)
    #
    #             # Controlling sequences while allowing 1 ED
    #             if Levenshtein.hamming(acg, "ACG") <=1 and Levenshtein.hamming(gac, "GAC") <= 1:
    #                 #print(matchLinker)
    #
    #                 # Comparing barcodes with witelist and extracting UMI sequence
    #                 bc1 = validBarcodes(matchLinker[0:6], barcodeList)
    #                 #print bc1
    #                 bc2 = validBarcodes(matchLinker[21:27], barcodeList)
    #                 #print bc2
    #                 bc3 = validBarcodes(matchLinker[42:48], barcodeList)
    #                 #print bc3
    #                 validBc = bc1 + bc2 + bc3
    #                 #print validBc
    #                 readUmi = matchLinker[51:59]
    #                 #print(validBc + "\t" + readUmi)
    #
    #                 # Composed barcodes have to have a length of exactly 18 bases
    #                 if len(validBc) == 18 and "Assigned" in str(read):
    #
    #                     bcRead = str(read)
    #                     spRead = bcRead.split("\t")
    #
    #                     a = pysam.AlignedSegment()
    #                     a.query_name = str(readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[4] + ":" + readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
    #                     a.mapping_quality = int(spRead[4])
    #                     a.query_sequence = spRead[9]
    #                     a.flag = int(spRead[1])
    #                     a.reference_start = int(spRead[3])
    #                     a.reference_id = int(spRead[2])
    #                     a.tags = (("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")), ("XT", read.get_tag("XT")))
    #
    #                     headerFile.write(a)
    #
    #                 elif len(validBc) == 18:
    #
    #                     bcRead = str(read)
    #                     spRead = bcRead.split("\t")
    #
    #                     a = pysam.AlignedSegment()
    #                     a.query_name = str(readSeq[0] + ":" + readSeq[1] + ":" + readSeq[2] + ":" + readSeq[3] + ":" + readSeq[4] + ":" + readSeq[5] + ":" + readSeq[6] + "_" + validBc + "_" + readUmi)
    #                     a.mapping_quality = int(spRead[4])
    #                     a.query_sequence = spRead[9]
    #                     a.flag = int(spRead[1])
    #                     a.reference_start = int(spRead[3])
    #                     a.reference_id = int(spRead[2])
    #                     a.tags = (("NH", read.get_tag("NH")), ("XS", read.get_tag("XS")))
    #
    #                     headerFile.write(a)


# bamFiles.close()


