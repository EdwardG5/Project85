from Bio import Align
from Bio import SeqIO
from Bio.Align import substitution_matrices
import os
from collections import Counter

###########################################################################

aligner = Align.PairwiseAligner()
aligner.open_gap_score = -10
aligner.extend_gap_score = -2 # Needs to be this high so that the program favors mismatch against N vs opening a long gap
aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")

# reference: str, dna : str
# Output: The aligned reference string, and the aligned DNA as a tuple. 
def align(reference, dna):
    output = format(aligner.align(reference, dna)[0])
    l = len(output)//3
    referenceAligned = output[:l-1]
    dnaAligned = output[2*l:3*l-1]
    return (referenceAligned, dnaAligned)

# Accepts an aligned reference and dna e.g. "AAA---GGG", "AAACCCGGG"
# Returns the contracted dna, and a list of indel locations, ("AAAGGG, [(3, 'CCC')]") meaning insert 'CCC' at index 3.
def getInsertions(refAligned, dnaAligned):
    indels = []
    while "-" in refAligned:
        start = refAligned.index("-")
        end = start
        while end < len(refAligned) and refAligned[end] == "-":
            end += 1
        refAligned = refAligned[:start]+refAligned[end:]
        insert = dnaAligned[start:end]
        dnaAligned = dnaAligned[:start]+dnaAligned[end:]
        indels.append((start, insert))
    return (refAligned, dnaAligned, indels)

# Test cases
assert(getInsertions("AAA---GGG", "AAACCCTTT") == ("AAAGGG", "AAATTT", [(3, "CCC")]))

# Accepts an aligned reference and dna e.g. "AAATTTGGG", "AAA---GGG"
# Returns the dna with "-"s replaced by the 'correct' dna from the reference, 
# and a list of delete locations, ("AAATTTGGG", [(3, 3)]) meaning delete 3 starting at 3
def getDeletions(refAligned, dnaAligned):
    indels = []
    while "-" in dnaAligned:
        start = dnaAligned.index("-")
        end = start
        while end < len(dnaAligned) and dnaAligned[end] == "-":
            end += 1
        insert = refAligned[start:end]
        dnaAligned = dnaAligned[:start]+insert+dnaAligned[end:]
        indels.append((start, len(insert)))
    return (dnaAligned, indels)

# Test cases
assert(getDeletions("AAATTTGGG", "AAA---GGG") == ("AAATTTGGG", [(3, 3)]))


# Replace Ns an aligned reference and dna e.g. "AAATTTGGG", "AAA---GGG"
# Returns the dna with "-"s replaced by the 'correct' dna from the reference, 
# and a list of delete locations, ("AAATTTGGG", [(3, 3)]) meaning delete 3 starting at 3
def replaceNs(refAligned, dnaAligned):
    indels = []
    while "N" in dnaAligned:
        start = dnaAligned.index("N")
        end = start
        while end < len(dnaAligned) and dnaAligned[end] == "N":
            end += 1
        insert = refAligned[start:end]
        dnaAligned = dnaAligned[:start]+insert+dnaAligned[end:]
        indels.append((start, len(insert)))
    return (dnaAligned, indels)

# Test cases
assert(replaceNs("AAATTTGGG", "AAANNNGGG") == ("AAATTTGGG", [(3, 3)]))


# Replaces all abnormal letters in dnaAligned with the 'correct' counterparts in refAligned.
# Creates a list of locations and letters i.e. "AAACCC", "AYACMC" => ("AAACCC", [(1, "Y"), (4, "M")])
def replaceOthers(refAligned, dnaAligned):
    indels = []
    dnaAligned = list(dnaAligned)
    for i in range(len(dnaAligned)):
        if dnaAligned[i] not in ["A", "C", "T", "G"]:
            indels.append((i, dnaAligned[i]))
            dnaAligned[i] = refAligned[i]
    dnaAligned = "".join(dnaAligned)
    return (dnaAligned, indels)

assert(replaceOthers("AAACCC", "AYACMC") == ("AAACCC", [(1, "Y"), (4, "M")]))


# Replaces all abnormal letters in dnaAligned with the 'correct' counterparts in refAligned.
# Creates a list of locations and letters i.e. "AAACCC", "ATACTC" => [(1, "T"), (4, "T")]
def getMismatches(refAligned, dnaAligned):
    mismatches = []
    for i in range(len(dnaAligned)):
        if dnaAligned[i] != refAligned[i]:
            mismatches.append((i, dnaAligned[i]))
    return mismatches

assert(getMismatches("AAACCC", "ATACTC") == [(1, "T"), (4, "T")])


# Accepts a reference and another dna molecule. Returns a list of positions
# where there were insertions, deletions, Ns, other characters, and mismatches.
# This function raises an error if applying the getStringFromInfo inverse
# function does not return the original string. This could be because the
# encoding is broken of that the inverse is broken. Either way, something has
# gone wrong and this ensures that the user is alerted to the issue and won't
# use the incorrect results assuming that they are correct.
def getInfoFromString(reference, dna):
    refAligned, dnaAligned = align(reference, dna)
    refAligned, dnaAligned, insertions = getInsertions(refAligned, dnaAligned)
    dnaAligned, deletions = getDeletions(refAligned, dnaAligned)
    dnaAligned, ns = replaceNs(refAligned, dnaAligned)
    dnaAligned, others = replaceOthers(refAligned, dnaAligned)
    mismatches = getMismatches(refAligned, dnaAligned)
    info = (insertions, deletions, ns, others, mismatches)
    assert(getStringFromInfo(ref, *info) == dna)
    return (insertions, deletions, ns, others, mismatches)

# Helper for getStringFromInfo
# dna : List[str]
def redoInsertions(dna, insertions):
    if insertions == []:
        return dna
    else:
        (position, s) = insertions.pop()
        before, after = dna[:position], dna[position:]
        rec = redoInsertions(before, insertions)
        return rec+list(s)+after

assert(redoInsertions(list("ACAACA---"), [(3, "AAA")]) == list("ACAAAAACA---"))

# The inverse for getInfoFromString.
def getStringFromInfo(reference, insertions, deletions, ns, others, mismatches):
    reference = list(reference)
    dna = list(reference)
    # Insert the mismatches
    for (i, l) in mismatches:
        dna[i] = l
    # Replace the strange characters
    for (i, l) in others:
        dna[i] = l
    # Replace the Ns
    for (startingPosition, length) in ns:
        dna[startingPosition:startingPosition+length] = ["N" for _ in range(length)]
    # Delete the inserted parts
    for (startingPosition, length) in deletions:
        dna[startingPosition:startingPosition+length] = ["-" for _ in range(length)]
    # Fix the insertions
    dna = redoInsertions(dna, insertions)
    # Now remove the "-"s. 
    dna = filter(lambda x: x != "-", dna)
    # Done : ). 
    return "".join(dna)


###########################################################################

charDict = {"A": "00", "C": "01", "T":"10", "G":"11"}
inverseDict = {v: k for k, v in charDict.items()}

# char -> binaryString
def encodeChar(c):
    return charDict[c]

# binaryString -> (char, rest of the string)
def decodeChar(binaryString):
    b = binaryString[:2]
    return (inverseDict[b], binaryString[2:])

# Test cases
# Unecessary

# Input: (number, how many bits to use), 
# Output: The number in binary format, with leading 0s such that length == bits
def encodeNumber(n, bits):
    return bin(n)[2:].rjust(bits, "0")

# Input: binary string * the number of bits used for number, output (position, rest)
def decodeNumber(binaryString, bits):
    b = binaryString[:bits]
    return (int(b, 2), binaryString[bits:])

###########################################################################

# Translation functions to and from binary

def encodeInsertions(insertions):
    # TODO FIX this is a hyperparameter based on an upper length of 2^15 = 32000 (suitable for COVID)
    bits = 15
    eN = lambda n: encodeNumber(n, bits)
    s = ""
    if len(insertions) > 0:
        s += "1"
        s += eN(len(insertions)) # Number of insertions
        for (position, insert) in insertions:
            s += eN(position)
            s += eN(len(insert))
            # Note: This next step assumes that the insert is composed of ACTG. I just realised that because I called getInsertions before I did getOthers,
            # this may not necessarily be true.
            # FIX / TODO
            # insertInBinary = "".join(map(lambda x: encodeChar(x), list(insert))) 
            # This was the original efficient version, using 2 bits per character. I changed it below to standard ASCII for ease. 
            insertInBinary = "".join(map(lambda x: encodeNumber(ord(x), 8), list(insert))) # 
            s += insertInBinary
    else:
        s += "0"
    return s

def encodeDeletions(deletions):
    # TODO FIX this is a hyperparameter based on an upper length of 2^15 = 32000 (suitable for COVID)
    bits = 15
    eN = lambda n: encodeNumber(n, bits)
    s = ""
    if len(deletions) > 0:
        s += "1"
        s += eN(len(deletions)) # Number of deletions
        for (position, length) in deletions:
            s += eN(position)
            s += eN(length)
    else:
        s += "0"
    return s

def encodeNs(ns):
    # TODO FIX this is a hyperparameter based on an upper length of 2^15 = 32000 (suitable for COVID)
    bits = 15
    eN = lambda n: encodeNumber(n, bits)
    s = ""
    if len(ns) > 0:
        s += "1"
        s += eN(len(ns)) # Number of insertions
        for (position, length) in ns:
            s += eN(position)
            s += eN(length)
    else:
        s += "0"
    return s

def encodeOthers(others):
    # TODO FIX this is a hyperparameter based on an upper length of 2^15 = 32000 (suitable for COVID)
    bits = 15
    eN = lambda n: encodeNumber(n, bits)
    s = ""
    if len(others) > 0:
        s += "1"
        s += eN(len(others)) # Number of insertions
        for (position, c) in others:
            s += eN(position)
            s += encodeNumber(ord(c), 8) # i.e. encode Y in ascii binary (can do this more efficiently but don't think it matters)
    else:
        s += "0"
    return s

def encodeMismatches(mismatches):
    # TODO FIX this is a hyperparameter based on an upper length of 2^15 = 32000 (suitable for COVID)
    bits = 15
    eN = lambda n: encodeNumber(n, bits)
    s = ""
    if len(mismatches) > 0:
        s += "1"
        s += eN(len(mismatches)) # Number of insertions
        for (position, c) in mismatches:
            s += eN(position)
            s += encodeChar(c)
    else:
        s += "0"
    return s

# Variables used to profile encode, i.e. see what is actually response for the
# majority of the bits.
# TODO: Haven't implemented this yet. Add in the encodeInfo section to keep a
# global counter of how many bits we have used for insertions, deletions, etc. 
INSERTION_BITS = 0
DELETION_BITS = 0
N_BITS = 0
OTHER_BITS = 0
MISMATCHE_BITS = 0

# We are going to do this super simply to establish a baseline. 
def encodeInfo(info):
    # TODO FIX this is a hyperparameter based on an upper length of 2^15 = 32000 (suitable for COVID)
    bits = 15
    eN = lambda n: encodeNumber(n, bits)
    # Begin encoding. 
    s = ""
    insertions, deletions, ns, others, mismatches = info
    # Insertions
    s += encodeInsertions(insertions)
    # Deletions
    s += encodeDeletions(deletions)
    # Ns
    s += encodeNs(ns)
    # Others
    s += encodeOthers(others)
    # Mismatches
    s += encodeMismatches(mismatches)
    # Done
    return s

# TODO: Should code the inverse for encodeInfo i.e. a function to take a binary input: str, and convert it
# into an info object. 

###########################################################################

# "101010111" => bytes
# First we pad the string with a single byte to indicate how many padding bits we have, and then 000s, and then the number. 
# (you can't have a partial number of bytes). 
# E.g. "11111" => "00000011 000 11111"
def convertToBinary(s):
    # Add the padding at the front
    paddingBits = (8-len(s)%8)%8
    padding = encodeNumber(paddingBits, 8)+"0"*paddingBits
    s = padding+s
    # Convert to bytes
    starts = range(len(s)//8)
    starts = map(lambda x: x*8, starts)
    bites = map(lambda x: s[x:x+8], starts)
    bites = map(lambda x: int(x, 2), bites)
    bites = bytes(bites)
    return bites

assert(convertToBinary("11111") == bytes([int("00000011", 2), int("00011111", 2)]))

# Inverse of convertToBinary. Removes the initial padding, and 
# then produces the original string. 
def convertFromBinary(b: bytes):
    s = list(b)
    s = map(lambda x: encodeNumber(x, 8), s)
    s = "".join(s)
    # Remove the initial padding
    paddingBits = int(s[:8], 2)
    s = s[8+paddingBits:]
    return s

assert(convertFromBinary(convertToBinary("11111")) == "11111")

###########################################################################

# We use the return convention the 0 = fine, -1 is error. 

# Accepts a reference, a dna string, and a file identifier. First it checks if
# the file has already been created. If it has, do nothing. Otherwise, create
# the file, and store the info representation of the binary string. 
def writeToFile(path, onezeroString):
    if os.path.exists(path):
        return 0
    else:
        bites = convertToBinary(onezeroString)
        try:
            with open(path, "wb") as binary_file:
                binary_file.write(bites)
            return 0
        except:
            return -1

###########################################################################

# We use the return convention the 0 = fine, -1 is error. 

# Accepts a reference, a dna string, and a file identifier. First it checks if
# the file has already been created. If it has, do nothing. Otherwise, create
# the file, and store the info representation of the binary string. 
def storeAsBinary(reference, dna, path):
    if os.path.exists(path):
        return 0
    else:
        info = getInfoFromString(reference, dna)
        onezeroString = encodeInfo(info)
        return writeToFile(path, onezeroString)

###########################################################################











# # Applied to the coronavirus multiple alignment
# data = ""
# with open("HCOV19-ENGLAND-081220-A.fasta", "r") as f:
#     data = f.read()
# originalSize = len(data)

# print(f"\n\n\nOriginal coronavirus file size (bytes): {originalSize}\n")

# data = data.split(">")
# data = data[1:]

# for i in range(len(data)):
#     data[i] = data[i].split("8")
#     data[i] = data[i][-1]
#     data[i] = "".join(data[i].split())
# print(f"There are {len(data)} genomes.\n")

# POSITION_BITS = len(bin(len(data))[2:])

# # TODO: There were 5 strange letters, K: 10, M: 5, Y: 22, R: 12, S:4, W:1, which appears in the multiple alignment. 
# # I don't know what these are. For now, I just replaced all of these by "A"s. 

# for i in range(len(data)):
#     s = list(data[i])
#     for j in range(len(s)):
#         if s[j] not in ["A", "C", "T", "G", "-", "N"]:
#             s[j] = "A"
#     s = "".join(s)
#     data[i] = s

# coronaInfo = createInfo(data)
# coronaCompressed = encodeInfo(coronaInfo)
# newSize = int(len(coronaCompressed)/8) # We didn't compress our "0101010" chars into binary bytes, but this is easy

# print(f"Compressed coronavirus file size (bytes): {newSize}\n") # We didn't compress our chars into binary, but this is easy

# print(f"Percentage reduction: {100*(1-round(newSize/originalSize, 3))}\n")


# # TODO REMOVE Tempoarary experiments

# # lengths = []
# # for i in range(len(data)):
# #     lengths.append(len(encodeInfo(createInfo(data[:i+1])))/8)
# # print(lengths)

# # """
# # lengths = [22429.5, 27858.0, 29450.5, 30187.25, 31138.375, 31761.375, 32394.0, 32806.125, 32842.0, 32875.25, 33280.375, 33657.5, 33816.75, 33868.375, 33907.75, 33942.75, 33998.0, 34033.0, 34082.0, 34094.25, 34103.875, 34109.125, 34116.125, 34157.25, 34258.875, 34264.125, 34266.75, 34272.0, 34281.625, 34290.0, 34296.125, 34341.875, 34380.375, 34383.0, 34384.75, 34390.875, 34396.125, 34534.375, 34538.75, 34647.5, 34737.5, 34770.75, 34909.375, 35002.125, 35066.0, 35130.75, 35182.375, 35242.75, 35300.375, 35363.375, 35538.0, 35599.375, 35661.5, 35720.125, 35923.125, 36054.75, 36107.75, 36164.625, 36322.25, 36380.125, 36517.0, 36571.25, 36625.5, 36791.625, 36855.5, 36915.125, 36960.25, 36986.5, 37018.0, 37048.625, 37084.5, 37106.0, 37183.875, 37213.625, 38708.125, 39446.625, 40182.125, 41851.625, 42835.125, 43996.25, 44498.5, 44801.75, 45212.125, 46653.5, 48104.75, 49169.625, 50200.375, 50752.625, 51242.25, 51896.375, 52383.875, 52788.75, 53070.5, 53488.5, 53839.875, 54190.5, 54564.125, 54919.375]
# # differences = [22429.5, 5428.5, 1592.5, 736.75, 951.125, 623.0, 632.625, 412.125, 35.875, 33.25, 405.125, 377.125, 159.25, 51.625, 39.375, 35.0, 55.25, 35.0, 49.0, 12.25, 9.625, 5.25, 7.0, 41.125, 101.625, 5.25, 2.625, 5.25, 9.625, 8.375, 6.125, 45.75, 38.5, 2.625, 1.75, 6.125, 5.25, 138.25, 4.375, 108.75, 90.0, 33.25, 138.625, 92.75, 63.875, 64.75, 51.625, 60.375, 57.625, 63.0, 174.625, 61.375, 62.125, 58.625, 203.0, 131.625, 53.0, 56.875, 157.625, 57.875, 136.875, 54.25, 54.25, 166.125, 63.875, 59.625, 45.125, 26.25, 31.5, 30.625, 35.875, 21.5, 77.875, 29.75, 1494.5, 738.5, 735.5, 1669.5, 983.5, 1161.125, 502.25, 303.25, 410.375, 1441.375, 1451.25, 1064.875, 1030.75, 552.25, 489.625, 654.125, 487.5, 404.875, 281.75, 418.0, 351.375, 350.625, 373.625, 355.25]
# # """

# print("All tests passed\n\n")