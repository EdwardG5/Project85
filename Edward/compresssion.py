from math import ceil
from collections import Counter

# This is a parameter used to determine how many bits to use to encode positions. 
# This will be determined by how many genomes you have. 
POSITION_BITS = 7

# Type 1
def variableLengthNumberEncoding(n):
    # This is a hyper parameter - adjust as needed
    # Most columns have <= 8 errors, so we set t = 4 (1 bit to indicate whether to look at the next one, and 3 to encode a number)
    t = 4
    b = bin(n)[2:]
    b = b.rjust(ceil(len(b)/3)*3, "0")
    out = ""
    i = 0
    for i in range(len(b)//3-1):
        out += "1"+b[0:3]
        b = b[3:]
    out += "0"+b
    return out
VLNE = variableLengthNumberEncoding

# Test cases
assert(VLNE(0) == "0000")
assert(VLNE(1) == "0001")
assert(VLNE(2) == "0010")
assert(VLNE(3) == "0011")
assert(VLNE(4) == "0100")
assert(VLNE(5) == "0101")
assert(VLNE(6) == "0110")
assert(VLNE(7) == "0111")
assert(VLNE(8) == "10010000")
assert(VLNE(9) == "10010001")
assert(VLNE(10) == "10010010")
assert(VLNE(11) == "10010011")
assert(VLNE(12) == "10010100")
assert(VLNE(13) == "10010101")
assert(VLNE(14) == "10010110")
assert(VLNE(15) == "10010111")
assert(VLNE(16) == "10100000")

# Assumes that the first 4 bits are of type 1
# Returns the number represented, as well as the binary string left over
def decodeVLNE(binaryString):
    n = ""
    while True:
        new = binaryString[:4]
        n += new[1:]
        binaryString = binaryString[4:]
        if new[0] == "0":
            break
    n = int(n, 2)
    return (n, binaryString)

# Test cases
assert(decodeVLNE(VLNE(0))[0] == 0)
assert(decodeVLNE(VLNE(1))[0] == 1)
assert(decodeVLNE(VLNE(2))[0] == 2)
assert(decodeVLNE(VLNE(3))[0] == 3)
assert(decodeVLNE(VLNE(4))[0] == 4)
assert(decodeVLNE(VLNE(5))[0] == 5)
assert(decodeVLNE(VLNE(6))[0] == 6)
assert(decodeVLNE(VLNE(7))[0] == 7)
assert(decodeVLNE(VLNE(8))[0] == 8)
assert(decodeVLNE(VLNE(9))[0] == 9)
assert(decodeVLNE(VLNE(10))[0] == 10)
assert(decodeVLNE(VLNE(11))[0] == 11)
assert(decodeVLNE(VLNE(12))[0] == 12)
assert(decodeVLNE(VLNE(13))[0] == 13)
assert(decodeVLNE(VLNE(14))[0] == 14)
assert(decodeVLNE(VLNE(15))[0] == 15)
assert(decodeVLNE(VLNE(16))[0] == 16)


charDict = {"A": "000", "C": "001", "T":"010", "G":"011", "N":"100", "-":"101"}
inverseDict = {v: k for k, v in charDict.items()}

# char -> binaryString
def encodeChar(c):
    return charDict[c]

# binaryString -> (char, rest of the string)
def decodeChar(binaryString):
    b = binaryString[:3]
    return (inverseDict[b], binaryString[3:])

# Test cases
# Unecessary

# Input: number * how many bits to use, output in binary format
def encodeNumber(n, bits):
    return bin(n)[2:].rjust(bits, "0")

# Input: binary string * the number of bits used for number, output (position, rest)
def decodeNumber(binaryString, bits):
    b = binaryString[:bits]
    return (int(b, 2), binaryString[bits:])

"""
We encode a multiple alignment in the following way:

position 0: (
                - consensusCharacter, 
                - # of different error characters appear (3 bits), 
                    - char1: # of positions (VLNE), position1, position2, ...,
                    - char2: ....
                    ...
            )
position 1: ...
...
position n: ...

Example: 

[
    ("A", 2, [("C", 3, [5, 12, 79]), ("G", 5, [2, 3, 80, 91, 92])] ), 
    ...
]

"""

def encodeInfo(info):
    binaryString = ""
    for entry in info:
        new = ""
        new += encodeChar(entry[0]) # The consensus character
        new += encodeNumber(entry[1], 3) # The number of different errors e.g. len("A", "-", "G")
        for triplet in entry[2]:
            new += encodeChar(triplet[0]) # The error character e.g. "G"
            new += VLNE(triplet[1]) # How many positions
            for j in range(triplet[1]):
                new += encodeNumber(triplet[2][j], POSITION_BITS) # Each position
        binaryString += new
    return binaryString

# Test Case
"""Consensus:
 AAC""" # Not really, but it lets us test multiple faulty characters

"AAC"
"AGA"
"AGA"
"AGA"
"AG-"
"AG-"
"AG-"
"AG-"
"AG-"
infoExample = [
        ("A", 0, []),
        ("A", 1, [("G", 8, [1, 2, 3, 4, 5, 6, 7, 8])]),
        ("C", 2, [("A", 3, [1, 2, 3]), ("-", 5, [4, 5, 6, 7, 8])])
]
"000 000"
"000 001 011 10010000 0000001 0000010 0000011 0000100 0000101 0000110 0000111 0001000"
"001 010 000 0011 0000001 0000010 0000011 101 0101 0000100 0000101 0000110 0000111 0001000"

assert(encodeInfo(infoExample) == "00000000000101110010000000000100000100000011000010000001010000110000011100010000010100000011000000100000100000011101010100001000000101000011000001110001000")

def decodeInfoDict(binaryString):
    info = []
    while len(binaryString) > 0:
        c, binaryString = decodeChar(binaryString)
        numberOfErrors, binaryString = decodeNumber(binaryString, 3)
        errors = []
        for i in range(numberOfErrors):
            errorChar, binaryString = decodeChar(binaryString)
            numberOfPositions, binaryString = decodeVLNE(binaryString)
            positions = []
            for i in range(numberOfPositions):
                position, binaryString = decodeNumber(binaryString, POSITION_BITS)
                positions.append(position)
            errors.append((errorChar, numberOfPositions, positions))
        info.append((c, numberOfErrors, errors))
    return info

# Test cases
assert(decodeInfoDict(encodeInfo(infoExample)) == infoExample)


###########################################################################

# Finds the consensus string of a multiple alignment
# list[str] -> str
# TODO: For whatever reason this function takes a while.
# I'm relatively sure it's linear time. I think the issue is that I loop through the rows -> columns, giving it bad cache locality. 
# Try doing column -> row instead. 
def consensus(data):
    consensus = ""
    for i in range(len(data[0])):
        characters = list(map(lambda x: x[i], data))
        c = Counter(characters)
        mode = c.most_common(1)[0][0]
        consensus += mode
    return consensus

# list[str] -> info (as above)
def createInfo(data):
    cS = consensus(data)
    info = []
    for i in range(len(cS)):
        c = cS[i]
        disagrees = {}
        for x in range(len(data)):
            if data[x][i] != c:
                if data[x][i] in disagrees:
                    disagrees[data[x][i]].append(x)
                else:
                    disagrees[data[x][i]] = [x]
        errors = []
        for k, v in disagrees.items():
            errors.append( (k, len(v), v ) )
        info.append( (c, len(errors), errors) )
    return info

# Test cases
data = [   
    "AAC",
    "AGA",
    "AGA",
    "AGA",
    "AG-",
    "AG-",
    "AG-",
    "AG-",
    "AG-"
]
infoCorrect = [
    ("A", 0, []),
    ("G", 1, [("A", 1, [0])]),
    ("-", 2, [("C", 1, [0]), ("A", 3, [1, 2, 3])])
]
assert(createInfo(data) == infoCorrect)


###########################################################################

# Applied to the coronavirus multiple alignment
data = ""
with open("HCOV19-ENGLAND-081220-A.fasta", "r") as f:
    data = f.read()
originalSize = len(data)

print(f"\n\n\nOriginal coronavirus file size (bytes): {originalSize}\n")

data = data.split(">")
data = data[1:]

for i in range(len(data)):
    data[i] = data[i].split("8")
    data[i] = data[i][-1]
    data[i] = "".join(data[i].split())
print(f"There are {len(data)} genomes.\n")

POSITION_BITS = len(bin(len(data))[2:])

# TODO: There were 5 strange letters, K: 10, M: 5, Y: 22, R: 12, S:4, W:1, which appears in the multiple alignment. 
# I don't know what these are. For now, I just replaced all of these by "A"s. 

for i in range(len(data)):
    s = list(data[i])
    for j in range(len(s)):
        if s[j] not in ["A", "C", "T", "G", "-", "N"]:
            s[j] = "A"
    s = "".join(s)
    data[i] = s

coronaInfo = createInfo(data)
coronaCompressed = encodeInfo(coronaInfo)
newSize = int(len(coronaCompressed)/8) # We didn't compress our "0101010" chars into binary bytes, but this is easy

print(f"Compressed coronavirus file size (bytes): {newSize}\n") # We didn't compress our chars into binary, but this is easy

print(f"Percentage reduction: {100*(1-round(newSize/originalSize, 3))}\n")





# TODO REMOVE Tempoarary experiments

# lengths = []
# for i in range(len(data)):
#     lengths.append(len(encodeInfo(createInfo(data[:i+1])))/8)
# print(lengths)

# """
# lengths = [22429.5, 27858.0, 29450.5, 30187.25, 31138.375, 31761.375, 32394.0, 32806.125, 32842.0, 32875.25, 33280.375, 33657.5, 33816.75, 33868.375, 33907.75, 33942.75, 33998.0, 34033.0, 34082.0, 34094.25, 34103.875, 34109.125, 34116.125, 34157.25, 34258.875, 34264.125, 34266.75, 34272.0, 34281.625, 34290.0, 34296.125, 34341.875, 34380.375, 34383.0, 34384.75, 34390.875, 34396.125, 34534.375, 34538.75, 34647.5, 34737.5, 34770.75, 34909.375, 35002.125, 35066.0, 35130.75, 35182.375, 35242.75, 35300.375, 35363.375, 35538.0, 35599.375, 35661.5, 35720.125, 35923.125, 36054.75, 36107.75, 36164.625, 36322.25, 36380.125, 36517.0, 36571.25, 36625.5, 36791.625, 36855.5, 36915.125, 36960.25, 36986.5, 37018.0, 37048.625, 37084.5, 37106.0, 37183.875, 37213.625, 38708.125, 39446.625, 40182.125, 41851.625, 42835.125, 43996.25, 44498.5, 44801.75, 45212.125, 46653.5, 48104.75, 49169.625, 50200.375, 50752.625, 51242.25, 51896.375, 52383.875, 52788.75, 53070.5, 53488.5, 53839.875, 54190.5, 54564.125, 54919.375]
# differences = [22429.5, 5428.5, 1592.5, 736.75, 951.125, 623.0, 632.625, 412.125, 35.875, 33.25, 405.125, 377.125, 159.25, 51.625, 39.375, 35.0, 55.25, 35.0, 49.0, 12.25, 9.625, 5.25, 7.0, 41.125, 101.625, 5.25, 2.625, 5.25, 9.625, 8.375, 6.125, 45.75, 38.5, 2.625, 1.75, 6.125, 5.25, 138.25, 4.375, 108.75, 90.0, 33.25, 138.625, 92.75, 63.875, 64.75, 51.625, 60.375, 57.625, 63.0, 174.625, 61.375, 62.125, 58.625, 203.0, 131.625, 53.0, 56.875, 157.625, 57.875, 136.875, 54.25, 54.25, 166.125, 63.875, 59.625, 45.125, 26.25, 31.5, 30.625, 35.875, 21.5, 77.875, 29.75, 1494.5, 738.5, 735.5, 1669.5, 983.5, 1161.125, 502.25, 303.25, 410.375, 1441.375, 1451.25, 1064.875, 1030.75, 552.25, 489.625, 654.125, 487.5, 404.875, 281.75, 418.0, 351.375, 350.625, 373.625, 355.25]
# """










print("All tests passed\n\n")

