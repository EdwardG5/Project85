
ENDIANNESS = 'big'  # default endianness for storing integers
NUM_LEN_BITS = 1    # number of bits used to store run lengths
K = 5               # checkpointed suffix array parameter
C = 10              # checkpointed counts array parameter

# encoding helpers
itoc = dict(zip(range(4), 'ACGT'))
ctoi = dict(zip('ACGT', range(4)))

# requires: 0 <= byteInt < 256, 0 <= start <= end <= 8
# start and end indexed by zero from the left of the byte
# represented by byteInt
def getDigits(byteInt, start, end):
    byteInt >>= 8 - end
    return byteInt & ~(-1 << (end - start))

# gets smallest type size that can fit i bytes
def getTypeSize(i):
    res = 1
    while res < i: res *= 2
    if res >= 16: raise OverflowError
    return res

def getBWT(s):
    s = s + '$'
    rotations = []
    for i in range(len(s)):
        rotations.append(s[i:] + s[:i])
    rotations.sort()
    bwt = ''.join([rot[-1] for rot in rotations])

    return bwt

def getBWTSA(s):
    s = s + '$'
    rotations = []
    for i in range(len(s)):
        rotations.append(s[i:] + s[:i])
    rotations.sort()
    bwt = ''.join([rot[-1] for rot in rotations])

    s = s + s
    sufarr = dict()
    for i, rot in enumerate(rotations):
        loc = s.index(rot)
        if loc % K == 0:
            sufarr[i] = loc

    return bwt, sufarr

def printBWProperties(s):
    s = s + '$'
    rotations = []
    for i in range(len(s)):
        rotations.append(s[i:] + s[:i])
    rotations.sort()
    bwt = ''.join([rot[-1] for rot in rotations])
    for rotation in rotations: print(rotation)
    print()
    print('test: ' + s)
    print('bwt:  ' + bwt)
    print()

    return bwt

def printBWSAProperties(s):
    s = s + '$'
    rotations = []
    for i in range(len(s)):
        rotations.append(s[i:] + s[:i])
    rotations.sort()
    bwt = ''.join([rot[-1] for rot in rotations])
    for rotation in rotations: print(rotation)
    print()
    print('test: ' + s)
    print('bwt:  ' + bwt)
    print()

    s = s + s
    sufarr = dict()
    for i, rot in enumerate(rotations):
        loc = s.index(rot)
        if loc % K == 0:
            sufarr[i] = loc
    
    print(sufarr)

    return bwt, sufarr

def printStats(srcPath, resPath):
    with open(srcPath, 'rb') as f:
        sl = len(f.read())
    with open(resPath, 'rb') as f:
        rl = len(f.read())
    
    print('Source length (bytes):\t' + str(sl))
    print('Compressed length:\t' + str(rl))
    print('Compression ratio:\t' + ('%.2f' % (100 * (1 - rl/sl))) + '%')

def getFirstOccurrence(totCounts):
    firstOccurrence = {}
    sofar = 1
    firstOccurrence['A'] = sofar
    sofar += totCounts[0]
    for charcode in range(1,4):
        firstOccurrence[itoc[charcode]] = sofar
        sofar += totCounts[charcode]
    return firstOccurrence