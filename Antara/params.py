
ENDIANNESS = 'big'  # default endianness for storing integers
NUM_LEN_BITS = 1    # number of bits used to store run lengths
K = 10              # checkpointed suffix array parameter
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

def getFirstOccurrence(totCounts):
    firstOccurrence = {}
    sofar = 1
    firstOccurrence['A'] = sofar
    sofar += totCounts[0]
    for charcode in range(1,4):
        firstOccurrence[itoc[charcode]] = sofar
        sofar += totCounts[charcode]
    return firstOccurrence

def removeNonACGT(s):
    nlCount = 0
    lowerCount = 0
    otherCount = 0
    toRemove = []; lowercase = []

    # collect info in one pass
    for i in range(len(s)):
        if s[i] not in 'ACGT':
            if s[i] in 'acgt':
                lowerCount += 1
                lowercase.append(i)
            elif s[i] == '\n':
                nlCount += 1
                toRemove.append(i)
            else:
                otherCount += 1
                toRemove.append(i)

    # construct result string in one pass
    res = ''
    start = 0
    while toRemove and lowercase:
        if toRemove[0] < lowercase[0]:
            ri = toRemove.pop(0)
            res += s[start:ri]
            start = ri + 1
        else:
            lc = lowercase.pop(0)
            res += s[start:lc] + s[lc].upper()
            start = lc + 1
    for ri in toRemove:
        res += s[start:ri]
        start = ri + 1
    for lc in lowercase:
        res += s[start:lc] + str(s[lc].upper())
        start = lc + 1

    # add the rest
    res += s[start:]

    return res, nlCount, lowerCount, otherCount



###############################################
# from here onward, don't include in paper
###############################################

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

if __name__ == '__main__':
    import bwrleFF as b
    import matches as m
    import time
    import numpy as np
    b.compress('sars.txt', 'scmpsd')
    m.matchRunner('scmpsd', 'pat.txt', 'mtchs.txt')

    # times = np.zeros((4,4))
    # testvals = [10,100,1000,5000]

    # for i,k in enumerate(testvals):
    #     for j,c in enumerate(testvals):
    #         m.setk(k)
    #         b.setk(k)
    #         m.setc(c)

    #         print()
    #         print('k =',k,'c =',c)
    #         b.compress('sars.txt', 'scmpsd')
    #         start = time.time()
    #         m.matchRunner('scmpsd', 'pat.txt', 'mtchs.txt')
    #         end = time.time()

    #         times[i,j] = end - start

    # with open('testingresults.txt', 'w') as f:
    #     f.write('Times to match. Rows are K, cols are C.\n')
    #     f.write('Test values are ' + str(testvals) + '.\n\n')
    #     f.write(str(times))