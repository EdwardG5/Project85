from params import NUM_LEN_BITS, ENDIANNESS, K
from params import getTypeSize, removeNonACGT
import rle
import numpy as np
from sys import argv

# bwrleCompress : string -> bytes
# s may only contain A, C, G, T.
def bwrleCompress(s):
    class Loc(object):
        def __init__(self, loc):
            self.loc = loc
        def __eq__(self, other):        # for finding zeroloc
            return self.loc == other
        def __lt__(self, other):        # for sorting
            if self.loc == len(s): return True
            if other.loc == len(s): return False
            offset = 0
            while s[self.loc + offset] == s[other.loc + offset]:
                offset += 1
                if self.loc + offset == len(s):   # self has $ first
                    return True
                if other.loc + offset == len(s):  # other has $ first
                    return False
            return s[self.loc + offset] < s[other.loc + offset]
        # def __repr__(self):
        #     return 'L' + str(self.loc)

    locs = list(map(Loc, list(range(len(s) + 1))))
    locs = sorted(locs)     # locs can be thought of as a representation
                            # of the sorted M array.

    # sufarr backwards maps locations in the original
    # to locations in the bwt.
    sufarrBackwards = dict.fromkeys(range(0, len(locs) - 1, K))
    for i, loc in enumerate(locs):
        if loc.loc % K == 0:
            sufarrBackwards[loc.loc] = i


    del locs[locs.index(0)]
    
    bwt = ''.join([s[locs[i].loc - 1] for i in range(len(s))])
    preFooter = rle.rleEncode(bwt)

    # add the footer
    bpv = ((len(locs) - 1).bit_length() + 7) // 8    # bytes per value
                                                     # for some footer values, 
                                                     # based on largest
                                                     # number to store

    sufarrbytes = bytearray()
    for i in range(0, len(sufarrBackwards)):
        sufarrbytes += int.to_bytes(sufarrBackwards[i * 5], bpv, ENDIANNESS)
    salenbytes = int.to_bytes(len(sufarrBackwards), bpv, ENDIANNESS)
    bpvbyte = int.to_bytes(bpv, 1, ENDIANNESS)
    kbyte = int.to_bytes(K, 1, ENDIANNESS)

    postFooter = preFooter + sufarrbytes + salenbytes + \
                 bpvbyte + kbyte
    return postFooter

# bwrleDecompress : bytes -> string
# b is a bytes object that was returned by bwrleCompress,
# using the same value for NUM_LEN_BITS.
def bwrleDecompress(b):
    # extract necessary information from footer
    # and isolate raw BWT encoding
    assert(K == b[-1])  # K should be the same as global K
    bpv = b[-2]
    salen = int.from_bytes(b[-2 - bpv:-2], ENDIANNESS)
    sastart = -2 - bpv - (bpv * salen)
    zeroloc = int.from_bytes(b[sastart:sastart + bpv], ENDIANNESS)
    withoutFooter = b[:sastart]
    bwt = rle.rleDecode(withoutFooter)

    # assemble array representing the BWT
    # collect total count information, occurrences information,
    # and insert $ character when necessary
    lastCol = np.empty(len(bwt) + 1, dtype='<U1')

    reqBytes = getTypeSize(((len(lastCol) - 2).bit_length() + 7) // 8)
    typeStr = '<u' + str(reqBytes)
    occurrences = np.empty(len(lastCol), dtype=typeStr)
    
    counts = dict.fromkeys('ACTG', 0)
    i = 0
    for char in bwt:
        if i == zeroloc:
            lastCol[i] = '$'
            occurrences[i] = 0
            i += 1
        lastCol[i] = char
        occurrences[i] = counts[char]
        counts[char] += 1
        i += 1
    if i != len(lastCol):
        lastCol[i] = '$'
        occurrences[i] = 0

    prevCounts = 0
    firstOccurrence = {}
    for char in 'ACG':
        firstOccurrence[char] = prevCounts + 1
        prevCounts += counts[char]
    firstOccurrence['T'] = prevCounts + 1

    i = 0
    currChar = lastCol[i]
    res = ''
    while currChar != '$':
        res = currChar + res
        i = firstOccurrence[currChar] + occurrences[i]
        currChar = lastCol[i]
    
    return res

# tests on n random strings
def testRandom(n):
    for _ in range(n):
        test = ''.join(np.random.choice(list('ACTG'), np.random.randint(3, 71)))
        try:
            if (test != bwrleDecompress(bwrleCompress(test))):
                print(test)
        except Exception:
            print(test)
            assert(False)

if __name__ == '__main__':
    if len(argv) != 4 or (argv[1] != '-c' and argv[1] != '-d'):
        print('Usage: python bwrleFF.py -c/-d ' +
              '[source filename] [destination filename]')
    elif argv[1] == '-c':
        print('Compressing ' + argv[2] + '...')
        with open(argv[2]) as f:
            s = f.read()
        newS, nlct, loct, othct = removeNonACGT(s)
        if nlct or loct or othct:
            print('Made '+str(loct)+' lowercase nucleotides uppercase, '+
                  'and removed '+str(nlct + othct)+' non-ACGT characters, '+
                  str(nlct)+' of which are newlines and '+str(othct)+' of '+
                  'which are not.')
        with open(argv[3], 'wb') as f:
            f.write(bwrleCompress(newS))
        print('Successful.')
    else:
        print('Decompressing ' + argv[2] + '...')
        with open(argv[2], 'rb') as f:
            b = f.read()
        with open(argv[3], 'w') as f:
            f.write(bwrleDecompress(b))
        print('Successful.')