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

    locs = list(map(Loc, list(range(len(s) + 1))))
    locs = sorted(locs)     # locs can be thought of as a representation
                            # of the sorted M array.

    zeroloc = locs.index(0)
    del locs[zeroloc]
    
    bwt = ''.join([s[locs[i].loc - 1] for i in range(len(s))])
    preFooter = rle.rleEncode(bwt)

    zls = hex(zeroloc).lstrip('0x')         # !! strips 0 and x sep
    if len(zls) % 2 != 0: zls = '0' + zls   # has to be bytes
    zerolocbytes = bytes.fromhex(zls)
    numbytes = bytes([len(zerolocbytes)])   # this will throw a ValueError
                                            # if too many bytes to store
                                            # in a byte (>256)
                                            # numbytes might be zero, if
                                            # zeroloc is zero

    postFooter = preFooter + zerolocbytes + numbytes
    return postFooter

# bwrleDecompress : bytes -> string
# b is a bytes object that was returned by bwrleCompress,
# using the same value for NUM_LEN_BITS.
def bwrleDecompress(b):
    zeroloc = int.from_bytes(b[-1 - (b[-1]):-1], 'big')
    withoutFooter = b[:-1 - (b[-1])]
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

def randomTests(n):
    # tests on a bunch of random strings
    for _ in range(n):
        test = ''.join(np.random.choice(list('ACTG'), np.random.randint(3, 71)))
        try:
            if (test != bwrleDecompress(bwrleCompress(test))): print(test)
        except Exception:
            print(test); assert(False)

if __name__ == '__main__':
    if len(argv) != 4 or (argv[1] != '-c' and argv[1] != '-d'):
        print('Usage: python bwrleCF.py -c/-d ' +
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