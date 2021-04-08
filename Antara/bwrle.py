import rle
import numpy as np

NUM_LEN_BITS = 1    # number of bits used to store run lengths

# bwrleCompress : string -> bytes
# s may only contain A, C, G, T.
def bwrleCompress(s):
    class Loc(object):
        def __init__(self, loc):
            self.loc = loc
        def __eq__(self, other):        # for finding zeroloc
            return self.loc == other
        def __lt__(self, other):        # for sorting
            # if self.loc == other.loc: return False  # never happens? TODO
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
        # def __repr__(self):                     # for debugging purposes
        #     return 'L' + str(self.loc)

    locs = list(map(Loc, list(range(len(s) + 1))))
    locs = sorted(locs)     # locs can be thought of as a representation
                            # of the sorted M array.
    zeroloc = locs.index(0)
    del locs[zeroloc]
    
    bwt = ''.join([s[locs[i].loc - 1] for i in range(len(s))])
    preFooter = rle.rleEncode(bwt, NUM_LEN_BITS)

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

# helper function
def annotate(L):
    counts = dict.fromkeys('ACGT', 0)
    result = []
    for item in L:
        result.append((item, counts[item]))
        counts[item] += 1
    return result

# bwrleDecompress : bytes -> string
# b is a bytes object that was returned by bwrleCompress,
# using the same value for NUM_LEN_BITS.
def bwrleDecompress(b):
    zeroloc = int.from_bytes(b[-1 - (b[-1]):-1], 'big')
    withoutFooter = b[:-1 - (b[-1])]
    bwtNoEOF = rle.rleDecode(withoutFooter, NUM_LEN_BITS)

    # assemble array representing the BWT
    # collect total count information, occurrences information,
    # and insert $ character when necessary
    lastCol = np.empty(len(bwtNoEOF) + 1, dtype='<U1')

    reqBytes = ((len(lastCol) - 2).bit_length() + 7) // 8
    typeStr = '<u' + str(reqBytes)
    occurrences = np.empty(len(lastCol), dtype=typeStr)
    
    counts = dict.fromkeys('ACTG', 0)
    i = 0
    for char in bwtNoEOF:
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

# for debugging purposes
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

def printStats(srcPath, resPath):
    with open(srcPath, 'rb') as f:
        sl = len(f.read())
    with open(resPath, 'rb') as f:
        rl = len(f.read())
    
    print('Source length (bytes):\t' + str(sl))
    print('Compressed length:\t' + str(rl))
    print('Compression ratio:\t' + ('%.2f' % (100 * (1 - rl/sl))) + '%')

if __name__ == '__main__':
    # tests on a bunch of random strings
    for i in range(100):
        test = ''.join(np.random.choice(list('ACTG'), np.random.randint(3, 71)))
        try:
            if (test != bwrleDecompress(bwrleCompress(test))):
                print(test)
        except Exception:
            print(test)
            assert(False)

    # test = 'AAAAA'
    # printBWProperties(test)
    # cmpsd = bwrleCompress(test)
    # print(cmpsd.hex())
    # print(bwrleDecompress(cmpsd))

    # with open('SC2genome.txt') as f:
    #     content = f.read()
    #     cmpsd = bwrleCompress(content)
    # with open('SC2cmpsd', 'wb') as f:
    #     f.write(cmpsd)
    # with open('SC2decmpsd.txt', 'w') as f:
    #     f.write(bwrleDecompress(cmpsd))
    # printStats('SC2genome.txt', 'SC2cmpsd')