import rle

NUM_LEN_BITS = 1    # number of bits used to store run lengths

# s may only contain A, C, G, T
def bwrleCompress(s):
    class Loc(object):
        def __init__(self, loc):
            self.loc = loc
        def __eq__(self, other):
            return self.loc == other
        def __lt__(self, other):                # only comparison that
                                                # is required for sorting
            if self.loc == other.loc: return False
            offset = 0
            while s[(self.loc + offset) % len(s)] == \
                  s[(other.loc + offset) % len(s)]:
                  offset += 1
                  if offset == len(s): return False
            return s[(self.loc + offset) % len(s)] < \
                   s[(other.loc + offset) % len(s)]
        # def __repr__(self):                     # for debugging purposes
        #     return 'L' + str(self.loc)

    locs = list(map(Loc, list(range(len(s)))))
    locs = sorted(locs)     # locs can be thought of as a representation
                            # of the sorted M array.
    
    bwt = ''.join([s[locs[i].loc - 1] for i in range(len(s))])
    preFooter = rle.rleEncode(bwt, NUM_LEN_BITS)

    zeroloc = hex(locs.index(0)).lstrip('0x')       # !! strips 0 and x sep
    if len(zeroloc) % 2 != 0: zeroloc = '0' + zeroloc   # has to be bytes
    zerolocbytes = bytes.fromhex(zeroloc)
    numbytes = bytes([len(zerolocbytes)])   # this will throw a ValueError
                                            # if too many bytes to store
                                            # in a byte (>256)
                                            # This might be zero, if
                                            # zeroloc is zero
    
    postFooter = preFooter + zerolocbytes + numbytes
    return postFooter


def annotate(L):
    counts = dict.fromkeys('ACGT', 0)
    result = []
    for item in L:
        result.append((item, counts[item]))
        counts[item] += 1
    return result


def bwrleDecompress(b):
    zeroloc = int.from_bytes(b[-1 - (b[-1]):-1], 'big')
    
    withoutFooter = b[:-1 - (b[-1])]
    bwt = rle.rleDecode(withoutFooter, NUM_LEN_BITS)
    lastCol = [char for char in bwt]

    # assemble helper hash dictionary
    firstCol = lastCol.copy(); firstCol.sort()
    firstToIndex = dict()

    currChar = 'A'
    charCount = 0
    for i, char in enumerate(firstCol):
        if char == currChar:
            firstToIndex[(currChar, charCount)] = i
            charCount += 1
        else:
            currChar = char
            charCount = 0
            firstToIndex[(currChar, charCount)] = i
            charCount += 1
    
    lastCol = annotate(lastCol)

    original = ''
    i = zeroloc
    original = lastCol[i][0] + original
    i = firstToIndex[lastCol[i]]
    while i != zeroloc:
        original = lastCol[i][0] + original
        i = firstToIndex[lastCol[i]]
    
    return original

# for debugging purposes
def printBWProperties(s):
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

def demonstrateIssues():
    # # tests on a bunch of random strings, shows ones that don't work
    doesntwork = []
    import random as r
    for i in range(10000):
        test = ''.join([r.choice('ACGT') for j in range(r.randrange(3, 70))])
        if (test != bwrleDecompress(bwrleCompress(test))):
            doesntwork.append(test)
    print(doesntwork)

    # This is because of an idiosyncracy with the first-last property
    # when a string is just multiple repeats of the same substring, and
    # there is no explicit end character.
    # There will be some cyclic rotations that are the same.
    # They are not guaranteed to be sorted in any particular order.
    # Because of this, the first-last property is not guaranteed to work.

    # It is unlikely that very long strings will have this property,
    # since 

    # TODO it wouldn't be too hard to update this so that
    # - zeroloc represents position to insert $ (you'd need to change
    #                                            sorting mechanism to
    #                                            take that into account)
    # - decoder inserts that, and decoding proceeds as normal with a few changes 

def printStats(srcPath, resPath):
    with open(srcPath, 'rb') as f:
        sl = len(f.read())
    with open(resPath, 'rb') as f:
        rl = len(f.read())
    
    print('Source length (bytes):\t' + str(sl))
    print('Compressed length:\t' + str(rl))
    print('Compression ratio:\t' + ('%.2f' % (100 * (1 - rl/sl))) + '%')

if __name__ == '__main__':
    with open('SARS-CoV-2genome.txt') as f:
        content = f.read().replace('\n', '')
    with open('SARScompressed', 'wb') as f:
        f.write(bwrleCompress(content))

    printStats('SARS-CoV-2genome.txt', 'SARScompressed')