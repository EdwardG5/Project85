from params import ENDIANNESS, NUM_LEN_BITS, K, C, itoc, ctoi
from params import getDigits, getFirstOccurrence
from params import printBWProperties
import numpy as np
import bwrleFF

# extracts info from compressed data.
# takes bytearray, returns info tuple
def extractInfo(b):
    # extract necessary information from footer
    # and isolate raw BWT encoding
    assert(K == b[-1])  # K should be the same as global K
    bpv = b[-2]
    salen = int.from_bytes(b[-2 - bpv:-2], ENDIANNESS)
    sastart = -2 - bpv - (bpv * salen)

    sufarr = dict()
    for i, start in enumerate(range(sastart, -2 - bpv, bpv)):
        sufarr[int.from_bytes(b[start:start + bpv], ENDIANNESS)] = i * K

    zeroloc = int.from_bytes(b[sastart:sastart+bpv], ENDIANNESS)
    raw = b[:sastart]

    end =  8 * len(raw)
    finalByte = int.from_bytes(raw[-1:], 'big')
    p = 1
    while (finalByte & (~(-1 << p)) == 0):
        p += 1
    end -= p

    # p is start location, length is number of bits
    def getCodeAt(p, length):
        i = (p + length) // 8
        j = (p + length) % 8
        if j == 0: i -= 1; j = 8
        if j < length:  # might be split across two bytes
            firstSeg = getDigits(raw[i-1], 8 + j - length, 8)
            secondSeg = getDigits(raw[i], 0, j)
            code = (firstSeg << j) + secondSeg
        else:
            code = getDigits(raw[i], j - length, j)
        return code

    return zeroloc, end, getCodeAt, sufarr

def assembleCounts(zeroloc, end, getCodeAt):
    counts = np.zeros((1, 4), dtype=int)    # TODO can tailor size to str len
    runningCts = np.zeros((1, 4), dtype=int)
    cOffsets = dict()
    countPtrs = []                          # integer pointers into data;
                                            # length should be len(counts) - 1

    loc = 0     # loc of char that starts at p
    p = 0       # points to the beginning of ch

    chcode = getCodeAt(0, 2)
    ch = itoc[chcode]
    nxcode = getCodeAt(2, 2)
    nx = itoc[nxcode]
    
    while True:
        # when we get here, p points to either
        # - beginning of singleton character, or
        # - beginning of run

        if ch != nx:    # ch is singleton
            if loc == zeroloc:
                if loc % C == 0:
                    counts = np.vstack((counts, runningCts))
                    countPtrs.append(p)
                loc += 1
            if loc % C == 0:
                counts = np.vstack((counts, runningCts))
                countPtrs.append(p)

            runningCts[0,chcode] += 1
            loc += 1
            p += 2
            if p + 2 == end:    # ends with singleton, singleton
                # update for final character; break
                if loc == zeroloc:
                    if loc % C == 0:
                        counts = np.vstack((counts, runningCts))
                        countPtrs.append(p)
                    loc += 1
                if loc % C == 0:
                    counts = np.vstack((counts, runningCts))
                    countPtrs.append(p)

                chcode = getCodeAt(p, 2)
                runningCts[0,chcode] += 1
                loc += 1    # loc is len
                break
            chcode = getCodeAt(p, 2)
            nxcode = getCodeAt(p + 2, 2)
            ch = itoc[chcode]
            nx = itoc[nxcode]
            continue
        else:           # ch at start of run
            runLength = getCodeAt(p + 4, NUM_LEN_BITS) + 2

            startLoc = loc
            while loc < startLoc + runLength:
                if loc == zeroloc:
                    if loc % C == 0:
                        counts = np.vstack((counts, runningCts))
                        countPtrs.append(p)
                        cOffsets[loc] = loc - startLoc
                    loc += 1
                    startLoc += 1   # incrementing loc messes loop guard up
                if loc % C == 0:
                    counts = np.vstack((counts, runningCts))
                    countPtrs.append(p)
                    cOffsets[loc] = loc - startLoc

                runningCts[0,chcode] += 1
                loc += 1

            # want to update ch, nx, p, loc to point to char that
            # immediately follows the run
            p += 4 + NUM_LEN_BITS
            if p == end: break  # ends with run, loc is already len
            if p + 2 == end:    # ends with run, singleton
                # update for final character; break
                if loc == zeroloc:
                    if loc % C == 0:
                        counts = np.vstack((counts, runningCts))
                        countPtrs.append(p)
                    loc += 1
                if loc % C == 0:
                    counts = np.vstack((counts, runningCts))
                    countPtrs.append(p)

                chcode = getCodeAt(p, 2)
                runningCts[0,chcode] += 1
                loc += 1    # loc is len
                break
            chcode = getCodeAt(p, 2)
            nxcode = getCodeAt(p + 2, 2)
            ch = itoc[chcode]
            nx = itoc[nxcode]
            continue

    # case where zeroloc is at the end
    if loc == zeroloc:
        if loc % C == 0:
            counts = np.vstack((counts, runningCts))
            countPtrs.append(end)
        loc += 1    # now loc is total len, including $

    # final row of counts table
    counts = np.vstack((counts[1:], runningCts))

    return counts, cOffsets, countPtrs, loc

# match : (bytes * string) -> int list
# returns list of starting locations of pattern
# in the string represented by b.
def match(b, pattern):
    zeroloc, end, getCodeAt, sufarr = extractInfo(bytearray(b))
    
    ###########################################
    # assemble auxiliary data
    ###########################################

    counts, cOffsets, countPtrs, bwlen = \
        assembleCounts(zeroloc, end, getCodeAt)

    firstOccurrence = getFirstOccurrence(counts[-1])

    # helper function to get the counts at a specific index
    # REQUIRES: index <= bwlen
    def getCounts(index):
        if index == bwlen: return counts[-1]

        # either index corresponds to a row of counts, or there's an offset
        if index % C == 0: return counts[index // C]

        baseRow = index // C
        currLoc = baseRow * C
        indexOffset = index % C     # always nonzero
        res = counts[baseRow].copy()

        p = countPtrs[baseRow]

        # either p starts at singleton, or p starts at run, or p is final char
        if p + 2 == end:
            # edge case; p points to the final character
            if currLoc == zeroloc:
                currLoc += 1
                indexOffset -= 1
                if indexOffset == 0: return res
            return counts[-1]
        
        # either p starts at singleton, or p starts at run
        chcode = getCodeAt(p, 2)
        nxcode = getCodeAt(p + 2, 2)
        ch = itoc[chcode]
        nx = itoc[nxcode]

        if ch == nx:    # p starts at run; currLoc, res, and indexOffset are wrong
                        # need to make them correct, using cof
            cof = cOffsets[currLoc]
            res[chcode] -= cof

            if currLoc == zeroloc:
                currLoc += 1
                indexOffset -= 1
            currLoc -= cof
            indexOffset += cof
            if currLoc <= zeroloc <= currLoc + cof:
                currLoc -= 1
                # res[chcode] += 1
                indexOffset += 1

        # just count up indexOffset times, then return res
        while indexOffset > 0:
            if ch != nx:
                if currLoc == zeroloc:
                    currLoc += 1
                    indexOffset -= 1
                    if indexOffset == 0: return res

                res[chcode] += 1

                p += 2
                if p + 2 == end:    # ends with singleton, singleton
                    indexOffset -= 1
                    currLoc += 1
                    chcode = getCodeAt(p, 2)
                    ch = itoc[chcode]

                    if currLoc == zeroloc:
                        currLoc += 1
                        indexOffset -= 1

                    if indexOffset <= 0: return res
                    else:
                        # $ is at end
                        res[chcode] += 1
                        return res

                chcode = getCodeAt(p, 2)
                nxcode = getCodeAt(p + 2, 2)
                ch = itoc[chcode]
                nx = itoc[nxcode]

                currLoc += 1
                indexOffset -= 1
                continue
            else:   # a run is starting
                runLength = getCodeAt(p + 4, NUM_LEN_BITS) + 2

                zlhere = False
                if currLoc <= zeroloc < currLoc + runLength:
                    # zeroloc is in the run (this will count before, but not after)
                    # one of the positions in the run corresponds to zeroloc, rather
                    zlhere = True
                    if indexOffset < runLength + 1:
                        zerolocOffset = zeroloc - currLoc
                        if indexOffset <= zerolocOffset:
                            res[chcode] += indexOffset
                            return res
                        else:
                            res[chcode] += indexOffset - 1
                            return res
                    
                # either index is within run, or it's after it
                if indexOffset < runLength:
                    res[chcode] += indexOffset
                    return res
                else:
                    res[chcode] += runLength
                    indexOffset -= runLength + int(zlhere)
                    currLoc += runLength + int(zlhere)
                    p += 4 + NUM_LEN_BITS

                    if p + 2 == end:    # ends with run, singleton
                        if indexOffset == 0:
                            return res

                        chcode = getCodeAt(p, 2)
                        ch = itoc[chcode]

                        if currLoc == zeroloc:
                            currLoc += 1
                            indexOffset -= 1

                        if indexOffset <= 0: return res
                        else:
                            # $ is at end
                            res[chcode] += 1
                            return res

                    if p == end:    # ends with run
                        return counts[-1]

                    chcode = getCodeAt(p, 2)
                    nxcode = getCodeAt(p + 2, 2)
                    ch = itoc[chcode]
                    nx = itoc[nxcode]
                    continue
        
        return res

    # helper function to get the character of BWT at a specific index
    # REQUIRES: index <= bwlen
    def getCharAt(index):
        if index == zeroloc: return '$'
        
        if index % C == 0:
            p = countPtrs[index // C]
            return itoc[getCodeAt(p, 2)]
        
        # otherwise...
        p = countPtrs[index // C]
        indexOffset = index % C         # guaranteed nonzero
        currLoc = index - indexOffset

        chcode = getCodeAt(p, 2)
        ch = itoc[chcode]

        if p + 2 == end: return ch

        nxcode = getCodeAt(p + 2, 2)
        nx = itoc[nxcode]

        if ch == nx:        # if p starts at run, need to adjust based on cof
            cof = cOffsets[currLoc]
            currLoc -= cof
            indexOffset += cof

            if currLoc <= zeroloc < currLoc + cof:
                currLoc -= 1
                indexOffset += 1

        # just count up indexOffset times
        while indexOffset > 0:
            if ch != nx:    # p at singleton
                if currLoc == zeroloc:
                    currLoc += 1
                    indexOffset -= 1
                    if indexOffset == 0: return ch

                indexOffset -= 1
                currLoc += 1

                p += 2
                chcode = nxcode
                ch = nx
                if p + 2 == end: return ch
                nxcode = getCodeAt(p + 2, 2)
                nx = itoc[nxcode]
                continue
            else:           # p at run
                runLength = getCodeAt(p + 4, NUM_LEN_BITS) + 2

                zlhere = False
                if currLoc <= zeroloc < currLoc + runLength:   # includes before the run, not after
                    zlhere = True

                if indexOffset < runLength + int(zlhere):
                    return ch
                
                # otherwise
                indexOffset -= (runLength + int(zlhere))
                currLoc += runLength + int(zlhere)

                p += 4 + NUM_LEN_BITS
                if p == end: # impossible
                    raise IndexError
                chcode = getCodeAt(p, 2)
                ch = itoc[chcode]
                if p + 2 == end:
                    return ch
                nxcode = getCodeAt(p + 2, 2)
                nx = itoc[nxcode]
                continue
        return ch

    ###########################################
    # search for pattern occurrence #
    ###########################################

    # init top & bottom pointers
    char = pattern[-1]
    pattern = pattern[:-1]

    top = firstOccurrence[char]
    bottom = top + counts[-1,ctoi[char]] - 1

    while len(pattern) > 0:
        char = pattern[-1]
        pattern = pattern[:-1]

        firstoc = getCounts(top)[ctoi[char]]
        lastoc = getCounts(bottom + 1)[ctoi[char]] - 1

        if firstoc > lastoc:    # no matches
            return []
        
        top = firstOccurrence[char] + firstoc
        bottom = firstOccurrence[char] + lastoc

    ###########################################
    # search for pattern locations
    ###########################################

    matches = []

    for bi in range(top, bottom + 1):
        if bi in sufarr:
            matches.append(sufarr[bi])
        else:
            offset = 1
            done = False
            while not done:
                ch = getCharAt(bi)
                occurrence = getCounts(bi)[ctoi[ch]]
                bi = firstOccurrence[ch] + occurrence

                if bi in sufarr:
                    matches.append(sufarr[bi] + offset); done = True

                offset += 1

    matches.sort()
    return matches

# test with n random strings and patterns
def testRandom(n):
    testStr = ''.join([np.random.choice(list('ACGT')) for _ in range(np.random.randint(10,20))])
    testPattern = ''.join([np.random.choice(list('ACGT')) for _ in range(np.random.randint(3,(len(testStr) - 1) // 2))])

    for _ in range(n):
    # if True:
        try:
            # bwt, sufarr = bwrleFF.getBWT(testStr)
            # bwt = printBWProperties(testStr)

            # print('test:', testStr, '\nbwt: ', bwt, '\npattern:', testPattern, '\nC:')
            # print()
            # print('returned: ', match(bwrleFF.bwrleCompress(testStr), testPattern))
            matches = match(bwrleFF.bwrleCompress(testStr), testPattern)

            for m in matches:
                if testStr[m:m + len(testPattern)] != testPattern: assert(False)

            testStr = ''.join([np.random.choice(list('ACGT')) for _ in range(np.random.randint(100,200))])
            testPattern = ''.join([np.random.choice(list('ACGT')) for _ in range(np.random.randint(3,(len(testStr) - 1) // 2))])

        except Exception:
            print(testStr, testPattern)
            assert(False)

if __name__ == '__main__':
    testRandom(100)