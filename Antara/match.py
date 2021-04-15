import numpy as np
import bwrle

NUM_LEN_BITS = 1    # number of bits used to store run lengths

# quick helper
# requires: 0 <= byteInt < 256, 0 <= start <= end <= 8
# start and end indexed by zero from the left of the byte
# represented by byteInt
def getDigits(byteInt, start, end):
    byteInt >>= 8 - end
    return byteInt & ~(-1 << (end - start))

# match : (bytes * string) -> int
#                             int list
# start with returning an int representing number of patterns,
# upgrade by storing checkpointed suffix array data and
# using results to find locations of patterns

def match(b, pattern, C, bwt):
    ###########################################
    # assemble auxiliary data
    ###########################################

    zeroloc = int.from_bytes(b[-1 - b[-1]:-1], 'big')
    raw = b[:-1 - b[-1]] # ; print(raw.hex())
    itoc = dict(zip(range(4), 'ACGT'))
    ctoi = dict(zip('ACGT', range(4)))

    end =  8 * len(raw)
    finalByte = int.from_bytes(raw[-1:], 'big')
    p = 1
    while (finalByte & (~(-1 << p)) == 0):
        p += 1
    end -= p

    counts = np.zeros((1, 4), dtype=int)    # TODO can tailor size to str len
    runningCts = np.zeros((1, 4), dtype=int)
    cOffsets = dict()
    countPtrs = []                          # integer pointers into data; length should be len(counts) - 1

    # helper
    # p is start location, length is number of bits
    def getCodeAt(p, length, raw=raw):
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

    loc = 0     # loc of char that starts at p
    p = 0       # points to the beginning of ch

    chcode = getDigits(raw[0], 0, 2)
    ch = itoc[chcode]
    nxcode = getDigits(raw[0], 2, 4)
    nx = itoc[nxcode]
    
    # assemble counts
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

    # assemble firstOccurrence
    firstOccurrence = {}
    sofar = 1
    firstOccurrence['A'] = sofar
    sofar += counts[-1,0]
    for charcode in range(1,4):
        firstOccurrence[itoc[charcode]] = sofar
        sofar += counts[-1,charcode]

    # helper function to get the counts at a specific index
    # REQUIRES: index <= loc
    def getCounts(index, counts=counts, C=C):
        if index == loc: return counts[-1]

        # index is compatible with counts

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

    ###########################################
    # search for pattern occurrences
    ###########################################

    # init top & bottom pointers
    char = pattern[-1]
    pattern = pattern[:-1]

    top = firstOccurrence[char]
    bottom = top + counts[-1,ctoi[char]] - 1

    while len(pattern) > 0:
        char = pattern[-1]
        pattern = pattern[:-1]

        firstoc = getCounts(top)[ctoi[char]] + 1
        lastoc = getCounts(bottom + 1)[ctoi[char]]

        if firstoc > lastoc:    # no matches
            return 0
        
        top = firstOccurrence[char] + firstoc
        bottom = firstOccurrence[char] + lastoc

    numMatches = bottom - top + 1
    return numMatches

pattern = 'CA'
C = 4

test = 'TCATAATCAGG'
bwt = bwrle.getBWT(test)
# print('test:', test, '\n bwt:', bwt, '\npattern:', pattern, '\nC:', C)
print('test:', test, '\npattern:', pattern, '\nC:', C)
print('returned: ', match(bwrle.bwrleCompress(test), pattern, C, bwt))