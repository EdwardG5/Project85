# Modifiable byte array object, for help in rleEncode
class ModByteArray(object):
    def __init__(self):
        self.bytesArr = bytearray()     # holds completed bytes
        self.currByte = 0               # represents next byte to add
        self.offset = 0                 # can be 0-7 inclusive

    # requires 0 < numBits <= 8
    def pushInt(self, toAddInt, numBits):
        if numBits < 8 - self.offset:
            mask = ~(-1 << numBits)
            self.currByte |= (toAddInt & mask) << (8 - numBits - self.offset)
            self.offset += numBits
        elif numBits == 8 - self.offset:
            mask = ~(-1 << numBits)
            self.currByte |= toAddInt & mask

            # reset
            self.bytesArr += bytes([self.currByte])
            self.currByte = 0
            self.offset = 0
        elif numBits <= 8:
            initOffset = self.offset
            self.currByte |= toAddInt >> (numBits - 8 + initOffset)

            # reset
            self.bytesArr += bytes([self.currByte])
            self.currByte = 0
            self.offset = 0

            mask = ~(-1 << (numBits - 8 + initOffset))
            self.currByte |= (toAddInt & mask) << (16 - numBits - initOffset)
            self.offset += numBits - 8 + initOffset
        else:
            # untested
            self.pushBytes(toAddInt.to_bytes((toAddInt.bit_length() + 7) // 8, 'big'), numBits)

    # numBits is from right of toAdd; can be anything,
    # so long as toAdd is long enough.
    def pushBytes(self, toAdd, numBits):
        if numBits == 0: return
        elif numBits < 8 - self.offset:
            asInt = int.from_bytes(toAdd[-1:], 'big')
            mask = ~(-1 << numBits)
            self.currByte |= (asInt & mask) << (8 - numBits - self.offset)
            self.offset += numBits
        elif numBits == 8 - self.offset:
            asInt = int.from_bytes(toAdd[-1:], 'big')
            mask = ~(-1 << numBits)
            self.currByte |= asInt & mask

            # reset
            self.bytesArr += bytes([self.currByte])
            self.currByte = 0
            self.offset = 0
        elif numBits <= 8:
            initOffset = self.offset
            asInt = int.from_bytes(toAdd[-1:], 'big')
            self.currByte |= asInt >> (numBits - 8 + initOffset)

            # reset
            self.bytesArr += bytes([self.currByte])
            self.currByte = 0
            self.offset = 0

            mask = ~(-1 << (numBits - 8 + initOffset))
            self.currByte |= (asInt & mask) << (16 - numBits - initOffset)
            self.offset += numBits - 8 + initOffset
        else:
            start = len(toAdd) - 1 - (numBits // 8)
            startByte = toAdd[start:start+1]
            self.pushBytes(startByte, numBits % 8)
            restBytes = toAdd[start+1:]
            for byte in restBytes:
                self.pushBytes(bytes([byte]), 8)

    # get a bytearray object representing content.
    def get(self):
        if self.offset == 0: return self.bytesArr
        else: return self.bytesArr + bytes([self.currByte])

# quick helper
# requires: 0 <= byteInt < 256, 0 <= start <= end <= 8
# start and end indexed by zero from the left of the byte
# represented by byteInt
def getDigits(byteInt, start, end):
    byteInt >>= 8 - end
    return byteInt & ~(-1 << (end - start))

# rleEncode : string * int -> bytes
# 
# Inputs:
#   s            The input string. ***Must contain only A, C, G, and T, and
#                must be of length greater than 2.***
#   numLenBits   The number of bits used to represent run lengths.
#                ***Must range from 1 to 8 (inclusive).***
def rleEncode(s, numLenBits):
    assert(len(s) > 2)
    assert(1 <= numLenBits <= 8)
    
    encoding = dict(zip('ACGT', range(4)))
    maxRunLength = 2**numLenBits + 1
    
    res = ModByteArray()
    prev = s[0]
    res.pushInt(encoding[prev], 2)
    i = 1
    curr = s[i]
    while i < len(s) - 1:
        res.pushInt(encoding[curr], 2)

        if curr == prev: # found a run
            runLength = 1
            while curr == prev and i < len(s) - 1 and runLength < maxRunLength:
                i += 1
                prev = curr
                curr = s[i]
                runLength += 1

            if i == len(s) - 1:     # finish it off
                if curr != prev:
                    res.pushInt(runLength - 2, numLenBits)
                    res.pushInt(encoding[curr], 2)
                    res.pushInt(1, 1)
                    return bytes(res.get())
                elif runLength < maxRunLength:
                    runLength += 1
                    res.pushInt(runLength - 2, numLenBits)
                    res.pushInt(1, 1)
                    return bytes(res.get())
                else:
                    # assert(runLength == maxRunLength)
                    res.pushInt(runLength - 2, numLenBits)
                    res.pushInt(encoding[curr], 2)
                    res.pushInt(1, 1)
                    return bytes(res.get())
            elif runLength == maxRunLength:
                res.pushInt(runLength - 2, numLenBits)
                res.pushInt(encoding[curr], 2)
                if i == len(s) - 2 and s[-1] == s[-2]:
                    # maxed out run, but there's two
                    # consecutive characters at the end
                    res.pushInt(encoding[s[-1]], 2)
                    res.pushInt(0, numLenBits)
                    res.pushInt(1, 1)
                    return bytes(res.get())
            else:   # curr != prev
                res.pushInt(runLength - 2, numLenBits)
                res.pushInt(encoding[curr], 2)

        i += 1
        prev = curr
        curr = s[i]

    res.pushInt(encoding[curr], 2)
    if curr == prev and prev != s[i - 2]:
        # ends with run of length 2
        # never entered the inner run loop.
        res.pushInt(0, numLenBits)
    res.pushInt(1, 1)       # always ends with a meaningless 1
    
    return bytes(res.get())

# rleDecode : bytes * int -> string
#  
# Inputs:
#   b            The input bytes object. Must have no zero padding at end,
#                in other words, b[-1] != 0. Should be returned from rleEncode.
#   numLenBits   The number of bits used to represent run lengths. ***Must be
#                the same as what was passed into rleEncode, or this will not
#                work. Must range from 1 to 8 (inclusive).***
# b : bytes
# requires that input was returned from rle
# has no zero padding at the end (final byte != 0)
def rleDecode(b, numLenBits):
    encoding = dict(zip(range(4), 'ACGT'))
    end =  8 * len(b)
    finalByte = int.from_bytes(b[-1:], 'big')
    p = 1
    while (finalByte & (~(-1 << p)) == 0):
        p += 1
    end -= p

    res = ''
    p = 2; prev = '\0'; ch = '\0'; justRan = False
    while p <= end:
        i = p // 8
        j = p % 8
        if j == 0: i -= 1; j = 8
        if j == 1:  # p might be split across two bytes
            firstBit = getDigits(b[i-1], 7, 8)
            secondBit = getDigits(b[i], 0, 1)
            ch = encoding[(firstBit << 1) + secondBit]
        else:
            ch = encoding[getDigits(b[i], j - 2, j)]
        
        res += ch

        if ch == prev: # run
            p += numLenBits
            i = p // 8
            j = p % 8
            if j == 0: i -= 1; j = 8
            if j < numLenBits:  # might be split across two bytes
                firstSeg = getDigits(b[i-1], 8 + j - numLenBits, 8)
                secondSeg = getDigits(b[i], 0, j)
                rl = (firstSeg << j) + secondSeg
            else:
                rl = getDigits(b[i], j - numLenBits, j)
            res += rl * ch
            justRan = True
        
        prev = ch if not justRan else '\0'
        p += 2
        justRan = False
    return res

if __name__ == '__main__':
    # tests on a bunch of random strings
    for lb in range(1, 9):
        import random as r
        for i in range(1000):
            test = ''.join([r.choice('ACGT') for j in range(r.randrange(3,20))])
            if not (test == rleDecode(rleEncode(test, lb), lb)):
                print(test)

    # test = 'CCCGG'
    # cmpsd = rleEncode(test, 1)
    # print(cmpsd.hex())
    # print(rleDecode(cmpsd, 1))