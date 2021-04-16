from sys import argv

# This function implements Burrows-Wheeler with run length encoding.
#
# Detailed explanation of compression mechanism:
# Characters are represented by two bits using the encoding A = 00, C = 01,
# G = 10, T = 11. A single character (that is not part of a run) is stored
# as its two-bit code, followed by a 1. A run of a character is stored as the
# character's two-bit code, followed by a 0, followed by three bits
# representing the run length minus 2. Note that the largest run that can be
# stored in six bits is 9, so any runs longer than that will get split up.
# At the end of the file is a footer, which aids in decompression. The final
# byte of the footer contains the number of bytes in the rest of the footer.
# The rest of the footer contains the number representing the location of the
# character in the file (once all the runs are expanded, indexed by zero) that
# corresponds to the last character in the original data string.
#
# Inputs:
# - srcPath: Path of the source file to compress. ***The file
#            must contain only A, C, T, G, or \n as of now.***
# - resPath: The path of the file in which to write the compressed data.
# Outputs: none
def bwrleCompress(srcPath, resPath):
    with open(srcPath) as f:
        everything = f.read().replace('\n', '')
    
    n = len(everything)
    
    class Loc(object):
        def __init__(self, loc):
            self.loc = loc
        def __eq__(self, other):
            return self.loc == other
        def __lt__(self, other):                # only comparison that
                                                # is required for sorting
            if self.loc == other.loc: return False
            offset = 0
            while everything[self.loc + offset] == \
                  everything[other.loc + offset]:
                offset += 1
                if (self.loc + offset == n):    # self contains $ first
                    return True
                if (other.loc + offset == n):   # other contains $ first
                    return False
            return everything[self.loc + offset] < \
                   everything[other.loc + offset]
        def __repr__(self):                     # for debugging purposes
            return 'L' + str(self.loc)
    
    locs = list(map(Loc, list(range(n))))
    locs.sort()     # locs can be thought of as a representation
                    # of the sorted M array.
    
    # loop through locs to write the main data
    # in increments of three bits
    with open(resPath, 'wb') as f:
        i = 0
        offset = 0
        buffer = bytearray(b'\x00\x00\x00')
        extra = bytearray(b'\x80')     # leftmost bit determines if invalid
                                       # can't use zero, since that represents
                                       # a run of length 2
        encoding = dict(zip('ACGT', range(4)))

        # quick helper for writing to the buffer in increments of 3 bytes
        # REQUIRES: offset < 24
        def writeBuffer(toWrite, offs, buffer):
            # print(toWrite, offs, buffer)
            if offs != 6 and offs != 15:
                    buffer[offs // 8] |= toWrite << (5 - offs % 8)
            else:
                buffer[offs // 8] |= toWrite >> (offs % 8 - 5)
                buffer[(offs // 8) + 1] |= (toWrite << (13 - offs % 8)) & 255
            # print(buffer)

        cur = everything[locs[0].loc - 1]
        nex = everything[locs[1].loc - 1]

        while (i < len(locs) - 2):
            if offset == 24: # buffer is full; clear buffer and
                             # write items in buffer to file
                f.write(buffer)
                buffer = bytearray(b'\x00\x00\x00')
                offset = 0

            if extra[0] != 128: # there is extra; write three bits in extra to buffer
                writeBuffer(extra[0], offset, buffer)
                offset += 3
                extra[0] = 128
                continue

            if cur != nex:  # single character
                writeBuffer((encoding[cur] << 1) | 1, offset, buffer)
                offset += 3

                i += 1
                cur = nex
                nex = everything[locs[i + 1].loc - 1]
            else:
                writeBuffer(encoding[cur] << 1, offset, buffer)
                offset += 3

                # found a run; loop here, incrementing i
                runLength = 1
                while (cur == nex and runLength < 9 and i < len(locs) - 2):
                    i += 1; cur = nex; nex = everything[locs[i + 1].loc - 1]
                    runLength += 1

                assert(runLength >= 2)  # loop must have gone at least once
                extra = bytearray([runLength - 2])

                # makes cur location of next to be read; could go over end
                i += 1; cur = nex; nex = everything[locs[i + 1].loc - 1] if i < len(locs) - 1 else None

        if (offset == 24):
            f.write(buffer)
            buffer = bytearray(b'\x00\x00\x00')
            offset = 0

        if (extra[0] != 128):
            writeBuffer(extra[0], offset, buffer)
            offset += 3
            if (offset == 24):
                f.write(buffer)
                buffer = bytearray(b'\x00\x00\x00')
                offset = 0

            # where we are depends on the value of i
            if i == len(locs) - 1:
                # final character is in cur
                writeBuffer((encoding[cur] << 1) | 1, offset, buffer)
                offset += 3
            else:
                # final characters in cur and nex
                writeBuffer((encoding[cur] << 1) | 1, offset, buffer)
                offset += 3
                if (offset == 24):
                    f.write(buffer)
                    buffer = bytearray(b'\x00\x00\x00')
                    offset = 0
                writeBuffer((encoding[nex] << 1) | 1, offset, buffer)
                offset += 3
        else:
            # final characters in cur and nex
            writeBuffer((encoding[cur] << 1) | 1, offset, buffer)
            offset += 3
            if (offset == 24):
                f.write(buffer)
                buffer = bytearray(b'\x00\x00\x00')
                offset = 0
            writeBuffer((encoding[nex] << 1) | 1, offset, buffer)
            offset += 3

        # write remaining in buffer to file; clear extra bits
        if offset < 8:
            f.write(bytes([buffer[0] & (-1 << (8 - offset))]))
        elif offset < 16:
            f.write(bytes([buffer[0]]))
            f.write(bytes([buffer[1] & (-1 << (16 - offset))]))
        else:
            f.write(bytes([buffer[0]]))
            f.write(bytes([buffer[1]]))
            f.write(bytes([buffer[2] & (-1 << (24 - offset))]))

        # write footer
        zeroloc = hex(locs.index(0)).lstrip('0x')       # !! strips 0 and x sep
        if len(zeroloc) % 2 != 0: zeroloc = '0' + zeroloc   # has to be bytes
        zerolocbytes = bytes.fromhex(zeroloc)
        numbytes = bytes([len(zerolocbytes)])   # this will throw a ValueError
                                                # if too many bytes to store
                                                # in a byte (>256)
                                                # This might be zero, if
                                                # zeroloc is zero
        f.write(zerolocbytes)
        f.write(numbytes)

def bwrleDecompress(srcPath, resPath):
    pass    # unimplemented
    # notes:
    # - when decompressing, read stuff into a hash dictionary
    # - index: (first, last) - can store in four bytes - can change later

def printProperties(srcPath):  # for debugging purposes
    with open(srcPath) as f:
        s = f.read()
    rotations = []
    for i in range(len(s)):
        rotations.append(s[i:] + s[:i])
    rotations.sort()
    for rotation in rotations: print(rotation)
    bwt = ''.join([rot[-1] for rot in rotations])
    print('bwt: ' + bwt)
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
    if len(argv) != 3:
        print('Usage: python bwrle.py [source file path] [destination file path]')
    else:
        src = argv[1]
        des = argv[2]
        bwrleCompress(src, des)
        printStats(src, des)