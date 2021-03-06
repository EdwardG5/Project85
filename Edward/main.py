import Bio
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import time
import timeit
from os import listdir, cpu_count
import argparse
import itertools
import os
from math import ceil
import notify

from compressionV2 import storeAsBinary, getInfoFromString, encodeInfo, writeToFile, getStringFromInfo
from multiprocessing import Pool
from multiprocessing import get_context

# Relative paths
genomesPath = "GenomicData/0-10k_genomes.fna"
referencePath = "GenomicData/COVID_Reference_Genome.fna"
compressed_file_location = "GenomicData/CompressedData/0-10/"

# Get reference
reference = next(SeqIO.parse(referencePath, "fasta"))

# toCompress = [ (0, ref, dna, path), (1, ref, dna, path), ...  ]
def transform(x):
    ID, ref, dna, path = x
    print(f"Starting compression of genome {ID}, {path[len(compressed_file_location):]}")
    try:
        info = getInfoFromString(ref, dna)
        onezeroString = encodeInfo(info)
        print(f"Finishing compression of {ID}, {path[len(compressed_file_location):]}")
        return (ID, path, onezeroString)
    except:
        print(f"Something went wrong with genome {ID}, {path[len(compressed_file_location):]}. Not storing.")
        return (ID, path, None)


p = print


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    # --N 5 --TYPE SEQUENTIAL
    parser.add_argument("--N", required=True, type=int, help="The number of additional genomes to compress.")
    parser.add_argument("--PARALLEL", action='store_const', const=True, default=False, help="A flag to indicate that we should run the program in parallel (multiprocesssing).")

    cmdline_args = parser.parse_args()
    N = cmdline_args.N
    PARALLEL = cmdline_args.PARALLEL

    p()

    # Count already compressed
    alreadyCompressed = listdir(compressed_file_location)
    alreadyCompressed = len(list(filter(lambda x: x[-4:] == ".bin", alreadyCompressed)))

    # Fetch genomes
    data = SeqIO.parse(genomesPath, "fasta")
    toCompress = list(itertools.islice(data, N+alreadyCompressed)) # This ensures safety if you give N > genomes available
    toCompress = map(lambda x: (x.seq, compressed_file_location+x.id+".bin"), toCompress)
    toCompress = filter(lambda x: not os.path.exists(x[1]), toCompress)
    toCompress = enumerate(toCompress)
    toCompress = map(lambda x: (x[0], str(reference.seq), x[1][0], x[1][1]), toCompress)
    toCompress = list(toCompress)
    # toCompress = [ (0, ref, dna, path), (1, ref, dna, path), ...  ]
    
    N = len(toCompress)

    if N == 1:
        print(f"Will compress: {N} genome. ", end="")
    else:
        print(f"Will compress: {N} genomes. ", end="")

    if PARALLEL:
        # This was actually changed, see later note. 
        print(f"In parallel, using {cpu_count()} processors.")
    else:
        print(f"Sequentially.")

    p("\nInitializing state...\n")

    # Begin compression
    p("Beginning compression...\n")
    start = time.time()
    if PARALLEL:
        # We compress and write the files in rounds of 1000. This means that if
        # we need to stop the script, we don't throw away massive amounts of
        # work. We also send an email notification to inform the user of
        # progress.
        batchSize = 1000
        for iteration in range(ceil(len(toCompress)/batchSize)):
            part = toCompress[(iteration)*batchSize:(iteration+1)*batchSize]
            # For some reason a deadlock seems to occur when you just let this
            # be the default ".Pool()" using all 40 processors. Unclear why. It
            # also only seem to have appeared on the second day of running the
            # program. Never managed to figure out what the issue is -
            # restricting the andrew machines to 30/40 of their cores seemed to
            # work.  
            with get_context("spawn").Pool(30) as pool:
                returnValues = pool.map(transform, part)
            print("\nBeginning storage of genomes, writing files...\n")
            for (ID, path, val) in returnValues:
                print(f"Writing {ID}, {path[len(compressed_file_location):]}...")
                if val != None:
                    r = writeToFile(path, val)
                if r == 0:
                    pass
                else:
                    print(f"Something went wrong while writing file {path[len(compressed_file_location):]}.")
            print("\nFinishing storage of genomes....\n")
            # Send an email notification.
            msg = notify.message(
                subject="Compression progress report",
                text=(f'Another {batchSize} genomes compressed.\n'
                    f'Only another {max(0, len(toCompress)-batchSize*(iteration+1))} to go!')
            )
            notify.send(msg)        
    else:
        # print("Sequential code needs to be fixed, nothing will be done.")
        returnValues = list(map(lambda x: storeAsBinary(x[1], x[2], x[3]), toCompress))
    end = time.time()
    elapsedTime = end - start

    # Send an email notification.
    msg = notify.message(
        subject="Compression progress report: COMPLETE",
        text=(f'All genomes compressed. Time for the next batch!\n')
    )
    notify.send(msg)

    # Report back
    if N == 1:
        print(f"{N} genome was succesfully compressed.")
    else:
        print(f"{N} genomes were succesfully compressed.")
    print(f"Total time: {round(elapsedTime, 2)}. Time per genome: {round(elapsedTime/N, 2)}")
    
    p()


