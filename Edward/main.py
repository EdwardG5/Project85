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

    if N == 1:
        print(f"Will compress: {N} genome. ", end="")
    else:
        print(f"Will compress: {N} genomes. ", end="")

    if PARALLEL:
        print(f"In parallel, using {cpu_count()} processors.")
    else:
        print(f"Sequentially.")

    p("\nInitializing state...\n")

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

    # Begin compression
    p("Beginning compression...\n")
    start = time.time()
    if PARALLEL:
        with get_context("spawn").Pool() as pool:
            returnValues = pool.map(transform, toCompress)
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
        
    else:
        # print("Sequential code needs to be fixed, nothing will be done.")
        returnValues = list(map(lambda x: storeAsBinary(reference.seq, x[0], compressed_file_location+x[1]+".bin"), toCompress))
    end = time.time()
    elapsedTime = end - start

    # Report back
    if N == 1:
        print(f"{N} genome was succesfully compressed.")
    else:
        print(f"{N} genomes were succesfully compressed.")
    print(f"Total time: {round(elapsedTime, 2)}. Time per genome: {round(elapsedTime/N, 2)}")
    
    p()


