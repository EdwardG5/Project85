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

from compressionV2 import storeAsBinary, getInfoFromString
from multiprocessing import Pool

# Relative paths
genomesPath = "GenomicData/85k_COVID_Genomes.fna"
referencePath = "GenomicData/COVID_Reference_Genome.fna"
compressed_file_location = "GenomicData/CompressedData/"

# Get reference
reference = next(SeqIO.parse(referencePath, "fasta"))

def transform(x):
    return storeAsBinary(reference.seq, x[0], compressed_file_location+x[1]+".bin")

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
    alreadyCompressed = listdir("GenomicData/CompressedData")
    alreadyCompressed = len(list(filter(lambda x: x[-4:] == ".bin", alreadyCompressed)))
    
    # Fetch genomes
    data = SeqIO.parse(genomesPath, "fasta")
    toCompress = list(itertools.islice(data, N+alreadyCompressed)) # This ensure safety if you give N > genomes available
    N = len(toCompress)-alreadyCompressed
    toCompress = map(lambda x: (x.seq, x.id), toCompress)
    
    # Begin compression
    p("Beginning compression...\n")
    start = time.time()
    if PARALLEL:
        returnValues = Pool().map(transform, toCompress)
    else:
        returnValues = list(map(lambda x: storeAsBinary(reference.seq, x[0], compressed_file_location+x[1]+".bin"), toCompress))
    end = time.time()
    elapsedTime = end - start

    # Report back
    if sum(returnValues) == 0:
        if N == 1:
            print(f"{N} genome was succesfully compressed.")
        else:
            print(f"{N} genomes were succesfully compressed.")
        print(f"Total time: {round(elapsedTime, 2)}. Time per genome: {round(elapsedTime/N, 2)}")
    else:
        failed = enumerate(returnValues)
        failed = filter(lambda x: x[1] != 0, failed)
        failed = list(map(lambda x: x[0], failed))
        print(f"Something went wrong. Genomes {failed} failed.")
    p()


