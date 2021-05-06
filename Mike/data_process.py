FOLDER_NAME = "genome_compression_datasets"

#####
#
# DATA BENCHMARKING
#
#####

def benchmark(function, dataset=None, sample=None, data_dir=None, check_atcg=False):
    import os
    import gzip
    from Bio import SeqIO
    import pandas as pd

    if dataset is None and data_dir is None:
        raise ValueError('Please specify data to compress.')
    elif dataset is not None and data_dir is not None:
        raise ValueError('Please specify either directory or dataset to compress, not both.')
    else:
        benchmark_results = []
        # Compress standard data
        from os.path import expanduser
        home = expanduser("~")
        dataset_fp = os.path.join(home, FOLDER_NAME, "datasets", f"dataset_{str(dataset)}") if data_dir is None else data_dir
        break_counter = 0
        for fasta in os.listdir(dataset_fp):
            if sample is not None and break_counter >= sample:
                break
            # Ignore DS_Store and other hidden files
            if fasta.endswith('fna.gz'):
                fasta_fp = os.path.join(dataset_fp, fasta)
                with gzip.open(fasta_fp, "rt") as handle:
                    # If want to consider filesize: print(os.fstat(handle.fileno()).st_size)
                    total_genome = ""
                    for record in SeqIO.parse(handle, "fasta"):
                        # sometimes genomes stored across records
                        total_genome += record.seq
                    # Only consider ACTG
                    if not check_atcg or len(set(list(total_genome)).difference(set(['A','C','T','G']))) == 0:
                        original = len(total_genome)
                        try:
                            compressed_str = function(total_genome)
                        except Exception as e:
                            raise Exception(f"Failed to run compression function on input genome {fasta}.")
                        new = len(compressed_str)
                        benchmark_results.append([fasta, original, new])
                    else:
                        extra_chars = list(set(list(total_genome)).difference(set(['A','C','T','G'])))
                        print(f"Did not run {fasta} as contains characters {str()}")
                break_counter += 1
        messy_df = pd.DataFrame(benchmark_results, columns=['genome', 'original_len', 'new_len'])
        agg_df = messy_df.groupby('genome').sum()
        agg_df['pct_compression'] = (agg_df['original_len'] - agg_df['new_len'])/agg_df['original_len'] * 100
    
        try:
            df = pd.read_csv(os.path.join(dataset_fp, f"dataset_{str(dataset)}.tsv"), header=None, delimiter='\t')
            def find_group(fn):
                return df[df[2] == fn.replace("_genomic.fna.gz", "")].iloc[0][0]

            agg_df['group'] = agg_df.index.map(find_group)
        except:
            raise ValueError("Can't find dataset? Did you run the dataset function to completion?")
        
        return agg_df

def benchmark_functions(func_list, dataset=None, sample=None, data_dir=None, check_atcg=False, verbose=True):
    import os
    import gzip
    from Bio import SeqIO
    import pandas as pd
    # order preserved in python list
    df_list = []
    for f in func_list:
        function_name = f.__name__
        if verbose:
            print(f"Running function {function_name}")
        tdf = benchmark(f, dataset=dataset, sample=sample, data_dir=data_dir, check_atcg=check_atcg)
        tdf[f"{function_name}_pct_compression"] = tdf['pct_compression']
        df_list.append(tdf)
    if verbose:
        print(f"Collating results...")
    df = pd.concat(df_list, axis=1).reset_index()
    df = df.loc[:,~df.columns.duplicated()]
    return df[[c for c in df.columns if c.endswith('_pct_compression') or c == 'genome' or c == 'group']].set_index('genome')

#####
#
# DATA HANDLING
#
#####

def download_org(data_payload, retries=0):
    import requests
    import os
    org = data_payload[0]
    retry = data_payload[1]
    dataset_folder = data_payload[2]
    try:
        group, name, genome = org.split("\t")

        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{group}/{name}/all_assembly_versions/{genome}/{genome}_genomic.fna.gz"

        r = requests.get(url, allow_redirects=True)

        fn = os.path.join(dataset_folder, f"{genome}_genomic.fna.gz")
        open(fn, 'wb').write(r.content)
    except Exception as e:
        print(e)
        if retries < retry:
            print(f"Failed to download {org}. Retrying...")
            download_org(org, retries+1)
    return True

import gzip 
import shutil
import os
import fnmatch

# gunzip and gunzip recurse taken and modified from 
# https://stackoverflow.com/questions/45596742/unzip-gz-files-within-folders-in-a-main-folder-using-python

# Accepts a file path to a .gz file
# Decompresses the file, and replaces it with the decompressed file, named the same (minus .gz)
def gunzip(compressedFilePath):
    assert(compressedFilePath[-3:] == ".gz")
    outPath = compressedFilePath[:-3]
    with gzip.open(compressedFilePath,"rb") as f_in, open(outPath,"wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
        os.remove(compressedFilePath) # Delete the file

# Accepts a path to a directory
# Decompresses all .gz files in the directory.
def recurse_and_gunzip(root):
    walker = os.walk(root)
    for root,dirs,files in walker:
        for f in files:
            if fnmatch.fnmatch(f,"*.gz"):
                print(f"Decompressing {f}...")
                gunzip(root+"/"+f)

def download_dataset(dataset_dir, dest=None, procs=1, retry=1):
    import requests
    import os
    read_data = {}

    if dest is None:

        home = os.path.expanduser("~")
        download_folder = os.path.join(home, FOLDER_NAME)

        if not os.path.exists(download_folder):
            os.makedirs(download_folder)

        dataset_folder = os.path.join(
            download_folder, dataset_dir.replace(".tsv", ""))

        if not os.path.exists(dataset_folder):
            os.makedirs(dataset_folder)

    else:
        dataset_folder = dest

    with open(dataset_dir, "r") as f:

        print("Reading dataset...\n")
        organisms = [o for o in f.read().split("\n") if o != '']
        retry_generator = [retry for _ in range(len(organisms))]
        dataset_generator = [dataset_folder for _ in range(len(organisms))]
        from multiprocessing import Pool
        print("Beginning downloads...\n")
        pool = Pool(procs)
        pool.map(download_org, zip(
            organisms, retry_generator, dataset_generator))
        print("Downloading finished.\n")

    print("Beginning file decompression...\n")
    recurse_and_gunzip(dataset_folder)
    print("\nFile decompression complete.\n")

    from shutil import copyfile
    copyfile(dataset_dir, os.path.join(
        dataset_folder, os.path.basename(dataset_dir)))

    print("Finished.\n")