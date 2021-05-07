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
    check_chars = set(['A','C','T','G','N','-'])
    
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
            try:
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
                        if not check_atcg or len(set(list(total_genome)).difference(check_chars)) == 0:
                            original = len(total_genome)
                            try:
                                compressed_str = function(total_genome)
                            except Exception as e:
                                raise Exception(f"Failed to run compression function on input genome {fasta}.")
                            new = len(compressed_str)
                            benchmark_results.append([fasta, original, new])
                        else:
                            break_counter -= 1
                            print(f"Did not run {fasta} as contains bad characters")
                    break_counter += 1
            except:
                print(f"Did not run {fasta} because bad.")
                pass
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
    function_names = []
    df_list = []
    for f in func_list:
        function_name = f.__name__
        function_names.append(function_name)
        if verbose:
            print(f"Running function {function_name}")
        tdf = benchmark(f, dataset=dataset, sample=sample, data_dir=data_dir, check_atcg=check_atcg)
        for c in tdf.columns:
            if c not in ['genome', 'group']:
                tdf[f"{function_name}_{c}"] = tdf[c]
        df_list.append(tdf)
    if verbose:
        print(f"Collating results...")
    df = pd.concat(df_list, axis=1).reset_index()
    df = df.loc[:,~df.columns.duplicated()]
    return df[[c for c in df.columns if any(c.startswith(fn) for fn in function_names) or c == 'genome' or c == 'group']].set_index('genome')

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

        print("Reading dataset...")
        organisms = [o for o in f.read().split("\n") if o != '']
        retry_generator = [retry for _ in range(len(organisms))]
        dataset_generator = [dataset_folder for _ in range(len(organisms))]
        from multiprocessing import Pool
        pool = Pool(procs)
        pool.map(download_org, zip(
            organisms, retry_generator, dataset_generator))

    from shutil import copyfile
    copyfile(dataset_dir, os.path.join(
        dataset_folder, os.path.basename(dataset_dir)))

###
#
# BEGIN MIKE'S TOOLKIT FOR PARALELLIZED GENOME ANALYSIS
# (many functions included in here in order to paralellize effectively
# using Jupyter notebooks)
#
###
    
def substring_finder(enter):
    string1 = enter[0]
    string2 = enter[1]
    answer = ""
    anslist=[]
    len1, len2 = len(string1), len(string2)
    for i in range(0, len1):
        match = ""
        for j in range(0, len2):
            if (i + j < len1 and string1[i + j] == string2[j]):
                match += string2[j]
            else:
                answer = match
                if answer != '' and len(answer) > 5:
                    anslist.append(answer)
                match = ""

        if match != '':
            anslist.append(match)
        # break
    return anslist

def substring_counter(genome):
    k = 14 # SET k to be length of k-mer
    
    count_list = {}
    i = 0
    while i < len(genome):
        kmer = genome[i:i+k]
        if kmer not in count_list:
            count_list[kmer] = 0
        count_list[kmer] += 1
        i+=1
    
    return count_list

def get_genomes(dataset, sample=None):
    import os
    import gzip
    from Bio import SeqIO
    genomes = []
    # Compress standard data
    from os.path import expanduser
    home = expanduser("~")
    dataset_fp = os.path.join(home, FOLDER_NAME, "datasets", f"dataset_{str(dataset)}")
    break_counter = 0
    for fasta in os.listdir(dataset_fp):
        print(fasta)
        if sample is not None and break_counter >= sample:
            break
        # Ignore DS_Store and other hidden files
        if fasta.endswith('fna.gz'):
            fasta_fp = os.path.join(dataset_fp, fasta)
            try:
                with gzip.open(fasta_fp, "rt") as handle:
                    # If want to consider filesize: print(os.fstat(handle.fileno()).st_size)
                    total_genome = ""
                    for record in SeqIO.parse(handle, "fasta"):
                        # sometimes genomes stored across records
                        total_genome += record.seq
                    genomes.append(total_genome)
                break_counter += 1
            except:
                print(f"Couldn't do {fasta}")
    return genomes

def get_tree(d, keys):
    for k in keys:
        d = d[k]
    return d

def set_tree(d, keys, value):
    d = get_tree(d, keys)
    d[value] = {}
    
def get_indices(d):
    for k in d.keys():
        if type(k) == tuple:
            return k
    return ()

def build_tree(reference):
    tree = {}
    # Convert reference into tree
    for i in range(len(reference)):
        if (i != 0 and i % 10000 == 0):
            print(i)
        for j in range(i, len(reference)):
            parents = [reference[x] for x in range(j-i,j)]
            set_tree(tree, parents, reference[j])
            set_tree(tree, parents + [reference[j]], (j-len(parents),j))
    return tree

def linear_ref_compress(reference, compress, tree):        
    compressed_ref = []
    temp_tree = tree
    for i in range(len(compress)):
        if compress[i] in temp_tree:
            temp_tree = temp_tree[compress[i]]
            if len(temp_tree.keys()) > 1:
                next
            else:
                #print(get_indices(temp_tree))
                compressed_ref.append(get_indices(temp_tree))
                temp_tree = tree
        else:
            #print(get_indices(temp_tree))
            compressed_ref.append(get_indices(temp_tree))
            temp_tree = tree[compress[i]]
    final_index = get_indices(temp_tree)
    if final_index != ():
        compressed_ref.append(final_index)
    return compressed_ref

def binarize(char):
    char_dict = {"A": "000", "C": "001", "T":"010", "G":"011", "N":"100", "-":"101", "_ref": "111"}
    return char_dict[char]

def compress_regular(genome):
    out_genome = ""
    for g in genome:
        out_genome += binarize(g)
    return out_genome

def linear_ref_compress_genome(reference):
    tree = build_tree(reference)
    def compress_genome(genome):
        char_dict = {"A": "000", "C": "001", "T":"010", "G":"011", "N":"100", "-":"101", "_ref": "111"}
        ref_indices = linear_ref_compress(reference, genome, tree)
        compressed = ""
        for (start, end) in ref_indices:
            if (end-start) > 6:
                compressed += binarize("_ref")
                compressed += str('{:08b}'.format(start))
                compressed += str('{:08b}'.format(end))
            else:
                for x in range(start, end+1):
                    compressed += binarize(reference[x])
        return compressed
    return compress_genome

def linear_ref_decompress_genome(reference):
    def decompress_genome(compressed):
        char_dict = {"A": "000", "C": "001", "T":"010", "G":"011", "N":"100", "-":"101", "_ref": "111"}
        rev_dict = {char_dict[k]:k for k in char_dict.keys()}
        i = 0
        decompressed = ""
        while i < len(compressed):
            tag = rev_dict[compressed[i:i+3]]
            if tag == "_ref":
                i += 3
                start = int(compressed[i:i+8], 2)
                i += 8
                end = int(compressed[i:i+8], 2)
                i += 8
                decompressed += reference[start:end+1]
            else:
                i += 3
                decompressed += tag
        return decompressed
    return decompress_genome


def test_ref(par_in):
    char_dict = {"A": "000", "C": "001", "T":"010", "G":"011", "N":"100", "-":"101"}
    r, reg_df = par_in
    for c in char_dict.keys():
        if c not in r:
            r += c
    compress_func = linear_ref_compress_genome(r)
    comp_df = benchmark_functions([compress_func], dataset="100bacteria", sample=20, check_atcg=True)
    return ((reg_df['compress_regular_new_len'] - comp_df['compress_genome_new_len']) / reg_df['compress_regular_new_len'] * 100).to_frame('pct_compression').mean()['pct_compression']