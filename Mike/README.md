# Using data_process

Raise an issue if anything doesn't make sense! :)

## Data Benchmarking:

### benchmark
Use data_process.benchmark to benchmark a single function.

Output: dataframe with relevant metrics

Example usage:
```
from data_process import benchmark

str -> str
def x(s):
    # Expected whole genome
    # Returns compressed genome
    return s

benchmark(x, dataset=1, sample=3, check_atcg=False)
```
i.e Benchmark function x on the first 3 elements of `dataset_1` (see datasets below). Do not skip genomes with bases other than ATCG.

Can also pass in `data_dir` instead of `dataset` if specify particular directory.

### benchmark_functions
Use data_process.benchmark_functions to benchmark a single function. 

Output: dataframe with relevant metrics

Example usage:
```
from data_process import benchmark_functions

str -> str
def x(s):
    # Expected whole genome
    # Returns compressed genome
    return s

str -> str
def y(s):
    # Expected whole genome
    # Returns compressed genome
    return s

benchmark_functions([x,y], dataset=1, sample=3, check_atcg=False)
```
Same as `benchmark` but runs on multiple functions. 


## Dataset Downloading

`download_dataset(datasets/dataset_1.tsv, procs=10, retry=2)`

Download `dataset_1` with 10 processors. Retry twice if a download fails (can happen with ncbi).

## Datasets

- `datasets/dataset_1.tsv`
- `datasets/dataset_10.tsv`
- `datasets/dataset_100.tsv`
- `datasets/dataset_500.tsv`

`dataset_{n}.tsv` samples `n` organism genomes from the following groups:

- archaea
- bacteria 
- fungi 
- invertebrate 
- plant 
- protozoa 
- vertebrate_mammalian 
- vertebrate_other 
- viral