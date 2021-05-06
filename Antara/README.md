
# Compressing/Decompressing Data

The following command compresses the data in `src.txt` and stores the result in `dest.txt`:

```
$ python bwrleFF.py -c src.txt dest.txt
```

Similarly, the following command decompresses the data, storing the decompressed version in `src_cpy.txt`:

```
$ python bwrleFF.py -d dest.txt src_cpy.txt
```

Replace `bwrleFF.py` with `bwrleCF.py` to compress data using a condensed footer rather than a full footer. Only files that are compressed with a full footer can be used with `matches.py`, while only files compressed with a condensed footer can be used with `occurrences.py`.


# Pattern Matching

The following command searches for the string stored in `pat.txt` in the file `cmpsd.txt` (which must have been produced by `bwrleFF.py`), and stores the results in `match_list.txt`:

```
python matches.py cmpsd.txt pat.txt match_list.txt
```