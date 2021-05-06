compresssion.py is legacy code. It was an attempt (succesful) at compressing a multiple alignment.

compressionV2.py is the file containing the code for the differential compression algorithm contained in our paper.

main.py is code to run this algorithm and compress COVID19 genomes on the Andrew cluster in parallel.

notify.py is code to send email notifications, used in main.py to provide status updates.

parallelExperiments.py is a tiny file where I figured out how to do parallel processing in python.

Random experiments contain a few jupyter notebooks. Not sure what they contain. Similarly for Comp\ Bio\ Project.ipyb.

analysis.py is old code - don't know what it does.

In addition to the files hosted on GitHub, locally I have a lot of genomic files. The most important are
- 2.57GB 85k COVID19 genomes downloaded from the NCBI
- This file partitioned into smaller pieces, in 10k increments. (2.57GB exceeds the space I have access to on the Andrew machines. To compress the genomes I had to do them bit by bit.
- The compressed versions of the COVID19 genomes.
These can't be uploaded because they're too large.

I also have some genomic files in a folder titled MikeData which contains genomes which I analyzed - just for fun / to get a sense for frequency distributions of Khmers etc - separately from the main algorithm.

I also have a file titled emailANDpassword.txt containing email details, needed for the notify.py script. See notify.py for details on required formatting.
