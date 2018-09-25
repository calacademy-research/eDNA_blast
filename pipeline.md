# Hybrid dbcAmplicons/jrussack eDNA pipeline instructions

Inserts a section that make it possible to run this pipeline against NCBI databases using BLAST.

[Github for BLAST python file](https://github.com/calacademy-research/eDNA_blastGithub) 

##Original UC Davis courseware:

[Full eDNA class by Matt Settles](https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/)

[Command line summary](https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/thursday/dbcAmplicons_commands.html)
​        
[Using phyloseq](https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html)

Pipeline, starting with fastq files from Illumina, containing primers. Derived from command line summary.

###Prerequisites:

Identify your per-sample barcodes. This pipeline is set up for dual ended bardcodes. 
Do not strip them prior to running this pipeline. In this example, we have three samples:

| Sample IDs|
| :------ |
|stream_May          |
| SouthShore_May     |
| SouthShore_April   |


Each sample has its attendant barcodes:


 Sample ID         | 5' Barcode | 3' Barcode
 :------:          | :--------: | :---------:
stream_May         | TATACAAG   | ACTGCATA
SouthShore_May     | ATATTGGC   | ACTGCATA 
SouthShore_April   | ATATTGGC   | GTAAGGAG 

________________
###1: Install the dbcAmplicons pipeline.

Install prerequisites (per readme: flash2, RDPClassifier, bowtie2). Ideally, do this using bioconda.
```
$ cd ~/
$ git clone https://github.com/msettles/dbcAmplicons
$ cd dbcAmplicons
```

Now, we need to patch dbcAmplicons:

Edit the file dbcAmplicons/sequenceReads.py. Change line 351 (may be slightly off if there are new versions) from:
```
​        elif len(self.barcode) is (bc1_length + bc2_length):
```
To
```
        elif len(self.barcode) is (bc1_length + bc2_length +1):
```

Edit the file dbcAmplicons misc.py and change line 131 from:
```
        basecomplement = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a', 'N': 'N', 'n': 'n'}
```
To
```
​        basecomplement = {'Y' : 'R','S':'S','K':'M','A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a', 'N': 'N', 'n': 'n'}
```

Edit the file dbcAmplicons/primers.py and change line 51 from:
```
READ, PAIR, ID, SEQ = row.split('\t')[0:4]
```
To:
```
READ, PAIR, ID, SEQ = row.split('')[0:4]
```



Now, finish installing:
```
$ python setup.py install
$ export PYTHONPATH=~/dbcAmplicons:$PYTHONPATH 
```
________________
### 2: Convert the 2-read samples to 4-read samples.
​        python ~/dbcAmplicons/scripts/python/convert2ReadTo4Read.py -v --debug -1 ./ML_SouthShore_April_S193_L001_R1_001.fastq.gz -2 ./ML_SouthShore_April_S193_L001_R2_001.fastq.gz


Do this for each sample pair. I recommend doing each in its own subdirectory, as it creates files that look like this:
DBCreads_R1.fastq.gz
And they will step on each other if you run with new input.


### 3: Compose input files
Personally, I place these in the “metadata” subdirectory of my working directory.
Primer table
#### This TSV table maps primers to loci. 

​    

| READ | Pair_ID   | Primer_ID      | Sequence                       |
| :-----: | :---------: | :--------------: | :------------------------------: |
| P5    | 16SrDNA   | MS-27F         | AGAGTTTGATCCTGGCTCAG           |
| P5    | 28SrDNA   | 28s_3F6        | TTTTGGTAAGCAGAACTGGYG          |
| P5    | 18SrDNAV1 | SSU_FO45       | GCTTGTCTCAAAGATTAAGCC          |
| P5    | COI-A     | mCOintF1       | GGWACWGGWTGAACWGTWTAYCCYCC     |
| P5    | COI-B     | ZBJ-ArtF1c     | AGATATTGGAACWTTATATTTTATTTTTGG |
| P7    | 16SrDNA   | MS-338R        | TGCTGCCTCCCGTAGGAGT            |
| P7    | 28SrDNA   | 28s_4R6        | ABTYGCTACTRCCACYRAGATC         |
| P7    | 18SrDNAV1 | SSU_R225       | GCCTGCTGCCTTCCTTGGA            |
| P7    | COI-A     | Fol-degen-rev2 | TANACYTCNGGRTGNCCRAARAAYCA     |
| P7    | COI-B     | ZBJ-ArtR2c     | WACTAATCAATTWCCAAATCCTCC       |

             

Ambiguity; some have three columns, some have 4. TBD
Sample sheet




Preprocess
```
dbcAmplicons preprocess -B ./metadata/dbcBarcodeLookupTable.txt -P ./metadata/PrimerTable.txt --debug -S ./metadata/workshopSamplesheet.txt -O intermediate -1 ./DBCreads_R1.fastq.gz
```