# Hybrid dbcAmplicons/jrussack eDNA pipeline instructions

Inserts a section that make it possible to run this pipeline against NCBI databases using BLAST.

[Github for BLAST python file](https://github.com/calacademy-research/eDNA_blastGithub) 

##Original UC Davis courseware:

[Full eDNA class by Matt Settles](https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/)

[Command line summary](https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/thursday/dbcAmplicons_commands.html)
​        
[Using phyloseq](https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html)

[github for above](https://github.com/ucdavis-bioinformatics-training/2017-September-Microbial-Community-Analysis-Workshop)

Pipeline, starting with fastq files from Illumina, containing primers. Derived from command line summary.

###Prerequisites:
Python scripts are written in python 2.7.c

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


And we are working with four loci (COI, 16S, 18S, 28S)
________________
###0: Install the dbcAmplicons pipeline.

Install prerequisites (per readme: flash2, RDPClassifier, bowtie2). Ideally, do this using bioconda.

installing DBCAmplicons:
```
$ cd ~/
$ git clone https://github.com/msettles/dbcAmplicons
$ cd dbcAmplicons
```

Now, we need to patch dbcAmplicons (to deal with bc1 '+' bc2 format). Make these changes manually
or use the files in ./dbcAmplicons_changes. Note that in either case, the source code may
have since evolved, so check to ensure that there aren't any conflicts.

Note that in the this repository, these files exist with the edits already made. Look in dbcAmpliconsChanges at the 
top level of this repository.

Edit the file dbcAmplicons/sequenceReads.py. Change line 351 (may be slightly off if there are new versions) from:
```
​        elif len(self.barcode) is (bc1_length + bc2_length):
```
To
```
        elif len(self.barcode) is (bc1_length + bc2_length +1):
```

and line 353 from:
```
                bc2 = self.barcode[bc2_length:]

```
To
```                
                bc2 = self.barcode[-bc2_length:]

```
Edit the file dbcAmplicons/sequenceReads.py. Change line 351 (may be slightly off if there are new versions) from:
```
​        elif len(self.barcode) is (bc1_length + bc2_length):
```
To
```
        elif len(self.barcode) is (bc1_length + bc2_length +1):
```
Edit the file dbcAmplicons/illumaRun.py. Change line 159 (may be slightly off if there are new versions) from:
```
                bc_2 = self.BC2.next().rstrip()  # read
```
To
```
                bc_2 = self.BC2.next().rstrip()[-8:]  # read
```



Edit the file dbcAmplicons misc.py to accommodate ambiguity codes in the barcode and change line 131 from:
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
### 1: Convert the 2-read samples to 4-read samples.
Do this for each of the raw read sample files. In this example, that's a total of three runs. 

Note that this outputs files with the names 'DBCreads_R1.fastq.gz' (one through four); 
these files will stomp each other if you do the entire run in one directory without renaming.
```
python ~/dbcAmplicons/scripts/python/convert2ReadTo4Read.py -v --debug -1 ./ML_SouthShore_April_S193_L001_R1_001.fastq.gz -2 ./ML_SouthShore_April_S193_L001_R2_001.fastq.gz
```

Do this for each sample pair. I recommend doing each in its own subdirectory, as it creates files that look like this:
DBCreads_R1.fastq.gz
And they will step on each other if you run with new input.


### 2: Compose three input files (aka: generate metadata)

Trailing newlines will cause an error in step 4.

Personally, I place these in the “metadata” subdirectory of my working directory.

#### PrimerTable.txt TSV table maps primers to loci. Note that the filename is arbitrary.

Note that this should be a proper TSV table; the formatting below is for ease of reading.
​                                 
Lists forward and reverse primers for each loci.

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

#### dbcBarcodeLookupTable.txt Maps samples to the corresponding dual barcodes.
We expect 8 base barcodes. 

| Sample name                      | leading barcode | trailing barcode 
| :------------------------------: | :-----------: | :-----------: 
|stream_May                        |  TATACAAG     |    ACTGCATA   
|SouthShore_May                    |  ATATTGGC     |    ACTGCATA   
|SouthShore_April                  |  ATATTGGC     |    GTAAGGAG   
             
#### sampleSheet.txt Defines the relationship between sample ID, barcodes, primer pair ID, and loci.

A sample ID is a unique identifier that applies to each combination of barcode pairs (Pair_ID) x Sample.
"ProjectID" is the loci. In the example below, there are two primers that map to the COI loci.

| SampleID | PairID   | BarcodeID      | ProjectID                       |
| :----------: | :---------: | :-------------------: | :----------: |                               
|SouthShore_April_0  | COI-A     | SouthShore_April | COI
|stream_May_1        | 16SrDNA   | stream_May       | 16S
|SouthShore_April_2  | 16SrDNA   | SouthShore_April | 16S
|SouthShore_May_3    | 16SrDNA   | SouthShore_May   | 16S
|stream_May_4        | COI-B     | stream_May       | COI
|SouthShore_April_5  | COI-B     | SouthShore_April | COI
|SouthShore_May_6    | COI-B     | SouthShore_May   | COI
|stream_May_7        | 18SrDNAV1 | stream_May       | 18S
|SouthShore_April_8  | 18SrDNAV1 | SouthShore_April | 18S
|SouthShore_May_9    | 18SrDNAV1 | SouthShore_May   | 18S
|stream_May_10       | 28SrDNA   | stream_May       | 28S
|SouthShore_April_11 | 28SrDNA   | SouthShore_April | 28S
|SouthShore_May_12   | 28SrDNA   | SouthShore_May   | 28S             






### 3: Preprocess
Run the preprocess step on the split files.
```
python ~/dbcAmplicons/bin/dbcAmplicons preprocess -B ../2.metadata/dbcBarcodeLookupTable.txt -P ../2.metadata/PrimerTable.txt -S ../2.metadata/samplesheet.txt -O SouthShore_April.intermediate -1 ../1.convert/SouthShore_April/DBCreads_R1.fastq.gz --debug
python ~/dbcAmplicons/bin/dbcAmplicons preprocess -B ../2.metadata/dbcBarcodeLookupTable.txt -P ../2.metadata/PrimerTable.txt -S ../2.metadata/samplesheet.txt -O SouthShore_May.intermediate -1 ../1.convert/SouthShore_May/DBCreads_R1.fastq.gz --debug
python ~/dbcAmplicons/bin/dbcAmplicons preprocess -B ../2.metadata/dbcBarcodeLookupTable.txt -P ../2.metadata/PrimerTable.txt -S ../2.metadata/samplesheet.txt -O stream_May.intermediate -1 ../1.convert/stream_May/DBCreads_R1.fastq.gz --debug
```
Check the outputs to ensure that the barcodes are matching.

### 4: Join
Wrapper for [flash2](https://github.com/dstreett/FLASH2)

Add it to your path, and make sure it shows up when you do a "which flash2".

```
PATH="~/FLASH2/flash2:$PATH"
```



Single example:
```
python ~/dbcAmplicons/bin/dbcAmplicons join -t 4 -O Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3 -1 Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3_R1.fastq.gz > join-16sV1V3.log

```

Handy script runs this on directories named in "samples" so you don't have to run join manually 12 times:
```
loci=( 16S 18S COI 28S )
samples=( 'SouthShore_April' 'SouthShore_May' 'stream_May')

for sample in "${samples[@]}"
do

    for locus in "${loci[@]}"
    do
        echo "Processing ${sample}_${locus}"
        python ~/dbcAmplicons/bin/dbcAmplicons join -t 4 -O ./${sample}_${locus} -1 ../3.preprocess/${sample}.intermediate/${locus}_R1.fastq.gz -2 ../3.preprocess/${sample}.intermediate/${locus}_R2.fastq.gz > $sample\_$locus.log
    done
done
```

# NCBI steps
At this point we diverge from the msettles pipeline. The 'extendedFrags' files in step 4 are the sequences, tagged by sample, that we want to match against our database.
```
    4224 -rw-r--r--  1 joe  wheel   1.8M Oct 24 19:42 stream_May_28S.extendedFrags.fastq.gz
    1024 -rw-r--r--  1 joe  wheel   507K Oct 24 19:42 stream_May_COI.extendedFrags.fastq.gz
    4224 -rw-r--r--  1 joe  wheel   1.7M Oct 24 19:42 stream_May_18S.extendedFrags.fastq.gz
    6272 -rw-r--r--  1 joe  wheel   2.4M Oct 24 19:42 stream_May_16S.extendedFrags.fastq.gz
    4224 -rw-r--r--  1 joe  wheel   1.7M Oct 24 19:42 SouthShore_May_28S.extendedFrags.fastq.gz
    1664 -rw-r--r--  1 joe  wheel   779K Oct 24 19:42 SouthShore_May_COI.extendedFrags.fastq.gz
    6272 -rw-r--r--  1 joe  wheel   2.5M Oct 24 19:42 SouthShore_May_18S.extendedFrags.fastq.gz
    4224 -rw-r--r--  1 joe  wheel   1.6M Oct 24 19:42 SouthShore_May_16S.extendedFrags.fastq.gz
   10368 -rw-r--r--  1 joe  wheel   4.6M Oct 24 19:42 SouthShore_April_28S.extendedFrags.fastq.gz
    8320 -rw-r--r--  1 joe  wheel   3.8M Oct 24 19:42 SouthShore_April_COI.extendedFrags.fastq.gz
   10368 -rw-r--r--  1 joe  wheel   4.2M Oct 24 19:42 SouthShore_April_18S.extendedFrags.fastq.gz
    6272 -rw-r--r--  1 joe  wheel   2.3M Oct 24 19:41 SouthShore_April_16S.extendedFrags.fastq.gz
```

### 5: Convert to fasta and run BLAST.

Fasta is what blast likes; convert "extendedfrags" files to fasta format.

Note that per ['Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows'](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty833/5106166?redirectedFrom=fulltext) 
we're not using max_target_seqs. You will want to post-process to remove low quality matches.
```
blastn -db /media/BigAssRAID/NCBI_database/nt/nt -query H3.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qlen” -out H3.out
```

### 6: Convert blast output to a "fixrank" file. 
This is similar to the output produced by 
rdpClassifier, and following this chain will make it possible rejoin the msettles pipeline
and run phyloseq. Note that abundance measurements are not calibrated when using NCBI, proceed
at your own risk.

get the most current lineage file from [here](https://gitlab.com/zyxue/ncbitax2lin-lineages/tree/master). 
The format converter, "blast_filter.py", expects this to be unzipped and named "lineages.csv". Using a softlink will
help keep track of whether you're using the most recent file.

```ln -s ./lineages.csv ./lineages-2019-01-07.csv```

usage:
```python blast_filter.py blast-input-file fastq-file```

example:

```python blast_filter.py ../5.blast/COI.blast_results ../4.join/COI.joined.extendedFrags.fastq > COI.fixrank```

The script uses the input fastq file from step 5 to get the names that the fasta format doesn't support.

The lineage generator that generates this file is [here](https://github.com/zyxue/ncbitax2lin).

Note (from fixrank documentation pipeline:)
"The "fixrank" format only outputs the results for a list of selected ranks in the following order: domain, phylum, class, order, 
family and genus. In case of missing ranks in the lineage, the bootstrap value and the taxon name from the immediate lower rank 
will be reported. This eliminates the gaps in the lineage, but also introduces non-existing taxon name and rank. Interpret the "fixrank" 
results with caution."

Note: 
ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz

 
 ### 7. Create abundance (from original pipeline):
 
 ```
 ~/dbcAmplicons/bin/dbcAmplicons abundance -S samplesheet.txt -O abundance_output -F 16S.fixrank  --biom > 16S.abundance.log
```
From the [docs](https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/thursday/dbcAmplicons_commands.html):
"dbcAmplicons abundance command above (with –biom) produces four files: abundance, proportions, tax_info and biom files. " 


Abundance tables are supported by phyloseq and qiime for analysis. The script "create_abundance.sh" is hardcoded
to perform these steps.

 ### 8. Optional: Check primer performance by generating a taxonomy based tree based on staxids. 
 
 run treemerge.py. Check comments for inputs. To generate staxids, modify blast_filter.py (step 6). Change
 generate_fixrank_output funciton to output only the taxid, one per line. Redirect to the input files for treemerge.py.
 