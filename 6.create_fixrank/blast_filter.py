#! /usr/bin/python
import itertools
import sys
# import biom

# input: lineages.csv (supplied by https://github.com/zyxue/ncbitax2lin)
# input: original fasta file, so we can recover full tag
# input: BLASTed file (format: -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qlen")

# Output: Simulated "fixrank" file; looks like RDPClassifier.

# Examples:
#M02326:25:000000000-BDDWL:1:1101:17065:1439|SouthShore_April_2:16SrDNA:304              Bacteria        domain  1.0     Firmicutes      phylum  0.79    Bacilli class   0.6     Bacillales      order   0.51    Paenibacillaceae 2      family  0.15    Oxalophagus     genus   0.15
#M02326:25:000000000-BDDWL:1:1101:11289:1559|SouthShore_April_2:16SrDNA:306              Bacteria        domain  1.0     "Actinobacteria"        phylum  0.96    Actinobacteria  class   0.96    Actinomycetales order   0.95    Intrasporangiaceae      family  0.55    Ornithinicoccus genus   0.11
#M02326:25:000000000-BDDWL:1:1112:2375:14789|SouthShore_May_3:16SrDNA:303	Bacteria	domain	1.0	Bacteroidetes	phylum	1.0	Flavobacteriia	class	1.0	Flavobacteriales	order	1.0	Flavobacteriaceae	family	1.0	Flavobacterium	genus	1.0	uncultured Flavobacterium sp.	species	1.0

#  starts from a blast results file and the taxonomy file
# filters blast results (only shows unique results)
# outputs staxids and/or lineages.

# M02326:25:000000000-BDDWL:1:1101:12251:1495	gi|257799000|gb|GQ844276.1|	80.62	289	54	2	20	307	141	428	4e-54	 222	265564
### Identify your blast results file #####
#print "What is your blast result file titled?"
#blast_input_filename = raw_input("Datafile Path: ")
#blast_input_filename= '/Users/joe/eDNA_Mountain_Lake/to_blast/good_blast/eDNA_Mountain_Lake/intermediate/28S/28S.blast_results'
if len(sys.argv) < 3:
    print "blast_filter.py blast-input-file fastq-file"
    sys.exit(1)
blast_input_filename = sys.argv[1]
fastq_file = sys.argv[2]
# print "Scanning file:",blast_input_filename
blast_results = open(blast_input_filename, 'rU')
lineages={}
kingdoms={}
phyla={}
all=[]
fasta_index={}

def load_fasta_tags(filename):
    lineages_file = open(filename, 'r')
    line_count=0
    for cur_line in lineages_file:
        if line_count % 4 == 0:
            # @M02326:25:000000000-BDDWL:1:1101:17065:1439 1:N:0:SouthShore_April_2:16SrDNA GCCAATAT|0|GTAAGGAG|0 MS-27F|1|20|
            array = cur_line.split(" ")
            id = array[0].strip('@')
            map_elements = array[1].split(":")
            map = map_elements[3] + ":" + map_elements[4]
            fasta_index[id] = map
        line_count +=1


def load_lineage():
    # print "Loading lineages.",
    lineages_file = open('lineages.csv', 'r')
    header=None
    line_count=0
    for lineage_line in lineages_file:
        line_count += 1
        # if (line_count % 181836/2 == 0):
        #     print ".",
        parts = lineage_line.strip().split(',')

        if header is None:
            header = parts
        else:
            staxid = parts[0]
            vals={}
            tuples = list(itertools.izip_longest(header,parts))


            for col_header,col_value in tuples:
                
                vals[col_header] = col_value
            lineages[staxid] = vals
                

def blast_filter(line):
    data = line.strip().split('\t')
    bitscore = float(data[11])
    e_value = float(data[10])
    s_end = int(data[9])
    s_start = int(data[8])
    q_start = int(data[6])
    q_end = int(data[7])
    gap_open =int(data[5])
    mmismatch=int(data[4])
    length=int(data[3])
    p_ident=float(data[2])
    qlen= float(data[13])
    if e_value < 10e-50 and bitscore > 175 and p_ident >= 100.0:


        # print "good blast"
        return True
    else:
        return False

def get_lineage(staxid):
    if staxid in lineages:
        return lineages[staxid]
    else:
        sys.stderr.write("Missing taxid:" + str(staxid))
        return None

def add_rank(lineage,i,translated_array, rank_array,prev_score):
    return "\t" + get_lowest(lineage,i,translated_array) +\
           "\t" +\
           rank_array[i]+\
           "\t"

def get_lowest(lineage,i,translated_array):
    if i > len(translated_array):
        print "failed lookup:", lineage
        return "failed_lookup"
    
    if len(lineage[translated_array[i]]) <= 1:
        return get_lowest(lineage,i+1,translated_array)
    else:
        return lineage[translated_array[i]]

def get_score(lineage,i, translated_array):
    if len(lineage[translated_array[i]]) <= 1:
        return float(0)
    else:
        return float(0.99)

def generate_fixrank_output(lineage,taxid,qlen,id):
    #Example fixrank output:
    #M02326:25:000000000-BDDWL:1:1101:17065:1439|SouthShore_April_2:16SrDNA:304              Bacteria        domain  1.0     Firmicutes      phylum  0.79    Bacilli class   0.6     Bacillales      order   0.51    Paenibacillaceae 2      family  0.15    Oxalophagus     genus   0.15

    # print "generating lineage for taxid", taxid, "with id:",id
    fixrank=id + "|" + fasta_index[id] + ':'+str(qlen) +'\t'
    rank_array = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    translated_array = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    # print fixrank
    prev_score = 1.0
    for i in enumerate(rank_array):
        fixrank = fixrank + add_rank(lineage,i[0],translated_array, rank_array,prev_score)
        score = get_score(lineage,i[0],translated_array)
        if prev_score == 0:
            score = 0
        fixrank = fixrank + str(score)
        prev_score = score


    return fixrank

def generate_staxid_output(taxid):
    print taxid


def main():


    taxid_output = open(blast_input_filename + '_taxids.txt', 'w')
    last_query = ''
    for line in blast_results:

        columns = line.split("\t")
        query = columns[0]
        staxid = columns[12].strip()
        qlen= int(columns[13].strip())

        if last_query != query and staxid.isdigit() and blast_filter(line):
            # print "new query:",query
            last_query = query
            lineage = get_lineage(staxid)
            # print (staxid)
            # print lineage
            if staxid not in all:
                all.append(staxid)
       # print line
            if lineage is not None:
                print generate_fixrank_output(lineage,staxid,qlen,line.strip().split()[0])
                #generate_fixrank_output(staxid)
                # print lineage
                # print line

                # if (len(lineage['kingdom']) > 0):
                #     kingdoms[lineage['kingdom']] = True
                # if (len(lineage['phylum']) > 0 and lineage['kingdom'] == 'Metazoa'):
                #     phyla[lineage['phylum']] = True
                # if lineage['phylum'] == 'Chordata':
                #     print lineage
                #     print line
            # exit(1)

    # for key in kingdoms:
    #     print key
    #
    # print "\nPhyla:"
    # for key in phyla:
    #     print key
    # for staxid in all:
    #     # lineage = s
    #     print staxid
    #     # lineage =  get_lineage(staxid)
    #     # if 'sp.' not in lineage['species']:
    #     # #     print "*", lineage['genus']
    #     # # else:
    #     #     print lineage['species']


load_lineage()
load_fasta_tags(fastq_file)
main()