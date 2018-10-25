python ~/dbcAmplicons/bin/dbcAmplicons preprocess -B ../2.metadata/dbcBarcodeLookupTable.txt -P ../2.metadata/PrimerTable.txt -S ../2.metadata/workshopSamplesheet.txt -O SouthShore_April.intermediate -1 ../1.convert/SouthShore_April/DBCreads_R1.fastq.gz --debug
python ~/dbcAmplicons/bin/dbcAmplicons preprocess -B ../2.metadata/dbcBarcodeLookupTable.txt -P ../2.metadata/PrimerTable.txt -S ../2.metadata/workshopSamplesheet.txt -O SouthShore_May.intermediate -1 ../1.convert/SouthShore_May/DBCreads_R1.fastq.gz --debug
python ~/dbcAmplicons/bin/dbcAmplicons preprocess -B ../2.metadata/dbcBarcodeLookupTable.txt -P ../2.metadata/PrimerTable.txt -S ../2.metadata/workshopSamplesheet.txt -O stream_May.intermediate -1 ../1.convert/stream_May/DBCreads_R1.fastq.gz --debug

