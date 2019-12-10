hybridPipeline:

takes illumina short reads and ONT long reads as input

how it works: construct a long read only assembly, then polish using the higher
nucleotide accuracy short reads

tools used:
flye
minimap2
racon
bwa
samtools
pilon
