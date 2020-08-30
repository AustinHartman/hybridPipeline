# $1 = assembly file
# $2 = nanopore fastq.gz
# $3 = illumina r1 fastq.gz
# $4 = illumina r2 fastq.gz
# $5 = output prefix

bwa index $1

bwa mem $1 $3 $4 > $5_illumina_alignment.sam
minimap2 -ax map-ont $1 $2 > $5_nanopore_alignment.sam

samtools view -Sb $5_illumina_alignment.sam > $5_illumina_alignment.bam
samtools view -Sb $5_nanopore_alignment.sam > $5_nanopore_alignment.bam

# rm large alignment files
rm $5_illumina_alignment.sam
rm $5_nanopore_alignment.sam

samtools sort $5_illumina_alignment.bam > $5_illumina_alignment_sorted.bam
samtools sort $5_nanopore_alignment.bam > $5_nanopore_alignment_sorted.bam

samtools index $5_illumina_alignment_sorted.bam
samtools index $5_nanopore_alignment_sorted.bam

samtools depth -a $5_illumina_alignment_sorted.bam $5_nanopore_alignment_sorted.bam > $5_output.txt
awk '{sum+=$3+$4;} END{print sum/2500000;}' $5_output.txt > $5_rd.txt



