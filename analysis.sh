NANOPORE_PATH=''
ILLUMINA_R1_PATH=''
ILLUMINA_R2_PATH=''
ILLUMINA_SINGLEEND_PATH=''

USE_PILON="false"
USE_RACON="false"
SINGLE_END="false" # will assume paired end illumina data unless -s flag used
GENOME_SIZE="5m"


# -n : nanopore reads path
# -1 : illumina read 1 path
# -2 : illumina read 2 path
# -s : illumina single end path
# -r : polish using racon
# -p : polish using pilon
# -g : genome size


while getopts "hrpn:1:2:s:g:o:" opt; do
  case ${opt} in
    h )
      echo "Help:"
      echo -e "\tInput:"
      echo -e "\t\t-1 : R1 for illumina paired end path"
      echo -e "\t\t-2 : R2 for illumina paired end path"
      echo -e "\t\t-s : single end illumina path"
      echo -e "\t\t-n : nanopore read path"
      echo -e "\tOptions:"
      echo -e "\t\t-r : polish flye contigs using racon + nanopore reads"
      echo -e "\t\t-p : polish racon / flye contigs using pilon + illumina reads"
      echo -e "\t\t-g : genome size (default: 5m)"
      echo -e "\t\t-o : output directory prefix"
      exit 1
      ;;
    n ) NANOPORE_PATH=$OPTARG
      ;;
    1 ) ILLUMINA_R1_PATH=$OPTARG
      ;;
    2 ) ILLUMINA_R2_PATH=$OPTARG
      ;;
    s ) ILLUMINA_SINGLEEND_PATH=$OPTARG; SINGLE_END="true"
      ;;
    g ) GENOME_SIZE=$OPTARG
      ;;
    o ) OUTPUT_PREFIX=$OPTARG
      ;;
    r ) USE_RACON="true"
      ;;
    p ) USE_PILON="true"
      ;;
    \? ) echo "Usage: analysis.sh [-rpn12s]"; exit 1
      ;;
  esac
done

OUTDIR="${OUTPUT_PREFIX}_hybrid_output"
mkdir $OUTDIR
FLYE_OUTDIR="${OUTDIR}/${OUTPUT_PREFIX}_flye_output"
MINIMAP_OUTPUT="${OUTDIR}/${OUTPUT_PREFIX}_nanopore_alignment.paf"
RACON_OUTPUT="${OUTDIR}/${OUTPUT_PREFIX}_racon_contigs.fasta"
BWA_OUTPUT="${OUTDIR}/${OUTPUT_PREFIX}_illumina_alignment.bam"
PILON_OUTPUT="${OUTDIR}/${OUTPUT_PREFIX}_pilon_output"

echo -e "\n\n### Beginning pipeline ###\n\n"
echo -e "\n\n### Assembling nanopore reads using flye long read assembler ###\n\n"
flye --meta --plasmids --nano-raw $NANOPORE_PATH --genome-size $GENOME_SIZE --out-dir $FLYE_OUTDIR

echo -e "\n\n### Run minimap to align nanopore reads to draft assembly ###\n\n"
minimap2 -x map-ont $FLYE_OUTDIR/assembly.fasta $NANOPORE_PATH > $MINIMAP_OUTPUT

echo -e "\n\n### Polish using racon ###\n\n"
racon $NANOPORE_PATH $MINIMAP_OUTPUT $FLYE_OUTDIR/assembly.fasta > $RACON_OUTPUT

echo -e "\n\n### Index racon contigs ###\n\n"
bwa index $RACON_OUTPUT

echo -e "\n\n### Align illumina reads: data is single-end? ${SINGLE_END}  ###\n\n"
if [ "$SINGLE_END" = "false" ]
then
  bwa mem $RACON_OUTPUT $ILLUMINA_R1_PATH $ILLUMINA_R2_PATH | samtools view -S -b -u - | samtools sort - -o $BWA_OUTPUT
fi

if [ "$SINGLE_END" = "true" ]
then
  bwa mem $RACON_OUTPUT $ILLUMINA_SINGLEEND_PATH | samtools view -S -b -u - | samtools sort - -o $BWA_OUTPUT
fi

samtools index $BWA_OUTPUT

echo -e "\n\n### Running pilon ###\n\n"
java -Xmx16G -jar /media/beastadmin/SeagateExpansion/programs/pilon-1.22.jar --genome $RACON_OUTPUT --bam $BWA_OUTPUT --outdir $PILON_OUTPUT --output pilon.contigs

