# script used to run custom pipeline on each of the sets of fastq files
# austin hartman 11/1/2019

# Austin Hartman 1/23/2020
# first command line arg is path to folder containing nanopore long reads
# second command line arg is the path to folder containing illumina reads


LONGREADS=$(ls $1*fastq.gz)
ID=()

for FILE in $LONGREADS
do
	ID+=(${FILE: -31:-25})
done

echo -e "Created list of bacteria IDs\n\n\n"

for PREFIX in ${ID[@]}
do
	NUM_ILLUMINA_FILES=$(find $2 -maxdepth 1 -name "*${PREFIX}*" | wc -l)

	if [ "${NUM_ILLUMINA_FILES}" -eq 2 ]; then
		echo "Running hybrid pipeline with paired-end illumina data"
		/media/beastadmin/SeagateExpansion/assemblies/helper_functions/hybridPipeline/analysis.sh -n ${1}${PREFIX}*.fastq.gz -1 ${2}${PREFIX}*R1*.fastq.gz -2 ${2}${PREFIX}*R2*.fastq.gz -r -p -g 5m -o ${PREFIX}


	else
		echo "Running hybrid pipeline with single-end illumina data"
		/media/beastadmin/SeagateExpansion/assemblies/helper_functions/hybridPipeline/analysis.sh -n ${1}${PREFIX}*.fastq.gz -s ${2}${PREFIX}*R1*.fastq.gz -r -p -g 5m -o ${PREFIX}
	fi
done
