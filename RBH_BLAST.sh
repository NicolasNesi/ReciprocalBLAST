#!/bin/bash
#SBATCH -J RBHBlast                         # Job name 
#SBATCH -o RBHBlast.%a.%A.out               # File to which stdout will be written
#SBATCH -e RBHBlast.%a.%A.err               # File to which stderr will be written
#SBATCH -p court                            # Partition
#SBATCH -N 1                                # Number of node
#SBATCH -c=96                               # Number of cores/cpus
#SBATCH -t 02-00:00                         # Runtime in DD-HH:MM
#SBATCH --mem 400G                          # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                     # BEGIN,END,ALL
#SBATCH --mail-user=mymail                  # Email address

# ---------------------------------
# Version v1.0
# ---------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date:27/06/2025
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
SCRIPTS="/home/2019013/PARTAGE/Various_Scripts_for_various_purposes"
PATH1="${FOLDER}/Input"
PATH2="/dlocal/home/2019013/Databases/RBH"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" $FOLDER/Scripts/List_samples_$RUN.txt)"
REFERENCE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" $FOLDER/Scripts/List_ReferenceGenome_$RUN.txt)"

RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Load Conda environment
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate NN_env
# include:
# blast+=2.16.0+
# ---------------------------------

# ---------------------------------
# Make subdirectories
# ---------------------------------

for dir in \
	"${FOLDER}/BLAST/${RUN}" \
	"${FOLDER}/Filter_RBH/${RUN}" \
	"${FOLDER}/Candidate_Orthologs/${RUN}" \
    "${FOLDER}/Log/${RUN}_run${RUNDATE}" \
    "${FOLDER}/Scripts/Sbatch_log"
do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${PATH1}"/"${SAMPLE}.fasta" "${LOCAL_WORK_DIR}"
cp "${PATH2}"/"${REFERENCE}.fna" "${LOCAL_WORK_DIR}"
cp "${PATH2}"/"${REFERENCE}.faa" "${LOCAL_WORK_DIR}"

cd "${LOCAL_WORK_DIR}"/

echo Analysis run the: "${RUNDATE}"
echo Working directory: "${PWD}"
echo ".fna Reference in nucleotide"
echo ".faa Reference in amino acids"
echo ".fasta sample query in nucleotide"
# ---------------------------------

# ---------------------------------
# STEP1 Make database
# ---------------------------------
echo "Starting to make database for Query sample and Database"

makeblastdb -in "${SAMPLE}.fasta" -dbtype 'nucl'

makeblastdb -in "${REFERENCE}.fna" -dbtype 'nucl'

makeblastdb -in "${REFERENCE}.faa" -dbtype 'prot'
# ---------------------------------

# ---------------------------------
# STEP2 Perform first blast query vs database
# ---------------------------------
echo "Starting blast: sample in nucleotide (query) versus reference in aa (database) using blastx"

blastx \
-query "${SAMPLE}.fasta" \
-db "${REFERENCE}.faa" \
-out "${SAMPLE}"_nuc_vs_"${REFERENCE}"_aa.blastx \
-evalue 1e-06 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads "${THREADS}"
# ---------------------------------

# ---------------------------------
# STEP3 Perform reverse blast database vs query
# ---------------------------------
echo "Starting blast: reference in aa (query) versus sample in nucleotide (database) using tblastn"

tblastn \
-db "${SAMPLE}.fasta" \
-query "${REFERENCE}.faa" \
-out "${REFERENCE}"_aa_vs_"${SAMPLE}"_nuc.tblastn \
-evalue 1e-06 \
-max_target_seqs 1 \
-outfmt 6 \
-num_threads "${THREADS}"
# ---------------------------------

# ---------------------------------
# STEP 4 Filtering of blastx and blastn results
# ---------------------------------
echo "Running trinity perl script to create an output file .w_pct_hit_length"

perl "${SCRIPTS}"/analyze_blastPlus_topHit_coverage.pl "${SAMPLE}"_nuc_vs_"${REFERENCE}"_aa.blastx "${SAMPLE}.fasta" "${REFERENCE}.faa"

echo "Running perl script filter reciprocal blast matches"

perl "${SCRIPTS}"/filter_reciprocal_blast_matches.pl \
-i "${SAMPLE}"_nuc_vs_"${REFERENCE}"_aa.blastx.w_pct_hit_length \
-r "${REFERENCE}"_aa_vs_"${SAMPLE}"_nuc.tblastn \
-d "${SAMPLE}.fasta" \
-s "${SAMPLE}" \
-o "${SAMPLE}"_candidate_orthologs_nuc
# ---------------------------------

# ---------------------------------
# Move results
# ---------------------------------
echo "Move output files"

mv "${SAMPLE}"_nuc_vs_"${REFERENCE}"_aa.blastx "${FOLDER}/BLAST/${RUN}"

mv "${REFERENCE}"_aa_vs_"${SAMPLE}"_nuc.tblastn "${FOLDER}/BLAST/${RUN}"

mv *w_pct_hit_length "${FOLDER}/Filter_RBH/${RUN}"

mv *.hist* "${FOLDER}/Filter_RBH/${RUN}"

mv "${SAMPLE}"_candidate_orthologs_nuc "${FOLDER}/Candidate_Orthologs/${RUN}"
# ---------------------------------

# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/RBHBlast.* "${FOLDER}/Scripts/Sbatch_log"
# ---------------------------------
