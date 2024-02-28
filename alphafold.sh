#!/bin/bash

#SBATCH --cpus-per-task 20
#SBATCH --gpus=1
#SBATCH --mem=300G
#SBATCH --tmp=300G
#SBATCH --time=2-0
#SBATCH --partition=special_features
#SBATCH --reservation=test_new_tcag_gpu
#SBATCH --job-name alphafold

module load alphafold/2.3.2
ALPHAFOLD_DATA="/hpf/tools/centos7/alphafold/2.3.2/data"

if [ $# -eq 0 ]
then
    echo "Error: No argument is supplied."
    echo "usage: <script> $input_fasta $output_directory"
    exit
fi

fasta=$1
outdir=$2
## AlphaFold2 flags
data_dir="$ALPHAFOLD_DATA"
uniref90="$ALPHAFOLD_DATA/uniref90/uniref90.fasta"
mgnify="$ALPHAFOLD_DATA/mgnify/mgy_clusters_2022_05.fa"
pdb70="$ALPHAFOLD_DATA/pdb70/pdb70"
template_mmcif="$ALPHAFOLD_DATA/pdb_mmcif/mmcif_files"
bfd="$ALPHAFOLD_DATA/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
max_date="2021-01-01"
obsolete_pdbs="$ALPHAFOLD_DATA/pdb_mmcif/obsolete.dat"
uniref30="$ALPHAFOLD_DATA/uniref30/UniRef30_2021_03"
model_preset="monomer"
pdb_seqres="$ALPHAFOLD_DATA/pdb_seqres/pdb_seqres.txt"
uniprot="$ALPHAFOLD_DATA/uniprot/uniprot.fasta"

alphafold \
--fasta_paths ${fasta} \
--output_dir ${outdir} \
--data_dir $ALPHAFOLD_DATA \
--uniref90_database_path $uniref90 \
--mgnify_database_path $mgnify \
--template_mmcif_dir $template_mmcif \
--max_template_date $max_date \
--obsolete_pdbs_path $obsolete_pdbs \
--bfd_database_path $bfd \
--uniref30_database_path $uniref30 \
--pdb70_database_path $pdb70 \
--model_preset $model_preset \
--use_gpu_relax=false
