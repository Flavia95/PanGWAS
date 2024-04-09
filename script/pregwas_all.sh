##Run the pipeline before run GEMMA on 160 linked-read folders, to run the job use:  sbatch -p workers --job-name extractlinked -c 1 --wrap "hostname; cd /scratch; /lizardfs/flaviav/mouse/assembly_D/GWAS/pregwas_all.sh"

#!/bin/bash
# Set the paths
BWA_PATH="/lizardfs/flaviav/tools/bwa-0.7.17/bwa"
GFAINJECT_PATH="/lizardfs/flaviav/tools/gfainject/target/release/gfainject"
REFERENCE="/lizardfs/flaviav/mouse/assembly_D/8_pggb/D_C.fa.gz"
REFERENCE_pggb="/lizardfs/flaviav/mouse/assembly_D/GWAS/pggb_out_ref/D_C_mm10.fa.gz"
DATA_DIR="/lizardfs/flaviav/data/linked_mouse_fastq"
OUTPUT_BWA="/lizardfs/flaviav/mouse/assembly_D/GWAS/out_bwa"
OUTPUT_GAF="/lizardfs/flaviav/mouse/assembly_D/GWAS/out_gfainj"

# Get a list of all folders in DATA_DIR
folders=($(find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d))

# Loop through each folder in the list
for folder_path in "${folders[@]}"; do
    # Get the folder name
    folder_name=$(basename "$folder_path")

    # Set the read file paths
    READ1="${folder_path}/${folder_name}_S*_R1_*.fastq.gz"
    READ2="${folder_path}/${folder_name}_S*_R2_*.fastq.gz"

    # Create a job script for the current folder
    job_script="job_${folder_name}.sh"
    cat > "$job_script" << EOF
#!/bin/bash

#SBATCH --job-name=${folder_name}
#SBATCH --output=${folder_name}.out
#SBATCH --error=${folder_name}.err

# Run the first command
$BWA_PATH mem $REFERENCE $READ1 $READ2 -t 40 | samtools view -hb | samtools sort -@20 > "${OUTPUT_BWA}/${folder_name}.bam"

# Run the second command
$GFAINJECT_PATH --gfa "${REFERENCE_pggb}.bf3285f.eb0f3d3.867196c.smooth.final.gfa" --bam "${OUTPUT_BWA}/${folder_name}.bam" > "${OUTPUT_GAF}/${folder_name}_inject.gaf"

# Run the third command
/lizardfs/guarracino/tools/gafpack/target/release/gafpack --graph "${REFERENCE_pggb}.bf3285f.eb0f3d3.867196c.smooth.final.gfa" --alignments "${OUTPUT_GAF}/${folder_name}_inject.gaf" > "${OUTPUT_GAF}/${folder_name}_pack.tsv"

rm "${OUTPUT_GAF}/${folder_name}_inject.gaf"
EOF

    # Submit the job script
    sbatch "$job_script"
done
