#!/bin/bash

# Set the paths
BWA_PATH="/lizardfs/flaviav/tools/bwa-0.7.17/bwa"
GFAINJECT_PATH="/lizardfs/flaviav/tools/gfainject/target/release/gfainject"
REFERENCE="/lizardfs/flaviav/mouse/assembly_D/8_pggb/D_C.fa.gz"
REFERENCE_pggb="/lizardfs/flaviav/mouse/assembly_D/GWAS/pggb_out_ref/D_C_mm10.fa.gz"
DATA_DIR="/lizardfs/flaviav/data/linked_mouse_fastq"
OUTPUT_BWA="/lizardfs/flaviav/mouse/assembly_D/GWAS/out_bwa"
OUTPUT_GAF="/lizardfs/flaviav/mouse/assembly_D/GWAS/out_gfainj"
JOB_SCRIPTS_DIR="/lizardfs/flaviav/mouse/assembly_D/GWAS/job_scripts"

# Create the job scripts directory if it doesn't exist
mkdir -p "$JOB_SCRIPTS_DIR"

# Find all directories within DATA_DIR and create a job script for each one
find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r folder_path; do
    folder_name=$(basename "$folder_path")
    job_script="${JOB_SCRIPTS_DIR}/${folder_name}.sh"

    # Create the job script
    cat << EOF > "$job_script"
#!/bin/bash

# Set the read file paths
READ1="${folder_path}/${folder_name}_S*_R1_*.fastq.gz"
READ2="${folder_path}/${folder_name}_S*_R2_*.fastq.gz"

# Run the first command
$BWA_PATH mem $REFERENCE \$READ1 \$READ2 -t 40 | samtools view -hb | samtools sort -@20 > "${OUTPUT_BWA}/${folder_name}.bam"

# Run the second command
$GFAINJECT_PATH --gfa "${REFERENCE_pggb}.bf3285f.eb0f3d3.867196c.smooth.final.gfa" --bam "${OUTPUT_BWA}/${folder_name}.bam" > "${OUTPUT_GAF}/${folder_name}_inject.gaf"

# Run the third command
/lizardfs/guarracino/tools/gafpack/target/release/gafpack --graph "${REFERENCE_pggb}.bf3285f.eb0f3d3.867196c.smooth.final.gfa" --alignments "${OUTPUT_GAF}/${folder_name}_inject.gaf" > "${OUTPUT_GAF}/${folder_name}_pack.tsv"

rm "${OUTPUT_GAF}/${folder_name}_inject.gaf"
EOF

    # Submit the job script using sbatch
    sbatch "$job_script"
done
