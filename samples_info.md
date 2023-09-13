# BXD1 and BXD40

1) 10x linked reads:
```
/lizardfs/flaviav/data/linked_mouse_fastq
```
Sample names on here https://docs.google.com/spreadsheets/d/16WzQc1qM-ehDar8UPmVVQArr41QTI5i54aMVVsDm8Kg/edit?usp=sharing

2) Illumina reads, only bam files (ENA/SRA has the bam files but the raw fastqs are not online): 

https://www.ebi.ac.uk/ena/browser/view/ERX9402074

https://www.ebi.ac.uk/ena/browser/view/ERX9393743

3) Nanopore data
   
```
/lizardfs/flaviav/mouse_ont/fastq
```
Informations here: https://www.dropbox.com/scl/fi/sm1afmyyqb7libqnf7fni/Long_read_BXD_stats.xlsx?activeCell=%27Sheet1%27!C8&cloud_editor=excel&dl=0&rlkey=0jsv3pjok1edozggcn1sh8dmv&wdinitialsession=f9c606b0-6a94-4817-8b79-6c8b41cb9ce0&wdrldc=1&wdrldr=AccessTokenExpiredWarning,AccessTokenExpiringWarni&wdrldsc=3 

4) PacBio data

```
mkdir -p /lizardfs/guarracino/mouse/data/hifi
cd /lizardfs/guarracino/mouse/data/hifi
wget -c https://palmerlab.s3.sdsc.edu/BXD1_40_HiFi/m84137_230816_004807_s3.hifi_reads.bc2021.bam
wget -c https://palmerlab.s3.sdsc.edu/BXD1_40_HiFi/m84137_230818_232354_s2.hifi_reads.bc2036.bam
wget -c https://palmerlab.s3.sdsc.edu/BXD1_40_HiFi/README

#Sample Name  File Name                                  HiFi Coverage
#BXD1         m84137_230816_004807_s3.hifi_reads.bc2021  11.9
#BXD40        m84137_230818_232354_s2.hifi_reads.bc2036  20.9
#Sequenced with PacBio Revio
```

5a) Assembly with verkko

```
mkdir -p /lizardfs/guarracino/mouse/assemblies/verkko
cd /lizardfs/guarracino/mouse/assemblies/verkko

conda activate /lizardfs/guarracino/condatools/verkko/1.4.1/
sbatch -c 48 -p workers --job-name BXD1-verkko --wrap "cd /scratch; \time -v verkko --local-cpus 4 -d BXD1.asm --hifi /lizardfs/guarracino/mouse/data/hifi/m84137_230816_004807_s3.hifi_reads.bc2021.fq.gz --nano /lizardfs/flaviav/mouse_ont/fastq/BXD1.fastq.gz"
sbatch -c 48 -p workers --job-name BXD40-verkko --wrap "cd /scratch; \time -v verkko --local-cpus 48 -d BXD40.asm --hifi /lizardfs/guarracino/mouse/data/hifi/m84137_230818_232354_s2.hifi_reads.bc2036.fq.gz --nano /lizardfs/flaviav/mouse_ont/fastq/BXD40.fastq.gz"
```

5b) Assembly with hifiasm

```
mkdir -p /lizardfs/guarracino/mouse/assemblies/hifiasm
cd /lizardfs/guarracino/mouse/assemblies/hifiasm

HIFIASM=/home/guarracino/tools/hifiasm/hifiasm-94a284b4309837417dd9951a5f72a13d513d826e

mkdir -p /lizardfs/guarracino/mouse/assemblies/hifiasm/BXD1
sbatch -c 48 -p workers -p workers --job-name BXD1-hifiasm --wrap "hostname; cd /scratch; \time -v $HIFIASM -o BXD1.asm -t 48 --ul /lizardfs/flaviav/mouse_ont/fastq/BXD1.fastq.gz /lizardfs/guarracino/mouse/data/hifi/m84137_230816_004807_s3.hifi_reads.bc2021.fq.gz; mv BXD1.asm* /lizardfs/guarracino/mouse/assemblies/hifiasm/BXD1"

mkdir -p /lizardfs/guarracino/mouse/assemblies/hifiasm/BXD40
sbatch -c 48 -p workers -p workers --job-name BXD40-hifiasm --wrap "hostname; cd /scratch; \time -v $HIFIASM -o BXD40.asm -t 48 --ul /lizardfs/flaviav/mouse_ont/fastq/BXD40.fastq.gz /lizardfs/guarracino/mouse/data/hifi/m84137_230818_232354_s2.hifi_reads.bc2036.fq.gz --nano /lizardfs/flaviav/mouse_ont/fastq/BXD40.fastq.gz; mv BXD40.asm* /lizardfs/guarracino/mouse/assemblies/hifiasm/BXD40"
```
6) Assembly and correction on ONT data
- First we generate the assembly from the ONT data using wtdbg2.

```
input_dir=/lizardfs/flaviav/mouse_ont/fastq/
wtdg_dir=/lizardfs/flaviav/tools/conda/wtdbg/bin/                                                                
output_dir=/lizardfs/flaviav/mouse_ont/assembly/

$wtdg_dir/wtdbg2 -x ont -g 2.5g -i $input_dir/BXD1.fastq.gz -t 30 -fo $output_dir/mouse_assembly.ont.wtdbg2.asm1
$wtdg_dir/wtpoa-cns -t 30 -i $output_dir/mouse_assembly.ont.wtdbg2.asm1.ctg.lay.gz -fo $output_dir/mouse_assembly.ont.wtdbg2.asm1.raw.fa
minimap2 -t30 -ax map-ont -r2k $output_dir/mouse_assembly.ont.wtdbg2.asm1.raw.fa $input_dir/BXD1.fastq.gz | samtools sort -@30 -o $output_dir/mouse_reads.asm1.ont.bam && samtools view -F0x900 $output_dir/mouse_reads.asm1.ont.bam 
$wtdg_dir/wtpoa-cns -t 30 -d $output_dir/mouse_assembly.ont.wtdbg2.asm1.raw.fa -i - -fo $output_dir/mouse_reads.ont.wtdbg2.asm1.cns.fa
```
