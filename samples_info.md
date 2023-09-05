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

5) Assembly with verkko

```
mkdir -p /lizardfs/guarracino/mouse/assemblies
cd /lizardfs/guarracino/mouse/assemblies

conda activate /lizardfs/guarracino/condatools/verkko/1.4.1/
sbatch -c 4 -p headnode --job-name verkko-BXD1 --wrap "cd /scratch; \time -v verkko --local-cpus 4 -d BXD1.asm --hifi /lizardfs/guarracino/mouse/data/hifi/m84137_230816_004807_s3.hifi_reads.bc2021.fq.gz --nano /lizardfs/flaviav/mouse_ont/fastq/BXD1.fastq.gz"
```
