# Assembly DBA/2J using Hi-C and Revio

### 0) Reports about data:
- Report for DBA/2J Hi-C [here](https://www.dropbox.com/home/97_Dovetail_Sequencing?di=left_nav_browse)
- Report for DBA/2J Revio [here](https://drive.google.com/file/d/1y35q3QDpy7X2XD-eef6Wuq3y3-PlxStf/view?usp=drive_link)

### 1) Assembly for DBA/2J using Revio data:

Data used:
-  Revio on octopus server, folder in which there are both runs: *flaviav@octopus02:/lizardfs/flaviav/data/DBA2J_revio* 

### 1.1). Convert bam file into fastq

I received a bam format, the first step is convert it into a fastq file.
`samtools fastq MouseStrainD2_TBG_4829_1.hifi_reads.bam > MouseStrainD2_TBG_4829_1.hifi_reads.fastq`

#### 1.2) Quality controls on raw data

- Check the quality, the length of the reads for both runs.
- Check if there are present adapters and if we need to trim the reads, with FASTQC for each run:
```
fastqc --threads=40 -o /scratch TBG-4829_m84078_231116_130625_s1.hifi_reads.default.fastq.gz
```

- Meryl and Genome scope. GenomeScope is an open-source web tool to rapidly estimate the overall characteristics of a genome, including genome size, heterozygosity rate and repeat content from unprocessed short reads.
```
#!/bin/bash
##To decide the best number of k-mers 
#$MERQURY/best_k.sh 2728222451                                                                                                 
#genome: 2728222451  #tolerable collision rate: 0.001   #20.6548

$meryl count k=21 /lizardfs/flaviav/data/DBA2J_revio/run2/MouseStrainD2_TBG_4829_1.hifi_reads.fastq.gz output $hifi.run2.meryl
$meryl count k=21 /lizardfs/flaviav/data/DBA2J_revio/TBG-4829_m84078_231116_130625_s1.hifi_reads.default.fastq.gz output $hifi.run1.meryl
```
Then upload both on GenomeScope.

### 2) Assembly using Hifiasm and using only Revio data
I used [Hifiasm](https://hifiasm.readthedocs.io/en/latest/index.html) tool, that works better with HiFi data.  I specified: -z 20 (to trim reads both ends) and -l0 (for haploid assembly).
```
$hifiasm/hifiasm-1ac574adc78fbdaed2d2dcd49d5ea3deed7478de -o /scratch/DBA2J_onlyhifi -z 20 -l0 -t40 /lizardfs/flaviav/data/DBA2J_revio/TBG-4829_m84078_231116_130625_s1.hifi_reads.default.fastq.gz /lizardfs/flaviav/data/DBA2J_revio/run2/MouseStrainD2_TBG_4829_1.hifi_reads.fastq.gz
```
### 2.1) Stats on the assembly:

Stats results: *flaviav@octopus02:/lizardfs/flaviav/mouse/assembly_D/stats/*

- Quast:

```
sbatch -p workers -x octopus07,octopus10 -c 48 --wrap "hostname; cd /scratch; python /lizardfs/flaviav/tools/quast-5.0.2/quast.py -t 30 -r /lizardfs/flaviav/mouse/148strains/UCSC_mm10.fa /lizardfs/flaviav/mouse/assembly_D/DBA2J.asm.hic.p_ctg.fa -o /lizardfs/flaviav/mouse/assembly_D/stats/"
```
- Compleasm:

Download mammalia’s BUSCO genes and run Compleasm:
```
mkdir -p busco/databases
/lizardfs/flaviav/tools/conda/compleasm/bin/compleasm download -L /lizardfs/flaviav/mouse/assembly_D/busco/databases mammalia
sbatch -p workers -x octopus07,octopus10 -c 48 --wrap "hostname; cd /scratch; /lizardfs/flaviav/tools/conda/compleasm/bin/compleasm run -a /lizardfs/flaviav/mouse/assembly_D/DBA2J.asm.hic.p_ctg.fa -o /lizardfs/flaviav/mouse/assembly_D/DBA2J -t 48 -l mammalia -L /lizardfs/flaviav/mouse/assembly_D/busco/databases; mv /lizardfs/flaviav/mouse/assembly_D/DBA2J/summary.txt /lizardfs/flaviav/mouse/assembly_D/DBA2J/summary.mammalia.DBA2J.txt"
```
- gfastats
```
$gfastats -j 35 --stats --seq-report "$input_path_fasta/DBA2J.asm.hic.p_ctg.fa"
```

### 3) Evaluation on the assembly:

- Map the assembly vs the reference genome. This step is highly recommended if a chro-mosome-level reference genome is available.
```
sbatch -p workers --job-name asvsref -c 48 --wrap "hostname; cd /scratch; /lizardfs/flaviav/tools/winnowmap -W /lizardfs/flaviav/mouse/assembly_BXD_ont/BXD101/repetitive_k15.txt --MD -Ha -x map-pb /lizardfs/erikg/mouse/assemblies/Mus_musculus.GRCm39.dna.toplevel.fa.gz /lizardfs/flaviav/mouse/assembly_D/3_hifiasm_hifi/DBA2J.asm.hic.p_ctg.fa | samtools view -hb | samtools sort -@20 > /scratch/Dvsref.bam"
```

- Scaffold step. If a chromosome-level genome assembly is already available, you can connect your contigs into larger chunks using homology between your contigs and the chromosomes. Such larger chunks with unidentified gaps are referred to as ‘‘scaffolds.’
```
ragtag.py scaffold -o /scratch/out_1 /lizardfs/erikg/mouse/assemblies/Mus_musculus.GRCm39.dna.toplevel.fa.gz /lizardfs/flaviav/mouse/assembly_D/3_hifiasm_hifi/DBA2J.asm.hic.p_ctg.fa

```

- Raw read vs the assembly
```
#!/bin/bash
query=$1
target=$2
minimap2 -x map-pb -k15 -a -t 36 $target $query | samtools view -hb | samtools sort -@20 >  /scratch/raw_as_minimap.bam
```

