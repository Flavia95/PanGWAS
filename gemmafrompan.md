# GEMMA from pangenome

### 0. Build the pangenome graph whole-genome with pggb
```
sbatch -p workers --job-name filterasvsref -c 48 --wrap "hostname; cd /scratch; /home/guarracino/tools/pggb/pggb-a814fa58b370eddd5a099053f39d792bceab9258 -i /lizardfs/flaviav/mouse/assembly_D/8_pggb/D_C.fa.gz -n 2  -o /lizardfs/flaviav/mouse/assembly_D/8_pggb/parts/pggb_out"
```
### 1. Mapping linked-read vs D+B assembly
### 2. Inject the BAM (alignment) file into the graph
### 3. Generated a matrix of coverage using GAFpack tool, merge all results in a unique format
```
/lizardfs/flaviav/mouse/assembly_D/GWAS/pre_gwas.sh
# Get the header from the first file of files 
head -n 1 4512-JFI-0346_pack.tsv > merged20.tsv
# Append the content from the second file (excluding the header)
for file in *_pack.tsv; do tail -n +2 "$file"; done >> merged20.tsv
```
### 4. Preparing files for GEMMA
```
python gfapacktogenotype.py out_gfainj/merged20.tsv in_gemma/genonodes.tsv
```

- Matrix coverage from the pangenome for two individuals: 

| #Sample   | node.1    | node.2|
| -------- | -------- | -------- |
|Ind1|100|0|

| #Sample    | node.1    | node.2|
| -------- | -------- | -------- |
|Ind2|40|0|

To starting, I put a threshold, if the coverage is >50 put 1 as genotype otherwise put 0. This is the result:

| #Sample    | node.1    | node.2|
| -------- | -------- | -------- |
|Ind1|1|0|

| #Sample    | node.1    | node.2|
| -------- | -------- | -------- |
|Ind2|0|0|

Then I converted above tables in the format that GEMMA likes, for each node (position), there is the genotype for each individual:

| Nodes | Allele1 | Allele2 |Ind1|Ind2|
| -------- | -------- | -------- |-------|-------|
| node.1     | A    | A   |   1 |   0   |
| node.2   | A    | A    |   0 |  0   |


Nodes in the pangenome graph: 
23,363,074
Absence (o) in the genotypes of nodes                                                   
23,363,075                                                                                                                             
Presence (1) in the genotypes of nodes                                                                                            
19,104,602                                                               


#### GEMMA wants as required inputs a genotype matrix and a phenotype file:

- Genotype matrix, GEMMA ignores allele types. Missing genotypes are as NA.

| RSID(position) | Allele1 | Allele2 |Ind1|Ind2|
| -------- | -------- | -------- |-------|-------|
| rs3144l     | X    | Y   |   1 |   0   |
| rs5173   | X    | Y    |   0 |  0   |

- Phenotype, each line is a number indicating the phenotype file value for each individual, same order of the individuals as above. 1.2 is for the Ind1 and NA for the Ind2.

| Values | 
| -------- | 
| 1.2     | 
| NA   | 

- Same format for the pangenome graph analyses

Optional files:
- SNP annotation file
  
| RSID | Position| Chr|
| -------- | -------- |--------|
| rs3144l     | 12000    | 1   | 
| rs5173   | 13000    | 1   |   

- Covariates file

| Individuals | Position|
| -------- | -------- |
| 1,1     | -1.5    | 
| 1,2   | 0.3    | 

### 5. Run GEMMA on 20 samples
- Check the location of nodes with low p-value.

- Using odgi position: check the location of each node, on which contigs it is
```
echo "#target.graph.pos target.path.pos dist.to.path strand.vs.ref" > results.txt; cat nodespvaluesignig.txt | { while read line; do output=$($ODGI position -i *.final.og -g "$line"); echo "$output" | sed '1d' >> results.txt; done; }
```
- Using mashmap to map the assembly vs the mm10 reference:

  while IFS=, read -r col1 col2 third_col rest; do     third_col=$(echo "$third_col" | sed 's/+//g');     echo "$col1 $third_col"; done < results.txt > poshighvalue.txt

  cut -f1,6 asvsmm10
  awk 'NR==FNR{a[$2];next} $1 in a' poshighvalue.txt asvsmm10_only > final_comparison.txt

```
/lizardfs/flaviav/tools/conda/mashmasp/bin/mashmap -r /lizardfs/flaviav/mouse/148strains/UCSC_mm10.fa -q /lizardfs/flaviav/mouse/assembly_D/3_hifiasm_hifi/DBA2J.asm.onlyhifi.rename.fa -t 35 -o mashmap/asvsmm10
``` 
- Write a script to see where the contigs are, on which chromosome, using odgi position and mashmap output.
- 
Using the reference:
/lizardfs/guarracino/git/odgi/bin/odgi position -t 40 -i /scratch/*.final.og -R /lizardfs/flaviav/mouse/assembly_D/GWAS/pggb_out_ref/refpos.txt -G /lizardfs/flaviav/mouse/assembly_D/GWAS/pggb_out_ref/nodespvaluesignig.txt > /scratch/results1.txt
cat results1.txt | cut -f 2 | cut -f 1 -d ',' | sort | uniq -c | sort -k 1,1nr | less -S
cat results1.txt | grep 'chr1,'
