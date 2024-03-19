# GEMMA from pangenome

### 0. Build the pangenome graph whole-genome with pggb
```
sbatch -p workers --job-name filterasvsref -c 48 --wrap "hostname; cd /scratch; /home/guarracino/tools/pggb/pggb-a814fa58b370eddd5a099053f39d792bceab9258 -i /lizardfs/flaviav/mouse/assembly_D/8_pggb/D_C.fa.gz -n 2  -o /lizardfs/flaviav/mouse/assembly_D/8_pggb/parts/pggb_out"
```
### 1. Mapping linked-read vs D+B assembly
### 2. Inject the BAM (alignment) file into the graph
### 3. Generated a matrix of coverage using GAFpack tool
```
/lizardfs/flaviav/mouse/assembly_D/GWAS/pre_gwas.sh
```
### 4. Preparing files for GEMMA

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
| node.1     | X    | Y   |   1 |   0   |
| node.2   | X    | Y    |   0 |  0   |


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


  