# Assembly DBA/2J using Hi-C and Revio

### 0) Reports about data:
- Report for DBA/2J Hi-C [here](https://www.dropbox.com/home/97_Dovetail_Sequencing?di=left_nav_browse)
- Report for DBA/2J Revio [here](https://drive.google.com/file/d/1y35q3QDpy7X2XD-eef6Wuq3y3-PlxStf/view?usp=drive_link)

### 1) Assembly for DBA/2J using Hi-C and Revio:

Data used:
-  Hi-C on octopus server: *flaviav@octopus02:/lizardfs/flaviav/data/DBA2J_HiC/Dovetail_DBA_2J*
-  Revio on octopus server: *flaviav@octopus02:/lizardfs/flaviav/data/DBA2J_revio* 

#### 1.1) I used [Hifiasm](https://hifiasm.readthedocs.io/en/latest/index.html) tool, that works better with HiFi data

- Script that I used to run the assembly: *flaviav@octopus02:/lizardfs/flaviav/mouse/assembly_D/script/as_hifiasm.sh*

### 2) Stats on the assembly:

- Quast:

```
sbatch -p workers -x octopus07,octopus10 -c 48 --wrap "hostname; cd /scratch; python /lizardfs/flaviav/tools/quast-5.0.2/quast.py -t 30 -r /lizardfs/flaviav/mouse/148strains/UCSC_mm10.fa /lizardfs/flaviav/mouse/assembly_D/DBA2J.asm.hic.p_ctg.fa -o /lizardfs/flaviav/mouse/assembly_D/stats/"
```
![Screenshot from 2024-01-02 14-06-15](https://github.com/Flavia95/BXD_ONT-PacBio/assets/52487106/e477e0ab-3b28-42b6-a77f-48afabdd8882)

- Compleasm:

Download mammalia’s BUSCO genes and run Compleasm:
```
mkdir -p busco/databases
/lizardfs/flaviav/tools/conda/compleasm/bin/compleasm download -L /lizardfs/flaviav/mouse/assembly_D/busco/databases mammalia
sbatch -p workers -x octopus07,octopus10 -c 48 --wrap "hostname; cd /scratch; /lizardfs/flaviav/tools/conda/compleasm/bin/compleasm run -a /lizardfs/flaviav/mouse/assembly_D/DBA2J.asm.hic.p_ctg.fa -o /lizardfs/flaviav/mouse/assembly_D/DBA2J -t 48 -l mammalia -L /lizardfs/flaviav/mouse/assembly_D/busco/databases; mv /lizardfs/flaviav/mouse/assembly_D/DBA2J/summary.txt /lizardfs/flaviav/mouse/assembly_D/DBA2J/summary.mammalia.DBA2J.txt"
```
