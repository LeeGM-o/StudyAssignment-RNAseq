# RNA-seq pipeline
## Purpose
* There are samples treated with various drugs. This analysis compares two samples to determine which genes are highly or lowly expressed in the sample treated with **Drug A**, and which pathways are enriched in the genes.  
* Since this analysis was conducted at the request of another laboratory, the specific sample names, genes, and pathways will not be disclosed.
* The drugs used are as follows:
> 1. **Drug A**
> 2. Drug B
> 3. Drug C
> 4. Control

## Processing
![image](https://github.com/user-attachments/assets/874b2e23-c994-43bd-9ddc-fea622fe965c)

## Result
* Drug A vs Drug B
![image](https://github.com/user-attachments/assets/ff379ce9-c444-4a04-8295-b2a2732efc5b)

* Drug A vs Drug C
![image](https://github.com/user-attachments/assets/0cfe914e-011c-4d84-9ad6-a8ec2b9f73d0)

* Drug A vs Control
![image](https://github.com/user-attachments/assets/fa66bf30-0ba2-4623-97f0-7758173a2035)

# RNA-seq Pipeline Expansion
## Snakemake
Snakemake is a workflow management system that efficiently automates data analysis pipelines and enhances reproducibility. It defines **rule** using Python-based syntax to identify data dependencies, automatically executes tasks using a **DAG(Directed Acyclic Graph)** approach, and is a universal tool scalable from single cores to clusters and the cloud. It is particularly widely used in bioinformatics, enabling complex analysis processes to be specified in scripts (Snakefiles) for easy management and sharing.
```
pip install snakemake
```
To process input files, a workflow file called a Snakefile is first required. Snakemake operates according to the contents of this workflow file.
The workflow consists of one or more task nodes. Snakemake refers to these nodes as `rule`, declared in the Snakefile using the keyword `rule`. A rule fundamentally has the following structure:
```
rule RULE_NAME:
    input:
    output:
    run, script, shell, notebook:
```
* shell: Describes shell commands
* run: Directly specifies Python scripts
* script: Specifies the path to a Python script file
* notebook: Specifies the path to a Jupyter notebook file   

Once a Snakefile is created, it is executed as follows:
```
$ snakemake -s Snakefile [bulid target] -j [number]
```

## Citation
[https://haje01.github.io/2020/04/21/snakemake-tutorial.html](https://haje01.github.io/2020/04/21/snakemake-tutorial.html)   
[https://snakemake.readthedocs.io/en/stable/](https://snakemake.readthedocs.io/en/stable/)
