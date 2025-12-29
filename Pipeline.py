import os
import subprocess
import glob

sample_list = glob.glob('*')
sample_list

paired_files = {}

for file in sample_list:
    base_name = file.replace('_1.fq.gz', '').replace('_2.fq.gz', '')
    if base_name not in paired_files:
        paired_files[base_name] = []
    paired_files[base_name].append(file)

for base_name in paired_files:
    paired_files[base_name] = sorted(paired_files[base_name], key=lambda x: '_2.fq.gz' in x)

## Trimming ##

trim_galore_path = os.path.expanduser("~/anaconda3/envs/covid/bin/trim_galore")
cutadapt_path = os.path.expanduser("~/anaconda3/envs/covid/bin/cutadapt")

# Trim Galore
for base, files in paired_files.items():
    if len(files) == 2:
        cmd = f"{trim_galore_path} --paired --quality 20 --length 20 --gzip --path_to_cutadapt {cutadapt_path} --cores 40 {files[0]} {files[1]}"
    else:
        cmd = f"{trim_galore_path} --quality 20 --length 20 --gzip --path_to_cutadapt {cutadapt_path} --cores 40 {files[0]}"

    print("Running:", cmd)
    os.system(cmd)

## Create HISAT2 index ##

cmd_hisat2_index = "~/anaconda3/envs/covid/bin/hisat2-build -p 20 ./HISAT2_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./HISAT2_index/reference_index_human"
print(cmd_hisat2_index)
os.system(cmd_hisat2_index)

## HISAT2 Execution ##

hisat2_path = "~/anaconda3/envs/NGS_2/bin/hisat2"
reference_index = "./HISAT2_index/reference_index_human"

output_dir = "HISAT2_result"
os.makedirs(output_dir, exist_ok=True)

paired_files_hisat2 = {}

val_1_files = sorted(glob.glob("*_val_1.fq.gz"))
val_2_files = sorted(glob.glob("*_val_2.fq.gz"))

for file in val_1_files:
    base_name = file.replace("_1_val_1.fq.gz", "").replace("_2_val_1.fq.gz", "")
    paired_files_hisat2[base_name] = [file]

for file in val_2_files:
    base_name = file.replace("_1_val_2.fq.gz", "").replace("_2_val_2.fq.gz", "")
    if base_name in paired_files_hisat2:  
        paired_files_hisat2[base_name].append(file)
    else:
        paired_files_hisat2[base_name] = [file] 

for base, files in paired_files_hisat2.items():
    output_sam = os.path.join(output_dir, f"{base}.sam")

    if len(files) == 2:
        cmd = f"{hisat2_path} -p 30 -x {reference_index} -1 {files[0]} -2 {files[1]} -S {output_sam}"
        print(f"Running HISAT2 for paired-end: {base}")
    else:
        cmd = f"{hisat2_path} -p 30 -x {reference_index} -U {files[0]} -S {output_sam}"
        print(f"Running HISAT2 for single-end (Check file pair!): {base}")

    subprocess.run(cmd, shell=True, check=True)

print("HISAT2 mapping completed!")

## SAM -> BAM ##

bam_output_dir = "bam_results"
os.makedirs(bam_output_dir, exist_ok=True)

sam_files = sorted(glob.glob(os.path.join(output_dir, "*.sam")))

samtools_path = os.path.expanduser('~/anaconda3/envs/NGS_2/bin/samtools')

for sam_file in sam_files:
    base_name = os.path.basename(sam_file).replace(".sam", "")
    bam_file = os.path.join(bam_output_dir, f"{base_name}_sorted.bam")

    cmd_bam = f"{samtools_path} view -@ 20 -bS {sam_file} | {samtools_path} sort -@ 20 -o {bam_file}"
    cmd_index = f"{samtools_path} index {bam_file}"

    print(f"Converting {sam_file} to sorted BAM...")
    subprocess.run(cmd_bam, shell=True, check=True)
    subprocess.run(cmd_index, shell=True, check=True)

print("BAM processing completed!")

## HTseq Execution ##

htseq_path = "~/anaconda3/envs/covid/bin/htseq-count"

gtf_file = "Homo_sapiens.GRCh38.113_clean.gtf"

bam_dir = "bam_results"

count_dir = "htseq_counts"
os.makedirs(count_dir, exist_ok=True)

bam_files = sorted(glob.glob(os.path.join(bam_dir, "*_sorted.bam")))

htseq_path = os.path.expanduser("~/anaconda3/envs/NGS_2/bin/htseq-count")

for bam_file in bam_files:
    sample_name = os.path.basename(bam_file).replace("_sorted.bam", "")
    output_file = os.path.join(count_dir, f"{sample_name}_counts.txt")

    cmd = f"{htseq_path} -f bam -r pos -s no -t exon -i gene_id --additional-attr gene_name -n 20 {bam_file} {gtf_file} > {output_file}"
    
    print(f"Processing {sample_name} ...")
    subprocess.run(cmd, shell=True, check=True)

print("All samples processed successfully!")
