# pksProfiler

**pksProfiler** is a scalable, Nextflow-based computational framework for the detection, quantification, and visualization of the *pks* (polyketide synthase) pathogenicity island from high-throughput sequencing data. It is designed to profile *pks* signals from cancer sequencing datasets (e.g., WGS, WES, RNA-seq) as well as microbiome-enriched metagenomic samples, with gene-level resolution across the 19-gene *clb* biosynthetic cluster.

---

## Key Features

1. Robust detection of *pks* and *clb* genes from host-depleted sequencing reads  
2. Gene-level quantification of the *pks* island using two complementary approaches:
   - Alignment-based profiling against a *pks*-positive *Escherichia coli* reference
   - Hidden Markov model (HMM)–based detection for divergent or fragmented sequences  
3. Visualization of *pks* island coverage and completeness  
4. Taxonomic attribution of *pks*-associated reads to infer putative microbial origins  

---

## Workflow

pksProfiler consists of three modular components:

1. Host read depletion and quality control  
2. *pks* island profiling using alignment-based and HMM-based approaches  
3. Visualization and quantitative summarization of gene-level coverage

---

## Workflow Schematic

<img src="https://github.com/ammalabbasi/pksProfiler/blob/main/workflow_logo/v1.png" width="800"/>

---

## 🚀 Quick Start

## How to run pksProfiler pipeline
1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment
2. Make sure you have the following folders and files in your working directory:
  - main.nf
  - nextflow.config
  - conf/base.config
  - sample.csv

  > You can find these shared files and folder in `/tscc/nfs/home/amabbasi/restricted/pksProfiler/`

3. Next, download the human reference genomes to be used for filtration. We recommend [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/), [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/), and all currently available pangenomes from the [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org). A download script is provided for convenience. Please update the reference paths in the main.nf
```bash
bash scripts/download_references.sh
```

4. Next, create Minimap2 indexes for the previously downloaded reference genomes. A script is provided for convenience to build minimap2 indexes. Please update the index paths in the main.nf
```bash
bash scripts/create_minimap2_indexes.sh
```

5. Prepare your sample.csv file (Example format below):
```
patient,bam
PD56137a,/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/test_data/PD56137a.unmapped.viral.bam
PD56137b,/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/test_data/PD56137b.unmapped.viral.bam
```

7. Request an interactive node and run Nextflow in your working directory under an interactive node:

```
# Node requesting
srun -N 1 -n 1 -c 8 --mem 125G -t 24:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash

# Activate your nextflow conda environment
conda activate env_nf

# Run nextflow
nextflow run main.nf

# If your pipeline terminates with external error, or the interactive node is killed, you can resume your task after setting up the previous steps again with the following command:
nextflow run main.nf -resume

# Optionally, you can recieve an notifiction email on completion with -N flag:
nextflow run main.nf -N your_email@gmail.com
```
8. All outputs, logs, and reports will be written to the RESULT/ directory.

---

## Outputs

- Per-sample *clb* gene abundance tables  
- Normalized *pks* island coverage profiles  
- Visualization-ready files for island completeness assessment  
- Taxonomic summaries of *pks*-associated reads  

---

## Citation


