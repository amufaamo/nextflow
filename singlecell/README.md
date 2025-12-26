```markdown
# Single-Cell RNA-Seq Pipeline (Pipseeker + Seurat)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![Seurat](https://img.shields.io/badge/R-Seurat-blue.svg)](https://satijalab.org/seurat/)

This is a Single-Cell RNA-Seq analysis pipeline built with Nextflow.
It processes 10x Genomics data using **Pipseeker** for alignment/quantification and **Seurat** for clustering and visualization.

## 📦 Requirements

* **Nextflow**: (v22.10.0 or higher recommended)
* **Apptainer** (formerly Singularity) or **Docker** (for Seurat)
* **Pipseeker Binary**: (Local installation required)
    * Since Pipseeker is proprietary software, it is executed locally (outside the container). Ensure the binary is accessible.

## 🚀 Usage

This pipeline supports two input modes: **Single Sample (Direct Input)** and **Batch Processing (Sample Sheet)**.

### 1. Single Sample Mode (Direct Input)

Use this mode to quickly analyze a single sample by specifying the FASTQ files directly.

**Command Example:**

```bash
nextflow run nextflow/singlecell/main.nf \
    --r1 /path/to/SampleA_R1.fastq.gz \
    --r2 /path/to/SampleA_R2.fastq.gz \
    --sample_name SampleA \
    --fasta /path/to/genome.fa \
    --gtf /path/to/genes.gtf \
    --outdir results_sc \
    --pipseeker_bin ./pipseeker/pipseeker

```

---

### 2. Batch Processing Mode (Sample Sheet)

Use this mode to analyze multiple samples at once.

**Step 1: Prepare Sample Sheet (`samples.csv`)**
Headers `sample`, `fastq_1`, `fastq_2` are mandatory.

```csv
sample,fastq_1,fastq_2
SampleA,SampleA_R1.fastq.gz,SampleA_R2.fastq.gz
SampleB,SampleB_R1.fastq.gz,SampleB_R2.fastq.gz

```

**Step 2: Run Command**

```bash
nextflow run nextflow/singlecell/main.nf \
    --samplesheet samples.csv \
    --fastq_dir /path/to/fastq_dir \
    --fasta /path/to/genome.fa \
    --gtf /path/to/genes.gtf \
    --outdir results_sc_batch

```

---

## 🛠️ Directory Structure

```text
nextflow/
 ├── module/                  # Shared Modules
 │    ├── pipseeker.nf
 │    └── seurat.nf
 │
 └── singlecell/              # Pipeline Logic
      ├── main.nf             # Main pipeline script
      ├── nextflow.config     # Configuration file
      └── README.md           # This document

```

## ⚠️ Notes

* **Pipseeker Binary**: The parameter `--pipseeker_bin` must point to the valid executable path on your machine.
* **File Naming**: Pipseeker requires specific file naming conventions (e.g., `_R1_`). The pipeline attempts to handle this automatically via symbolic links.

## ⚙️ Arguments

### Input Data (Choose Mode A or B)

**Mode A: Single Sample**
| Parameter | Description | Required |
|-----------|-------------|:--------:|
| `--r1` | Path to Read 1 FASTQ file. | **Yes** |
| `--r2` | Path to Read 2 FASTQ file. | **Yes** |
| `--sample_name` | Sample identifier for output files. | Default: `sample` |

**Mode B: Batch Processing**
| Parameter | Description | Required |
|-----------|-------------|:--------:|
| `--samplesheet` | Path to the CSV file (`sample`, `fastq_1`, `fastq_2`). | **Yes** |
| `--fastq_dir` | Directory containing the FASTQ files listed in CSV. | **Yes** |

### Reference Genomes

| Parameter | Description | Required |
| --- | --- | --- |
| `--fasta` | Path to the reference genome FASTA file. | **Yes** |
| `--gtf` | Path to the GTF annotation file. | **Yes** |
| `--pipseeker_index` | (Optional) Path to pre-built Pipseeker index directory. If omitted, index is built from FASTA/GTF. | No |

### Single Cell Parameters

| Parameter | Description | Default |
| --- | --- | --- |
| `--sc_chemistry` | 10x Genomics Chemistry version (e.g., `V`, `v3`). | `V` |
| `--sc_sensitivity` | Pipseeker sensitivity level for cell filtering (folder name in output). | `3` |
| `--pipseeker_bin` | Path to the local Pipseeker executable. | `./pipseeker/...` |

### System & Output

| Parameter | Description | Default |
| --- | --- | --- |
| `--outdir` | Directory where results will be saved. | `results_sc` |
| `--cpus` | Number of CPUs to use for local execution. | `16` |

```

