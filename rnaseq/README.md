# Universal RNA-Seq Pipeline (Nextflow)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/docker/automated/trinityrnaseq/trinityrnaseq.svg)](https://hub.docker.com/r/trinityrnaseq/trinityrnaseq)

This is an RNA-Seq analysis pipeline built with Nextflow.
It supports both **Reference-based analysis (STAR + featureCounts)** and **De novo analysis (Trinity + Salmon)** modes.

## 📦 Requirements

* **Nextflow**: (v22.10.0 or higher recommended)
* **Apptainer** (formerly Singularity) or **Docker**
    * The pipeline uses **Apptainer** by default.

## 🚀 Usage

### 1. Prepare Sample Sheet (`samples.csv`)
Create a CSV file containing information about your FASTQ files. The headers (`sample`, `fastq_1`, `fastq_2`) are mandatory.

**Example:**
```csv
sample,fastq_1,fastq_2
SampleA,SampleA_R1.fastq.gz,SampleA_R2.fastq.gz
SampleB,SampleB_R1.fastq.gz,SampleB_R2.fastq.gz
````

-----

### 2\. Reference-based Analysis

Use this mode when a reference genome is available (e.g., Human, Mouse).

**Command Example:**

```bash
nextflow run rnaseq/main.nf \
    --samplesheet samples.csv \
    --fastq_dir /path/to/fastq \
    --outdir results_ref \
    --ref_gtf /path/to/genes.gtf \
    --ref_fasta /path/to/genome.fa \
    --cpus 16
```

**Key Options:**

  * `--single_end`: Add this flag for single-end reads (default is paired-end).
  * `--fc_group_features`: The attribute type used by featureCounts (default: `gene_id`).

-----

### 3\. De novo Analysis (De novo Assembly)

Use this mode when no reference genome is available (e.g., non-model organisms).
This workflow performs assembly using **Trinity**, quality assessment with **BUSCO**, and quantification with **Salmon**.

**Command Example:**

```bash
nextflow run rnaseq/main.nf \
    --denovo true \
    --samplesheet samples.csv \
    --fastq_dir /path/to/fastq \
    --outdir results_denovo \
    --busco_lineage vertebrata_odb10 \
    --cpus 16
```

**Mandatory Options:**

  * `--denovo true`: Enables De novo mode.
  * `--busco_lineage`: Specifies the lineage dataset for BUSCO evaluation (e.g., `mammalia_odb10`, `vertebrata_odb10`, `eukaryota_odb10`).

-----

## 🛠️ Directory Structure

```text
rnaseq/
 ├── main.nf              # Main pipeline script
 ├── nextflow.config      # Configuration file
 └── modules/             # Process modules
      ├── fastp.nf
      ├── star.nf
      ├── featurecounts.nf
      ├── trinity.nf
      ├── busco.nf
      └── salmon.nf
```

## ⚠️ Notes

  * **Memory Usage**: De novo assembly (Trinity) requires a significant amount of RAM. We recommend using a machine with at least **200GB** of memory.
  * **Cache**: The first execution may take some time to download the necessary container images.

## ⚙️ Arguments

### Input / Output
| Parameter | Description | Default | Required |
|-----------|-------------|:-------:|:--------:|
| `--samplesheet` | Path to the input sample sheet (CSV). Columns `sample`, `fastq_1`, `fastq_2` are required. | `./samples.csv` | **Yes** |
| `--fastq_dir` | Directory path containing the raw FASTQ files. | `~/fastq` | No |
| `--outdir` | Directory path where results will be saved. | `results` | No |

### Mode & Reads
| Parameter | Description | Default |
|-----------|-------------|:-------:|
| `--denovo` | Set to `true` to run the **De novo assembly** pipeline (Trinity). If `false`, runs Reference-based (STAR). | `false` |
| `--single_end` | Set to `true` for single-end reads. Default is paired-end. | `false` |

### Reference-based Analysis (only when `--denovo false`)
| Parameter | Description | Default | Required |
|-----------|-------------|:-------:|:--------:|
| `--ref_gtf` | Path to the GTF annotation file. | `null` | **Yes** |
| `--ref_fasta` | Path to the reference genome FASTA file. | `null` | **Yes** |
| `--star_index` | Path to a pre-built STAR index directory. If provided, index building is skipped. | `null` | No |
| `--fc_group_features` | The feature attribute used by featureCounts (e.g., `gene_id`, `gene_name`, `gene`). | `gene_id` | No |

### De novo Analysis (only when `--denovo true`)
| Parameter | Description | Default | Required |
|-----------|-------------|:-------:|:--------:|
| `--busco_lineage` | Lineage dataset for BUSCO assessment (e.g., `vertebrata_odb10`, `mammalia_odb10`). | `eukaryota_odb10` | **Yes** |
| `--busco_download_path` | Directory to store downloaded BUSCO datasets. | `./busco_downloads` | No |

### System & Adapters
| Parameter | Description | Default |
|-----------|-------------|:-------:|
| `--cpus` | Number of CPUs to use for multi-threaded processes (STAR, Trinity, etc.). | `8` |
| `--adapter_fasta` | Path to the adapter FASTA file for trimming (if needed). | `~/fasta/adapter.fasta` |