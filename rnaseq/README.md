# Universal RNA-Seq Pipeline (Nextflow)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/docker/automated/trinityrnaseq/trinityrnaseq.svg)](https://hub.docker.com/r/trinityrnaseq/trinityrnaseq)

This is an RNA-Seq analysis pipeline built with Nextflow.
It supports both **Reference-based analysis (STAR + featureCounts)** and **De novo analysis (Trinity + Salmon)** modes.

## 📦 Requirements

* **Nextflow**: (v22.10.0 or higher recommended)
* **Apptainer** (formerly Singularity) or **Docker**
    * The pipeline uses **Apptainer** by default.

---

## 📘 Configuration Guide by Sample Type

Before running the pipeline, check your organism type and data size to select the best parameters.

### 🧬 Case 1: Model Organisms (Human, Mouse, etc.)
* **Genome Size**: Large (> 1GB)
* **Annotation**: Standard GTF (contains `exon` features)

| Parameter | Value | Reason |
| :--- | :--- | :--- |
| `--star_index_nbases` | `14` (Default) | Standard for large genomes. |
| `--gtf_feature_type` | `exon` (Default) | Standard annotations use 'exon'. |
| `--star_bam_sort_ram` | `30.GB` | Sufficient for standard depth. |

### 🦠 Case 2: Microbes (Bacteria, Fungi, Yeast)
* **Genome Size**: Small (< 100MB)
* **Annotation**: **Often lacks `exon` features** (uses `CDS` instead).

| Parameter | Value | Reason |
| :--- | :--- | :--- |
| `--star_index_nbases` | **`11`** | **Critical!** Prevents segfaults on small genomes. |
| `--gtf_feature_type` | **`CDS`** | **Critical!** If your featureCounts result is 0, check your GTF. Many microbial GTFs use 'CDS'. |

### 🧠 Case 3: Deep Sequencing / High Memory Machine
* **Data Size**: Large FASTQ files
* **Hardware**: High RAM availability (> 64GB)

| Parameter | Value | Reason |
| :--- | :--- | :--- |
| `--star_bam_sort_ram` | `60.GB` | Speeds up BAM sorting significantly. Ensure Total RAM > (Index ~30GB + Sort 60GB). |

---

## 🚀 Usage

### 1. Prepare Sample Sheet (`samples.csv`)
Create a CSV file containing information about your FASTQ files. The headers (`sample`, `fastq_1`, `fastq_2`) are mandatory.

**Example:**
```csv
sample,fastq_1,fastq_2
SampleA,SampleA_R1.fastq.gz,SampleA_R2.fastq.gz
SampleB,SampleB_R1.fastq.gz,SampleB_R2.fastq.gz

```

### 2. Reference-based Analysis

**Command Example (for Fungi/Bacteria):**

```bash
nextflow run rnaseq/main.nf \
    --samplesheet samples.csv \
    --fastq_dir /path/to/fastq \
    --outdir results_ref \
    --ref_gtf /path/to/genes.gtf \
    --ref_fasta /path/to/genome.fa \
    --gtf_feature_type CDS \
    --star_index_nbases 11 \
    --cpus 16

```

**Key Options:**

* `--single_end`: Add this flag for single-end reads (default is paired-end).
* `--fc_group_features`: The attribute type used by featureCounts (default: `gene_id`).

### 3. De novo Analysis (De novo Assembly)

Use this mode when no reference genome is available.

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
* `--busco_lineage`: Specifies the lineage dataset for BUSCO evaluation (e.g., `mammalia_odb10`, `eukaryota_odb10`).

---

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

* **Memory Usage (De novo)**: Trinity requires a significant amount of RAM. We recommend using a machine with at least **200GB** of memory.
* **Memory Usage (STAR)**: STAR alignment requires memory for both the **Genome Index** (~30GB for Human/Mouse) and **BAM Sorting** (configurable via `--star_bam_sort_ram`, default 30GB). Ensure your system has `Index Size + Sort RAM` available (e.g., 64GB+ recommended).

## ⚙️ Arguments

### Input / Output

| Parameter | Description | Default | Required |
| --- | --- | --- | --- |
| `--samplesheet` | Path to input sample sheet (CSV). | `./samples.csv` | **Yes** |
| `--fastq_dir` | Directory path containing raw FASTQ files. | `~/fastq` | No |
| `--outdir` | Directory path where results will be saved. | `results` | No |

### Mode & Reads

| Parameter | Description | Default |
| --- | --- | --- |
| `--denovo` | Set to `true` to run **De novo assembly** (Trinity). | `false` |
| `--single_end` | Set to `true` for single-end reads. | `false` |

### Reference-based Analysis (only when `--denovo false`)

| Parameter | Description | Default | Required |
| --- | --- | --- | --- |
| `--ref_gtf` | Path to GTF annotation file. | `null` | **Yes** |
| `--ref_fasta` | Path to reference genome FASTA file. | `null` | **Yes** |
| `--star_index` | Path to pre-built STAR index directory. | `null` | No |
| `--star_index_nbases` | STAR parameter. **Use 11 for small genomes (Bacteria/Yeast).** | `14` | No |
| `--star_bam_sort_ram` | Memory limit for BAM sorting (e.g., `10.GB`, `60.GB`). | `30.GB` | No |
| `--gtf_feature_type` | Feature type to count (e.g., `exon`, `CDS`). **Use `CDS` for microbes.** | `exon` | No |
| `--fc_group_features` | Attribute type for featureCounts (e.g., `gene_id`). | `gene_id` | No |

### De novo Analysis (only when `--denovo true`)

| Parameter | Description | Default | Required |
| --- | --- | --- | --- |
| `--busco_lineage` | Lineage for BUSCO (e.g., `vertebrata_odb10`). | `eukaryota_odb10` | **Yes** |