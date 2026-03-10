# Universal RNA-Seq Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-23aa62.svg)](https://www.nextflow.io/)
[![Apptainer](https://img.shields.io/badge/Apptainer-Supported-blue)](https://apptainer.org/)

本リポジトリは、次世代シーケンサー（NGS）によるRNA-Seqデータを解析するための、**Nextflowベースの自動化パイプライン**です。
初心者の方でも、わずかな設定ファイルとコマンド1つで「学会・論文レベル」の解析結果を得られるように設計されています。

## ✨ パイプラインの特長

このパイプラインには、目的に応じて使い分けられる **2種類** の主要な解析モードが搭載されています。

1. **リファレンスありモード (Reference-based)**
   ヒトやマウスなどの「参照ゲノム」が存在するモデル生物向けの標準的な解析。`fastp ➔ STAR ➔ featureCounts` という王道の流れに加え、インデックスの自動生成機能や、全段階（Pre/Post/マッピング/カウント）の品質管理を一目で確認できる `summary.tsv` 自動生成機能を備えています（`rnaseq_pipeline.nf`）。
2. **De novo 生物向け「真の完全版」モード (Reference-free)**
   参照ゲノムが存在しない非モデル生物向けの、最高峰のアーキテクチャ。
   `Trinity` によるアセンブリから始まり、**発現変動解析（定量ルート）**と**機能アノテーション（アノテーションルート）**を、サーバーリソースが許す限り**完全並列**で同時に走り抜ける超高速・高精度仕様です（`rnaseq/main.nf`）。

---

## 🚀 動作環境

* **OS**: Linux
* **Nextflow**: バージョン 22.10.0 以上
* **Apptainer (旧 Singularity)**: 必須 (各種解析ツールは自動でコンテナから呼び出されるため、ローカルへの面倒なインストールは不要です)

---

## 🔰 初心者向け：使い方ガイド

### Step 1: サンプルシート (`samples.csv`) の作成
解析したいデータの情報を記載したCSVファイルを準備します。一番上の行（ヘッダー）は絶対に文字を変更しないでください。

```csv
sample,fastq_1,fastq_2
SampleA,SampleA_R1.fastq.gz,SampleA_R2.fastq.gz
SampleB,SampleB_R1.fastq.gz,SampleB_R2.fastq.gz
```

### Step 2: 目的の解析を実行する

#### 💡 A. 参照ゲノムがある場合 (モデル生物など)

FASTQファイルが入っているディレクトリ（例: `/home/user/fastq`）、ゲノム配列（`.fa`）、アノテーションファイル（`.gtf`）を指定して実行します。

```bash
nextflow run rnaseq_pipeline.nf \
    --samplesheet ./samples.csv \
    --fastq_dir /path/to/fastq \
    --ref_fasta /path/to/genome.fa \
    --ref_gtf /path/to/genes.gtf \
    --outdir results
```
**🎉 実行後の見どころ**: `results/summary/summary.tsv` を見ると、全サンプルのトリミング前後、マッピング率、カウント数が一覧で確認できます！

#### 💡 B. 参照ゲノムがない場合 (De novo 解析 真の完全版)

De novoモードを使用する場合は `rnaseq/main.nf` を使用し `--denovo true` を指定します。
定量やアノテーションをフルに行いたい場合は、各種データベースへのパスを指定して実行します。これらのパスが指定されると、パイプラインが自動的に並列処理を開始します。

```bash
nextflow run rnaseq/main.nf \
    -profile apptainer \
    --denovo true \
    --samplesheet ./samples.csv \
    --fastq_dir /path/to/fastq \
    --busco_lineage eukaryota_odb10 \
    --eggnog_db /path/to/eggnog_data_dir \
    --nr_db /path/to/nr.dmnd \
    --interpro_db /path/to/interproscan_data \
    --outdir results_denovo
```

> **注意:** アノテーションDB（EggNOG, nr, InterProScan）へのパスを指定しなかった場合、その工程はスキップされ、定量（Salmon ➔ Corset ➔ DESeq2）のみが実行されます。

---

## 🧬 De novo 解析モードの内部フロー

「真の完全版」モードでは、以下の処理が完全自動で進行します。

1. **品質管理:** `fastp` によるアダプタートリミングとQC。
2. **アセンブリ:** `Trinity` によるRNAからの直接転写産物構築。
3. **品質評価:** `BUSCO` によるアセンブリの完全性評価。
4. **【分岐並列スタート！】**
   * **🔴 定量・クラスタリングルート:**
     * `Salmon`: アセンブルされた全配列に対するマッピングと定量。
     * `Corset`: Transcriptレベルの出力を生物学的な「Gene（遺伝子）レベル」に美しく束ねる（これをしないと偽陽性が爆発します）。
     * `DESeq2`: Corsetで得られたカウント結果から、サンプルの発現変動を統計モデリング。
   * **🔵 完璧な機能アノテーションルート (DB指定時のみ作動):**
     * `TransDecoder`: RNA配列から高精度なタンパク質配列(ORF)を予測。
     * `EggNOG-mapper` / `DIAMOND` / `InterProScan`: タンパク質配列を3つのツールに同時投下し、最も強力な意味付け・相同性検索・ドメイン検索を行う。

---

## 📁 出力構造（De novo の場合）

```text
results_denovo/
 ├── fastp/            # トリミングされたクリーンなFASTQファイル
 ├── trinity/          # 最終アセンブリ済み FASTAファイル
 ├── busco/            # BUSCOによる完全性評価スコア
 ├── salmon/           # リードカウント結果
 ├── corset/           # クラスタリング後の counts.txt と clusters.txt
 ├── deseq2/           # deseq2_results.csv (発現変動一覧)
 └── annotation/
      ├── transdecoder/    # 予測されたタンパク質配列(.pep)
      ├── eggnog/          # 遺伝子の機能・パスウェイ情報
      ├── diamond/         # BLASTp互換の類似タンパク検索結果
      └── interproscan/    # PFAM等の機能ドメイン結果
```
