# rafilter2: Rare-kmer based Alignment Filter using O(n log n) LIS Engine

`rafilter2` is a high-performance command-line tool designed to validate and filter genomic alignment results (SAM/PAF). It identifies "rare k-mers" from a reference genome and uses the **Longest Increasing Subsequence (LIS)** algorithm to evaluate the syntenic consistency of alignments, effectively removing false positives in complex or repetitive genomic regions.

## 🚀 Key Features

- **Dual Format Support**: Handles both `SAM` (via pipes) and `PAF` formats.
- **Rare k-mer Fingerprinting**: Filters out repetitive k-mers to focus on unique genomic anchors.
- **O(n log n) LIS Engine**: Uses an efficient mathematical approach to calculate the best collinear path of anchors.
- **Memory Efficiency**: Implements 2-bit sequence compression to minimize RAM usage for query reads in PAF mode.
- **Multithreaded Pipeline**: Leverages a producer-consumer model to maximize throughput on multi-core systems.
- **Zero Dependencies**: Written in standard C++11; easy to compile and deploy.

## 🧠 Core Algorithm

The validation process follows these steps:

1. **Indexing**: The reference genome is scanned to build a dictionary of k-mers that appear $\le n$ times (default: 4).
2. **Anchor Mapping**: For each alignment, the tool extracts shared rare k-mers between the reference interval and the query sequence.
3. **Collinearity Check**:

  - Shared k-mers are treated as 2D anchors $(ref\_pos, query\_pos)$.
  - The LIS algorithm finds the longest subset of anchors that are strictly increasing in both coordinates.
4. **Scoring**:
A Phred-scaled collinearity score is calculated:

$$Score = -10 \times \log_{10}(1 - \frac{LIS\_length}{Total\_Rare\_Kmers})$$

A higher score indicates a more "trustworthy" alignment (scaled 0–60).

## 🛠 Installation

Requires a C++11 compatible compiler (e.g., `g++` or `clang++`).

Bash

```
g++ -O3 -std=c++11 -pthread rafilter.cpp -o rafilter

```

## 📖 Usage

### Command Syntax

Bash

```
rafilter2 [options] <ref.fa> <in.aln> [query.fa]

```

### Options

| Option | Description | Default |
| :--- | :--- | :--- |
| `-k INT` | k-mer length | 21 |
| `-n INT` | Max frequency for a k-mer to be considered "rare" | 4 |
| `-c INT` | LIS score threshold for filtering | 15 |
| `-t INT` | Number of threads | 4 |
| `-m STR` | Output mode: `filter` (remove) or `tag` (append `km:i:score`) | filter |
| `-f STR` | Input format: `sam` or `paf` | auto-detect |
| `-o FILE` | Output file name | stdout |

导出到 Google 表格

### Examples

**1. Filtering SAM stream (with samtools)**

Bash

```
samtools view -h input.bam | rafilter2 -k 21 -n 4 -t 8 -o filtered.sam ref.fasta -

```

**2. Filtering PAF files**
*Note: PAF mode requires the query fasta as the third positional argument.*

Bash

```
./rafilter2 -c 20 -t 16 ref.fasta alignments.paf queries.fasta > high_quality.paf

```

**3. Tagging alignments without removing them**

Bash

```
./rafilter2 -m tag ref.fasta in.sam - > tagged.sam

```

## ⚠️ Important Notes

- **Memory Usage**: The reference genome and rare k-mer map are stored in RAM. For the human genome, expect ~8-16 GB of memory consumption.
- **N Bases**: Any k-mer containing 'N' (or non-ACGT bases) is ignored.
- **PAF Mode**: Entire query sequences are loaded into memory (compressed) to allow random access during validation.

**Would you like me to help you implement support for Gzip compressed inputs or add a detailed performance benchmarking section?**
