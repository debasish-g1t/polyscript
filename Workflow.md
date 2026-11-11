# Poly(script) Workflow

This document uses ASD-STE100 Simplified Technical English.

## Purpose

This document gives the steps to reproduce the PolyScript pipeline.
Each step has a purpose, an input, a script, and an output.
Run the steps in the sequence in this document.

## Terms

| Term | Meaning |
| --- | --- |
| BRICS | Breaking of Retrosynthetically Interesting Chemical Substructures. |
| Class split | A directory that contains the parquet batches for one monomer class. |
| MPI | Message Passing Interface. |
| Monomer class | A functional-group type that the classifier gives to a monomer. |
| SMARTS | SMILES Arbitrary Target Specification. |
| SMILES | Simplified Molecular-Input Line-Entry System. |
| VPF | Virtual polymerization framework. |

## Prerequisites

1. Install the package and the test dependencies:

   ```bash
   pip install -e ".[test]"
   ```

2. Run all commands from the repository root.

## Pipeline overview

| Stage | Steps | Output data |
| --- | --- | --- |
| 1. Monomer dataset creation | 1.1 to 1.5 | [1.5.0_cls_split_batch_5K.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.5.0_cls_split_batch_5K.zip) |
| 2. Virtual polymerization | 2.1 | [2.1.0_ps-run-500_10b_val.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/2.1.0_ps-run-500_10b_val.zip) |
| 3. Extender evaluation | 3.1 to 3.2 | [3.2.0_extender_benchmarks.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/3.2.0_extender_benchmarks.zip) |
| 4. VPF resource benchmarking | 4.1 | [4.1.0_VPF_benchmarking.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/4.1.0_VPF_benchmarking.zip) |

## Stage 1 – Monomer dataset creation

Stage 1 makes the class split data.
Do Steps 1.1 to 1.5 in the given sequence.

### Step 1.1 – Prepare the existing monomer dataset

**Purpose:** Supply the monomer SMILES data.

**Input:** None.

**Script:** None.

**Output:** The combined monomer dataset
[1.1.0_merged_smiles_set.csv](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.1.0_merged_smiles_set.csv).

**Note:** This step is a data collection step.
The released CSV file is the input for Step 1.2.

### Step 1.2 – Decompose the monomers with BRICS

**Purpose:** Split each monomer into BRICS fragments.

**Input:** The CSV file from Step 1.1.

**Script:** `runners/1_2_brics_decompose.py`

**Procedure:**

1. Open the script.
2. Set `input_path` to the CSV file from Step 1.1.
3. Set `column_name` to the SMILES column. The default is `SMILES`.
4. Set `output_path` to the output CSV file.
5. Run the script:

   ```bash
   python runners/1_2_brics_decompose.py
   ```

**Output:** The BRICS fragments in
[1.2.0_decomp.csv](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.2.0_decomp.csv).
The fragments are in the `decomp_SMILES` column.
The script removes duplicate fragments.

### Step 1.3 – Compose the BRICS fragments

**Purpose:** Make new molecules from the BRICS fragments.
Save the fragment, seed, and product SMILES.
Save the reaction SMARTS.

**Input:** The CSV file from Step 1.2.

**Script:** `runners/1_3_brics_compose.py`

**Procedure:**

1. Put the input CSV file at `./source.csv`.
2. Open the script.
3. Set the file path and the SMILES column if necessary.
4. Run the script:

   ```bash
   python runners/1_3_brics_compose.py
   ```

**Output:** The BRICS compose data in
[1.3.0_all_rxn_data_BRICS.parquet](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.3.0_all_rxn_data_BRICS.parquet).

**Note:** The script writes files only if you set `output_dir` in the
`apply_BRICS` call.

### Step 1.4 – Classify the monomers

**Purpose:** Give each monomer in the BRICS output a functional-group class.

**Input:** The parquet file from Step 1.3.

**Script:** `runners/1_4_monomer_classification.py`

**Procedure:**

1. Run the script:

   ```bash
   python runners/1_4_monomer_classification.py \
       1.3.0_all_rxn_data_BRICS.parquet \
       -c product_smiles -o 1.4.0_classified.parquet
   ```

2. The `-c` value must be the SMILES column name in the input file.
   If you omit `-c`, the script finds the column.

**Output:** `1.4.0_classified.parquet`.
The file has one boolean column for each monomer class.

**Notes:**

- The script writes the canonical SMILES to the input column.
- The script adds two reference rows for CO and HCHO.
  Use `--no-carbonate` to prevent this.

### Step 1.5 – Make the class split

**Purpose:** Count the monomers in each class.
Write 5000-row batches for each class.

**Input:** The parquet file from Step 1.4.

**Script:** `runners/1_5_compare_class_split.py`

**Procedure:**

1. Run the script:

   ```bash
   python runners/1_5_compare_class_split.py \
       1.4.0_classified.parquet \
       --smiles-column product_smiles \
       -o outputs/cls_split_batch_5K
   ```

2. The `--smiles-column` value must be the same as the `-c` value in Step 1.4.

**Output:** The classification report and the class split data:

- `classification_report.csv`
- `classification_report.png`
- `classification_report.txt`
- One directory for each class. Each directory contains parquet batches.

The released class split data is in
[1.5.0_cls_split_batch_5K.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.5.0_cls_split_batch_5K.zip).

**Notes:**

- The input can be one parquet file or a directory of parquet files.
- The `--batch-size` option changes the batch size. The default is 5000.

## Stage 2 – Virtual polymerization

### Step 2.1 – Run the virtual polymerization

**Purpose:** Generate polymer sequences for each monomer class pair.

**Input:** The class split directory from Step 1.5.

**Script:** `runners/mpi_runner.py`

**Procedure:**

1. Run the script:

   ```bash
   mpirun -np 8 python runners/mpi_runner.py \
       ps-run-500 outputs/cls_split_batch_5K \
       --col-name product_smiles -o outputs
   ```

2. The `-np` value must be the number of MPI processes.
3. The `--col-name` value must be the SMILES column name in the split files.

**Output:** The polymer run data in `outputs/ps-run-500`.

The released run data is in
[2.1.0_ps-run-500_10b_val.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/2.1.0_ps-run-500_10b_val.zip).

## Stage 3 – Extender mechanism evaluation

Stage 3 is an optional evaluation stage.
The extender is under active development.
The mechanism can change.

**Note:** Step 3.1 gives a static reference dataset for the paper.

### Step 3.1 – Extract the evaluation subset

**Purpose:** Select polymer rows for the extender evaluation.

**Input:** The existing polymer output in `tests/test_data/sample_outputs/`.
The script reads all parquet files in this directory.

**Script:** `tests/extract_extender_test_data.py`

**Procedure:**

1. Run the script:

   ```bash
   python tests/extract_extender_test_data.py
   ```

**Output:** `tests/test_data/extender_test_data_all.csv`.

### Step 3.2 – Benchmark the extender

**Purpose:** Measure the extender result for each reaction SMARTS key at
different cycle counts.

**Input:** The CSV file from Step 3.1.

**Script:** `tests/benchmark_extender.py`

**Procedure:**

1. Run the script:

   ```bash
   python tests/benchmark_extender.py
   ```

2. Use the options to change the run. For example:
   - `--cycles` sets the cycle counts.
   - `--max-products` sets the frontier cap.
   - `--limit` sets the maximum rows for each reaction.
   - `--class` sets the polymer classes.
   - `--output` sets the result directory.

**Output:** The result directory `benchmarks/results/<timestamp>/`.

This benchmark is external to the Python package.
It is part of the research paper.
The released benchmark data is in
[3.2.0_extender_benchmarks.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/3.2.0_extender_benchmarks.zip).

## Stage 4 – VPF resource benchmarking

### Step 4.1 – Compare the VPF resources

**Purpose:** Compare the run time and the memory of PolyScript with SMiPoly
and OMG.

**Input:** The benchmark data file.

**Script:** The external scripts in the benchmark data file.

**Procedure:**

1. Download
   [4.1.0_VPF_benchmarking.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/4.1.0_VPF_benchmarking.zip).
2. Extract the zip file.
3. Run the scripts in the extracted directory.

**Output:** The benchmark scripts, the output data, and the plots.

The comparison is external to the Python package.
It is part of the research paper.

## Released data

For a dataset overview, file sizes, download examples, and citation, see
[DATASET_README.md](DATASET_README.md).

| Step | File |
| --- | --- |
| 1.1 | [1.1.0_merged_smiles_set.csv](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.1.0_merged_smiles_set.csv) |
| 1.2 | [1.2.0_decomp.csv](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.2.0_decomp.csv) |
| 1.3 | [1.3.0_all_rxn_data_BRICS.parquet](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.3.0_all_rxn_data_BRICS.parquet) |
| 1.5 | [1.5.0_cls_split_batch_5K.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.5.0_cls_split_batch_5K.zip) |
| 2.1 | [2.1.0_ps-run-500_10b_val.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/2.1.0_ps-run-500_10b_val.zip) |
| 3.2 | [3.2.0_extender_benchmarks.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/3.2.0_extender_benchmarks.zip) |
| 4.1 | [4.1.0_VPF_benchmarking.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/4.1.0_VPF_benchmarking.zip) |
