# PolyScript

An optimized virtual polymerization framework that utilizes a custom polymer representation to induce better programmability for the chain formation.

## Installation

Install directly from PyPI:

```bash
pip install polyscript-pypi
```

To install from source for development:

```bash
pip install -e .
```

## Reproducible pipeline

The study pipeline is organized into four stages. The table below gives the
entry points and released artifacts; [`Workflow.md`](Workflow.md) documents each
distinct step, its input data, and the evaluation scripts in detail. File sizes,
download examples, and the citation are in
[`DATASET_README.md`](DATASET_README.md).

| Stage | Entry point | Output |
| --- | --- | --- |
| 1. Monomer dataset creation | `runners/1_2_brics_decompose.py` → `runners/1_3_brics_compose.py` → `runners/1_4_monomer_classification.py` → `runners/1_5_compare_class_split.py` | [1.5.0_cls_split_batch_5K.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.5.0_cls_split_batch_5K.zip) |
| 2. Virtual polymerization | `runners/mpi_runner.py` | [2.1.0_ps-run-500_10b_val.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/2.1.0_ps-run-500_10b_val.zip) |
| 3. Extender evaluation | `tests/extract_extender_test_data.py` → `tests/benchmark_extender.py` | [3.2.0_extender_benchmarks.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/3.2.0_extender_benchmarks.zip) |
| 4. VPF resource benchmarking | External scripts (see [`Workflow.md`](Workflow.md)) | [4.1.0_VPF_benchmarking.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/4.1.0_VPF_benchmarking.zip) |

The core path from BRICS output to polymerization:

```bash
# 1.4 Classify the BRICS output (CSV or parquet).
python runners/1_4_monomer_classification.py 1.3.0_all_rxn_data_BRICS.parquet \
    -c product_smiles -o 1.4.0_classified.parquet

# 1.5 Report class counts and split each class into 5K batches.
python runners/1_5_compare_class_split.py 1.4.0_classified.parquet \
    --smiles-column product_smiles -o outputs/cls_split_batch_5K

# 2.1 Run virtual polymerization over the class splits.
mpirun -np 8 python runners/mpi_runner.py ps-run-500 outputs/cls_split_batch_5K \
    --col-name product_smiles -o outputs
```

Steps 1.1–1.3 build the BRICS input from the merged monomer dataset, and the
evaluation stages are run separately. See [`Workflow.md`](Workflow.md) for the
exact files and datasets.

## Documentation

To build the documentation locally, install the docs dependencies and run mkdocs:
```bash
pip install .[docs]
mkdocs build
```

## Testing

```bash
pip install .[test]
pytest tests/
```

## Citation

If you use PolyScript in your research, please cite it using the following BibTeX entry:

```bibtex
@misc{polyscript2026,
  author = {Debasish Mohanty},
  title = {poly(script): The Programmable Virtual Polymerization
  Framework},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/debasish-g1t/polyscript}}
}
```
