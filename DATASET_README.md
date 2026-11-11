---
pretty_name: Poly(script) Data
tags:
  - chemistry
  - polymers
  - virtual-polymerization
  - smiles
  - brics
---

# Poly(script) Data

Data artifacts for the Poly(script) virtual polymerization pipeline. The pipeline
builds a monomer dataset, decomposes and recombines it with BRICS, classifies
the monomers, and runs virtual polymerization. Optional stages evaluate the
extender mechanism and benchmark VPF resources.

The exact commands, inputs, and outputs for each step are in
[Workflow.md](https://github.com/debasish-g1t/polyscript/blob/master/Workflow.md).
The accompanying research paper describes the
chemistry, benchmark design, and analysis.

## Files

| Step | File | Size | Contents |
| --- | --- | --- | --- |
| 1.1 | [1.1.0_merged_smiles_set.csv](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.1.0_merged_smiles_set.csv) | 1.97 MB | Combined monomer SMILES (78,247 rows) |
| 1.2 | [1.2.0_decomp.csv](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.2.0_decomp.csv) | 572 KB | Unique BRICS fragments (26,602 rows) |
| 1.3 | [1.3.0_all_rxn_data_BRICS.parquet](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.3.0_all_rxn_data_BRICS.parquet) | 498 MB | BRICS compose reactions and products (40.4M rows) |
| 1.5 | [1.5.0_cls_split_batch_5K.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/1.5.0_cls_split_batch_5K.zip) | 367 MB | Monomer class splits in 5,000-row batches |
| 2.1 | [2.1.0_ps-run-500_10b_val.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/2.1.0_ps-run-500_10b_val.zip) | 14.6 GB | Virtual polymerization output |
| 3.2 | [3.2.0_extender_benchmarks.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/3.2.0_extender_benchmarks.zip) | 1.17 MB | Extender benchmark results and traces |
| 4.1 | [4.1.0_VPF_benchmarking.zip](https://huggingface.co/datasets/debasishmohanty/poly-script/blob/main/4.1.0_VPF_benchmarking.zip) | 74.2 MB | VPF resource benchmark bundle |

Step 1.4 produces the classified table that Step 1.5 consumes; the workflow
links the split archive from Step 1.5.

## Pipeline

```text
1.1 SMILES -> 1.2 BRICS decompose -> 1.3 BRICS compose
    -> 1.4 classify -> 1.5 class split -> 2.1 virtual polymerization
    -> 3.x extender evaluation -> 4.x resource benchmarking
```

The 1.5 class directories are not mutually exclusive. A monomer can occur in
more than one class because the classifier tests each functional-group type
independently. The `unknown` rows are excluded from the split archive.

## Download

Large files use Git LFS and Xet. Use `huggingface_hub` to download them:

```python
from huggingface_hub import hf_hub_download

path = hf_hub_download(
    repo_id="debasishmohanty/poly-script",
    filename="1.5.0_cls_split_batch_5K.zip",
    repo_type="dataset",
)
```

You can also clone the dataset repository:

```bash
git lfs install
git clone https://huggingface.co/datasets/debasishmohanty/poly-script
```

The 2.1 archive is 14.6 GB. Download only the files that you need.

## Reproduction

The scripts, step order, and released outputs are documented in
[Workflow.md](https://github.com/debasish-g1t/polyscript/blob/master/Workflow.md).
The same workflow is reproduced in the
[PolyScript source repository](https://github.com/debasish-g1t/polyscript)
under `runners/` and `tests/`.

## License

The Poly(script) source repository uses the MIT license. The 4.1 archive includes
third-party code and data from SMiPoly and OpenMacromolecularGenome; those parts
keep their original licenses. Check the source datasets before redistribution.

## Citation

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
