# PolyScript

Polymer chemistry sequence generation, classification, validation, and depolymerization.

## Overview

PolyScript is a Python toolkit for computational polymer chemistry. It provides:

- **Polymerizer** — Generate polymer sequences from monomer inputs with reaction rules
- **Classifier** — Classify monomers by functional group types
- **Depolymerizer** — Depolymerize polymer sequences back to constituent monomers
- **BRICS** — Fragment molecules using the BRICS retrosynthetic rules

## Quick Start

```python
from polyscript.polymerizer import Polymerizer

# Initialize with default reaction rules
poly = Polymerizer()
```

```python
from polyscript.classifier import PolyScriptClassifier

# Classify monomers
classifier = PolyScriptClassifier(MinFG=2, MaxFG=4)
```

## Installation

```bash
pip install polyscript
```

## Dependencies

- RDKit — Cheminformatics toolkit
- Pandas — Data manipulation
- MPI (optional) — Parallel execution support
