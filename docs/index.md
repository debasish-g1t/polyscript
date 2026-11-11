# PolyScript

Polymer chemistry sequence generation, classification, validation, and depolymerization.

## Overview

PolyScript is a Python toolkit for computational polymer chemistry. It provides:

- **Polymerizer** — Generate polymer sequences from monomer inputs with reaction rules
- **Extender** — Grow polymer chains from monomers with addition, step-growth,
  ring-opening (ROP/ROCOP), and ROMP mechanisms
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

## Extender mechanisms

`polyscript.extender.map.rxn_map` registers one recipe per reaction SMARTS and
selects a driver by `type`:

| type | growth event |
| --- | --- |
| `addition` | wildcard chain growth; end groups from initiation/termination |
| `step_growth` | condensation of real monomer functional groups |
| `ring_opening` | open the cyclic monomer (or a balanced anhydride + epoxide + water unit), then propagate the opened unit |
| `romp` | ring-opening metathesis by programmatic ring cut and insertion |
| `ring_open` | legacy open-then-condense form with a generic coupling |

Ring-opening polymerisation is an insertion process, not generic step growth.
The `ring_opening` driver opens the ring into its real open form and propagates
the opened unit with a ROP-specific coupling.  Water is the default initiator;
its oxygen becomes the acid O-H, which is why the first water-initiated residue
has one more oxygen than the cyclic monomer (or the anhydride + epoxide pair).
This is a stoichiometric/topological model: it reproduces the terminated chain
directly without tracking an active-end dummy.  ROMP uses a programmatic ring
cut because a variable-size ring cannot be opened by one generic two-reactant
reaction SMARTS in RDKit. See [Extender](api/extender.md) for the mechanism
details and the terminal-group policy.

## Dependencies

- RDKit — Cheminformatics toolkit
- Pandas — Data manipulation
- MPI (optional) — Parallel execution support
