# Extender

The extender grows a polymer from one or more monomer structures by applying
the reaction recipe registered for its polymer family. The registry
(`polyscript.extender.map.rxn_map`) is pure data; each entry selects a driver
and supplies reaction SMARTS.

## Mechanism classes

| `type` | driver | growth event |
| --- | --- | --- |
| `addition` | `ChainGrowthReaction` | wildcard propagation; initiation and termination provide the end groups |
| `step_growth` | `StepGrowthReaction` | condensation of one chain functional group with one monomer functional group |
| `ring_opening` | `RingOpeningReaction` | open the cyclic monomer (or a balanced anhydride + epoxide + water unit), then propagate the opened unit with a ROP-specific coupling |
| `romp` | `RingOpeningMetathesisReaction` | programmatic ring cut and insertion of one cycloalkene per propagation |
| `ring_open` | `StepGrowthReaction` | legacy open-then-condense form with a generic coupling |

`reaction_config(key)` is the single entry point for building driver kwargs. It
injects the map key as `ring_open_smarts` for `ring_open`/`ring_opening` and as
`ring_bond_smarts` for `romp`.

## Ring-opening polymerisation

Ring-opening polymerisation (ROP), ring-opening copolymerisation (ROCOP), and
ring-opening metathesis polymerisation (ROMP) are insertion/chain-growth
processes. A cyclic monomer is not an `AB` monomer even though its opened form
can look like one.

The registry uses `ring_opening` to open the ring into its real open form and
then propagate the opened unit:

* **ROP.** Water is the initiator: opening the first lactone/lactam adds water
  across the scissile bond, so the first residue is a hydroxyacid/aminoacid and
  has one more oxygen than the cyclic monomer. Later residues add only the
  monomer's own atoms. `lactone_rop`/`lactam_rop` avoid the generic amide/ester
  guards that rejected the R11 and R12 opened end groups.
* **ROCOP.** The balanced `ROCOP_OPENING_SMARTS` consumes one water molecule:
  the anhydride ring O becomes the ester O, the water O becomes the acid O-H,
  and the epoxide O becomes the terminal alcohol. One cycle consumes one
  anhydride + one epoxide pair.

For PES, succinic anhydride (3 O) + ethylene oxide (1 O) + water (1 O) gives
the first opened residue with 5 O. The extra oxygen is the water initiator.

This is a **stoichiometric/topological** model: it reproduces the terminated
chain directly and does not track an active-end dummy or a separate termination
event. The insertion/active-end variant was prototyped and removed from the
active code path because the topological model produces the same final chains
for the known polymer examples while remaining faster and simpler. ROMP uses a
programmatic ring cut because a variable-size ring cannot be opened by one
generic two-reactant SMARTS in RDKit.

## ROMP

ROMP cuts the metathesis ``C=C`` double bond named by the map key and joins
repeat units through new double bonds at the junctions. `RingOpeningMetathesisReaction` performs the
cut programmatically: initiation opens one ring, caps one side with the
initiator fragment through a double bond and leaves a double-bonded dummy at
the growing end; propagation opens the next ring, forms a double bond to the
chain and moves the dummy; termination replaces the dummy with the terminator
fragment through a double bond. This is necessary because the ring path is
variable-length and cannot be expressed as one generic two-reactant SMARTS.

## API

::: polyscript.extender.chain_growth.ChainGrowthReaction

::: polyscript.extender.step_growth.StepGrowthReaction

::: polyscript.extender.ring_opening.RingOpeningReaction

::: polyscript.extender.romp.RingOpeningMetathesisReaction

::: polyscript.extender.map
