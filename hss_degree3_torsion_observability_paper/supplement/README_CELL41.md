# Cell 41 — minimal external-data bridge

Cell 40 proves that nine coefficients of the torsion-summed sector plus six
coefficients of one viable nontrivial character sector reconstruct the complete
15-dimensional observable degree-three ambiguity.

Cell 41 turns that theorem into a concrete data request. It reads the exact
Cell-40 atlas, chooses the preferred extra sector, and exports:

- the exact nine torsion-summed q-coefficient positions;
- the exact six extra-sector q-coefficient positions;
- the verified 15 by 15 system and inverse matrices;
- a TSV worksheet for matching published or independently computed data;
- a fillable Sage input template;
- an optional exact reconstruction when the completed input file is present.

## Inputs

Place this existing file in `results/`:

- `degree3_exact_reconstruction_atlas_cell40.sobj`

## Install

Copy:

- `41_degree3_minimal_external_data_bridge.sage` to `code/`
- `run_cell41_external_data_bridge.sh` to the workspace root

## Run

```bash
cd /home/will/research
mamba activate sage
chmod +x run_cell41_external_data_bridge.sh
bash run_cell41_external_data_bridge.sh
```

## Outputs

- `results/degree3_minimal_external_data_bridge_cell41.sobj`
- `results/degree3_minimal_external_data_bridge_cell41.txt`
- `results/degree3_minimal_external_data_request_cell41.tsv`
- `results/degree3_external_data_input_cell41_template.sage`

The template contains fifteen `None` placeholders. Replace each with an exact
Sage value, run the template, then rerun Cell 41. It will reconstruct the
fifteen observable coefficients and verify the residual exactly.

## Scope

This cell packages the exact reconstruction map. It does not identify BKOS
variables with the HSS character labels. That geometric dictionary must be
established independently before published coefficients are inserted.
