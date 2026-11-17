# RNArefine examples

This directory provides **ready-to-use input structures and reference outputs** for two RNArefine workflows:

- **Backbone-based interaction prediction and refinement** (`Backbone/`)
- **C4′-trace-based interaction prediction and refinement** (`C4'_trace/`)

Each example includes an RNA structure, predicted base-pairing and base-stacking interactions, intermediate structures, and a final refined structure. You can run either example and compare your output with the supplied reference files.

## 1. Prerequisites

First, install RNArefine and download its pretrained models and parameter files as described in the [main RNArefine README](https://github.com/ShoTsukiyama/RNArefine#1-installation).

Place the `Example/` directory in the **root of the RNArefine repository** (alongside `RNArefinebp/` and `environment.yaml`). Run the following commands from that repository root. Activate the installed Conda environment and make the RNArefine executable available:

```bash
conda activate rnarefine-predict
export PATH="$PWD/.local/bin:$PATH"
```

## 2. Example files

```text
Example/
├── Backbone/
│   ├── input.pdb
│   └── output/
│       ├── RNArefine_bp.txt
│       ├── RNArefine_stack.txt
│       ├── superimpose.pdb
│       ├── mc.pdb
│       └── rnaref3d.pdb
├── C4'_trace/
│   ├── input.pdb
│   └── output/
│       ├── RNArefine_bp.txt
│       ├── RNArefine_stack.txt
│       ├── superimpose.pdb
│       ├── mc.pdb
│       └── rnaref3d.pdb
└── README.md
```

| File | Description |
| --- | --- |
| `input.pdb` | Input RNA 3D structure. |
| `output/RNArefine_bp.txt` | Predicted base-pairing interactions for the example. |
| `output/RNArefine_stack.txt` | Predicted base-stacking interactions for the example. |
| `output/superimpose.pdb` | Reconstructed full-atom structure before optimization. |
| `output/mc.pdb` | Intermediate structure after Monte Carlo (MC) refinement. |
| `output/rnaref3d.pdb` | Final structure after refinement. |

The interaction files contain comma-separated records without a header, in the format:

```text
chain_1,chain_2,residue_1,residue_2,score
```

Here, `score` is the predicted interaction score. The provided input structures are different RNA targets: `Backbone/input.pdb` contains **392 nucleotides** and uses a reduced-atom representation; `C4'_trace/input.pdb` contains **119 nucleotides** and includes all atoms. In the latter example, the **`-r` option instructs the interaction-prediction model to use C4′-trace information**, even though the supplied PDB contains additional atoms.

## 3. Run the examples

The commands below first predict RNA interactions and then use them as restraints during structure refinement. **New outputs are saved under `run/`**, leaving the provided reference results in `output/` unchanged.

### A. Backbone-based refinement

**Step 1 — Predict base-pairing and base-stacking interactions**

```bash
python RNArefinebp/main.py \
  --pdb Example/Backbone/input.pdb \
  -o Example/Backbone/run
```

**Step 2 — Refine the RNA structure**

```bash
RNArefine \
  -i Example/Backbone/input.pdb \
  -o Example/Backbone/run \
  -bp Example/Backbone/run/RNArefine_bp.txt \
  -st Example/Backbone/run/RNArefine_stack.txt
```

The final structure is saved as `Example/Backbone/run/rnaref3d.pdb`.

### B. C4′-trace-based refinement

**Step 1 — Predict interactions using the C4′-trace mode (`-r`)**

```bash
python RNArefinebp/main.py \
  --pdb "Example/C4'_trace/input.pdb" \
  -o "Example/C4'_trace/run" \
  -r
```

**Step 2 — Refine the RNA structure**

```bash
RNArefine \
  -i "Example/C4'_trace/input.pdb" \
  -o "Example/C4'_trace/run" \
  -bp "Example/C4'_trace/run/RNArefine_bp.txt" \
  -st "Example/C4'_trace/run/RNArefine_stack.txt"
```

The final structure is saved as `Example/C4'_trace/run/rnaref3d.pdb`.

**Important:** The `-r` flag applies to **interaction prediction** (`RNArefinebp/main.py`). It is not an option of the `RNArefine` refinement executable. Quotation marks are needed around paths containing the apostrophe in `C4'_trace`.

## 4. Check the results

After running either example:

1. Open `run/rnaref3d.pdb` in a molecular viewer such as PyMOL or UCSF ChimeraX and compare it with the corresponding `input.pdb` and `output/rnaref3d.pdb`.
2. Inspect `run/RNArefine_bp.txt` and `run/RNArefine_stack.txt`, and compare them with the reference interaction files in `output/`.
3. Optionally inspect `run/superimpose.pdb` and `run/mc.pdb` to see the intermediate structures.

**Reproducibility:** MC refinement involves stochastic sampling, so a new run may not produce exactly the same atomic coordinates as the reference output. Results may also differ if runtime settings are changed.

By default, RNArefine allows up to **7,200 seconds for MC** (`-mtime`) and **3,600 seconds for L-BFGS** (`-ltime`). For a quick functionality test, you may append the following options to either `RNArefine` command:

```bash
-mtime 120 -ltime 60
```

Shortened runs are intended for testing the workflow, **not** for reproducing the supplied reference structures.

For installation details, all available options, and troubleshooting, see the [main RNArefine documentation](https://github.com/ShoTsukiyama/RNArefine).
