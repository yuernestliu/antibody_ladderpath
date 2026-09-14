# LPAB Workflow

This repository follows the workflow in the figure below:

1. Clustal Omega alignment and Ladderpath sequence generation
2. Rosetta antibody chain, framework, and six-CDR checks
3. MMseqs2 clustering of `heavy + light` 
4. Cluster-balanced Monte Carlo sampling with weight `sqrt(cluster size)`
5. AlphaFold 3 prediction
6. Official `ipsae.py` calculation and one final CSV

![LPAB workflow](lpab_tools/workflow.jpg)

## Files you need to touch

- `run_lpab.py`: the only program to run.
- `parameters.txt`: paths and adjustable parameters.
- `example_seeds.csv`: example antibody table.
- `example_antigen.fasta`: example antigen sequence.

All implementation details are kept in `lpab_tools/`.

## Install

Create a Python environment and install the small Python part:

```bash
conda create -n lpab python=3.11 -y
conda activate lpab
python -m pip install -r lpab_tools/requirements.txt
```

This installs Ladderpath directly from its official repository. Rosetta,
MMseqs2, Clustal Omega, and AlphaFold 3 are separate scientific programs.
See [lpab_tools/INSTALL.md](lpab_tools/INSTALL.md) for short installation guides.

Edit the three AlphaFold paths and any site-specific executable names in
`parameters.txt`, then check everything:

```bash
python run_lpab.py --check
```

## Run

Run the complete workflow:

```bash
python run_lpab.py --step all
```

Long calculations are resumable by running one stage at a time:

```bash
python run_lpab.py --step generate
python run_lpab.py --step rosetta
python run_lpab.py --step cluster
python run_lpab.py --step sample
python run_lpab.py --step af3-input
python run_lpab.py --step af3-run
python run_lpab.py --step collect
```

Each stage writes a numbered result inside `results/`, so it is clear where a
failed or interrupted run should restart.

## Input

The seed CSV accepts either of these column styles:

```text
name,heavy,light,score
name,VH,VL,Score
```

Higher scores increase the chance that motifs from that seed are reused.
Sequences should contain amino-acid letters only.

## Final output

`results/07_final_scores.csv` contains exactly:

```text
name,heavy,light,iptm,ipsae_min,plddt,ptm,pae,ipae
```

- `iptm` and `ptm` come from AlphaFold 3 summary confidences.
- `plddt` is the mean atom pLDDT.
- `pae` is the mean of the full PAE matrix.
- `ipae` is the mean PAE between different chains.
- `ipsae_min` is the lower of the H-antigen and L-antigen `max ipSAE` values
  produced by the bundled, unmodified `ipsae.py` using the cutoffs in
  `parameters.txt`.

The generated AlphaFold jobs use chain IDs `H`, `L`, and `A`.
