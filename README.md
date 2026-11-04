# cryo-em-res

Feature engineering scripts for relating sequence-derived properties of protein assemblies to the
resolution of their cryo-EM structures. Per-chain data (sequence, length, and predicted B-values,
solvent accessibility and disorder) are aggregated into per-assembly features, each feature is
scored against the native resolution, and a non-redundant subset is selected.

## Pipeline

1. **Data preparation** (`peripheral-code/`, run against `Prepared_Data.xlsx`):
   - `expand_assemblies.py` copies the chain rows of every training assembly into a `Train_Chains` sheet.
   - `export_seqs.py` writes each chain's sequence to `fastas/<row>.fasta`.
   - `batch_profbval_serial.py` / `batch_profbval_parallel.py` run `blastpgp` + PROF + PROFbval on
     each sequence; `batchASAquick.py` runs ASAquick.
   - `Import_ProfBVAL_res.py`, `import_ASAQ_res.py`, `import_iupred_res.py` (and `import_asaquick.py`
     at the top level) read the predictor outputs back into the workbook or a CSV as comma-separated
     per-residue profiles.
2. **Feature computation** — `gen_asm_features.py` reads `Train_Set_equal_dist.csv` (PDB IDs with a
   `Native_Resolution` column) and `Chain_Data.csv` (one row per chain), computes every feature in
   the `compute_feature` registry of `feature_implement.py` (chain counts and lengths, and
   max/average statistics of the IUPred, ASAquick and PROFbval profiles) plus amino-acid fraction
   features (`feature_901`–`feature_940`), and writes `Train_Set_equal_dist.csv_features.csv`.
   `feature_implement.py` reads amino-acid attributes from `amino_acid_data.csv`.
3. **Feature evaluation** — `eval_features.py` computes, for each feature, the Spearman correlation
   with `Native_Resolution` on many random subsets, compares it with a random feature by a paired
   t-test, and writes `..._features.csv_eval_res.csv` plus a histogram `feature_<name>.png` per
   feature. `eval_features_2.py` and `eval_features_3.py` are variants with a larger, fixed set of
   subsets; `eval_features_4.py` is empty.
4. **Feature selection** — `feature_selection.py` greedily picks up to 20 features by effect size and
   p-value while rejecting features whose Spearman correlation with an already selected one exceeds
   0.2, and prints the selection.

`Plottings.py` plots PROFbval profiles of low- vs high-resolution assemblies; `plotting 2.py` is a
matplotlib demo; `get_amino_acid_attributes.py` duplicates `get_amino_acid_data`; `create_fasta.py`
is a one-off converter from `sequences.csv`.

## Running

The code is Python 2 (print statements, `xrange`, `sys.maxint`, `'wb'` CSV output). The Python
libraries used are listed in `requirements.txt`:

```bash
pip install -r requirements.txt
python gen_asm_features.py
python eval_features.py
python feature_selection.py
```

Scripts take no arguments; the input and output file names are set at the top of each file and are
read from the current directory.

## Caveats

- No data is included: `Prepared_Data.xlsx`, `Chain_Data.csv`, `Train_Set_equal_dist.csv`,
  `amino_acid_data.csv` and the predictor outputs must be supplied.
- The `peripheral-code` scripts use Windows-style paths (`fastas\\`, `ASA_Quick_Res\\`) and
  hard-coded locations of `blastpgp`, the NR database and the PROF Perl tools; `batchASAquick.py`
  calls `os(command)` and does not run as written.
- `.idea/` is PyCharm project configuration.
