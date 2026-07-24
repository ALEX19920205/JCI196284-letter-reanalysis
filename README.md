# Verification code for a Letter concerning JCI196284

This repository contains the verification code accompanying a Letter concerning:

Bandaru R, et al. *Genome-wide variation in cell-free DNA end-motif entropy
predicts immunotherapy response in head and neck cancer.* J Clin Invest.
2026;136(13):e196284.

## Contents

- `verify_table1_and_baseline_v2_20260724.py`: reproduces the patient-level
  isolation-batch table and the exploratory isolation-date-only baseline.
- `requirements_verify.txt`: Python package versions used for verification.

No patient-level source data are stored in this repository. The script reads
public files from the authors' repository at the fixed commit below.

## Public source

- Repository: `https://github.com/epifluidlab/headneck`
- Commit: `dcc7a36a46f3d36e34992b08db8a57f0d24aef85`

## Reproduction

```bash
git clone https://github.com/epifluidlab/headneck.git
git clone <URL-OF-THIS-REPOSITORY>

cd headneck
git checkout dcc7a36a46f3d36e34992b08db8a57f0d24aef85

python3 -m venv .venv
source .venv/bin/activate
pip install -r ../JCI196284-letter-reanalysis/requirements_verify.txt
python ../JCI196284-letter-reanalysis/verify_table1_and_baseline_v2_20260724.py
```

## Expected output

- Four of nine isolation batches are single-class.
- Three multi-patient single-class batches cover 15 patients.
- Date-only baseline, 68 patients: AUC 0.76 ± 0.16.
- Date-only baseline, 58-patient analysis set: AUC 0.71 ± 0.18.

The date-only analysis is exploratory and unpaired with the reported rMDS
results because the exact published resampling indices were not available.
