# Analysis Code

This repository contains the analysis code submitted for peer review of
the manuscript:

**"Insulin response and protein intake are associated with distinct aspects
of weight loss during obesity treatment"**


---

## What is in this repository

| File | Description |
|---|---|
| `analysis_code_submission.R` | R script reproducing all tables and figures in the manuscript |
| `sample_data_for_review.csv` | De-identified sample dataset for verifying code execution |
| `column_dictionary.pdf` | Description of all variables in the dataset |
| `session_info.txt` | R session information (version and package details) |

> **Note:** The sample dataset is a de-identified subset provided for
> code verification only. Output values will differ from the manuscript;
> the full dataset is available upon request (see below).

---

## Requirements

- R version 4.5.1
- All required packages are loaded automatically via `pacman::p_load()`
  at the top of the script. Install `pacman` first if needed:
  ```r
  install.packages("pacman")
  ```

---

## How to run

1. Open `analysis_code_submission.R`
2. Edit the two path variables near the top of the script (lines ~49–50):
   ```r
   DATA_DIR <- ""   # folder containing sample_data_for_review.csv
   OUT_DIR  <- ""   # folder where outputs will be written
   ```
3. Run the script from top to bottom. All tables (XLSX) and figures (PDF)
   are written to `OUT_DIR`.

---

## How the code maps to the manuscript

| Manuscript element | Output file |
|---|---|
| Table 1 | `table1_baseline_clinical_characteristics.xlsx` |
| Fig. 2 (a–m) | `Fig2a–m_*.pdf` |
| Fig. 3 (a–c) | `Fig3a–c_*.pdf` |
| Fig. 4 (a–h) | `Fig4a–h_*.pdf` |
| Fig. 5 (b–e) | `Fig5b–e_*.pdf` |
| Extended Data Fig. 1 | `ExtDataFig1a–h_*.pdf` |
| Extended Data Fig. 2 (a–d) | `ExtDataFig2a–d_*.pdf` |
| Supplementary Tables 1–9 | `Sup_data_table*.xlsx` / `lollipop1_data.csv` |

> **Note:** Fig. 1 and Fig. 5a are schematic diagrams not generated
> by this script.

---

## Data availability

The full analytical dataset cannot be shared publicly due to patient
privacy regulations (Taiwan's Personal Data Protection Act) and
institutional ethics requirements (TSGHIRB No. B202505071).

Reviewers requiring access to the full dataset for verification
purposes may contact the corresponding authors directly.

---

## Corresponding authors

- Andy P. Tsai (andytsai@stanford.edu)
- Yen-Jen Sung (yjsungmstw@gmail.com)
