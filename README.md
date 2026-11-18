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
| `column_dictionary.pdf` | Description of all variables in the analysis dataset |
| `session_info.txt` | R session information (version and package details) |

> **Note:** This repository does not contain any participant-level data.
> See [Data availability](#data-availability) for how the data may be
> accessed.

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

1. Obtain the analysis dataset (see [Data availability](#data-availability)).
   The script reads it from `DATA_DIR/sample_data.csv`.
2. Open `analysis_code_submission.R`
3. Edit the two path variables near the top of the script (lines ~49–50):
   ```r
   DATA_DIR <- ""   # folder containing the analysis dataset
   OUT_DIR  <- ""   # folder where outputs will be written
   ```
4. Run the script from top to bottom. All tables (XLSX) and figures (PDF)
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

The participant-level data analysed in this study contain sensitive health
information and are not publicly available, in order to protect the privacy
of study participants.

De-identified data supporting the findings of this study may be made
available from the corresponding author upon reasonable request, subject to
applicable ethical approvals and a data use agreement. The structure and
definitions of all variables are documented in `column_dictionary.pdf`.
