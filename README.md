# Pathogenic Germline Variants, Family History, and Breast Cancer Risk by Menopausal Status

**Statistical analysis and figure-generation code | Stata 18 · R**

This repository accompanies the manuscript *Pathogenic Germline Variants, Family History, and Breast Cancer Risk by Menopausal Status*. It documents analyses of the joint associations of pathogenic germline variants (PGVs) and breast cancer family history (FH) with breast cancer risk among Japanese women, including analyses stratified by menopausal status.

> **Data access:** Individual-level study data are restricted and are **not included** in this repository. Access to the underlying data requires the appropriate institutional approvals.

## Contents

- [Study overview](#study-overview)
- [Repository files](#repository-files)
- [Analysis approach](#analysis-approach)
- [Software requirements](#software-requirements)
- [How to run the scripts](#how-to-run-the-scripts)
- [Data availability](#data-availability)
- [Reproducibility notes](#reproducibility-notes)

## Study overview

The study examines how PGV carrier status and family history of breast cancer are jointly associated with breast cancer risk, with particular attention to differences by menopausal status. The code covers the manuscript's main and supplementary analyses and the generation of selected figures.

## Repository files

| File | Output or purpose |
| --- | --- |
| [`01_Main_tables_and_supplementary_analyses.do`](01_Main_tables_and_supplementary_analyses.do) | Constructs the analytic dataset; defines PGV groups and reference categories; runs analyses for Tables 1–3 and eTables 2–9. Gene-group definitions corresponding to eTable 1 are included. |
| [`02_Figure_1_forest_plot.R`](02_Figure_1_forest_plot.R) | Produces Figure 1 from precomputed odds ratios, confidence intervals, and interaction P values. |
| [`03_Figure_2_PGV_distribution.do`](03_Figure_2_PGV_distribution.do) | Produces Figure 2 from the intermediate analytic dataset. |
| [`04_eFigure_2_PGV_prevalence_by_age.do`](04_eFigure_2_PGV_prevalence_by_age.do) | Performs the prevalence-by-age analysis and produces eFigure 2. |
| [`05_eFigure_3_TNBC_by_BRCA_menopause.do`](05_eFigure_3_TNBC_by_BRCA_menopause.do) | Produces eFigure 3 using aggregate counts from the subtype analyses. |

*eFigure 1 is a participant-selection flow diagram and does not have a separate generation script in this repository.*

## Analysis approach

The main Stata script includes joint-category association models and separate factorial models for interaction analyses. The PGV reference category is defined consistently across gene-group comparisons:

- **Noncarriers:** Participants with no PGV in any of the 27 evaluated genes.
- **Group-specific carriers:** Participants carrying a PGV in the gene group being evaluated.
- Participants with PGVs **only outside** the evaluated gene group are **not** assigned to that group's noncarrier category.

The figure scripts differ in their inputs: Figure 2 and eFigure 2 require restricted study data, whereas Figure 1 and eFigure 3 plot estimates or aggregate counts embedded in their scripts.

## Software requirements

- **Stata 18** for `.do` scripts. Firth penalized logistic regression sensitivity analyses use the user-written `firthlogit` command.
- **R** with `ggplot2`, `dplyr`, `tidyr`, `grid`, and `patchwork` for Figure 1.

R package versions are not pinned.

## How to run the scripts

The scripts use generic input filenames rather than internal study dataset names:

| Local file | Purpose |
| --- | --- |
| `study_input.dta` | Restricted primary study input |
| `subtype_input.dta` | Restricted subtype-analysis input |
| `analysis_intermediate.dta` | Intermediate participant-level dataset created by script 01 and used by script 03 |

Authorized users must provide approved input data with the required variable structure and configure local paths as needed.

**Suggested execution order:**

1. `01_Main_tables_and_supplementary_analyses.do` — main and supplementary analyses; generates `analysis_intermediate.dta`.
2. `02_Figure_1_forest_plot.R` — Figure 1, using embedded estimates.
3. `03_Figure_2_PGV_distribution.do` — Figure 2, using `analysis_intermediate.dta`.
4. `04_eFigure_2_PGV_prevalence_by_age.do` — eFigure 2, using restricted source data.
5. `05_eFigure_3_TNBC_by_BRCA_menopause.do` — eFigure 3, using embedded aggregate counts.

## Data availability

Individual-level genetic, family history, and clinical data are not publicly available because access is subject to participant consent and institutional ethics and data-governance requirements. Requests for de-identified data may be considered through the appropriate institutional review and approval processes. Public access to this code does not grant access to the underlying data.

Neither the original participant-level datasets nor the intermediate participant-level dataset are distributed here.

## Reproducibility notes

The repository documents the analysis specifications and figure-generation procedures. Analyses that require individual-level data cannot be run from this repository alone. Figure 1 and eFigure 3 can be generated from the precomputed values embedded in their scripts, but those scripts do not independently estimate or verify the underlying results.
