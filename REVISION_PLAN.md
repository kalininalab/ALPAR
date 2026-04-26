# ALPAR Revision & Implementation Plan

This document serves as the master roadmap for addressing reviewer comments and upgrading ALPAR (Automated Learning Pipeline for Antimicrobial Resistance).

---

## 1. High-Level Objectives
- **Substantiate Claims:** Transition from "demonstrated performance" to rigorous head-to-head benchmarking.
- **Mitigate Bias:** Improve handling of accessory genomes and reference-based limitations.
- **Enhance Rigor:** Formalize Quality Control (QC) and expand ML performance reporting.
- **Modernize Infrastructure:** Replace unmaintained dependencies (Prokka) and clarify architecture.

---

## 2. To-Do List by Category

### Phase A: Benchmarking & Validation
*Goal: Address the lack of head-to-head comparisons and limited validation sets.*

- [ ] **Standardized Benchmarking Suite:**
    - Compare ALPAR against **ARIBA**, **AMRFinderPlus**, **Mykrobe**, and **pyseer** pipelines.
    - Measure: Accuracy, MCC, ROC-AUC, F1, Precision, Recall, and Runtime.
- [ ] **External Validation:**
    - Test on at least one new clinical dataset from a different geographical region (e.g., public datasets from NCBI or ENA).
- [ ] **Rules-Based Comparison:**
    - Compare ML performance against a rules-based approach for well-characterized resistance (e.g., Ciprofloxacin/gyrA).

### Phase B: Methodological Enhancements
*Goal: Address reference bias and lineage associations.*

- [ ] **Baseline Model Integration:**
    - Add **Regularized Logistic Regression** (Lasso/Ridge) to the ML pipeline as a performance baseline.
- [ ] **Lineage Correction:**
    - Incorporate population structure (e.g., Mash distances or Pyseer similarity matrices) as features in the ML models to account for lineage bias.
- [ ] **Tool Modernization:**
    - Implement **Bakta** as an alternative/replacement for **Prokka**.
- [ ] **Reference Bias Assessment:**
    - Run the pipeline on the same dataset using 2-3 different reference genomes and quantify the variance in results.

### Phase C: Data QC & Resource Transparency
*Goal: Improve reproducibility and clarify computational requirements.*

- [ ] **Automated QC Module:**
    - Formalize genome exclusion criteria (e.g., length, contig count, N50).
    - Generate a standard `qc_report.txt` for every run.
- [ ] **Scalability Benchmarks:**
    - Profile runtime and RAM usage across varying dataset sizes (e.g., 10, 50, 100, 500, 1000 genomes).
    - Provide a "Minimum Requirements" guide for newcomers.

### Phase D: Architecture & Visualization
*Goal: Clarify modularity and data flow.*

- [ ] **Modular Refactoring:**
    - Break down `amr.py` into smaller sub-modules (Preprocessing, Feature Extraction, Modeling).
- [ ] **Pipeline Visualization:**
    - Redesign Figure 1 to include a clear legend and better distinguish between stages (Preprocessing vs. Analysis).
- [ ] **Automated Plotting:**
    - Add functions to generate Confusion Matrices, ROC curves, and Feature Importance plots automatically in the output folder.

---

## 3. Implementation Pipeline (Step-by-Step)

### Step 1: Metrics & Baselines
*Files: `sr_amr/ml.py`, `sr_amr/ml_lite.py`*
- Update `output_file_writer` to include ROC-AUC, Precision, and Recall.
- Integrate `sklearn.linear_model.LogisticRegressionCV` into the model selection.

### Step 2: Annotation Engine (Bakta)
*Files: `sr_amr/binary_tables.py`, `sr_amr/amr.py`*
- Add `bakta_runner` function.
- Add `--use_bakta` flag to CLI subcommands.

### Step 3: Formalized QC Logic
*Files: `sr_amr/qc.py` (New), `sr_amr/full_automatix.py`*
- Create a dedicated `qc.py` to handle genome filtering before feature extraction.
- Implement explicit thresholds for genome size (e.g., +/- 10% of median) and contig fragmentation.

### Step 4: GWAS vs. ML Comparison
*Files: `sr_amr/gwas.py`*
- Implement a comparison function that overlaps Pyseer hits with ML Feature Importance lists.
- Highlight features identified by both methods in the final report.

---

## 4. File-by-File Change Matrix

| File Path | Action | Reviewer Link |
| :--- | :--- | :--- |
| `sr_amr/ml.py` | Add Logistic Regression & expanded metrics. | Rev 2-1, Rev 1-4 |
| `sr_amr/binary_tables.py` | Integrate Bakta; improve Panaroo handling. | Rev 2-4, Rev 1-2 |
| `sr_amr/amr.py` | Update CLI args; modularize orchestration. | Rev 1-8 |
| `sr_amr/qc.py` | **NEW**: Formalize QC steps and reporting. | Rev 1-6 |
| `sr_amr/gwas.py` | Add GWAS hit vs. ML feature importance plots. | Rev 1-9 |
| `docs/` | Update documentation with resource scaling graphs. | Rev 1-5 |

---

## 5. Success Metrics for Final Submission
1. **Benchmarking Table:** Direct comparison with ARIBA/AMRFinderPlus on Accuracy/MCC.
2. **Resource Curve:** A plot showing ALPAR is usable on standard hardware for smaller datasets.
3. **Clinical Generalizability:** Performance metrics from at least one non-PATRIC dataset.
4. **Baseline Proof:** Statistical comparison showing ML outperforms simple Logistic Regression on complex phenotypes.
