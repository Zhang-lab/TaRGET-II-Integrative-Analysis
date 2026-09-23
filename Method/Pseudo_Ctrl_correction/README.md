## Normalization and Pseudo-control–based Batch Correction

RNA-seq and ATAC-seq datasets were processed using a unified normalization and batch-correction framework to minimize technical variation associated with data production centers while preserving exposure-related biological signals.

### Step 1. Within-center normalization and removal of unwanted variation

RNA-seq gene-count matrices and ATAC-seq peak-count matrices were normalized using relative log expression (RLE) normalization. Residual unwanted variation was subsequently estimated using `RUVr` from the `RUVSeq` package.

Three unwanted factors (`k = 3`) were included to capture technical variation unrelated to exposure, including variation associated with tissue dissection, library preparation, and other experimental factors.

**Script**

- `1. RLE_normalization_and_RUVr_function.R`

### Step 2. Assembly of exposure-specific datasets

RUVr-adjusted matrices and the corresponding raw count matrices were consolidated across exposure groups to generate harmonized datasets for pseudo-control construction and subsequent cross-center correction.

**Scripts**

- `2.1 combine_Exposure_RUVr_output.R`
- `2.2 combine_Exposure_raw_table.R`

### Step 3. Construction of pseudo-control references

Pseudo-control references were generated separately for matched biological strata defined by relevant covariates, including tissue, age, and sex.

Within each stratum, control samples from multiple production centers were used to estimate a common cross-center reference profile. Center-specific control profiles were then compared with the corresponding pseudo-control reference to quantify systematic production-center effects.

**Script**

- `3. build_pseudo_control.R`

### Step 4. Pseudo-control–based production-center correction

For each production center, the deviation between the center-specific control profile and the corresponding pseudo-control reference was used to estimate a center-specific correction factor.

This correction was applied consistently to both control and exposed samples generated at the same center. This procedure reduces systematic center-associated variation while retaining the relative molecular differences between exposed and control samples.

**Script**

- `4. Run_correction_based_on_pseudo_ctrl.R`

### Step 5. Differential analysis

The batch-corrected matrices were used as input for downstream differential analyses.

For RNA-seq, differentially expressed genes (DEGs) were identified using:

- [`Run_DEG_analysis.R`](https://github.com/Zhang-lab/TaRGET-II-Integrative-Analysis/blob/main/Method/Pseudo_Ctrl_correction/RNA/5.%20Run_DEG_analysis.R)

For ATAC-seq, differentially accessible regions (DARs) were identified using:

- [`Run_DAR_analysis.R`](https://github.com/Zhang-lab/TaRGET-II-Integrative-Analysis/blob/main/Method/Pseudo_Ctrl_correction/ATAC/5.%20Run_DAR_analysis.R)
