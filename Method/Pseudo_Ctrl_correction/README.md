## Normalization and Pseudo-control Batch Correction

RNA-seq and ATAC-seq data were normalized and corrected for technical variation and production-center effects using the following workflow.

### Step 1. RLE normalization and RUVr correction

RNA-seq gene counts and ATAC-seq peak counts were first normalized using the relative log expression (RLE) method. Residual unwanted technical variation was estimated using `RUVr` from the `RUVSeq` package.

Three unwanted factors (`k = 3`) were used to account for technical variation unrelated to exposure, including variation associated with tissue dissection, library preparation, and other technical sources.

Scripts:

- `1. RLE_normalization_and_RUVr_function.R`

### Step 2. Combine exposure datasets

Normalized/RUVr-corrected datasets and corresponding raw count tables from different exposure groups were combined for subsequent pseudo-control construction.

Scripts:

- `2.1 combine_Exposure_RUVr_output.R`
- `2.2 combine_Exposure_raw_table.R`

### Step 3. Build pseudo-controls

Pseudo-controls were constructed separately for matched biological conditions, including sex, tissue, and age.

For each condition, control samples from different production centers were combined to generate a common pseudo-control reference. Center-specific control profiles were then compared with the corresponding pseudo-control to estimate production-center-specific differences.

Script:

- `3. build_pseudo_control.R`

### Step 4. Pseudo-control-based batch correction

For each production center, the difference between the center-specific control profile and the corresponding pseudo-control was used to derive the batch-correction factor.

The same center-specific correction was applied to both control and exposed samples from that center, reducing systematic production-center effects while preserving exposure-associated differences.

Script:

- `4. Run_correction_based_on_pseudo_ctrl.R`

### Step 5. Differential analysis

The batch-corrected datasets were used for downstream differential analysis.

For RNA-seq:

- Differentially expressed genes (DEGs) were identified using [`Run_DEG_analysis.R`](https://github.com/Zhang-lab/TaRGET-II-Integrative-Analysis/blob/main/Method/Pseudo_Ctrl_correction/RNA/5.%20Run_DEG_analysis.R).

For ATAC-seq:

- Differentially accessible regions (DARs) were identified using [`Run_DAR_analysis.R`](https://github.com/Zhang-lab/TaRGET-II-Integrative-Analysis/blob/main/Method/Pseudo_Ctrl_correction/ATAC/5.%20Run_DAR_analysis.R).
