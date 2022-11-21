# Sequence processing instructions
Author: Winona Wijaya
Last updated: 20 Jul 2022

Steps in Linux terminal
- Download fastq files from SRA (Bioproject number PRJNA848014) into a folder named "FASTQ"
- Separate v1 & v2 files
- Run Step 1 on v1 & v2 separately
- Run Step 2 on both together (at the start of the file, there is a step to merge outputs from both v1 & v2 runs)
- Output file called "17_seqtab_final_noNAs.RDS" will be used in the analysis section


Session info:

```
R version 4.1.2 (2021-11-01)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /home/winona/miniconda3/envs/r4-base/lib/libopenblasp-r0.3.18.so

locale:
 [1] LC_CTYPE=C                 LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] ShortRead_1.52.0            GenomicAlignments_1.30.0
 [3] SummarizedExperiment_1.24.0 Biobase_2.54.0
 [5] MatrixGenerics_1.6.0        matrixStats_0.61.0
 [7] Rsamtools_2.10.0            GenomicRanges_1.46.0
 [9] BiocParallel_1.28.0         forcats_0.5.1
[11] stringr_1.4.0               dplyr_1.0.7
[13] purrr_0.3.4                 readr_2.1.1
[15] tidyr_1.1.4                 tibble_3.1.6
[17] ggplot2_3.3.5               tidyverse_1.3.1
[19] Biostrings_2.62.0           GenomeInfoDb_1.30.0
[21] XVector_0.34.0              IRanges_2.28.0
[23] S4Vectors_0.32.0            BiocGenerics_0.40.0
[25] dada2_1.22.0                Rcpp_1.0.7

loaded via a namespace (and not attached):
 [1] httr_1.4.2             bit64_4.0.5            vroom_1.5.7
 [4] jsonlite_1.7.2         modelr_0.1.8           RcppParallel_5.1.5
 [7] assertthat_0.2.1       latticeExtra_0.6-29    GenomeInfoDbData_1.2.7
[10] cellranger_1.1.0       pillar_1.6.4           backports_1.4.1
[13] lattice_0.20-45        glue_1.6.0             digest_0.6.29
[16] RColorBrewer_1.1-2     rvest_1.0.2            colorspace_2.0-2
[19] Matrix_1.4-0           plyr_1.8.6             pkgconfig_2.0.3
[22] broom_0.7.11           haven_2.4.3            zlibbioc_1.40.0
[25] scales_1.1.1           jpeg_0.1-9             tzdb_0.2.0
[28] farver_2.1.0           generics_0.1.1         ellipsis_0.3.2
[31] withr_2.4.3            cli_3.1.0              magrittr_2.0.1
[34] crayon_1.4.2           readxl_1.3.1           fs_1.5.2
[37] fansi_1.0.0            xml2_1.3.3             hwriter_1.3.2
[40] tools_4.1.2            hms_1.1.1              lifecycle_1.0.1
[43] munsell_0.5.0          reprex_2.0.1           DelayedArray_0.20.0
[46] compiler_4.1.2         rlang_0.4.12           grid_4.1.2
[49] RCurl_1.98-1.5         rstudioapi_0.13        labeling_0.4.2
[52] bitops_1.0-7           gtable_0.3.0           DBI_1.1.2
[55] reshape2_1.4.4         R6_2.5.1               lubridate_1.8.0
[58] bit_4.0.4              utf8_1.2.2             stringi_1.7.6
[61] parallel_4.1.2         vctrs_0.3.8            png_0.1-7
[64] dbplyr_2.1.1           tidyselect_1.1.1
```
