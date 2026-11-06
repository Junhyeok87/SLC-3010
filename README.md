# SLC-3010: Data Analysis and Reproducibility

R-based data analysis, statistical testing, and figure generation for the SLC-3010 study.

This repository accompanies the research article entitled “SLC-3010 enhances antitumor immunity with CD122-biased, CD25-accessible primary signaling and persistent antibody activity”.

It is intended to support the reproducibility of the research presented in the article.

## Data availability

### Source Data

Numerical data supporting the findings in this manuscript will be provided in the accompanying Source Data files.

### Sequencing data

Raw RNA-seq and TCR-seq data have been submitted to the NCBI Sequence Read Archive (SRA) and are currently not publicly available. BioProject accession details and confidential reviewer links for accessing the submission metadata will be provided to the journal during peer review. Public accession links will be added here upon release.

## Reproducing the analyses

Each analysis script reads its specified Source Data workbook, writes numerical results to CSV files and can generate the associated plots. The starting inputs differ between panels; this repository does not provide a single workflow beginning with FASTQ files.

## Run your first analysis in RStudio

The following example reproduces the pSTAT5 half-response concentrations in Fig. 1i. The steps use RStudio's **Source** button and **Console** pane.

1. Extract the code ZIP and open `SLC3010.Rproj` in RStudio.
2. Obtain **`Source_Data.xlsx`** from the Source Data files provided with the manuscript. Create a folder called `data` beside `SLC3010.Rproj` and put the workbook in it. If you downloaded only the individual R script, put the workbook beside that `.R` file instead. If RStudio opens a file chooser, select this workbook. Fig. 1i uses this file alone.
3. Open [Fig_01I_pSTAT5_CD8.R](analysis/Fig_01I_pSTAT5_CD8.R). To install its packages on first use, paste the following R command in the **Console** and press Enter.

```r
install.packages(c("digest", "dplyr", "ggplot2", "minpack.lm",
                   "readr", "readxl", "scales", "tidyr"))
```

4. With the script open, click **Source** at the top of its editor pane. The script calls `analyze_data()` to calculate the numerical results and save the CSVs. The Console prints the panel-specific estimates or statistical tests and saves the same report as `analysis_summary.txt`. It then calls `plot_results()` with those results. Look in **Plots** for the figure. For scripts with several plots, use the arrows in that pane.
5. To inspect the estimates, paste these lines into the **Console** after the script finishes.

```r
analysis_result$tables$ec50        # Half-response concentrations in ng/mL
analysis_result$tables$fit_status  # Convergence and parameter bounds
analysis_result$output_dir         # Full path to this run's CSVs and plots
```

Open `Fig_01I_gMFI_5PL_EC50.csv` in that directory for the full-precision estimates. With the journal workbook, the rounded values are:

| Drug | CD25-positive CD8 T cells (ng/mL) | CD25-negative CD8 T cells (ng/mL) |
|---|---:|---:|
| Beta IL-2 | 22.794 | 16.208 |
| SLC-3010 | 0.7804 | 84.532 |

For another panel, open its R file, read the comments at the top, and click **Source**. Those comments identify the workbook, packages, result table and corresponding CSV. Each Console report includes all rows of the selected summary or comparison tables. Complete returned tables are available in the analysis result object, and the numerical CSV outputs are saved in the results directory. In any script, `names(analysis_result$tables)` lists the tables returned by that analysis.

### Calculate only the numbers in RStudio

Install only the numerical packages listed at the top of your chosen script. In the **Console**, enter the line below and press Enter. Then click **Source** on the script. This run writes the CSVs and returns the numerical results without calling the plotting stage. Analysis packages can themselves depend on graphics packages.

```r
options(SLC3010.tables_only = TRUE)
```

For Fig. 1i, inspect `analysis_result$tables$ec50` just as above. The option applies to subsequent runs in this R session. To include plots again, enter the following line in the **Console**, press Enter, and click **Source** again.

```r
options(SLC3010.tables_only = NULL)
```

### If the numbers were calculated but a plot failed

If the Console reports `Analysis tables saved in ... Plotting failed: ...`, the numerical stage completed. In the same RStudio session, use `names(analysis_result$tables)` to see the results and `analysis_result$output_dir` to find the CSVs. If the message names a missing plotting package, install that package in RStudio's Console and click **Source** again. An error during workbook reading or numerical analysis occurs earlier and can prevent a complete set of results.

## Find an analysis

**Main**, **RNA**, and **TCR** below denote `Source_Data.xlsx`, `Source_Data_RNA_inputs.xlsx`, and `Source_Data_TCR_inputs.xlsx`, respectively. The workbooks are provided separately as Source Data accompanying the manuscript. Each link takes you directly to the script that computes the listed result.

| Analysis and script | Panel | Workbook | Example numerical output |
|---|---|---|---|
| [CD8 pSTAT5 concentration responses](analysis/Fig_01I_pSTAT5_CD8.R) | 1i | Main | `Fig_01I_gMFI_5PL_EC50.csv` |
| [pSTAT5 responses in unstimulated and activated CD8 cells](analysis/Fig_01J_pSTAT5_Unstimulated_Activated.R) | 1j | Main | `Fig_01J_5PL_EC50_summary.csv` |
| [pSTAT5 concentration responses in CD4 and CD8 cells](analysis/Suppl_Fig_02B_pSTAT5_CD25posCD4_CD8.R) | S2b | Main | `Suppl_Fig_02B_gMFI_5PL_EC50.csv` |
| [Transferred and host CD8 counts](analysis/Fig_02A_Suppl_Fig_03A_OT1_Host_CD8.R) | 2a, S3a | Main | `Fig_02A_Suppl_Fig_03A_summary.csv` |
| [CD25 blockade and OVA-tetramer CD8 responses](analysis/Fig_02D-G_OT1_CD25_Blockade_OVA_Tetramer.R) | 2d–g | Main | `Fig_02D_summary.csv` |
| [MO5 growth and day-21 transfer comparisons](analysis/Fig_03A_MO5_tumor_growth.R) | 3a | Main | `Fig_03A_day21_ttest_stats.csv` |
| [Intratumoral transferred and host CD8 counts](analysis/Fig_03C-D_intratumoral_CD8_counts.R) | 3c/d | Main | `Fig_03C-D_vehicle_ttest_stats.csv` |
| [AH1-reference TCR sequences detected per tumor](analysis/Fig_03E-G_AH1_reference_TCR.R) | 3e–g | TCR | `Fig_03E_Student_contrast.csv` |
| [Single-dose serum pharmacokinetics](analysis/Fig_04A-C__Fig_S06A_single_dose_PK.R) | 4a–c, S6a | Main | `pk_metrics_from_mean_profiles.csv` |
| [Sequential serum pharmacokinetics](analysis/Fig_04E-F_sequential_serum_PK.R) | 4e/f | Main | `Fig_04E-F_summary.csv` |
| [Splenic responses after sequential administration](analysis/Fig_04G-I_sequential_splenic_responses.R) | 4g–i | Main | `Fig_04G-I_summary.csv` |
| [Splenic responses by administration interval](analysis/Fig_04K-M_interval_splenic_responses.R) | 4k–m | Main | `Fig_04K-M_recomputed_Dunnett.csv` |
| [Vehicle-normalized blood Treg responses](analysis/Fig_05A-B_blood_Treg.R) | 5a/b | Main | `Fig_05A-B_summary.csv` |
| [Genotype comparisons of splenic responses](analysis/Fig_05E-G_splenic_responses.R) | 5e–g | Main | `Fig_05E-G_test_status.csv` |
| [MP-CD8 principal component analysis](analysis/Fig_05H_PCA.R) | 5h | RNA | `Fig_05H_PCA_coordinates.csv` |
| [Count-based C1D4 differential expression](analysis/Fig_05I_volcano.R) | 5i | RNA | `Fig_05I_results.csv` |
| [Reported pathway enrichment and independent comparison](analysis/Fig_05J_fgsea.R) | 5j | RNA | `Fig_05J_results.csv` |
| [Mean MC38 tumor growth](analysis/Fig_06A_MC38_mean_growth.R) | 6a | Main | `Fig_06A_MC38_mean_growth_summary.csv` |
| [Individual MC38 tumor trajectories](analysis/Fig_06B_MC38_individual_growth.R) | 6b | Main | `Fig_06B_MC38_individual_growth_individual_data.csv` |
| [Tumor-free MC38 mice at each assessment](analysis/Fig_06C_MC38_tumor_free.R) | 6c | Main | `Fig_06C_MC38_tumor_free_percentages.csv` |
| [Mean EO771 tumor growth](analysis/Fig_06D_EO771_mean_growth.R) | 6d | Main | `Fig_06D_EO771_mean_growth_summary.csv` |
| [Individual EO771 tumor trajectories](analysis/Fig_06E_EO771_individual_growth.R) | 6e | Main | `Fig_06E_EO771_individual_growth_individual_data.csv` |
| [EO771 tumor-free percentages across experiments](analysis/Fig_06F_EO771_tumor_free.R) | 6f | Main | `Fig_06F_EO771_tumor_free_percentages.csv` |
| [Intratumoral immune-cell comparisons](analysis/Fig_06G-I_intratumoral_responses.R) | 6g–i | Main | `Fig_06G-I_comparisons.csv` |
| [Bulk-TIL expression programs](analysis/Fig_S04A-B_bulk_TIL_transcriptome.R) | S4a/b | RNA | `Fig_S04B_individual_program_scores.csv` |
| [Tissue radiotracer uptake](analysis/Fig_S05B_biodistribution.R) | S5b | Main | `Fig_S05B_mean_SEM.csv` |
| [Il2rb-positive atlas cells by tissue](analysis/Fig_S05C_CD122_expression.R) | S5c | Main | `Fig_S05C_CD122_positive_cells_by_tissue_summary.csv` |
| [Tissue CD122 and tracer-retention correlation](analysis/Fig_S05D_CD122_retention_correlation.R) | S5d | Main | `Fig_S05D_CD122_correlation_statistics.csv` |
| [Tissue-level Pearson correlations](analysis/Fig_S05E_tissue_correlation_matrix.R) | S5e | Main | `Fig_S05E_correlation_matrix_statistics.csv` |
| [Splenic dose-response comparisons](analysis/Fig_S06C-H_splenic_dose_responses.R) | S6c–h | Main | `Fig_S06C-H_recomputed_Dunnett.csv` |
| [Individual-mouse sequential serum exposure](analysis/Fig_S06I_sequential_PK_AUC.R) | S6i | Main | `Fig_S06I_mouse_AUC.csv` |
| [Representative enrichment curves](analysis/Fig_S08B_representative_GSEA.R) | S8b | RNA | `Fig_S08B_running_scores.csv` |
| [Tumor volumes by treatment cycle](analysis/Fig_S09A-B_tumor_volume_cycles.R) | S9a/b | Main | `Fig_S09A-B_comparisons.csv` |

Fig. 4e/f reports concentration means, SD, SEM and observed sample sizes for each treatment. Additional descriptive PK tables give per-mouse linear-trapezoidal AUC from 0 to 24 h and observed Cmax/Tmax, with group summaries based on complete profiles only. Missing time points are not imputed. These PK summaries use 3–6 complete profiles per analyte/treatment, whereas plotted concentrations use all 5–6 available observations at each time. No inferential test is used for this panel. The figure uses 11-point axis labels, 10-point ticks and a two-row treatment legend.

For panels with compact publication typography, the RStudio Plots pane uses a separate, larger-text display object, available as `analysis_result$display_plots`. The saved publication plots retain their original typography. For the dense S5e matrix, the display text is enlarged before layout.

Supplementary Fig. 8b is saved as `plot_01.svg` at 171.979166667 × 89.164583333 mm. The vector figure uses Arial, with nominal 8-point pathway titles, 7-point axis labels and statistics, and 6.5-point tick labels.

The [input map](docs/analysis_inputs.csv) lists the worksheets and ranges used by each script. Fig. 5i and Supplementary Fig. 8b calculate differential expression from count matrices; S8b also calculates pathway enrichment and reads pathway names from `S8_pathway_annotations`; the supplied S8 result tables are not numerical inputs to that calculation. Fig. 6c calculates tumor-free percentages from individual MC38 tumor-volume records, whereas Fig. 6f calculates the mean of three experiment-level EO771 percentages. Fig. 5j visualizes supplied pathway statistics; its independently calculated enrichment estimates are exported separately and do not replace the plotted statistics.

### AH1-reference TCR comparison (Fig. 3e)

[Fig_03E-G_AH1_reference_TCR.R](analysis/Fig_03E-G_AH1_reference_TCR.R) compares the number of distinct observed AH1-reference CDR3 matches per tumor using an unpaired, two-sided Student's t-test with equal variances and no multiple-testing adjustment. Tumors are the independent observations. The results are available in `analysis_result$tables$student_comparison`, `Fig_03E_Student_statistics.csv` and `Fig_03E_Student_contrast.csv`. Common-depth expected richness is a separate descriptive result for Fig. 3f.

### RNA sample metadata (Fig. 5h)

`Source_Data_RNA_inputs.xlsx`, `Sample_information!A3:D21`, links `Count_matrix_column` to `Sample.ID`, `Genotype` and `Time_point`. The table identifies all 18 MP-CD8 samples. The PCA workflow uses all 18 samples for normalization and displays the 12 Baseline/C1D4 samples. The matching reader is [Fig_05H_PCA.R](analysis/Fig_05H_PCA.R).

### Reproduce the EO771 tumor-free percentages (Fig. 6f)

Open [Fig_06F_EO771_tumor_free.R](analysis/Fig_06F_EO771_tumor_free.R) and click **Source** using the matching `Source_Data.xlsx` provided with the manuscript. On the **Figure 6** worksheet, `A378:G402` contains the tumor-free count and the number of mice with a recorded volume for each experiment, strain and cycle. A mouse contributes to the tumor-free count when its recorded tumor volume is zero on that assessment day. In the RStudio Console, use `analysis_result$tables$experiment_percentages` for all 24 experiment-level records and `analysis_result$tables$reported_percentages` for the cycle means.

For Triple hIL-2 mice at cycle 3, the counts are 3/8, 2/8 and 0/7. The arithmetic mean of the three percentages is 20.8% when rounded to one decimal place. At cycle 4, 4/8, 2/8 and 1/7 yield 29.8% by the same calculation. Each experiment has equal weight.

### Trace a sample across the three splenic readouts (Supplementary Fig. 6f–h)

Open [Fig_S06C-H_splenic_dose_responses.R](analysis/Fig_S06C-H_splenic_dose_responses.R) and click **Source**. The **Supplementary Figure 6** worksheet in `Source_Data.xlsx` records the experiment year and original FCS filename for each observation in `A97:I295`. A filename can recur in both experiments, so use the year and filename together when locating a sample. The script uses those identifiers to match the two cell counts and the recorded MP CD8/Treg ratio for each of the 66 samples.

Enter `analysis_result$tables$paired_samples` in the RStudio Console, or open `Fig_S06F-H_paired_samples.csv` in `analysis_result$output_dir`. Each row shows the counts, the recorded ratio, the ratio calculated from the counts and their difference. Use `analysis_result$tables$comparisons` for the panel-specific Dunnett comparisons.

## Packages, results and reuse

- [Analysis and plotting dependencies](environment/analysis_packages.csv) are listed separately per script. Install CRAN packages with `install.packages()`. For Bioconductor packages such as DESeq2 and fgsea, use `BiocManager::install()` with the Bioconductor release appropriate for your R version. The [reference package versions](environment/package_versions.csv) and [R session information](environment/R_session_info.txt) record the reference environment.
- With the complete repository, each run writes to a new directory under `<repository-name>_results` beside the repository. For a standalone file, the default is `results` beside the script. If necessary the script tries the user's application data directory. Set `SLC_OUTPUT_DIR` to choose a writable output root; set `SLC_SOURCE_DATA_DIR`  to use a different input. [Workbook checksums](docs/source_data_checksums.csv) record the reference size and SHA-256 for each of the three Source Data workbooks. Match inputs by the required worksheets and the workbook checksum; identical filenames do not establish identical contents.
- To use the functions directly, set `options(SLC3010.autorun = FALSE)` before sourcing the file. Call `analyze_data(workbook, output_dir)` with an existing, writable output directory. Set `options(SLC3010.autorun = NULL)` to restore normal Source behavior.
- The optional `RUN_ALL.R` requires workbooks containing all worksheets listed in the input map, including `MP_CD8_GSEA_results` for Fig. 5j. With those inputs in `data/`, open `RUN_ALL.R` and click **Source**. The same `SLC3010.tables_only` R option described above selects numerical results alone. The independent Fig. 5j enrichment calculation explicitly uses one worker through `BiocParallel::SerialParam()`, including within `RUN_ALL.R` child processes; the figure uses the supplied pathway statistics.

## Statistical reporting

The `Statistical tests` and `Statistics details` sheets in the Source Data workbooks report the inferential analyses. Main contains 119 records and TCR contains one: 110 comparisons, nine omnibus ANOVA results and one reported test without a supplied numerical result. Missing or inapplicable fields remain blank and are explained in `Test_status` and `Missing_value_notes`. The [analysis overview](docs/analysis_overview.csv) links the analyses, experimental units, methods, inputs and result locations, including descriptive panels and supplied RNA result tables.

The two statistics sheets share `Figure`, `Panel` and the unique `Comparison_ID`. `Statistical tests` contains 32 principal fields; `Statistics details` contains the other 28 fields together with those three identifiers. For t-tests, `Estimate_direction` defines the direction of the mean difference, t statistic and interval. Fig. 6g-i displays Holm-adjusted P values with unadjusted mean-difference intervals. Pearson intervals are unadjusted Fisher-z intervals; unavailable Spearman degrees of freedom and intervals remain blank.

For Dunnett comparisons, `Adjusted_P_value` and `P_used_for_display` contain the reported Source Data probabilities. `Comparative_adjusted_P_value` contains the separately recomputed probability. Reported and recomputed probabilities remain distinct; unavailable raw P values and simultaneous intervals remain blank. The recorded discrepancy threshold of 1e-8 is a numerical diagnostic rather than a criterion for scientific equivalence. The recorded 1e-15 precision limit is retained where applicable.

Individual analysis scripts write their numerical results to the CSV locations listed in the analysis overview. Fig. 5i calculates gene-level differential expression from raw counts. Fig. 5j plots supplied pathway statistics and writes its independent enrichment estimates separately. Supplementary Fig. 8b calculates enrichment from raw counts and fixed references. Descriptive panels have no added hypothesis tests, and numerical results are not inferred from significance symbols.

## Questions

For questions about the code or to report a problem, please open a [GitHub issue](https://github.com/Junhyeok87/SLC-3010/issues). Include the relevant script, software versions, and any error message when reporting an execution problem.
