# run_dpo_analysis.R
library(limma)
library(reshape2)
library(plyr)


#' Remove Low-Quality Controls
#'
#' Identifies and removes phosphosites from the dataset that demonstrate poor quantification 
#' in control samples.
#'
#' @param data: DataFrame with phosphosite data, columns as samples.
#' @param mdat: Metadata DataFrame.
#' @param condition_col: Name of the column in mdat to use for determining control samples.
#' @param control_name: The value within condition_col that denotes control samples.
#' @param min_ctrs: Minimum number of control measurements required to keep a phosphosite.
#' 
#' @return Returns a list with two elements:
#'         - `filtered_data`: A DataFrame containing the subset of the original `data` where phosphosites with adequate quantification in control samples are retained. Rows represent phosphosites, and columns represent samples.
#'         - `removed_ctrs`: A DataFrame of the phosphosites removed from `data` due to insufficient quantification in control samples. This set includes phosphosites that did not meet the `min_ctrs` threshold. Similar to `filtered_data`, rows represent phosphosites, and columns represent samples. This DataFrame allows for further inspection of phosphosites excluded from the analysis and may be used for subsequent analysis or quality control steps.
#'
remove_low_quality_controls <- function(data, mdat, condition_col='condition', control_name='control', min_ctrs=4) {
  # Check if 'sampleID' column exists in mdat
  if(!'sampleID' %in% names(mdat)) {
    stop("Error: Metadata must contain a 'sampleID' column.")
  }
  
  # Extract control samples based on metadata
  ctr_samples <- mdat$sampleID[mdat[[condition_col]] == control_name]
  
  # Ensure control samples are present
  if(length(ctr_samples) == 0) {
    stop(paste("No control samples found with condition_col:", condition_col, "and control_name:", control_name))
  }
  
  ctr_data <- data[, ctr_samples, drop=FALSE]
  
  # Determine phosphosites to keep based on the minimum number of control measurements
  keep <- rowSums(!is.na(ctr_data)) >= min_ctrs
  removed_ctrs <- data[!keep, , drop=FALSE]  # Data to add back later
  
  return(list(filtered_data=data[keep, , drop=FALSE], removed_ctrs=removed_ctrs[, ctr_samples, drop=FALSE]))
}


#' Calculate Differential Phosphosite Occupancy with Limma
#'
#' Performs differential phosphosite occupancy analysis using the limma package. It constructs 
#' a linear model incorporating sample conditions and potential batch effects to identify 
#' statistically significant changes in phosphosite occupancy across conditions.
#' 
#' @param data: Data frame with rows as phosphosites and columns as samples.
#' @param mdat: Metadata data frame with sample information including conditions and batch data.
#' @param condition_col: Column name in 'mdat' that specifies the condition (e.g., perturbagen) of each sample.
#' @param batch_col: Column name in 'mdat' for batch information to account for in the model.
#' @param control_name: Name of the control group within the condition column.
#' 
#' @return A DataFrame containing the results of the differential phosphosite occupancy analysis. The DataFrame includes:
#'         - `phosphosite`: Identifiers for each analyzed phosphosite.
#'         - `condition`: The condition name for each comparison made against the control group.
#'         - `fold_change`: Computed fold change values for each phosphosite in each condition relative to the control, indicating the magnitude of occupancy change.
#'         - `pval_eb`: eBayes-adjusted p-values associated with each fold change, quantifying the statistical significance of observed changes.
#'         - `error`: An optional column indicating any errors or warnings associated with the calculation for each phosphosite, such as 'condition_missing' if data for the condition was insufficient. This column aids in identifying potential issues or limitations in the analysis results.
#'
calculate_dpo_limma <- function(data, mdat, condition_col='condition', batch_col='batch', control_name='control') {
  # Ensure the condition and batch columns exist in metadata
  if (!(condition_col %in% names(mdat))) {
    stop(paste("The condition column", condition_col, "does not exist in metadata."))
  }
  if (!(batch_col %in% names(mdat))) {
    stop(paste("The batch column", batch_col, "does not exist in metadata."))
  }
  
  # Prepare design matrix with conditions and batch as covariates
  conditions <- factor(mdat[[condition_col]])
  batch <- factor(mdat[[batch_col]])
  design <- model.matrix(~0 + conditions + batch)
  colnames(design) <- gsub('conditions', '', colnames(design))
  colnames(design) <- gsub('batch', '', colnames(design))
  
  # Prepare contrasts matrix
  contrasts <- lapply(levels(conditions)[levels(conditions) != 'control'], paste, control_name, sep='-')
  contrasts$levels <- design
  cont.matrix <- do.call(makeContrasts, contrasts)
  
  # Fit linear model and compute eBayes statistics
  fit <- lmFit(data, design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, trend=TRUE)
  
  # Format the results
  dpoa_res <- melt(fit2$coefficients)
  colnames(dpoa_res) <- c('phosphosite', condition_col, 'fold_change')
  dpoa_res$phosphosite <- as.character(dpoa_res$phosphosite)
  dpoa_res$perturbagen <- gsub(paste0('-', control_name), '', dpoa_res$perturbagen)
  pvals <- melt(fit2$p.value)
  dpoa_res$pval_eb <- pvals$value
  dpoa_res$error[is.na(dpoa_res$fold_change)] <- 'condition_missing'
  
  return(dpoa_res)
}


#' Add Back Low-Quality Controls
#' 
#' Add back phosphosites removed due to poor quantification in control samples. This step 
#' reintroduces phosphosites that were excluded from the initial analysis because of insufficient 
#' data quality in control samples, marking them accordingly for further review or analysis.
#'
#' @param dpoa_res: DataFrame containing the results of DPOA analysis.
#' @param removed_ctrs: DataFrame containing phosphosites removed due to poor control quantification.
#' @param mdat: Metadata DataFrame.
#' @param condition_col: The column in mdat that indicates the condition of the samples.
#' @param control_name: The name indicating control samples within the condition_col.
#'   
#' @return Returns an updated DataFrame combining the original DPOA results with the phosphosites added back. The output includes:
#'         - `phosphosite`: Identifiers for each phosphosite, including those added back.
#'         - condition column (specified by `condition_col`): Indicates the experimental condition associated with each phosphosite.
#'         - `error`: An error annotation indicating phosphosites that had insufficient quantification in control samples.
#' 
add_back_controls <- function(dpoa_res, removed_ctrs, mdat, condition_col='condition', control_name='control') {
  # Ensure the sample ID and condition columns exist in metadata
  if(!'sampleID' %in% names(mdat)) {
    stop("Error: Metadata must contain a 'sampleID' column.")
  }
  if (!(condition_col %in% names(mdat))) {
    stop(paste("The condition column", condition_col, "does not exist in metadata."))
  }
  
  # Identify unique conditions
  conditions <- unique(mdat[[condition_col]])
  to_add_back <- data.frame(phosphosite=character(), condition_col=character(), error=character(), stringsAsFactors=FALSE)
  
  # Iterate through conditions
  for (cnd in conditions) {
    if(cnd == control_name) next # Skip the control group itself
    # Extract control samples for the current condition
    ctr_samples <- mdat$sampleID[mdat[[condition_col]] == control_name]
    to_add_back <- rbind(to_add_back, data.frame(phosphosite=rownames(removed_ctrs), condition_col=cnd, error='control_missing', stringsAsFactors=FALSE))
  }
  colnames(to_add_back)[colnames(to_add_back) == 'condition_col'] <- condition_col
  
  # Merge back into dpoa results
  final_res <- rbind.fill(dpoa_res, to_add_back)
  
  return(final_res)
}


#' Add Metrics to DPOA Results
#'
#' This function enriches the Differential Phosphosite Occupancy Analysis (DPOA)
#' results with additional metrics including the total runs per phosphosite, counts of
#' control samples, counts and mean signals for each treatment condition.
#'
#' @param data Data frame with rows representing phosphosites and columns representing samples.
#'             It should contain log2-transformed area under the peak (AUP) values or similar 
#'             metrics indicating the presence or absence of phosphosites.
#'
#' @param mdat Metadata data frame corresponding to samples in 'data'. It should include
#'             at least sample IDs and condition information.
#'
#' @param dpoa Data frame containing the initial DPOA results, including phosphosite
#'             identifiers and initial metrics such as fold changes and p-values.
#'
#' @param condition_col Name of the column in 'mdat' that specifies the experimental condition
#'                      for each sample, distinguishing between control and various treatment groups.
#'
#' @param control_name Value within 'condition_col' that identifies control samples, used to
#'                     differentiate controls from treatments in the analysis.
#'
#' @return Updated 'dpoa' data frame, now including the additional columns: 'n_runs' (total
#'         successful runs per phosphosite), 'n_ctr' (number of control samples per phosphosite),
#'         'n_cnd' and 'meansig_cnd' (count and mean signal of each treatment condition per phosphosite,
#'         respectively).
#'
add_dpoa_metrics <- function(data, mdat, dpoa, condition_col='condition', control_name='control') {
  if(!'sampleID' %in% names(mdat)) {
    stop("Error: Metadata must contain a 'sampleID' column.")
  }
  
  # Total runs per phosphosite (non-NA counts)
  dpoa$n_runs <- rowSums(!is.na(data))
  
  # Control samples: counts and mean signal
  control_samples <- which(mdat[[condition_col]] == control_name)
  dpoa$n_ctr <- rowSums(!is.na(data[, control_samples]))[dpoa$phosphosite]
  dpoa$meansig_ctr <- rowMeans(data[, control_samples], na.rm=TRUE)[dpoa$phosphosite]
  
  # Unique conditions excluding control: counts and mean signal
  unique_conditions <- setdiff(unique(mdat[[condition_col]]), control_name)
  for (cond in unique_conditions) {
    cond_samples <- which(mdat[[condition_col]] == cond)
    n_cnd_values <- rowSums(!is.na(data[, cond_samples]))
    meansig_cnd_values <- rowMeans(data[, cond_samples], na.rm=TRUE)
    
    # Update dpoa with values for the current condition
    condition_rows <- dpoa[[condition_col]] == cond
    dpoa$n_cnd[condition_rows] <- n_cnd_values[dpoa$phosphosite[condition_rows]]
    dpoa$meansig_cnd[condition_rows] <- meansig_cnd_values[dpoa$phosphosite[condition_rows]]
  }
  
  return(dpoa)
}


#' Run Differential Phosphosite Occupancy Analysis (DPOA)
#'
#' This function serves as a comprehensive workflow for conducting DPOA on LC/MS-MS
#' phosphoproteomics data. It involves several key steps: quality control (QC) to
#' exclude low-quality phosphosites, differential analysis using the limma package,
#' and the integration of phosphosites removed during QC back into the analysis results.
#'
#' @param log2aup_path String path to the log2-transformed Area Under the Peak (AUP)
#'        data file. The file should be in tab-separated values (TSV) format with
#'        phosphosites as rows and samples as columns. Row names should be phosphosite identifiers.
#' @param metadata_path String path to the metadata file containing sample information.
#'        The file should be in TSV format and must include at least 'sampleID' and
#'        'condition' (or another condition column) columns.
#'
#' @return Returns a data frame with the results of the DPOA analysis, including
#'         phosphosite identifiers, conditions, fold changes, p-values, and annotations 
#'         providing insights into potential errors or limitations in fold change calculations 
#'         for each phosphosite ('condition_missing' or 'control_missing' or None).
#'
run_dpo_analysis <- function(log2aup_path, metadata_path) {
  message("Loading data and metadata...")
  data <- read.csv(log2aup_path, sep='\t', row.names=1)
  mdat <- read.csv(metadata_path, sep='\t')
  
  message("Aligning data and metadata...")
  mdat <- mdat[mdat$sampleID %in% colnames(data), ]
  
  message("Removing phosphosites with insufficient data in control samples...")
  qc_res <- remove_low_quality_controls(data, mdat, condition_col='perturbagen')
  
  message("Performing DPOA using limma...")
  dpoa_res <- calculate_dpo_limma(qc_res$filtered_data, mdat, condition_col='perturbagen')
  
  message("Adding back phosphosites with insufficient data in control samples...")
  dpoa_res <- add_back_controls(dpoa_res, qc_res$removed_ctrs, mdat, condition_col='perturbagen')
  
  message("Annotating DPOA results with additional metrics...")
  final_res <- add_dpoa_metrics(data, mdat, dpoa_res, condition_col='perturbagen', control_name='control')
  
  message("DPOA analysis completed.")
  return(final_res)
}


if (!interactive() && exists('snakemake')) {
  # Access Snakemake variables
  log2aup_path <- snakemake@input[[1]]
  metadata_path <- snakemake@input[[2]]
  output_path <- snakemake@output[[1]]
  
  # Run the DPOA analysis with progress messages
  message("Starting DPOA analysis...")
  phosphodata_dpoa <- run_dpo_analysis(log2aup_path, metadata_path)
  
  # Write the results to the output path
  message(paste("Saving results to", output_path, "..."))
  write.table(phosphodata_dpoa, output_path, sep='\t', quote=FALSE, row.names=FALSE)
  message("Results saved successfully.")
}

# log2aup_path <- 'workspace/normalised_data.tsv'
# metadata_path <- 'resources/raw_data/metadata.tsv'
# output_path <- 'workspace/phosphodata_dpoa.tsv'
