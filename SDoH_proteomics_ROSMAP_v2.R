library(dplyr)
library(limma)
library(ggplot2)
library(EnhancedVolcano)

###############################################################################
# Read Data
# dir <- "/Users/mamourie/Library/CloudStorage/Box-Box/"
# dir_ROSMAP <- "/Users/mamourie/Library/CloudStorage/Box-Box/Elizabeth-Project-Data/ROSMAP/"
dir <- "/Users/elizabethmamourian/Library/CloudStorage/Box-Box/"
dir_ROSMAP <- "/Users/elizabethmamourian/Library/CloudStorage/Box-Box/Elizabeth-Project-Data/ROSMAP/"

## SDoH genes of interest
sdoh.orig <- readxl::read_xlsx(paste0(dir,"mamourie/ShenLab/SDoH/potential_edges_gene_sdoh(TZedited).xlsx"), sheet= "unique genes")

## QC'ed proetomics data
prot.orig <- read.csv(paste0(dir_ROSMAP, "proteomics/syn21261728/output/C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv"), stringsAsFactors = FALSE)

## clinical data w/ Subject ID
clin.orig <- read.csv(paste0(dir_ROSMAP, "proteomics/syn21261728/input/rosmap_50batch_specimen_metadata_for_batch_correction.csv"), stringsAsFactors = FALSE)

## metadata
meta.orig <- read.csv(paste0(dir_ROSMAP, "ROSMAP_clinical.csv"), stringsAsFactors = FALSE)


################################################################################
# Prepare data format, primary keys
## create subject metadata file
table(clin.orig$EmoryStrictDx.2019, exclude=NULL)
clin <- merge(clin.orig, meta.orig) %>%
  dplyr::filter(EmoryStrictDx.2019 %in% c("AD", "Control")) %>%
  dplyr::rename(Diagnosis= EmoryStrictDx.2019)
clin[clin$age_death=="90+", "age_death"] <- 90

## determine demographics of each diagnosis group
clin_ad <- clin %>%
  dplyr::filter(Diagnosis== "AD")
clin_ctrl <- clin %>%
  dplyr::filter(Diagnosis== "Control")
# sex
table(clin_ad$msex) # 1=Male
table(clin_ctrl$msex) 
# age
mean(as.numeric(clin_ad$age_death))
sd(as.numeric(clin_ad$age_death))
mean(as.numeric(clin_ctrl$age_death))
sd(as.numeric(clin_ctrl$age_death))

## match genes list to SDoH genes of interest
prot_ids <- as.data.frame(stringr::str_extract(prot.orig$X, "[^|]+")) %>%
  dplyr::rename(Gene= "stringr::str_extract(prot.orig$X, \"[^|]+\")")
prot_ids$X <- prot.orig$X 
prot_ids <- prot_ids %>%
  dplyr::filter(Gene != 0)
sdoh <- merge(sdoh.orig, prot_ids)  # 225 / 758 match

## select genes of interest and subjects with metadata
prot_sel <- merge(prot.orig, sdoh) 
prot_sel_sbj <- prot_sel[, names(prot_sel) %in% clin$SampleID]
row.names(prot_sel_sbj) <- prot_sel$X  # multiple UniprotIDs per gene

## create order of subjects/samples to apply to metadata
prot_order <- as.data.frame(names(prot_sel_sbj))
names(prot_order) <- "SampleID"
prot_order$order <- 1:length(prot_order$SampleID)

## apply order of samples in "prot_sel_sbj" to subject metadata
clin_order <- merge(clin, prot_order) %>%
  dplyr::arrange(order)

## check order of metadata against order of proteomics data
clin_order$SampleID[1:10]
names(prot_sel_sbj)[1:10]

## Add pseudo count to avoid log(0)
prot_sel_sbj_orig <- prot_sel_sbj
prot_sel_sbj <- log2(prot_sel_sbj_orig + 1)  


###############################################################################
# We will use the `lmFit()` function from the `limma` package to test each gene 
# for differential abundance between the treated and untreated groups using a 
# linear model.
# After fitting our data to the linear model, in this example we apply empirical 
# Bayes smoothing with the `eBayes()` function.

## Create the design matrix from the metadata
des_mat <- model.matrix(~Diagnosis+0, data = clin_order)

# Clean design matrix column names
colnames(des_mat) <- stringr::str_remove(colnames(des_mat), "Diagnosis")

# Look at the design matrix
head(des_mat)

# Apply linear model to data
fit <- limma::lmFit(prot_sel_sbj, design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# It is also necessary to perform some multiple test corrections, which we will 
# do with the Benjamini-Hochberg method while making a table of results with `topTable()`.

# deisgn a contrast matrix
contrast_matrix <- makeContrasts("ADvsControl" = AD - Control, levels = des_mat)

# Fit the model according to the contrasts matrix
contrasts_fit <- contrasts.fit(fit, contrast_matrix)

# Re-smooth the Bayes
contrasts_fit <- eBayes(contrasts_fit)

# Apply multiple testing correction and obtain stats
stats_df.orig <- topTable(contrasts_fit, number = nrow(prot_sel_sbj)) %>%
  tibble::rownames_to_column("X")
stats_df_all <- merge(sdoh.orig, stats_df.orig, all=TRUE) # include missing genes
stats_df <- merge(sdoh, stats_df.orig, all=TRUE) # includes only results

stats_df_sig <- stats_df %>%
  dplyr::filter(adj.P.Val < .05)

## Write output
write.csv(stats_df, paste0(dir,"mamourie/ShenLab/SDoH/SDoH_differentialAbundance_ROSMAP_v3.csv"), row.names = FALSE)


###############################################################################

volcano_plot <- EnhancedVolcano::EnhancedVolcano(stats_df,
                                                 lab = stats_df$Gene, 
                                                 x = "logFC", 
                                                 y = "P.Value", 
                                                 subtitle ="AD vs. Control in ROSMAP",
                                                 xlab = bquote(~Log[2]~ 'fold change'),
                                                 pCutoff = 0.05,
                                                 FCcutoff = .3,
                                                 pointSize = 4.0,
                                                 labSize = 2.0,
                                                 colAlpha = 0.5,
                                                  xlim= c(-1,1),
                                                  ylim= c(0,10),
                                                 legendPosition = 'top',
                                                 legendLabSize = 12,
                                                 legendIconSize = 4.0,
                                                 drawConnectors = TRUE,
                                                 widthConnectors = 0.35,
                                                 arrowheads = FALSE
)
# Print out our plot
volcano_plot

# Save plot using `ggsave()` function
ggsave(
  file.path(paste0(dir,"mamourie/ShenLab/SDoH/SDoH_VolcanoPlot_ROSMAP_v3.png")),
  plot = volcano_plot, # The plot object that we want saved to file
  bg = 'white'
)


###############################################################################

volcano_plot_FC2 <- EnhancedVolcano::EnhancedVolcano(stats_df,
                                                 lab = stats_df$Gene, 
                                                 x = "logFC", 
                                                 y = "P.Value", 
                                                 subtitle ="AD vs. Control in ROSMAP",
                                                 xlab = bquote(~Log[2]~ 'fold change'),
                                                 pCutoff = 0.05,
                                                 FCcutoff = 2,
                                                 pointSize = 4.0,
                                                 labSize = 2.0,
                                                 colAlpha = 0.5,
                                                 xlim= c(-3,3),
                                                 ylim= c(0,10),
                                                 legendPosition = 'top',
                                                 legendLabSize = 12,
                                                 legendIconSize = 4.0,
                                                 drawConnectors = TRUE,
                                                 widthConnectors = 0.35,
                                                 arrowheads = FALSE
)
# Print out our plot
volcano_plot_FC2

# Save plot using `ggsave()` function
ggsave(
  file.path(paste0(dir,"mamourie/ShenLab/SDoH/SDoH_VolcanoPlot_ROSMAP_v3_FC2.png")),
  plot = volcano_plot_FC2, # The plot object that we want saved to file
  bg = 'white'
)

###############################################################################
