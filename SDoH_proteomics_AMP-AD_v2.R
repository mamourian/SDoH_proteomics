library(dplyr)
library(limma)
library(ggplot2)
library(EnhancedVolcano)

###############################################################################
# Read Data
# dir <- "/Users/mamourie/Library/CloudStorage/Box-Box/"
# dir_ampAD <- "/Users/mamourie/Library/CloudStorage/Box-Box/Elizabeth-Project-Data/AMP-AD/DiversityCohort/"
dir <- "/Users/elizabethmamourian/Library/CloudStorage/Box-Box/"
dir_ampAD <- "/Users/elizabethmamourian/Library/CloudStorage/Box-Box/Elizabeth-Project-Data/AMP-AD/DiversityCohort/"

## SDoH genes of interest
sdoh.orig <- readxl::read_xlsx(paste0(dir,"mamourie/ShenLab/SDoH/potential_edges_gene_sdoh(TZedited).xlsx"), sheet= "unique genes")

## QC'ed proetomics data
prot.orig <- read.csv(paste0(dir_ampAD, "syn59611693/QC/normAbundances_post-FP_frontal.csv"), stringsAsFactors = FALSE)

## clinical data w/ Subject ID
clin.orig <- read.csv(paste0(dir_ampAD, "syn51732482/Data/Metadata/AMP-AD_DiverseCohorts_individual_metadata.csv"), stringsAsFactors = FALSE)

## mapping between sample ID & Subject ID
id_map <- read.csv(paste0(dir_ampAD, "syn51732482/Data/Metadata/AMP-AD_DiverseCohorts_biospecimen_metadata.csv"), stringsAsFactors = FALSE)



################################################################################
# Prepare data format, primary keys
## create subject metadata file
clin_id <- merge(clin.orig, id_map, by="individualID") %>%
  dplyr::filter(specimenID %in% names(prot.orig)) 
table(clin_id$ADoutcome, exclude=NULL)
clin <- clin_id %>%
  dplyr::filter(ADoutcome %in% c("AD", "Control"))
clin[clin$ageDeath=="90+", "ageDeath"] <- 90

## determine demographics of each diagnosis group
clin_ad <- clin %>%
  dplyr::filter(ADoutcome== "AD")
clin_ctrl <- clin %>%
  dplyr::filter(ADoutcome== "Control")
# sex
table(clin_ad$sex) 
table(clin_ctrl$sex) 
# age
mean(as.numeric(clin_ad$ageDeath), na.rm = TRUE)
sd(as.numeric(clin_ad$ageDeath), na.rm = TRUE)
mean(as.numeric(clin_ctrl$ageDeath))
sd(as.numeric(clin_ctrl$ageDeath))

## match genes list to SDoH genes of interest
prot_ids <- as.data.frame(stringr::str_extract(prot.orig$X, "[^|]+")) %>%
  dplyr::rename(Gene= "stringr::str_extract(prot.orig$X, \"[^|]+\")")
prot_ids$X <- prot.orig$X
sdoh <- merge(sdoh.orig, prot_ids)  # 298 / 758 match

## select genes of interest and subjects with metadata
prot_sel <- merge(prot.orig, sdoh)
prot_sel_sbj <- prot_sel[, names(prot_sel) %in% clin$specimenID]
row.names(prot_sel_sbj) <- prot_sel$Gene

## create order of subjects/samples to apply to metadata
prot_order <- as.data.frame(names(prot_sel_sbj))
names(prot_order) <- "specimenID"
prot_order$order <- 1:length(prot_order$specimenID)

## apply order of samples in "prot_sel_sbj" to subject metadata
clin_order <- merge(clin, prot_order) %>%
  dplyr::arrange(order)

## check order of metadata against order of proteomics data
clin_order$specimenID[1:10]
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
des_mat <- model.matrix(~ADoutcome+0, data = clin_order)

# Clean design matrix column names
colnames(des_mat) <- stringr::str_remove(colnames(des_mat), "ADoutcome")

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
  tibble::rownames_to_column("Gene")
stats_df_all <- merge(sdoh.orig, stats_df.orig, all=TRUE) # include missing genes
stats_df <- merge(sdoh, stats_df.orig, all=TRUE) %>% # includes only results
  dplyr::mutate(FoldChange = 2^logFC - 1)

stats_df_sig <- stats_df %>%
  dplyr::filter(adj.P.Val < .05)

## Write output
write.csv(stats_df, paste0(dir,"mamourie/ShenLab/SDoH/SDoH_differentialAbundance_AMP-AD_v3.csv"), row.names = FALSE)


###############################################################################

volcano_plot <- EnhancedVolcano::EnhancedVolcano(stats_df,
                                                 lab = stats_df$Gene,
                                                 x = "logFC",
                                                 y = "P.Value", 
                                                 subtitle ="AD vs. Control in AMP-AD",
                                                 xlab = bquote(~Log[2]~ 'fold change'),
                                                 pCutoff = 0.05,
                                                 FCcutoff = .3,
                                                 pointSize = 4.0,
                                                 labSize = 2.0,
                                                 colAlpha = 0.5,
                                                 xlim= c(-1,1),
                                                  ylim= c(0,17),
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
  file.path(paste0(dir,"mamourie/ShenLab/SDoH/SDoH_VolcanoPlot_AMP-AD_v3.png")),
  plot = volcano_plot, # The plot object that we want saved to file
  bg = 'white'
)


###############################################################################

volcano_plot_FC2 <- EnhancedVolcano::EnhancedVolcano(stats_df,
                                                 lab = stats_df$Gene, 
                                                 x = "logFC", 
                                                 y = "P.Value", 
                                                 subtitle ="AD vs. Control in AMP-AD",
                                                 xlab = bquote(~Log[2]~ 'fold change'),
                                                 pCutoff = 0.05,
                                                 FCcutoff = 2,
                                                 pointSize = 4.0,
                                                 labSize = 2.0,
                                                 colAlpha = 0.5,
                                                 # xlim= c(-.3,.3),
                                                 ylim= c(0,17),
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
  file.path(paste0(dir,"mamourie/ShenLab/SDoH/SDoH_VolcanoPlot_AMP-AD_v3_FC2.png")),
  plot = volcano_plot_FC2, # The plot object that we want saved to file
  bg = 'white'
)

###############################################################################