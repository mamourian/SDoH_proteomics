library(dplyr)

# read data
stats_ampad <- read.csv("/Users/mamourie/Library/CloudStorage/Box-Box/mamourie/ShenLab/SDoH/SDoH_differentialAbundance_AMP-AD_v4.csv", 
                        stringsAsFactors = FALSE)

stats_rosmap <- read.csv("/Users/mamourie/Library/CloudStorage/Box-Box/mamourie/ShenLab/SDoH/SDoH_differentialAbundance_ROSMAP_v3.csv", 
                        stringsAsFactors = FALSE)

# select significant Genes (based on adjusted P)
stats_ampad_sig <- stats_ampad %>%
  dplyr::filter(adj.P.Val < .05)

stats_rosmap_sig <- stats_rosmap %>%
  dplyr::filter(adj.P.Val < .05)

# prepare lists before merging
ampad_sig <- stats_ampad_sig %>%
  dplyr::mutate(AMPAD= "AMP-AD") %>%
  dplyr::rename(AMPAD_P.adj= adj.P.Val) %>%
  dplyr::select(Gene, AMPAD, AMPAD_P.adj)

rosmap_sig <- stats_rosmap_sig %>%
  dplyr::mutate(ROSMAP= "ROSMAP") %>%
  dplyr::rename(ROSMAP_P.adj= adj.P.Val) %>%
  dplyr::select(Gene, ROSMAP, ROSMAP_P.adj)

# merge significant genes from two datasets
sig <- merge(ampad_sig, rosmap_sig, all=TRUE)
sig_f <- sig %>%
  dplyr::mutate(Study= ifelse(!is.na(ROSMAP) & is.na(AMPAD), "ROSMAP", NA)) %>%
  dplyr::mutate(Study= ifelse(is.na(ROSMAP) & !is.na(AMPAD), "AMPAD", Study)) %>%
  dplyr::mutate(Study= ifelse(!is.na(ROSMAP) & !is.na(AMPAD), "Both", Study)) %>%
  dplyr::select(-c(AMPAD, ROSMAP))

write.csv(sig_f, "/Users/mamourie/Library/CloudStorage/Box-Box/mamourie/ShenLab/SDoH/SDoH_differentialAbundance_combined_v4.csv",
          row.names = FALSE)  

###############

overlap <- sig_f %>%
  dplyr::filter(Study == "Both")



