
# DeepCellProfiler features without preprocessing and whitening step


library(htmlTable)
library(knitr)
library(magrittr)
library(stringr)
library(tidyverse)
library(cytominer)
library(readbulk)
library(dplyr)
library(tidyverse)

# Importing Normalized csv file to extract metadata
metadata.df <- data.frame(readr::read_csv("../input/metadata_CDRP.csv"),  stringsAsFactors =F)
path <- "/Users/habbasi/Desktop/stuff/Mohammad_paper/CDRP/DeepProfiler/"
profiles <- read_bulk(directory = path, extension = "_normalized.csv", stringsAsFactors=FALSE)

profiles$File <- NULL  # removing File column
profile.type <- "mean"
quant <- 0.99

mean.na <- function(x) {mean(x, na.rm = T)}

metadata_profiles <- c(
  stringr::str_subset(colnames(profiles), "^Meta")
)
pmeta <- profiles %>% 
  select(one_of(metadata_profiles)) %>% 
  dplyr::collect()


# Importing deepCellProfiler features

f.name <- "InceptionResnetV2-mean_profiles.csv"
data <- readr::read_csv(f.name)

data$Plate_Well <- NULL # Extra column generated
data$Metadata_Plate_1 <- NULL
# Renaming deepCellProfiler features
colnames <- colnames(data)

metadata_data  <- c(
  stringr::str_subset(colnames, "^Meta")
) 


common_features  <- setdiff(colnames, metadata_data)

# Adding Feature_ prefix to columns to be able to read properly

common_features  <- paste0("Feature_",common_features)

variables <- c(metadata_data, common_features)
#variables <- c(common_features, metadata)  #for TA-ORF dataset

# Defining new colnames
data <- setNames(data, variables)

# Combining metadata infomation to deepCellProfiler features

df <- merge(data, pmeta, by.x=c("Metadata_Plate", "Metadata_Well"), by.y=c("Metadata_Plate", "Metadata_Well")) 

# defining features 
colnames <- colnames(df)

metadata_df  <- c(
  stringr::str_subset(colnames, "^Meta")
)
variables <- c(
  stringr::str_subset(colnames, "^Feature_")
) 

profiles.nrm <- df
if (!"Metadata_mmoles_per_liter" %in% colnames(profiles.nrm)) {
  profiles.nrm <- profiles.nrm %>%
    mutate(Metadata_mmoles_per_liter = 10)
}

prf <- profiles.nrm %>% 
  group_by(Metadata_broad_sample, Metadata_mmoles_per_liter, Metadata_Plate_Map_Name) %>%
  summarise_at(.vars = variables, .funs = "mean.na")

prf %<>%
  arrange(abs(Metadata_mmoles_per_liter - 10)) %>%
  group_by(Metadata_broad_sample) %>%
  slice(1) %>%
  ungroup()

if (is.null(metadata.df)) {
  prf %<>% 
    left_join(profiles.nrm %>% select(Metadata_broad_sample, Metadata_moa) %>% unique, by = "Metadata_broad_sample")
} else {
  prf %<>% 
    left_join(metadata.df %>% select(Metadata_broad_sample, Metadata_moa) %>% unique, by = "Metadata_broad_sample")
}



metadata_prf  <- c(
  stringr::str_subset(colnames(prf), "^Meta")
)

feats  <- c(
  stringr::str_subset(colnames(prf), "^Feature")
)

profiles.meta <- prf %>% select(Metadata_broad_sample, Metadata_moa) %>% dplyr::collect()

evaluate.moa <- function(cr, profiles.meta, quant = 0.95, type.eval = "global", k = 1, skip.comp = F) {
  cr.melt <- cr %>% reshape2::melt()
  cr.melt <- merge(cr.melt, profiles.meta,
                   by.x = "Var1",
                   by.y = "Metadata_broad_sample")
  cr.melt <- merge(cr.melt, profiles.meta,
                   by.x = "Var2",
                   by.y = "Metadata_broad_sample")
  
  
  match.moas <- function(moa1, moa2) {
    if (is.na(moa1) | is.na(moa2) | moa1 == "" | moa2 == "") {
      return(FALSE)
    }
    x <- str_split(moa1, "\\|")[[1]]
    y <- str_split(moa2, "\\|")[[1]]
    return(any(x %in% y) | any(y %in% x))
  }
  match.moas <- Vectorize(match.moas)
  saveRDS(cr.melt, paste0("dp_test_cr_",  profile.type, ifelse(profile.type == "mix", paste0(mix1, mix2), ""), ".rds"))
  
  if (skip.comp) {
    return(NULL)
  }
  
  if (type.eval == "global") {
    cr.melt <- cr.melt %>% filter(Var1 < Var2) %>% View()
    
    thr <- quantile(cr.melt$value, quant, na.rm = T)
    
    v11 <- cr.melt %>% filter(value > thr & match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    v12 <- cr.melt %>% filter(value > thr & !match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    v21 <- cr.melt %>% filter(value <= thr & match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    v22 <- cr.melt %>% filter(value <= thr & !match.moas(Metadata_moa.x, Metadata_moa.y)) %>% NROW()
    
    V <- rbind(c(v11, v12), c(v21, v22))
    
    return(fisher.test(V, alternative = "greater")$estimate)
  } else if (type.eval == "classification") {
    res <- cr.melt %>%
      filter(Var1 != Var2) %>% 
      arrange(-value) %>%
      group_by(Var1, Metadata_moa.x) %>%
      slice(1:1) %>%
      summarise(good = any(match.moas(Metadata_moa.x, Metadata_moa.y))) %>%
      ungroup() %>%
      filter(good) 
    return(res)
  } else if (type.eval == "lift") {
    same.moa <- match.moas
    
    u <- cr.melt %>%
      filter(as.character(Var1) < as.character(Var2) &
               Var1 != "DMSO" & Var2 != "DMSO" & !is.na(Var1) & !is.na(Var2)) %>%
      arrange(-value) %>%
      mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))
    
    x <- u$same.moa %>% as.numeric() %>% cumsum()
    x <- x/x[length(x)]
    y <- (!u$same.moa) %>% as.numeric() %>% cumsum()
    y <- y/y[length(y)]
    
    u <- data.frame(x - y)
    
    colnames(u) <- profile.type
    u <- cbind(data.frame(n = 1:NROW(u)), u)
    
    return(u)
  }
}


cr <- cor(prf[, feats] %>% t)
rownames(cr) <- prf$Metadata_broad_sample
colnames(cr) <- prf$Metadata_broad_sample

res <- evaluate.moa(cr = cr, profiles.meta = profiles.meta, quant = quant, type.eval = type.eval, skip.comp = T)
















