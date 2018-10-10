#!/usr/bin/env Rscript
#setwd("~/Desktop/Mohammad_paper/TA/code/")
#setwd("~/Desktop/Mohammad_paper/Bioactives-BBBC022-Gustafsdottir/code/")
setwd("~/Desktop/Mohammad_paper/CDRP/code/")

rm(list = ls())
args = commandArgs(trailingOnly = T)

extends <- methods::extends



library(dplyr)
library(foreach)
library(stringr)
library(readr)
library(magrittr)
library(SNFtool)
library(ggplot2)
library(readbulk)

# batch.name = "SIGMA2_Pilot_2013_10_11"
# # Selecting only column features from the defined feature list
# feat.list <- read.table("../input/feature_list.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist() %>% unname()
# plate.list <- read.table("../input/processed_plates_TA.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist()
# metadata.df <- data.frame(readr::read_csv("../input/metadata_TA.csv"),  stringsAsFactors =F)
# variables <- feat.list
# 
# Metadata_broad_sample == "Untreated"
----------------------------------------------------------------------------------------------------------
# batch.name = "BBBC022_2013"
# # Selecting only column features from the defined feature list
# feat.list <- read.table("../input/feature_list_BBBC022.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist() %>% unname()
# plate.list <- read.table("../input/processed_plates_BBBC022.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist()
# metadata.df <- data.frame(readr::read_csv("../input/metadata_BBBC022.csv"),  stringsAsFactors =F)
# variables <- feat.list
# Metadata_broad_sample == "DMSO"
----------------------------------------------------------------------------------------------------------
batch.name = "CDRP"
# Selecting only column features from the defined feature list
feat.list <- read.table("../input/feature_list.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist() %>% unname()
plate.list <- read.table("../input/processed_plates_CDRP_bio.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist()
metadata.df <- data.frame(readr::read_csv("../input/metadata_CDRP.csv"),  stringsAsFactors =F)
variables <- feat.list
Metadata_broad_sample = "DMSO"


profile.type <- "median"

quant <- 0.99

mean.na <- function(x) {mean(x, na.rm = T)}

# reading all the PCA
path <- "../PCA/"

read.and.summarize <- function(profile.type) {
    if (profile.type == "mean") {
      x <- read_bulk(directory = path, extension = paste0(profile.type,".csv"), stringsAsFactors=FALSE)
    } else if (profile.type == "median") {
      x <- read_bulk(directory = path, extension = paste0(profile.type,".csv"), stringsAsFactors=FALSE)
    } else if (profile.type == "mad") {
      x <- read_bulk(directory = path, extension = paste0(profile.type,".csv"), stringsAsFactors=FALSE)
    } else {
      x <-NULL
    }
    x
}
profiles.nrm <- read.and.summarize(profile.type)
# Removing Filename column in the end
profiles.nrm$File <- NULL 
  variable.names <- colnames(profiles.nrm)
  meta.cols <- variable.names[which(str_detect(variable.names, "Metadata_"))] 
  variable.names <- setdiff(variable.names,  meta.cols)

  
  
  # Selecting variables with SD less than 10 
  
ids <-  apply(profiles.nrm[,variable.names], 2, function(x) !any(is.na(x) | is.nan(x) | is.infinite(x) | sd(x) > 10)) %>% which

  variable.names <- variable.names[ids]
  
  profiles.nrm <- profiles.nrm %>% select(one_of(c(meta.cols, variable.names)))
  


  
  if (!"Metadata_mmoles_per_liter" %in% colnames(profiles.nrm)) {
    profiles.nrm <- profiles.nrm %>%
      mutate(Metadata_mmoles_per_liter = 10)
  }
  
  prf <- profiles.nrm %>% 
    group_by(Metadata_broad_sample, Metadata_mmoles_per_liter, Metadata_Plate_Map_Name) %>%
    summarise_at(.vars = variable.names, .funs = "mean.na")
  
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
  

  feats <- colnames(prf)
  feats <- feats[which(!str_detect(feats, "Metadata_"))]
  profiles.meta <- colnames(prf)[which(str_detect(colnames(prf), "Metadata"))]
  profiles.meta <- prf %>% select(Metadata_broad_sample, Metadata_moa) %>% dplyr::collect()
  
  # pm <- profiles.nrm %>% select(Metadata_broad_sample, Metadata_Plate_Map_Name) %>% unique 
  # profiles.meta <- profiles.meta %>% left_join(pm, by = "Metadata_broad_sample")
  # 



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
  saveRDS(cr.melt, paste0("pc_cr_",  profile.type, ifelse(profile.type == "mix", paste0(mix1, mix2), ""), ".rds"))
  
  if (skip.comp) {
    return(NULL)
  }
  
  if (type.eval == "global") {
    cr.melt <- cr.melt %>% filter(Var1 < Var2)
    
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
      slice(1:k) %>%
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
  
 









