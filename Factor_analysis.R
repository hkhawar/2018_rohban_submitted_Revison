#setwd("~/Desktop/Mohammad_paper/TA/code/")
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)

extends <- methods::extends


library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(iterators)
library(doParallel)
library(reshape2)
library(psych)
library(ggplot2)
#install.packages("readbulk")
library(readbulk)

batch.name = "SIGMA2_Pilot_2013_10_11"
# Selecting only column features from the defined feature list
feat.list <- read.table("../input/feature_list.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist() %>% unname()
plate.list <- read.table("../input/processed_plates_TA.txt", header = F) %>% as.matrix() %>% as.vector() %>% unlist()
metadata.df <- data.frame(readr::read_csv("../input/metadata_TA.csv"),  stringsAsFactors =F)
variables <- feat.list


# Combining all the metadata information from normalized csv of all plates
f.path <- NULL
for (p in 1:length(plate.list)) {
  f.path[p]<- paste0("../input/", plate.list[p], "_normalized.csv")
}

# Reading the sqlite files path
sql.path <- NULL
for (pl in seq_along(plate.list)) {
  sql.path[pl] <- as.vector(paste0("../backend/", batch.name, "/", plate.list[pl], "/", plate.list[pl], ".sqlite"))
}

# reading sqlite
read_sql<- function(sql.path) {
  db <- DBI::dbConnect(RSQLite::SQLite(), sql.path)
  RSQLite::initExtension(db)
  
  image <- RSQLite::dbReadTable(conn = db, "Image")
  cells <- RSQLite::dbReadTable(conn = db, "Cells")
  nuclei <- RSQLite::dbReadTable(conn = db, "Nuclei")
  cytoplasm <- RSQLite::dbReadTable(conn = db, "Cytoplasm")
  
  dt <- cells %>%
    left_join(cytoplasm, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
    left_join(nuclei, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
    left_join(image, by = c("TableNumber", "ImageNumber"))
  
  return(dt)
  
}

dmso.all <- data.frame(readr::read_csv("../FA/dmso.all.csv"),  stringsAsFactors =F)

dmso.all$File <- NULL
dmeta <- colnames(dmso.all)[which(str_detect(colnames(dmso.all), "Metadata"))]
variables <- setdiff(colnames(dmso.all), dmeta)

# Standardization step
mn <- apply(dmso.all %>% select(one_of(variables)), 2, function(x) mean(x, na.rm = T))
sdv <- apply(dmso.all %>% select(one_of(variables)), 2, function(x) sd(x, na.rm = T))

dmso.all <- scale(dmso.all[, variables], center = mn, scale = sdv)
dmso.all <- as.data.frame(dmso.all)

# doing the factor analysis

sample_size <- 25000
set.seed(123)
train_indx <- sample(seq_len(nrow(dmso.all)), size = sample_size)
dmso_train <- dmso.all[train_indx, ]
dmso_test <- dmso.all[-train_indx, ]

fa <- factanal(dmso_train, factors = 50, rotation = "varimax", lower = 0.05, scores = "regression")

# Calculating factor-score coefficient for the "regression" method with orthogonally rotated factors and standardize variables

coef <- solve(fa$correlation) %*% fa$loadings


profile <- NULL
for (i in 1:length(sql.path)) {
  prf <- read.csv(f.path[i], stringsAsFactors = FALSE)
  # Extracting metadata
  meta <- colnames(prf)[which(str_detect(colnames(prf), "Metadata_"))]
  pmeta <- prf %>% select(one_of(meta)) %>% dplyr::collect()
  # reading sqlite file
  sql_data <- as.data.frame(lapply(sql.path[i], read_sql))
  image.col <- sql_data %>%
    select(Image_Metadata_Well, Image_Metadata_Plate) %>%
    colnames()
  sql_data <- sql_data %>%
    select(image.col, variables) %>%
    dplyr::collect()
  # removing NAs
  sql_data[is.na(sql_data)] <- 0
  sql_data <- merge(sql_data, pmeta,
                    by.x = c("Image_Metadata_Plate","Image_Metadata_Well"),
                    by.y = c("Metadata_Plate", "Metadata_Well"))
  
  # sql_data <- merge(sql_data, pmeta,
  #                   by.x = "Image_Metadata_Well",
  #                   by.y = "Metadata_Well")
  
  metadata <- colnames(sql_data)[which(str_detect(colnames( sql_data), "Metadata_"))]
  metadata <- sql_data %>% select(one_of(metadata)) %>% dplyr::collect()
  
  dmso <- sql_data %>%
    filter(Metadata_ASSAY_WELL_ROLE== "Untreated") %>%
    dplyr::collect()
  
  # calculating mean and standard_deviation of control condition
  mn <- apply(dmso %>% select(one_of(variables)), 2, function(x) mean(x, na.rm = T))
  sdv <- apply(dmso %>% select(one_of(variables)), 2, function(x) sd(x, na.rm = T))
  print("10")
  sql_data <- scale(sql_data[, variables], center = mn, scale = sdv)
  #sql_data <- as.data.frame(sql_data)
  
  # predict factors of data using control condition
  prediction <- sql_data[, variables] %*% coef 
  # converting matrix into dataframe
  prediction <- as.data.frame(prediction)
  
  profiles.nrm  <- cbind(metadata, prediction)
  
  # selecting variables
  fa_variables <- colnames(profiles.nrm)
  fa_metavariables <- fa_variables[which(str_detect(fa_variables, "Metadata"))]
  factor_variables <- setdiff(colnames(profiles.nrm), fa_metavariables)
  
  
  
  operation <- c("mean", "median", "mad")
  
  for (l in 1:length(operation)) {
    profiles.nrm[, factor_variables] <- apply(profiles.nrm[, factor_variables],
                                              2, function(x) as.numeric(x))
    profile <- cytominer::aggregate(population = profiles.nrm,
                                    strata = c("Image_Metadata_Plate", "Image_Metadata_Well"),
                                    variables = factor_variables,
                                    operation = operation[l])
    
    profile <- cbind(profile, pmeta)
    readr::write_csv(profile, paste0("../FA/", plate.list[i], "_FA_", operation[l], ".csv"))
    
    print("13")
    
  }
  
}


print("Successfully Executed")



