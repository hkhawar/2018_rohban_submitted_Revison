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
library(xtable)
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


dmso <- NULL

for (i in 1:length(sql.path)) {
  prf <- read.csv(f.path[i], stringsAsFactors = FALSE)
  # Extracting metadata
  meta <- colnames(prf)[which(str_detect(colnames(prf), "Metadata_"))]
  pmeta <- prf %>% select(one_of(meta)) %>% dplyr::collect()
  # reading sqlite file
  sql_data <- as.data.frame(lapply(sql.path[i], read_sql))
  image.col <- sql_data %>% select(Image_Metadata_Well, Image_Metadata_Plate) %>% colnames()
  sql_data <- sql_data %>% select(image.col, variables) %>% dplyr::collect()
  # removing NAs
  sql_data[is.na(sql_data)] <- 0
  sql_data <- merge(sql_data, pmeta,
                    by.x = c("Image_Metadata_Plate","Image_Metadata_Well"),
                    by.y = c("Metadata_Plate", "Metadata_Well"))
  # sql_data <- merge(sql_data, pmeta,
  #                   by.x = "Image_Metadata_Well",
  #                   by.y = "Metadata_Well")
  metadata <- colnames(sql_data)[which(str_detect(colnames(sql_data), "Metadata_"))]
  dmso <- sql_data %>%
    filter(Metadata_ASSAY_WELL_ROLE == "Untreated") %>%
    dplyr::collect()
  set.seed(123)
  sample_size <- 455
  train_indx <- sample(seq_len(nrow(dmso)), size = sample_size)
  dmso <- dmso[train_indx, ]
  readr::write_csv(dmso, paste0("../FA/dmso/", plate.list[i], "_dmso", ".csv"))
  
}
# combining all dmso cells collected from each replicate plate

path <- "../FA/dmso/"
dmso.all <- read_bulk(directory = path, extension = "_dmso.csv", stringsAsFactors=FALSE)
readr::write_csv(dmso.all, paste0("../FA/dmso/", "dmso.all", ".csv"))

print("Successfully Executed")