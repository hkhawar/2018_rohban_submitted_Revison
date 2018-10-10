
#setwd("/Users/habbasi/Desktop/Mohammad_paper/TA/code")
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


dmeta <- colnames(dmso.all)[which(str_detect(colnames(dmso.all), "Metadata"))]
variables <- setdiff(colnames(dmso.all), dmeta)

# Standardization step
mn <- apply(dmso.all %>% select(one_of(variables)), 2, function(x) mean(x, na.rm = T))
sdv <- apply(dmso.all %>% select(one_of(variables)), 2, function(x) sd(x, na.rm = T))

# dmso.all <- scale(dmso.all[, variables], center = mn, scale = sdv)
# dmso.all <- as.data.frame(dmso.all)

# doing the Principle component analysis
pca <- prcomp(dmso.all[, variables], center=mn, scale=sdv)
pc <- (pca$rotation[,1:50])

#screeplot(pca, type = "line")



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
    print("8")
     # sql_data <- merge(sql_data, pmeta,
     #                 by.x = "Image_Metadata_Well",
     #                 by.y = "Metadata_Well")
  

    metadata <- colnames(sql_data)[which(str_detect(colnames(sql_data), "Metadata_"))]
    metadata <- sql_data %>% select(one_of(metadata)) %>% dplyr::collect()
    
    dmso <- sql_data %>%
      filter(Metadata_ASSAY_WELL_ROLE == "Untreated") %>%
      dplyr::collect()
    
    # calculating mean and standard_deviation of control condition
    mn <- apply(dmso %>% select(one_of(variables)), 2, function(x) mean(x, na.rm = T))
    sdv <- apply(dmso %>% select(one_of(variables)), 2, function(x) sd(x, na.rm = T))
    
    sql_data <- scale(sql_data[, variables], center = mn, scale = sdv)

    # prediction of PCs for test dataset
    prediction <- sql_data %*% pc

    # converting matrix into dataframe
    prediction <- as.data.frame(prediction)
    
    profiles.nrm  <- cbind(metadata, prediction)
    
    # selecting variables
    pc_variables <- colnames(profiles.nrm)
    pc_metavariables <- pc_variables[which(str_detect(pc_variables, "Metadata"))]
    pc_variables <- setdiff(colnames(profiles.nrm), pc_metavariables)
    
    
    operation <- c("mean", "median", "mad")
    
    for (l in 1:length(operation)) {
      profiles.nrm[, pc_variables] <- apply(profiles.nrm[, pc_variables],
                                                2, function(x) as.numeric(x))
      profile <- cytominer::aggregate(population = profiles.nrm,
                                      strata = c("Image_Metadata_Plate", "Image_Metadata_Well"),
                                      variables = pc_variables,
                                      operation = operation[l])
                                    
      
      profile <- cbind(profile, pmeta)
      readr::write_csv(profile, paste0("../PCA/", plate.list[i], "_PC_", operation[l], ".csv"))
      
    }
}


print("Successfully Executed")






