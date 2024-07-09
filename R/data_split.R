#' getMetadata - It splits the metadata table based on different conditions (i.e. cell type cluster, diagnosis, brain region)
#' It is an auxiliary function that we will use in getExprDataAndMetadata function
#'
#' @param globalMetadata data frame containing metadata of all the cells (one row per cell, one column per variable)
#' @param condition1 the variable 1 (column name) use for metadata split
#' @param condition2 the variable 2 (column name) use for metadata split
#' @return a list containing all partitions of initial metadata table
#' based on conditions proposed
#' @export
#' @examples

getMetadata <- function(globalMetadata, condition1, condition2=NULL) {

  # Check arguments
  if (missing(globalMetadata)) { stop ("Please provide a global metadata table\n")}
  if (missing(condition1)) { stop ("Please provide at least a variable name in order to split global metadata table\n")}

  # Split the global metadata table based on the proposed conditions
  if(!is.null(condition2)) {
    conditions <- list(globalMetadata[, condition1], globalMetadata[, condition2])
  } else {
    conditions <- list(globalMetadata[, condition1])
  }

  sp <- split(globalMetadata, conditions, drop=TRUE)

  return(sp)
}




#' metadata2Covs - It combines single-cell metadata and corresponding individual covariates in the same data frame
#' It is an auxiliary function that we will use in getExprDataAndMetadata function
#'
#' @param metadata data frame containing single-cell metadata (one row per cell, one column per variable)
#' @param cellID the variable name representing the ID of the cells
#' @param samplesCovs the data frame containing individuals covariates
#' @param sampleID the variable name representing the ID of the individuals
#' @return a data frame containing both single-cell metadata and corresponding individual covariates
#' @export
#' @examples

metadata2Covs <- function(exprData,
                          metadata,
                          cellID,
                          samplesCovs,
                          sampleID,
                          condition1,
                          condition2=NULL) {

  # Check arguments
  if (missing(metadata)) { stop("Please provide metadata table\n")}
  if (missing(cellID)) { stop("Please provide the name of the column that represents cellID\n")}
  if (missing(samplesCovs)) { stop("Please provide samplesCovs table\n")}
  if (missing(sampleID)) { stop("Please provide a the name of the column that represents sampleID\n")}

  # First of all, we order cells per sampleID
  metadata <- metadata[order(metadata[, sampleID]), ]

  # The sample IDs in the individual's covariate table must be in the same order as the sample IDs in the cell metadata table
  samplesNames <- unique(metadata[, sampleID])
  index <- match(samplesNames, samplesCovs[, sampleID])
  newSamplesCovs <- samplesCovs[index, ]
  stopifnot(identical(samplesNames, newSamplesCovs[, sampleID]))

  # We write down the number of cells of each individual
  reps <- table(metadata[, sampleID])
  index <- match(samplesNames, names(reps))
  reps <- reps[index]
  stopifnot(identical(samplesNames, names(reps)))

  # We repeat the covariates of the individuals as many times as we have cells for each one
  cellsCovs <- newSamplesCovs[rep(seq_len(nrow(newSamplesCovs)), reps), ]
  stopifnot(identical(cellsCovs[, sampleID], metadata[, sampleID]))

  # Add cellID and remove sampleID from table
  cellsCovs <- cbind(metadata[, -match(sampleID, colnames(metadata))], cellsCovs[, -match(sampleID, colnames(cellsCovs))])

  cellsCovs[, condition1] <- NULL

  if(!is.null(condition2)) {
    cellsCovs[, condition2] <- NULL
  }

  index <- match(colnames(exprData), cellsCovs[, cellID])
  cellsCovs <- cellsCovs[index, ]
  stopifnot(identical(colnames(exprData), cellsCovs[, cellID]))

  return(cellsCovs)
}



#' getExprDataAndMetadata - It creates the necessary files to create the single-cell gene co-expression networks.
#' First, it split metadata table into as many combinations as condition1-condition2 level factors using getMetadata function.
#' For each combination, it combines single-cell metadata and individuals covariates using metadata2Covs function and it obtains
#' the specific matrix for this condition.
#'
#' @param globalExprData the gene expression matrix containing all genes (rows) and all cells (columns)
#' @param globalMetadata the data frame containing metadata from all the cells (one row per cell, one column per variable)
#' @param condition1 the variable 1 (column name) used to split the metadata table
#' @param condition2 the variable 2 (column name) used to split the metadata table
#' @param cellID the variable name (column name) representing cells ID
#' @param outputDir the directory where output files will be saved
#' @param initialCondition if we already have the global gene expression matrix split, which was the criteria used for this (usually, we have one gene expression matrix
#' per cell type)
#' @param initialExprDataDir if we already have the global gene expression matrix split, in which directory these expression matrices are located
#' @param samplesCovs the data frame containing the covariates of each individual
#' @param sampleID the variable name (column name) representing the samples IF
#' @return it returns a message indicating in which directory the matrix files and the metadata files have been saved.
#' @export
#' @examples

getExprDataAndMetadata <- function(globalExprData=NULL,
                                   globalMetadata,
                                   condition1,
                                   condition2=NULL,
                                   cellID,
                                   outputDir=getwd(),
                                   initialCondition=NULL,
                                   initialExprDataDir=NULL,
                                   samplesCovs=NULL,
                                   sampleID=NULL,
                                   cutoff=100) {

  # Check arguments
  if (missing(globalMetadata)) { stop("Please provide a global covs table\n")}
  if (missing(condition1)) { stop("Please provide at least a variable name in order to split globalMatrix table\n")}
  if (missing(globalExprData) & missing(initialExprDataDir)) { stop("Please provide a global counts matrix or specify the directory of the expression files\n")}
  if (missing(cellID)) { stop("Please provide the name of the column that represents cellID\n")}

  # Load libraries
  require(Matrix)

  # Check if the directories exists and, if not, create it
  exprDataDir <- paste0(outputDir, "/exprData/")
  metadataDir <- paste0(outputDir, "/metadata/")

  if(!(dir.exists(exprDataDir))) {
    cat("Directory for the expression data not found, creating one.\n")
    dir.create(exprDataDir)
  }

  if(!(dir.exists(metadataDir))) {
    cat("Directory for the metadata not found, creating one.\n")
    dir.create(metadataDir)
  }

  # If we have a global matrix
  if (!is.null(globalExprData)) {

    print("We start fom a globalExprData\n")
    # We are going to create a metadata table for each level of condition 1
    spMetadata <- getMetadata(globalMetadata=globalMetadata, condition1=condition1, condition2=condition2)

    # For each level of condition 1
    for (i in 1:length(spMetadata)) {
      mymetadata <- spMetadata[[i]]

      if(nrow(mymetadata)>cutoff) {
        index <- match(mymetadata[, cellID], colnames(globalExprData))
        exprData <- globalExprData[, index]
        stopifnot(identical(mymetadata[, cellID], colnames(exprData)))

        name <- gsub(" ", "_", names(spMetadata)[i])
        name <- gsub("-", "_", name)

        # If the metadata file already contains all the information we want
        if (is.null(samplesCovs)) {

          stopifnot(identical(colnames(exprData), mymetadata[, cellID]))
          saveRDS(exprData, paste0(exprDataDir, name, "_", "exprData", ".rds"))
          saveRDS(mymetadata, paste0(metadataDir, name, "_", "metadata", ".rds"))

          # If we want to create a metadata file from the covariates of the individuals
        } else {
          cellsCovs <- metadata2Covs(exprData, mymetadata, cellID, samplesCovs, sampleID, condition1, condition2)

          stopifnot(identical(colnames(exprData), mymetadata[, cellID]))
          saveRDS(exprData, paste0(exprDataDir, name, "_", "exprData", ".rds"))
          saveRDS(cellsCovs, paste0(metadataDir, name, "_", "metadata", ".rds"))
        }

      } else {
        cat("Skipping combination", i, "because it has less than", cutoff, "cells\n\n")
      }

    }
  }


  # If we do not have a global gene expression matrix but instead we have this matrix divided into a set of matrices

  if (!is.null(initialCondition)) {
    print("We start from many gene expression matrices\n")

    # We are going to create a metadata table for each array using the initial condition
    spMetadata <- getMetadata(globalMetadata=globalMetadata, condition1=initialCondition)

    # We are going to list the expression files
    myFiles <- list.files(path=initialExprDataDir, pattern=".rds")

    # Now we have a covariate file and expression file for each cell type
    stopifnot(setequal(names(spMetadata), gsub(".rds", "", myFiles)))

    # The covariate tables must be in the same order as the names of the expression arrays
    index <- match(names(spMetadata), gsub(".rds", "", myFiles))
    spMetadata <- spMetadata[index]
    stopifnot(identical(names(spMetadata), gsub(".rds", "", myFiles)))

    # For each combination table of covariates-matrix of expression
    for (i in 1:length(spMetadata)) {
      metadata <- spMetadata[[i]]

      if(nrow(metadata)>cutoff) {

        # We select the normalized expression array
        exprData <- as.matrix(readRDS(paste0(initialExprDataDir, myFiles[i])))

        # Take only the expression data of the cells for which we have the metadata
        index <- match(metadata[, cellID], colnames(exprData))
        exprData <- exprData[, index]
        stopifnot(identical(metadata[, cellID], colnames(exprData)))

        # Now split this table of covariates according to the desired final condition
        newSpMetadata <- getMetadata(metadata, condition1=condition1, condition2=condition2)

        # For each new partition
        for (j in 1:length(newSpMetadata)) {
          newMetadata <- newSpMetadata[[j]]
          index <- match(newMetadata[, cellID], colnames(exprData))
          newExprData <- exprData[, index]
          stopifnot(identical(newMetadata[, cellID], colnames(newExprData)))

          firstName <- gsub(" ", "_", names(spMetadata)[i])
          firstName <- gsub("-", "_", firstName)
          secondName <- names(newSpMetadata)[j]

          cat("Creating the files for", firstName, secondName, "\n")

          saveRDS(newExprData, paste0(exprDataDir, firstName, "_", secondName, "_", "exprData", ".rds"))

          # If the metadata file already contains all the information we want
          if (is.null(samplesCovs)) {
            saveRDS(newMetadata, paste0(metadataDir, firstName, "_", secondName, "_", "metadata", ".rds"))

            # If we want to create a metadata file from the covariates of the individuals
          } else {
            cellsCovs <- metadata2Covs(newMetadata, cellID, samplesCovs, sampleID, condition1, condition2)
            stopifnot(identical(colnames(newExprData), cellsCovs[, cellID]))

            saveRDS(cellsCovs, paste0(metadataDir, firstName, "_", secondName, "_", "metadata", ".rds"))
          }
        }
      }
    }
  }
  cat("Expression files are available in the", exprDataDir, "directory\n")
  cat("The metadata files are available in the", metadataDir, "directory\n")
}




