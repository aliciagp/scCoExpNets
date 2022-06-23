#' getSubclusterName - It extracts the subcluster name and the iteration from the whole name of the network. It is an auxiliary function.
#'
#' @param netName the path of the selected network
#' @return a list containing the subcluster name and the iteration
#' based on conditions proposed
#' @export
#' @examples

getSubclusterName <- function(netName) {

  require(stringr)

  subcluster <- stringr::str_split(netName, "/Net/")[[1]][1]
  subcluster <- stringr::str_split(subcluster, "/")[[1]]
  subcluster <- subcluster[length(subcluster)]
  index <- grep("iter", subcluster)

  if(length(index)>0) {
    sub <- stringr::str_split(subcluster, "_iter")[[1]][1]
    iter <- paste0("iter", stringr::str_split(subcluster, "_iter")[[1]][2])

  } else {
    sub <- stringr::str_split(subcluster, "_iter")[[1]][1]
    iter <- "iter0"
  }
  return(list(subcluster=sub, iter=iter))
}




#' getModulesPreservedPerIter - It identifies the modules from the original network that are preserved through all pseudo-cells iterations.
#'
#' @param dir the path where the networks from the same subcluster are located or the bulk network is located
#' @param plot if plot=TRUE, plots will be shown
#' @return a list containing: i) a data frame showing which modules are preserved through all pseudo-cells iterations and which not,
#' ii) the results of the tests for each iteration.
#' @export
#' @examples

getModulesPreservedPerIter <- function(dir=getwd(),
                                       plot=F) {

  dirs <- list.dirs(dir, recursive=FALSE)
  netIter0 <- list.files(paste0(dirs[-grep("iter", dirs)], "/Net/"), pattern=".50.rds$", full.names=T)
  subcluster <- getSubclusterName(netIter0)[["subcluster"]]
  netIter0 <- readRDS(netIter0)
  dirs <- dirs[grep("iter", dirs)]
  allResults <- data.frame()

  for (i in 1:length(dirs)) {
    netIterX <- readRDS(list.files(paste0(dirs[i], "/Net/"), pattern=".50.rds$", full.names=T))
    results <- myGenCrossTabPlot(netIter0$moduleColors,
                                 netIterX$moduleColors,
                                 tissue1="Original",
                                 tissue2=paste0("Iter", i),
                                 plot=plot,
                                 plot.file=NULL)

    for(j in 1:nrow(results)) {
      length(rep(subcluster, nrow(results)))
      length(rep(rownames(results)[j], ncol(results)))

      moduleResults <- cbind(rep(subcluster, ncol(results)), rep(rownames(results)[j], ncol(results)), as.numeric(results[j, ]), colnames(results), rep(paste0("iter", i), ncol(results)))
      allResults <- rbind(allResults, moduleResults)
    }

  }

  colnames(allResults) <- c("subcluster", "originalModule", "pvalue", "iterModule", "iter")

  allResults$pvalue <- as.numeric(as.character(allResults$pvalue))
  allResults$pvalue <- p.adjust(allResults$pvalue, method="BH")
  allResults <- allResults[allResults$pvalue<0.05, ]

  df <- table(allResults$originalModule, allResults$iter)
  index <- which(rowSums(df)>=ncol(df))

  toReturn <- as.data.frame(cbind(rep(subcluster, nrow(df)), rep("iter0", nrow(df)), rownames(df), rep(0, nrow(df))))
  toReturn[index, 4] <- "preserved"
  toReturn[-index, 4] <- "nonPreserved"

  colnames(toReturn) <- c("subcluster", "iter", "module", "criteria")

  return(list(summary=toReturn, allTests=allResults))
}




#' getModulesNonPreserved - It tests if modules from the iteration !=0 are not preserved in the T0 GCN.
#' To this end, it applies a Fisher's exact test using myGenCrossTabPlot function.
#'
#' @param dir the path where the networks from the same subcluster are located
#' @param plot if plot=TRUE, plots will be shown
#' @return a data frame containing one row per module from iteration !=0.
#' If the criteria=="exclusive", it means that this module is not preserved in the original network
#' @export
#' @examples

getModulesNonPreserved <- function(dir=getwd(),
                                   plot=F) {

  dirs <- list.dirs(dir, recursive=FALSE)
  netIter0 <- list.files(paste0(dirs[1], "/Net/"), pattern=".50.rds$", full.names=T)
  subcluster <- getSubclusterName(netIter0)[["subcluster"]]
  netIter0 <- readRDS(netIter0)
  dirs <- dirs[-1]
  allResults <- vector(mode="list", length=length(dirs))


  for (i in 1:length(dirs)) {
    netIterX <- readRDS(list.files(paste0(dirs[i], "/Net/"), pattern=".50.rds$", full.names=T))
    results <- myGenCrossTabPlot(netIter0$moduleColors,
                                 netIterX$moduleColors,
                                 tissue1="Original",
                                 tissue2=paste0("Iter", i),
                                 plot=plot,
                                 plot.file=NULL)

    allResults[[i]] <- results
  }

  toReturn <- data.frame(subcluster=as.character(),
                         ier=as.character(),
                         module=as.character(),
                         criteria=as.character())

  for (iter in 1:length(allResults)) {
    for (module in 1:ncol(allResults[[iter]])) {
      index <- which(allResults[[iter]][, module]<0.05)

      if(length(index)==0) {
        toReturn <- rbind(toReturn, c(subcluster, paste0("iter", iter), colnames(allResults[[iter]])[module], "exclusive"))
      }
    }
  }

  colnames(toReturn) <- c("subcluster", "iter", "module", "criteria")

  return(toReturn)
}




#' getMarkersList - It returns a list of markers.
#'
#' @param type it indicates the kind of markers you want to get
#' @return the list of markers selected
#' @export
#' @examples

getMarkersList <- function(type=c("GenomicsEngland", "PD", "cellTypeMarker", "DNs", "DEGs", "ferroptosis")) {

  markersList <- list()

  the.dir = system.file("", "markers", package = "scGCNs")
  files = list.files(path=the.dir, full.names = T)

  if(type=="cellTypeMarker") {
    markersList = readRDS(files[grep(pattern="CellType_Markers_CoExpNets.rds",files)])

  } else if (type=="GenomicsEngland") {
    markersList = readRDS(files[grep(pattern="neurologyPanel.rds",files)])

  } else if (type=="PD") {
    markersList = readRDS(files[grep(pattern="PDgenes.rds",files)])

  } else if (type=="DNs") {
    markersList = readRDS(files[grep(pattern="DNsmarkers.rds",files)])

  } else if (type=="DEGs") {
    markersList = readRDS(files[grep(pattern="DNs_DEGs.rds",files)])

  } else if (type=="ferroptosis") {
    markersList = readRDS(files[grep(pattern="ferroptosis_markers.rds",files)])
  } else {
    stop("We didn't find any coincidence, please check the name of the marker list and try again")
  }

  return(markersList)

}




#' getFisherTest - It applies a Fisher's exact to check the overlap between a list of markers and the genes that made up each module.
#'
#' @param dir the path where the networks from the same subcluster are located
#' @param markersList which markers list we want to select to check the enrichment
#' @return A data frame containing the results of the Fisher's exact test. We will find one row for each combination module-iter-markers sublist
#' @export
#' @examples

getFisherTest <- function(dir=getwd(),
                          markersList=c("GenomicsEngland", "cellTypeMarker", "PD", "DNs", "DEGs", "ferroptosis")) {


  genesList <- list()

  if(is.character(markersList)) {
    genesList <- getMarkersList(type=markersList)
  } else {
    genesList <- markersList
  }

  # List network if bulk or networks if single-cell
  files <- list.files(dir, recursive=T)
  nets <- paste0(dir, "/", files[grep(".50.rds$", files)])
  results <- data.frame()
  adjustPvalue <- c()

  # For each set of genes
  for (set in 1:length(genesList)) {
    genes <- genesList[[set]]
    pvalue <- c()

    # For each net if single-cell
    for (net in nets) {
      subcluster <- getSubclusterName(net)[["subcluster"]]
      iter <- getSubclusterName(net)[["iter"]]

      net <- readRDS(net)
      myGenes <- genes[which(genes %in% names(net$moduleColors))]
      total.specific <- length(myGenes)
      total.net <- length(net$moduleColors)
      modules <- unique(net$moduleColors)

      # For each module
      for (module in modules) {

        n.module <- length(net$moduleColors[net$moduleColors==module])
        n.module.and.specific <- length(intersect(names(net$moduleColors[net$moduleColors==module]), myGenes))

        test <- CoExpNets::testGeneSet(n.module=n.module,
                                       n.module.and.specific=n.module.and.specific,
                                       total.specific=total.specific,
                                       total.net=total.net,
                                       test="fisher",
                                       oldform=F)

        netResult <- c(subcluster, iter, module, "PDgenesClustering", names(genesList)[set], total.specific, total.net, n.module, n.module.and.specific, paste0(intersect(names(net$moduleColors[net$moduleColors==module]), myGenes), collapse=", "), test$p.value)
        results <- rbind(results, netResult)
        pvalue <- c(pvalue, test$p.value)
      }
    }
    pvalue <- p.adjust(pvalue, method="BH")
    adjustPvalue <- c(adjustPvalue, pvalue)
  }

  results <- cbind(results, adjustPvalue)
  colnames(results) <- c("subcluster", "iter", "module", "criteria", "set", "total.specific", "total.net", "n.module", "n.module.and.specific", "overlapGenes", "pvalue", "adjust_pvalue")
  # results <- results[results$adjust_pvalue<0.05, ]

  return(results)
}




#' getNewAnnotations - It identifies the annotations that we only find when creating pseudo-cells
#'
#' @param dir the path where the networks from the same subcluster are located
#' @param types if we want to find new annotations from functional enrichment analysis (FEA) or phenotype enrichment analysis (PEA)
#' @return
#' @export
#' @examples

getNewAnnotations <- function(dir=getwd(),
                              types=c("FEA", "PEA")) {

  allFiles <- list.files(dir, recursive=T, full.names=T)
  toReturn <- list()

  for (type in types) {
    files <- allFiles[grep(type, allFiles)]
    iterFiles <- files[grep("iter", files)]
    originalFile <- setdiff(files, iterFiles)

    subcluster <- getSubclusterName(originalFile)[["subcluster"]]
    iter <- getSubclusterName(originalFile)[["iter"]]

    originalTerms <- read_csv(originalFile)
    originalTermsNames <- unique(originalTerms$term_name)
    newTerms <- data.frame()

    for (iter in 1:length(iterFiles)) {
      iterNewTerms <- read_csv(iterFiles[iter])
      iterNewTerms <- iterNewTerms[-which(iterNewTerms$term_name %in% originalTermsNames), ]

      if(type=="FEA") {
        cols <- c("query", "term_name", "p_value", "IC", "source")
      } else {
        cols <- c("module", "term_name", "adjust_pvalue", "source")
      }

      iterNewTerms <- cbind(rep(subcluster, nrow(iterNewTerms)), rep(paste0("iter", iter), nrow(iterNewTerms)), iterNewTerms[, cols])
      newTerms <- rbind(newTerms, iterNewTerms)
    }
    newTerms <- cbind(newTerms[, 1:3], "newAnnotations", newTerms[, 4:ncol(newTerms)])
    colnames(newTerms)[1:4] <- c("subcluster", "iter", "module", "criteria")
    toReturn[[type]] <- newTerms
  }

  return(toReturn)
}



#' getAnnotations -  It returns both the functional and phenotype enrichment annotations of the modules selected.
#'
#' @param df the data frame obtained from any of the previous functions that shows some interesting modules.
#' @param dir the path where the networks from the same subcluster are located
#' @return a list containing both the functional and the phenotype enrichment analysis annotations for the modules selected from each iteration
#' @export
#' @examples

getAnnotations <- function(df,
                           dir=getwd()) {

  dirs <- list.dirs(dir, recursive=T)
  pegTerms <- data.frame()
  gprofTerms <- data.frame()

  for (i in 1:nrow(df)) {
    if(df$iter[i]!="iter0") {
      name <- paste0(df$subcluster[i], "_", df$iter[i])
      nameDir <- dirs[grep(name, dirs)][1]
    } else {
      name <- df$subcluster[i]
      nameDir <- paste0(dirs[grep(name, dirs)][1], "/", name)
    }

    gprof <- read_csv(list.files(paste0(nameDir, "/Net"), pattern=".gprof.csv", full.names=T))
    gprof <- gprof[gprof$query==df$module[i], ]
    gprof <- cbind(rep(df$subcluster[i], nrow(gprof)), rep(df$iter[i], nrow(gprof)), gprof[, c("query", "term_name", "p_value", "IC", "source")])
    gprofTerms <- rbind(gprofTerms, gprof)

    peg <- read_csv(list.files(paste0(nameDir, "/Net"), pattern=".PEG.csv", full.names=T))
    peg <- peg[peg$module==df$module[i], ]
    peg <- cbind(rep(df$subcluster[i], nrow(peg)), rep(df$iter[i], nrow(peg)), peg[, c("module", "term_name", "adjust_pvalue", "source")])
    pegTerms <- rbind(pegTerms, peg)
  }

  colnames(gprofTerms)[1:2] <- c("subcluster", "iter")
  colnames(pegTerms)[1:2] <- c("subcluster", "iter")

  return(list(gprofTerms=gprofTerms, pegTerms=pegTerms))
}




#' myfgsea - It applies a gene set enrichment analysis on one module based on a list of markers selected
#'
#' @param net the path of the selected network on which we want to apply GSEA
#' @param markersList which markers list we want to use
#' @param stats a matrix containing the mean expression of each gene
#' @return a data frame containing the results of the GSEA. We will find one row per module-markers sublist enrichment.
#' @export
#' @examples

myfgsea <- function(net,
                    markersList,
                    stats) {

  require(fgsea)
  set.seed(1234)

  subcluster <- getSubclusterName(net)[["subcluster"]]
  iter <- getSubclusterName(net)[["iter"]]

  net <- readRDS(net)
  modules <- unique(net$moduleColors)
  df <- data.frame()

  for (module in modules) {
    genes <- names(net$moduleColors[net$moduleColors==module])
    mystats <- stats[match(genes, names(stats))]

    result <- fgsea::fgsea(pathways = markersList,
                           stats    = mystats,
                           eps = 0.0,
                           minSize  = 15,
                           maxSize  = 500)

    df <- rbind(df, cbind(rep(subcluster, nrow(result)), rep(iter, nrow(result)), rep(module, nrow(result)), result))
  }
  colnames(df)[1:3] <- c("subcluster", "iter", "module")

  df$padj2 <- p.adjust(df$padj, method="BH")
  df <- df[df$padj2<0.05, ]

  return(df)
}




#' getfgsea - It applies a gene set enrichment analysis on each module of each iteration of the selected cell type for a list markers
#'
#' @param dir the path where the networks from the same cell type are located
#' @param exprDataDir the path where the gene expression matrix for each cluster are located
#' @param metadataDir the path where the metadata tables for each cluster are located
#' @param markersType which markers list we want to use
#' @param type specify if we are dealing with bulk or single-cell RNA-seq data
#' @return a data frame containing the results of the GSEA. We will find one row per module-iter-markers sublist enrichment.
#' @export
#' @examples

getfgsea <- function(dir=getwd(),
                     exprData=NULL,
                     metadata=NULL,
                     markersType=c("GenomicsEngland", "cellTypeMarker", "PD", "DNs", "DEGs", "ferroptosis"),
                     type=c("bulk", "sc"),
                     myiter=F,
                     cutoff=200) {


  require(fgsea)

  finaldf <- data.frame()

  if(is.character(markersType)) {
    markersList <- getMarkersList(markersType)
    prefix <- markersType
  } else {
    markersList <- markersType
    prefix="myList"
  }
  cat("Number of markers list to check", length(markersList), "\n")

  nets <- list.files(dir, pattern=".50.rds$", recursive=T, full.names=T)
  nets <- nets[!duplicated(nets)]

  subcluster <- getSubclusterName(nets[1])[["subcluster"]]
  iter <- getSubclusterName(nets[1])[["iter"]]
  exprData <- readRDS(exprData)

  if(length(nets)==1) {
    net <- nets
  } else {
    net <- nets[-grep("iter", nets)]
  }

  if(type=="bulk") {
    exprData <- as.data.frame(scale(exprData, center = TRUE, scale = TRUE))
  }

  stats <- apply(exprData, 1, mean)
  iter <- 0

  df <- myfgsea(net, markersList, stats)
  finaldf <- rbind(finaldf, df)


  if(type=="sc" & myiter==T) {
    print(metadataDir)
    metaData <- readRDS(metadataDir)

    while(ncol(exprData) > cutoff) {
      iter <- iter + 1
      myPseudoCells <- pseudoCells(exprData, metaData)
      newExprData <- myPseudoCells[["pseudoCellsExprData"]]
      newMetaData <- myPseudoCells[["pseudoCellsCovs"]]
      stats1 <- apply(newExprData, 1, mean)
      net <- nets[grep(paste0("iter", iter), nets)]

      df <- myfgsea(net, markersList, stats1)
      finaldf <- rbind(finaldf, df)

      exprData <- newExprData
      metaData <- newMetaData
    }
  }

  finaldf <- cbind(finaldf[, 1:3], rep(prefix, nrow(finaldf)), finaldf[, 4:ncol(finaldf)])
  colnames(finaldf)[4] <- "criteria"

  return(finaldf=finaldf)
}



#' getModuleTraitCorr - It applies a linear regression approach to evaluate the association between the eigengenes of each module and the covariates (age, diagnosis or sex)
#'
#' @param nets_dir the path where the networks from the same cell type are located
#' @param metadata the data frame containing the covariates information for each cell
#' @param cellID the name of the column that contains the cell id
#' @param covariates the possible covariates to be evaluated
#' @return a data frame containing the results of the module-trait association. You will find as many rows as combinations module*covariate
#' @export
#' @examples

getModuleTraitCorr <- function(nets_dir,
                               metadata,
                               cellID="cell_id",
                               covariates=c("Clinical_Dx", "Donor_age_years", "Sex")) {

  require(stringr)

  results <- data.frame()

  # Load DNs modules eigengenes
  nets <- list.files(nets_dir, pattern=".50.rds$", full.names=T, recursive=T)

  # Load metadata
  if(is.character(metadata)) {
    metadata <- readRDS(metadata)
  }

  metadata$sampleID <- gsub("\\_.*", "", metadata$cell_id)
  metadata <- metadata[, which(colnames(metadata) %in% c("sampleID", cellID, covariates))]

  # Remove donors with less than 100 cells
  tab <- table(metadata$sampleID)
  IDs <- names(tab)[which(tab>100)]
  metadata <- metadata[which(metadata$sampleID %in% IDs), ]

  for (net in nets) {
    iter <- substring(str_split(net, "iter")[[1]][3],1,1)
    if(is.na(iter)) {iter=0}
    net <- readRDS(net)
    MEs <- net$MEs

    cells_selected <- intersect(rownames(MEs), metadata$cell_id)

    MEs <- MEs[which(rownames(MEs) %in% cells_selected), ]
    mymetadata <- metadata[which(metadata$cell_id %in% cells_selected), ]

    stopifnot(identical(mymetadata$cell_id, rownames(MEs)))

    mymetadata <- mymetadata[, -1]

    mymodel <- getModuleTraitCorrModel(MEs, mymetadata)
    iter <- rep(paste0("iter", iter), nrow(mymodel))
    results <- rbind(results, cbind(iter, mymodel))

  }

  return(results)
}


#' getModuleTraitCorrModel - It applies a linear regression approach to evaluate the association between the eigengenes of one module and the covariates (age, diagnosis or sex)
#'
#' @param MEs the data frame containing the module eigengenes
#' @param metadata the data frame containing the covariates information for each cell
#' @return a data frame containing the results of the module-trait association. This data frame contains as many rows as covariates.
#' @export
#' @examples

getModuleTraitCorrModel <- function(MEs, metadata) {

  require(caret)
  results <- data.frame()

  for (i in 1:ncol(MEs)) {
    mydf <- as.data.frame(cbind(MEs[, i], metadata))
    colnames(mydf)[1] <- "eigengenes"
    mymodel <- lm(eigengenes~., data=mydf)
    sum <- summary(mymodel)

    results <- rbind(results, c(gsub("ME", "", colnames(MEs)[i]),
                                format(as.numeric(sum$coefficients[, 4][2]),sci=T,digits=4),
                                format(as.numeric(sum$coefficients[, 4][3]), sci=T, digits=4),
                                format(as.numeric(sum$coefficients[, 4][4]), sci=T, digits=4),
                                round(sum$adj.r.squared,3),
                                round(RMSE(as.vector(unlist(mymodel$fitted.values)), mydf$eigengenes),3),
                                round(cor(as.vector(unlist(mymodel$fitted.values)), mydf$eigengenes),3),
                                format(mean(as.vector(unlist(mymodel$residuals))), sci=T, digits=4)))
  }

  colnames(results) <- c("module", "dx_pval", "age_pval", "sex_pval", "adj.r.squared", "rmse", "cor", "mean_resids")

  return(results=results)

}

