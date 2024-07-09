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


#' getMarkersList - It returns a list of markers.
#'
#' @param type it indicates the kind of markers you want to get
#' @return the list of markers selected
#' @export
#' @examples

getMarkersList <- function(type=c("GenomicsEngland", "PD", "cellTypeMarker", "DNs", "DEGs", "ferroptosis")) {

  markersList <- list()

  # the.dir = system.file("", "markers", package = "scCoExpNets")
  the.dir <- "C:/Users/alici/OneDrive/Documentos/collabBackUp/scCoExpNets/inst/markers/"
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
    markersList <- "newList"
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

        netResult <- c(subcluster, iter, module, markersList, names(genesList)[set], total.specific, total.net, n.module, n.module.and.specific, paste0(intersect(names(net$moduleColors[net$moduleColors==module]), myGenes), collapse=", "), test$p.value)
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
                           minSize  = 10,
                           maxSize  = 500)

    df <- rbind(df, cbind(rep(subcluster, nrow(result)), rep(iter, nrow(result)), rep(module, nrow(result)), result))
  }
  colnames(df)[1:3] <- c("subcluster", "iter", "module")

  df$padj2 <- p.adjust(df$padj, method="BH")
  # df <- df[df$padj2<0.05, ]

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
                     pseudoCell=F,
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
  cat("Number of lists to check", length(markersList), "\n")

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


  if(type=="sc" & pseudoCell==T) {
    print(metadata)
    metaData <- readRDS(metadata)

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


#' myPreservationOneWay
#'
#' @param
#' @return
#' @export
#' @examples

myPreservationOneWay <- function (network, expr.data.files = NULL, tissues = c("snig",
                                                                               "putm"), permutations = 200, maxModuleSize = 5000, maxGoldModuleSize = 400,
                                  randomSeed = 1)
{
  cat("Entering preservation\n")
  n1.shortname = tissues[1]
  n2.shortname = tissues[2]
  if (typeof(network) == "character") {
    print(paste0("Reading network ", network))
    network1 <- readRDS(network)
  }
  else {
    network1 = network
  }
  print(tissue1 <- tissues[1])
  print(tissue2 <- tissues[2])
  print(paste0(tissue1, " vs. ", tissue2))
  if (tissues[1] == tissues[2])
    stop("Can't do a preservation against self")
  options(stringsAsFactors = FALSE)
  print(paste0("Reading expression data for tissue ", tissues[1]))
  if (typeof(expr.data.files[[1]]) == "character") {
    expression.data1 <- readRDS(expr.data.files[[1]])

    if(nrow(expression.data1)>1000) {
      expression.data1 <- expression.data1[sample(1:nrow(expression.data1), 1000), ]
    }

    cat(expr.data.files[[1]], "\n")
    cat(nrow(expression.data1), "\n")
  }
  else {
    expression.data1 = expr.data.files[[1]]
  }
  print(expression.data1[1:5, 1:5])
  print(paste0("Reading expression data for tissue ", tissues[2]))
  if (typeof(expr.data.files[[2]]) == "character") {
    expression.data2 <- readRDS(expr.data.files[[2]])

    if(nrow(expression.data2)>1000) {
      expression.data2 <- expression.data2[sample(1:nrow(expression.data2), 1000), ]
    }

    cat(expr.data.files[[2]], "\n")
    cat(nrow(expression.data2), "\n")
  }
  else expression.data2 = expr.data.files[[2]]
  print(expression.data2[1:5, 1:5])
  intersect.g = intersect(colnames(expression.data1), colnames(expression.data2))
  expression.data1 = expression.data1[, match(intersect.g,
                                              colnames(expression.data1))]
  expression.data2 = expression.data2[, match(intersect.g,
                                              colnames(expression.data2))]
  network1 = network1$moduleColors[match(intersect.g, names(network1$moduleColors))]
  network2 = network1
  cat("We'll use", ncol(expression.data1), "genes for the preservation analysis\n")
  multiExpr <- list()
  multiExpr[[1]] <- list(data = expression.data1)
  multiExpr[[2]] <- list(data = expression.data2)
  names(multiExpr) <- c(tissues[1], tissues[2])
  checkSets(multiExpr, checkStructure = FALSE, useSets = NULL)
  multiColor <- list(network1)
  names(multiColor) <- c(tissues[1])
  enableWGCNAThreads(4)
  print(WGCNAnThreads())
  system.time({
    mp <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1,
                             nPermutations = permutations, networkType = "signed",
                             maxGoldModuleSize = maxGoldModuleSize, randomSeed = randomSeed,
                             verbose = 3, maxModuleSize = maxModuleSize, parallelCalculation = TRUE)
  })
  return(list(stats=getPreservationStatisticsOneWay(tissues = tissues, presRes = mp), mp=mp))
}




#' getNewTerms - It returns the new terms obtained only when creating pseudo-cells iterations
#'
#' @param nets_dir the path where the networks from the same cell type are located
#' @param source the source that will be used to compared the enrichment
#' @return a data frame containing the results of the module-trait association. You will find as many rows as combinations module*covariate
#' @export a list containing A) the new terms obtained in each pseudo-cell iteration compared to T0; B) a data frame with the number of new terms
#' obtained in each pseudo-cell iteration compared to T0; c) the enrichment analysis results for all iterations.
#' @examples

getNewTerms <- function(nets_dir, source) {

  # Load library
  require(GOSim)
  require(stringr)

  # List nets
  nets <- list.files(nets_dir, recursive=T, pattern=".50.rds$", full.names=T)
  results <- data.frame()

  # Get enrichment per module per net
  for(net in nets) {
    name <- gsub("\\..*", "", gsub("net", "", basename(net)))
    n <- str_split(name, "_iter")[[1]]
    if(length(n)==1) {
      iter <- 0
      ct <- n[1]
    } else {
      iter <- n[2]
      ct <- n[1]
    }
    net <- readRDS(net)
    all.genes <- list()

    for(m in unique(net$moduleColors)) {
      all.genes[[m]] <- names(net$moduleColors[net$moduleColors==m])
    }

    enrich <- gprofiler2::gost(all.genes,
                               correction_method="g_SCS",
                               sources=source,
                               organism = "hsapiens",
                               exclude_iea = T,
                               highlight = T)$result

    enrich$ct <- rep(ct, nrow(enrich))
    enrich$iter <- rep(iter, nrow(enrich))
    results <- rbind(results, enrich)
  }

  # Calculate information content
  # setOntology("CC", loadIC=FALSE)
  # calcICs(getwd())
  # setOntology("MF", loadIC=FALSE)
  # calcICs(getwd())
  # setOntology("BP", loadIC=FALSE)
  # calcICs(getwd())

  # Load information content
  ic <- c()
  load(file="C:/Users/alici/Downloads/paperSeb/aarmdv2/files/ICsCChumanall.rda")
  ic <- c(ic, IC)
  load(file="C:/Users/alici/Downloads/paperSeb/aarmdv2/files/ICsBPhumanall.rda")
  ic <- c(ic, IC)
  load(file="C:/Users/alici/Downloads/paperSeb/aarmdv2/files/ICsMFhumanall.rda")
  ic <- c(ic, IC)
  length(ic)

  # Get IC for my GO terms
  myic <- ic[match(results$term_id, names(ic))]
  results$IC <- as.vector(myic)

  # Get statistics
  stats <- data.frame()
  iters <- unique(results$iter)
  T0 <- results[results$iter==0, ]
  new_highlighted_terms <- list()

  for (iter in iters) {
    # New terms in Tn compared to T0
    Tn <- results[results$iter==iter, ]
    Tn_terms <- Tn$term_name[Tn$highlighted==T]
    T0_terms <- T0$term_name[T0$highlighted==T]
    newTerms_h <- sort(setdiff(Tn_terms, T0_terms))
    toRemove <- c()

    if(length(newTerms_h)>0) {
      for (t in 1:length(newTerms_h)) {
        term <- newTerms_h[t]
        len <- lengths(str_split(term, " "))
        if(len==1) {
          parent <- T0_terms[grep(term, T0_terms)]
          if(length(parent)>0) {
            toRemove <- c(toRemove, t)
          }
        }
      }
    }
    if(length(toRemove)>0) {
      newTerms_h <- newTerms_h[-toRemove]
    }
    new_highlighted_terms[[iter]] <- newTerms_h
    IC_h <- round(mean(Tn$IC[Tn$highlighted==T]), 2)
    stats <- rbind(stats, c(ct, iter, IC_h, length(newTerms_h)))
  }

  colnames(stats) <- c("ct", "iter", "IC_highlighted", "new_highlighted")

  return(list(newTerms=new_highlighted_terms, stats=stats, enrich=enrich))
}

