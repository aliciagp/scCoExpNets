#' getMarkersList - It returns a list of markers.
#'
#' @param type it indicates the kind of markers you want to get
#' @return the list of markers selected
#' @export

getMarkers <- function(type=c("cellTypeMarkers", "dopaminergicMarkers", "ferroptosisMarkers", "PDgenes")) {

  markersList <- list()

  # the.dir = system.file("", "inst/data/markers/", package = "scCoExpNets")
  the.dir <- "C:/Users/alici/OneDrive/Documentos/collabBackUp/scCoExpNets/inst/data/markers/"
  files = list.files(path=the.dir, full.names = T)

  if(type=="cellTypeMarkers") {
    markersList = readRDS(files[grep(pattern="cellType_markers_CoExpNets.rds",files)])

  } else if (type=="PDgenes") {
    markersList = readRDS(files[grep(pattern="PD_GWAS_gene_markers.rds",files)])

  } else if (type=="dopaminergicMarkers") {
    markersList = readRDS(files[grep(pattern="DNs_scRNAseq_markers_curated_ownlist.rds",files)])

  } else if (type=="ferroptosisMarkers") {
    markersList = readRDS(files[grep(pattern="ferroptosis_markers_FerrDbV2.rds",files)])
  } else {
    stop("We didn't find any coincidence, please check the name of the marker list and try again")
  }

  return(markersList)

}



#' getMarkersListEnrich - It applies a gene set enrichment analysis on each module of the network for a list of markers.
#' @param net the name of the network
#' @param markersList the list of markers to be used for the enrichment. It can be one of the list markers available in scCoExpNets or a list of markers provided by the user.
#' @param exprData the gene expression matrix used to create the network
#' @return a data frame containing the results of the GSEA. We will find one row per module-markers.
#' @export

getMarkersListEnrich <- function(net,
                                 markersList="cellTypeMarkers",
                                 exprData) {

  if(is.character(net)) {
    net <- readRDS(net)
  }

  if(is.character(exprData)) {
    exprData <- readRDS(exprData)
  }

  if(is.character(markersList)) {
    markersList <- getMarkers(markersList)
  }

  stats <- apply(exprData, 1, mean)
  modules <- unique(net$moduleColors)
  results <- data.frame()

  for (module in modules[1:2]) {
    genes <- names(net$moduleColors[net$moduleColors==module])
    mystats <- stats[match(genes, names(stats))]

    res <- fgsea::fgsea(pathways = markersList,
                        stats    = mystats,
                        eps = 0.0,
                        minSize  = 10,
                        maxSize  = 500)

    res <- data.frame(module=rep(module, nrow(res)), res)
    results <- rbind(results, res)
  }

  results$padj2 <- p.adjust(results$padj, method="BH")
  results$leadingEdge <- unlist(lapply(results$leadingEdge, function(x) paste0(x, collapse=", ")))

  return(results)
}



#' getTraitCorr - It applies a linear regression approach to evaluate the association between the eigengenes of one module and donors' covariates
#'
#' @param net the network with the modules to be annotated
#' @param metadata the data frame containing the covariates information for each cell
#' @param covsToCorr covariates to evaluate
#' @return a data frame containing the results of the module-trait association. This data frame contains as many rows as covariates.
#' @export

getTraitCorr <- function(net, metadata, covsToCorr=NULL) {

  require(caret)

  if(is.character(net)) {
    net <- readRDS(net)
  }

  results <- data.frame()
  MEs <- net$MEs

  if(!is.null(covsToCorr)) {
    metadata <- metadata[, match(covsToCorr, colnames(metadata))]
  }

  for (i in 1:ncol(MEs)) {

    df <- data.frame(MEs=MEs[, i], metadata)
    model <- lm(MEs~., data=df)
    sum <- summary(model)

    results <- rbind(results, c(gsub("ME", "", colnames(MEs)[i]),
                                format(as.numeric(sum$coefficients[-1, 4]), digits=4, sci=T),
                                round(sum$adj.r.squared,3),
                                round(RMSE(as.vector(unlist(model$fitted.values)), df$MEs),3),
                                round(cor(as.vector(unlist(model$fitted.values)), df$MEs),3),
                                format(mean(as.vector(unlist(model$residuals))), sci=T, digits=4)))

  }

  colnames(results) <- c("module", paste0(colnames(metadata), "_pval"), "adj.r.squared", "rmse", "cor", "mean_resids")

  return(results)

}


#' getFunctionalEnrichment - It returns the results of a functional enrichment analysis applied on each module of the selected network.
#'
#' @param net the pathway of the network selected for the functional enrichment analysis
#' @param correction.method the correction method used for functional enrichment analysis (deafult method is gSCS)
#' @param sources the databases we want to use for functional enrichment analysis (default databases are GO, KEGG and REAC)
#' @param organism the specie of the organism
#' @param exclude.iea if exclude.iea=TRUE, we remove automatically annotations (not manually curated)
#' @return it returns the data frame containing the results of gprofiler2, where each row represents
#' one annotation associated with a specific module
#' @export

getFunctionalEnrichment <- function(net,
                                    correction_method="gSCS",
                                    sources=c("GO", "KEGG", "REAC"),
                                    organism="hsapiens",
                                    exclude_iea=T) {

  if(is.character(net)) {
    net <- readRDS(net)
  }

  modules <- unique(net$moduleColors)
  all.genes <- list()

  for(module in modules){
    genes <- names(net$moduleColors)[net$moduleColors == module]
    all.genes[[module]] <- genes
  }

  go <- gprofiler2::gost(query=all.genes,
                         correction_method=correction_method,
                         sources=sources,
                         organism=organism,
                         exclude_iea=exclude_iea)

  go <- as.data.frame(go$result)
  go$parents <- unlist(lapply(go$parents, function(x) paste0(x, collapse=", ")))

  return(go)
}




#' getPhenotypeTerms - It returns the results of a phenotype enrichment analysis applied on each module of the selected network.
#'
#' @param net.file the pathway of the network selected for the functional enrichment analysis
#' @return it returns the data frame containing the results of the phenotype enrichment analysis, where each row represents
#' one annotation associated with a specific module
#' @export

getPhenotypeTerms <- function(net) {

  require(PhenoExam)

  if(typeof(net) == "character")
    net <- readRDS(net)

  modules = unique(net$moduleColors)

  databases <- getdbnames()

  results <- data.frame()

  for (module in modules) {
    genes <- names(net$moduleColors[net$moduleColors==module])
    enrich <- PhenoExam::PhenoEnrichGenes(genes= genes, database=databases, url=F)
    enrich <- enrich$alldata
    enrich <- data.frame(module=rep(module, nrow(enrich)), enrich)
    results <- rbind(results, enrich)
  }

  # Removing not interesting terms
  index <- which(results$term_name %in% c("No diseases associated", "No CRB phenotype", "No HPO phenotype"))

  if(length(index)>0) {
    results <- results[-index, ]
  }

  return(results)

}







