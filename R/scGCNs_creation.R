#' pair - It combines single-cell metadata and corresponding individual covariates in the same data frame
#' It is an auxiliary function for pseudoCells function
#'
#' @param indexes it represents a range of numbers from 1 to the number of columns of the gene expression matrix
#' @param pair it represents the ID of the cells. To group the cells per individual, it is necessary that individual ID is included in cell ID
#' @return a data frame containing pairs of cells grouped per individual ID. The data frame contains three columns: cell1, cell2, sampleID
#' @export
#' @examples

pair <- function(indexes, paired){

  set.seed(1234)

  # Check arguments
  if (missing(indexes)) { stop("Please provide indexes\n")}
  if (missing(paired)) { stop("Please provide ids\n")}

  stopifnot(length(indexes) == length(paired))

  # Deal with ids with one sample
  iout = NULL
  mytab = table(paired)
  keys = names(mytab)[mytab == 1]

  for(key in keys){
    positions = as.character(paired) == key
    lindexes = indexes[positions]
    ids = which(positions)
    pairs = do.call(rbind,lapply(seq(1,length(lindexes),2),
                                 function(x){
                                   c(lindexes[x],lindexes[x])
                                 }))

    iout = rbind(iout,cbind(pairs,rep(key,nrow(pairs))))
  }

  # Deal with ids with a set of samples
  keys = names(mytab)[mytab > 1]
  for(key in keys){
    positions = as.character(paired) == key
    lindexes = indexes[positions]
    ids = which(positions)
    pairs = do.call(rbind,lapply(seq(1,length(lindexes)-1,2),
                                 function(x){
                                   c(lindexes[x],lindexes[x+1])
                                 }))

    iout = rbind(iout,cbind(pairs,rep(key,nrow(pairs))))
  }

  return(iout)
}




#' pseudoCells - It creates a pseudo-cells matrix and the corresponding pseudo-cells metadata table.
#'
#' @param exprData it represents the gene expression matrix of a specific cell type cluster
#' @param covs it represents the metadata table of a specific cell type cluster
#' @param save if we want to save the created files
#' @param outputDir if save=TRUE, in which directory we would like to save the files created
#' @return it returns a list where the first element represents the pseudo-cells matrix and the second one represents the pseudo-cells metadata table
#' @export
#' @examples

pseudoCells <- function(exprData, covs, save=F, outputDir=paste0(getwd(), "/pseudoCells/")) {

  require(stringr)

  # Check arguments
  if (missing(exprData)) { stop("Please provide an gene expression matrix\n")}
  if (missing(covs)) { stop("Please provide a table of covariates\n")}

  cat("Creating pseudo-cells from", ncol(exprData), "cells \n")

  # Check the name of the cells to get the indexes correctly
  colnames(exprData) <- stringr::str_extract_all(colnames(exprData), "[0-9]+")

  # Get the indices using the pair function
  round1 = pair(1:ncol(exprData), colnames(exprData))

  # Create pseudo-cells gene expression matrix from these indices
  r1expr <- apply(round1, 1, function(x) {exprData[, as.numeric(x[1])] + exprData[, as.numeric(x[2])]})
  colnames(r1expr) <- round1[, 3]
  rownames(r1expr) <- rownames(exprData)

  # Update the covariate table to match with the new expression array
  r1covs <- covs[as.numeric(round1[,1]), ]

  # Save the results
  myResults <- NULL
  myResults[[1]] <- r1expr
  myResults[[2]] <- r1covs
  names(myResults) <- c("pseudoCellsExprData", "pseudoCellsCovs")

  cat(ncol(r1expr), "pseudo-cells created\n")

  # Save the results as files
  if (save==T) {
    name = deparse(substitute(exprData))

    if(!(dir.exists(outputDir))){
      cat("Directory for the pseudo-data not found, creating one.")
      dir.create(outputDir)
    }

    saveRDS(myResults, outputDir, paste0(name, "PseudoCells.rds"))

    cat("The pseudo-cells file has been stored in the", outputDir, "directory.\n")
  }
  return(myResults)
}




#' fromEnsemblIDtoExternalName - It changes from ensemble genes nomenclature to HGNC nomenclature if necessary.
#' It is an auxiliary function.
#'
#' @param genes a list of ensembl gene names
#' @return it returns a list where the first element represents the external gene names and the second one represents the ensembl gene names
#' @export
#' @examples

fromEnsemblIDtoExternalName <- function(genes) {

  # Check arguments
  if (missing(genes)) { stop("Please provide a list of genes\n")}

  # Load library
  require(biomaRt)

  match <- grep("ENSG", genes)

  if (length(match)>0) {
    cat("A nomenclature of ensemble_gene_id has been detected, changing to external gene name...\n")
    cat("Number of genes before changing the nomenclature:", length(genes), "\n")

    # Get the external gene names
    mart=biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
    genesNames <- biomaRt::getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                                filters="ensembl_gene_id",
                                values=genes,
                                mart=mart,
                                useCache=FALSE)

    # Remove duplicated genes
    duplicatedGenes <- duplicated(genesNames$external_gene_name)
    genesNames <- genesNames[grep("FALSE", duplicatedGenes), ]

    # Change the nomenclature of genes in the original matrix
    index <- match(genes, genesNames$ensembl_gene_id)
    genesNames2 <- genesNames[index[!is.na(index)], ]
    genes_HGCN <- genesNames2$external_gene_name
    genes_ensembl <- genesNames2$ensembl_gene_id

    cat("Number of genes after changing the nomenclature:", length(genes_HGCN), "\n")
    cat("Showing some examples of gene names with the new nomenclature", genes_HGCN[1:5], "\n")

    return(list(genes_HGCN=genes_HGCN, genes_ensembl=genes_ensembl))

  } else {
    cat("An external gene name nomenclature has been detected, nothing needs to be changed.\n")
    return(genes)
  }
}




#' prepareData - It checks if the gene expression matrix and the metadata files are ready for batch effect correction and/or surrogate variable analysis.
#'
#' @param exprData a normalized cell type cluster matrix
#' @param covs a cell type cluster metadata table
#' @param batchCov the name of the variable that represents the batch
#' @param idCov the name of the variable that represents the individuals ID
#' @param filteringGenes if true, genes are filtered based on their minimum expression (filteringCutoff) in a percentage of cells (filteringProportion)
#' @param filteringCuroff if filteringGenes=TRUE, filteringCutoff represents the minimum value expression
#' @param filteringProportion if filteringGenes=TRUE, filteringProportion represents the proportion of cells
#' @return it returns a list containing our data ready. The first element represents the gene expression matrix and the second one, the metadata table.
#' @export
#' @examples

prepareData <- function(exprData, covs, batchCov=NULL, idCov, filteringGenes=F, filteringCutoff=NULL, filteringProportion=NULL) {

  # Check arguments
  if (missing(exprData)) { stop("Please provide an gene expression matrix\n")}
  if (missing(covs)) { stop("Please provide a table of covariates\n")}
  if (!is.null(batchCov) && is.null(idCov)) { stop("Please provide the name of the variable that represents the id of the samples\n")}

  # Load library
  require(preprocessCore)
  require(stringr)
  require(CoExpNets)
  require(biomaRt)

  cat("Preparing data\n")
  cat("Checking the name of cells in the gene expression matrix\n")
  # Check the name of the cells (columns)
  covs <- as.data.frame(covs)

  cat("Checking that all variables are numeric except ID (and batch if there is one)\n")
  # Convert all covariates to numeric type except ID and batch effect
  covs <- as.data.frame(covs)
  covs[, -match(idCov, colnames(covs))] <- apply(covs[, -match(idCov, colnames(covs))], 2, function(x) as.numeric(as.character(x)))
  covs[, idCov] <- as.character(covs[, idCov])

  if (!is.null(batchCov)) {

    cat("Check batch effect\n")

    # Check batch effect as factor
    covs[, batchCov] <- as.factor(as.character(covs[, batchCov]))

    # Check complete cases
    covs <- covs[complete.cases(covs[ , batchCov]), ]

    # Check number of samples of each batch effect
    batchDistr <- table(covs[, batchCov])
    batchEffectsToSelect <- names(batchDistr[batchDistr>1])
    covs <- covs[covs[, batchCov] %in% batchEffectsToSelect, ]
    exprData <- exprData[, colnames(exprData) %in% covs[, idCov]]

    colnames(exprData) <- covs[ , idCov]
  }

  colnames(exprData) <- unlist(lapply(stringr::str_extract_all(colnames(exprData), "[0-9]+"), function(x) x[1]))

  cat("Number of genes before filtering:", dim(exprData)[1], "\n")
  if (filteringGenes) {
    expressedGenes <- rowSums(exprData > filteringCutoff)  > (filteringProportion * ncol(exprData))
    filteredData <- exprData[expressedGenes,]
    cat("- Number of genes with expression >", filteringCutoff, "in at least", filteringProportion, "of the samples:", nrow(filteredData), "\n")

  } else {
    filteredData <- exprData
  }

  # Check the name of the genes (rows)
  cat("Checking the name of genes in the gene expression matrix\n")

  match <- grep("ENSG", rownames(filteredData))

  if (length(match)>0) {

    genesNames <- fromEnsemblIDtoExternalName(rownames(filteredData))

    if (length(genesNames[["genes_ensembl"]])<nrow(filteredData)) {
      filteredData <- filteredData[match(genesNames[["genes_ensembl"]], rownames(filteredData)), ]
      rownames(filteredData) <- genesNames[["genes_HGCN"]]
    } else {
      rownames(filteredData) <- genesNames[["genes_HGCN"]]
    }
    stopifnot(identical(rownames(filteredData), genesNames[["genes_HGCN"]]))
  }

  cat("Quantile-normalisation\n")
  genesNames = rownames(filteredData)
  samplesNames = colnames(filteredData)
  exprDataQN = preprocessCore::normalize.quantiles(as.matrix(filteredData))
  rownames(exprDataQN) <- genesNames
  colnames(exprDataQN) <- samplesNames
  cat("All done\n")

  return(list(normalized_expression=as.matrix(exprDataQN), covariates=covs))
}




#' plotMDS_Batch - It plots the cells based on their expression using a MDS plot in order to visualize a possible batch effect
#'
#' @param exprData a normalized cell type cluster matrix
#' @param covs a cell type cluster metadata table
#' @param batchCov the name of the variable that represents the batch
#' @param name the name of the network
#' @param save if save=TRUE, the plot will be saved
#' @param outputDir if save=TRUE, in which directory the plot will be saved
#' @return it shows a MDs plot with a small set of cells
#' @export
#' @examples

plotMDS_Batch <- function(exprData, covs, batchCov, name,
                          save=F, outputDir=getwd()) {

  # Check arguments
  if (missing(exprData)) { stop("Please provide an gene expression matrix\n")}
  if (missing(covs)) { stop("Please provide a table of covariates\n")}

  # Load library
  require(limma)

  # Represent batch effect using limma::plotMDs function
  name= deparse(substitute(exprData))
  covs[, batchCov] <- as.factor(covs[, batchCov])

  # Take a sample to make the representation
  mask = sample(1:ncol(exprData), 100)

  # Color each sample according to the batchCov level to which it belongs
  colors = rainbow(length(unique(as.numeric(covs[,batchCov]))))
  finalcolors = colors[as.numeric(covs[, batchCov])][mask]

  # If save==TRUE, open a pdf file
  if (save == TRUE){
    cat("Starting MDS plots\n")
    pdf(file = paste0(outputDir,"/", gsub(" ", "_", name),
                      "_MDS_using_", batchCov,".pdf"))
  }
  limma::plotMDS(exprData[, mask], col=finalcolors,
                 main=paste0(name, "_MDS_using_", batchCov))
  legend("topright", fill=colors,
         legend=levels(covs[,batchCov]))

  if (save==TRUE) {
    dev.off()
    cat("Done\n")
  }

}




#' combatCorrect - It corrects the batch effect from the gene expression matrix using comBat function
#'
#' @param preparedExprData the gene expression matrix of a specific cell type cluster obtained from prepareData function
#' @param preparedCovs the metadata table of a specific cell type cluster obtained from prepareData function
#' @param batchCov the name of the variable that represents the batch
#' @param bioCovsToCorrect other covariates besides batch
#' @param PCA_plots if true, PCA plots will be saved
#' @param name the name of the cell type cluster (i.e.microglia_cluster1)
#' @param outputDir if save=TRUE, in which directory the plot will be saved
#' @return it returns the gene expression matrix corrected by batch effect
#' @export
#' @examples

combatCorrect <- function(preparedExprData,
                          preparedCovs,
                          batchCov,
                          bioCovsToCorrect,
                          PCA_plots = F,
                          name=deparse(substitute(preparedExprData)),
                          outputDir = paste0(getwd(), "pcaPlots/")) {

  # Check arguments
  if (missing(preparedExprData)) { stop("Please provide an gene expression matrix\n")}
  if (missing(preparedCovs)) { stop("Please provide a table of covariates\n")}
  if (missing(batchCov)) { stop("Please provide the name of the variable that represents the batch effect\n")}
  if (missing(bioCovsToCorrect)) { stop("Please provide the name of the biological variables that you want not to correct\n")}

  # Load libraries
  require(sva)
  require(swamp)
  require(CoExpNets)

  # Prepare data
  preparedCovs[, bioCovsToCorrect] <- as.numeric(as.character(preparedCovs[, bioCovsToCorrect]))

  # Correct data
  cat("Correcting batch effect with combat.\n")
  combatExprData = sva::ComBat(dat=as.matrix(preparedExprData), batch=as.factor(preparedCovs[, batchCov]),
                          mod=model.matrix(~1, data=as.data.frame(preparedCovs[, bioCovsToCorrect])))
  combatExprData = combatExprData - min(combatExprData)

  cat("Done\n")

  if (PCA_plots ==T){
    if(!(dir.exists(outputDir))){
      cat("Directory for the plots not found, creating one.\n")
      dir.create(outputDir)
    }
    cat("PCA uncorrected and plotting in process\n")
    pdf(file = paste0(outputDir, name, "_PCA_uncorrected.pdf"))
    pcres = swamp::prince(preparedExprData, preparedCovs[, c(batchCov, bioCovsToCorrect)], top=20)
    CoExpNets::princePlot(prince=pcres,
                          main=paste0("PCA for" ,
                                      name, "\n uncorrected"))
    dev.off()
    cat("PCA uncorrected done.\nPCA corrected and plotting in process\n")
    pdf(file = paste0(outputDir, name, "_PCA_corrected.pdf"))
    pcres = swamp::prince(combatExprData, preparedCovs[, c(batchCov, bioCovsToCorrect)],top=20)
    CoExpNets::princePlot(prince=pcres,
                          main=paste0("PCA for" ,
                                      name, "\n corrected"))
    dev.off()
    cat("PCA corrected done\n")
  }
  return(combatExprData)
}



#' svaCorrection - It applies a surrogate variable analysis in order to remove biological and/or technical effects that we are not interested in
#'
#' @param combatExprData the gene expression matrix of a specific cell type cluster obtained from combatCorrect function
#' @param preparedCovs the metadata table of a specific cell type cluster obtained from prepareData function
#' @param bioCovsToCorrect the covariates (biological, technical or both) that we are not interested in and we want to remove their effect on our dataset
#' @param idCov the name of the variable that represents the individual ID
#' @param batchCov the name of the variable that represents the batch
#' @param plots if plots=TRUE, plots will be saved
#' @param outputDir if save=TRUE, in which directory the plot will be saved
#' @param name the name of the cell type cluster (i.e.microglia_cluster1)
#' @param save if save=TRUE, the residuals file will be saved
#' @param residsOutout if save=TRUE, the directory where the residuals file will be saved
#' @return it returns the residuals matrix (genes at columns and cells at rows)
#' @export
#' @examples

svaCorrection <- function(combatExprData,
                         preparedCovs,
                         bioCovsToCorrect,
                         idCov,
                         batchCov=NULL,
                         plots = T,
                         outputDir = paste0(getwd(),"/svas_plots/"),
                         name=deparse(substitute(combatExprData)),
                         save = T,
                         residsOutput=paste0(getwd(), "/residuals/")) {

  # Check arguments
  if (missing(combatExprData)) { stop("Please provide an expression matrix\n")}
  if (missing(preparedCovs)) { stop("Please provide a table of covariates\n")}
  if (missing(bioCovsToCorrect)) { stop("Please provide the name of the biological variables that you want not to correct\n")}

  # Load libraries
  require(CoExpNets)
  require(swamp)
  require(sva)
  require(dplyr)
  require(gplots)

  # Prepare data
  preparedCovs <- as.data.frame(preparedCovs)
  preparedCovs[, idCov] <- as.character(preparedCovs[, idCov])

  if(!is.null(batchCov)) {
    preparedCovs[, batchCov] <- as.factor(as.character(preparedCovs[, batchCov]))
  }

  if(length(bioCovsToCorrect)>1) {
    preparedCovs[, bioCovsToCorrect] <- apply(preparedCovs[, bioCovsToCorrect], 2, function(x) as.numeric(as.character(x)))
  } else {
    preparedCovs[, bioCovsToCorrect] <- as.numeric(as.character(preparedCovs[, bioCovsToCorrect]))
  }

  # Creating the models
  myCovs <- paste0(bioCovsToCorrect, collapse=" + ")
  mm = model.matrix(~ eval(parse(text=myCovs)), data=as.data.frame(preparedCovs))
  nullmm = model.matrix(~ 1,data=as.data.frame(preparedCovs))

  cat("Launching svaseq, this will take some time\n")
  svas=sva::svaseq(dat=as.matrix(combatExprData),mod=mm,mod0=nullmm)
  colnames(svas$sv) <- paste0("SV", seq(1:ncol(svas$sv)))

  if (svas$n.sv != 0){
    ## LM correction.
    cat("\nStarting linear model correction\n")
    numericCovs = dplyr::select_if(preparedCovs, is.numeric) # In this way, we avoid selecting the ID and batch
    index <- match(bioCovsToCorrect, colnames(numericCovs))
    numericCovs = numericCovs[, -index] # We remove bioCovsToCorrect because we want to study the association between the SVs and the covariates of interest

    cat("Studying the correlation between svs and the covariates:", paste0(colnames(numericCovs), collapse=", "), "\n")

    linp = matrix(ncol=svas$n.sv,nrow=ncol(numericCovs))
    rownames(linp) = colnames(numericCovs)
    colnames(linp) = paste0("SV",1:svas$n.sv)
    linp[] = 0

    for(cov in 1:ncol(numericCovs)){
      for(sva in 1:svas$n.sv){ if(svas$n.sv == 1)
        axis = svas$sv else
          axis = svas$sv[,sva]
        linp[cov,sva] = cor.test(as.numeric(unlist(numericCovs[,cov])),axis)$p.value
      }
    }
    smallest = -10
    linp10 <- log10(linp)
    linp10 <- replace(linp10, linp10 <= smallest, smallest)
    tonote <- signif(linp, 1)

    # We remove the SVs that are significantly associated with some covariate of interest.
    linp <- as.data.frame(linp)
    mins <- sapply(linp, min)
    svToRemove <- names(mins[which(mins < 0.05)])

    if (length(svToRemove>0)) {
      index <- which(colnames(svas$sv) %in% svToRemove)
      svas$sv <- svas$sv[, -index]
    }
  }

  if (svas$n.sv != 0) {
    covs.rs = as.data.frame(preparedCovs[,match(bioCovsToCorrect, colnames(preparedCovs))])
    cat("Creating residuals taking into account", paste0(colnames(svas$sv), collapse=", "), "\n")
    resids <- apply(combatExprData, 1, function(y){
      lm( y ~ . , data=cbind(covs.rs,svas$sv))$residuals
    })

  } else{
    covs.rs = preparedCovs[,match(bioCovsToCorrect,colnames(preparedCovs))]
    cat("Creating residuals taking into account all svs \n")
    resids <- apply(combatExprData, 1, function(y){
      lm( y ~ . , data=cbind(covs.rs))$residuals
    })
  }


  if ((plots == T) & (svas$n.sv != 0) & (ncol(linp10)>=2) & (nrow(linp10)>=2) ){
    cat("Creating plots\n")

    if(!(dir.exists(outputDir))){
      cat("Directory for the plots not found, creating one.")
      dir.create(outputDir)
    }

    pdf(file = paste0(outputDir, name, "_Corr_of_SVs_and_covs.pdf"))
    gplots::heatmap.2(linp10, Colv = F, Rowv = F, dendrogram = "none",
                      trace = "none", symbreaks = F, symkey = F,
                      breaks = seq(-20, 0, length.out = 100),
                      key = T, colsep = NULL, rowsep = NULL, sepcolor = "black",
                      sepwidth = c(0.05, 0.05), main = "Corr. of SVs and covs.\n Combat-FPKM QN",
                      labCol = colnames(linp10),
                      labRow = colnames(numericCovs),
                      xlab = "Surrogate variables",
                      cexRow=0.5,
                      margins=c(5, 10))
    dev.off()
    pdf(file = paste0(outputDir, name, "_residuals_corr.pdf"))
    pcres = swamp::prince(as.matrix(CoExpNets::trasposeDataFrame(resids,F)),covs.rs,top=20)
    CoExpNets::princePlot(prince=pcres, margins=c(5, 10),
                          main=paste0(name,": residual cor.\n with ", paste0(bioCovsToCorrect, collapse=", ")))
    dev.off()
    cat("Plots done\n")

  } else if ((plots == T) & (svas$n.sv == 0)){

    if(!(dir.exists(outputDir))){
      cat("Directory for the plots not found, creating one.")
      dir.create(outputDir)
    }

    pdf(file = paste0(outputDir, name, "_residuals_corr.pdf"))
    pcres = swamp::prince(as.matrix(CoExpNets::trasposeDataFrame(resids,F)),covs.rs,top=20)
    CoExpNets::princePlot(prince=pcres, margins=c(5, 10), cexRow=0.5,
                          main=paste0(name,": \n residual cor.\n with ", paste0(bioCovsToCorrect, collapse=", ")))
    dev.off()
  }
  if(save==T){
    cat("Saving residuals\n")

    if(!(dir.exists(residsOutput))){
      cat("Directory for the residuals not found, creating one.")
      dir.create(residsOutput)
    }
    rownames(resids) <- preparedCovs[, idCov]
    saveRDS(resids,paste0(residsOutput, name, "_resids", ".rds"))
    cat(paste0("Residuals saved at: ", residsOutput, "\n"))
  }
  cat("Done\n")
  return(resids)
}




#' myGetProfilerOnNet - It returns the results of a functional enrichment analysis applied on each module of the selected network.
#' It is based on CoExpNets function getGprofilerOnNet, but this one is based on gprofiler2 package.
#'
#' @param net.file the pathway of the network selected for the functional enrichment analysis
#' @param filter the databases we want to use for functional enrichment analysis (default databases are GO, KEGG and REAC)
#' @param ensembl if ensembl=TRUE, it means that the names of the genes are the ensembl ones
#' @param exclude.iea if exclude.iea=TRUE, we remove automatically annotations (not manually curated)
#' @param organism the specie of the organism
#' @param correction.method the correction method used for functional enrichment analysis (deafult method is gSCS)
#' @param out.file the name of the file where functional enrichment analysis results will be saved
#' @return it returns the data frame containing the results of the functional enrichment analysis, where each row represents
#' one annotation associated with a specific module
#' @export
#' @examples

myGetGProfilerOnNet <- function(net.file,
                                filter=c("GO","KEGG","REAC"),
                                ensembl=FALSE,
                                exclude.iea=T,
                                organism = "hsapiens",
                                correction.method="gSCS",
                                out.file=NULL) {

  # require(gprofiler2)
  # require(CoExpNets)

  if(typeof(net.file) == "character")
    net <- readRDS(net.file)
  else
    net = net.file

  modules = unique(net$moduleColors)

  background <- names(net$moduleColors)
  if(ensembl)
    background <- CoExpNets::fromEnsembl2GeneName(background)

  all.genes <- NULL
  for(module in modules){

    genes <- names(net$moduleColors)[net$moduleColors == module]
    if(ensembl)
      genes <- fromEnsembl2GeneName(genes)
    all.genes[[module]] <- genes
  }

  go <- gprofiler2::gost(all.genes,
                         correction_method=correction.method,
                         #custom_bg=background,
                         sources=filter,
                         organism=organism,
                         exclude_iea=exclude.iea)

  go <- as.data.frame(go$result)

  if(nrow(go) == 0)
    return(go)
  #png_fn = paste0(out.file,".png"),
  #no_isects=T)
  go.out = cbind(go,rep(NA,nrow(go)))
  colnames(go.out) = c(colnames(go),"IC")

  go.out$parents <- gsub("\\(|\\)", "", go.out$parents)
  go.out$parents <- gsub("c", "", go.out$parents)
  go.out$parents <- gsub("\"", "", go.out$parents)
  go.out$parents <- as.character(go.out$parents)

  #	#Now we will add our own column with the information content of each term
  #	#So it is possible to evaluate them
  cat("Generating IC-BP")
  loadICsGOSim("GO:BP")
  mask = go$source == "GO:BP"
  bp.terms.go = go$term_id[mask]
  bp.ic.go = IC[bp.terms.go]
  go.out$IC[mask] = bp.ic.go

  cat("Generating IC-MF")
  loadICsGOSim("GO:MF")
  mask = go$source == "GO:MF"
  mf.terms.go = go$term_id[mask]
  mf.ic.go = IC[mf.terms.go]
  go.out$IC[mask] = mf.ic.go

  cat("Generating IC-CC")
  loadICsGOSim("GO:CC")
  mask = go$source == "GO:CC"
  cc.terms.go = go$term_id[mask]
  cc.ic.go = IC[cc.terms.go]
  go.out$IC[mask] = cc.ic.go
  #
  #	data("ICsMFhumanall")
  #	data("ICsCChumanall")

  if(!is.null(out.file)){
    write.csv(go.out,out.file)
    cat("Obtained",nrow(go),"terms saved at",out.file,"\n")
  }else{
    cat("Obtained",nrow(go),"terms\n")
  }
  return(go.out)
}

loadICsGOSim = function(onto){
  stopifnot(onto %in% c("GO:BP","GO:MF","GO:CC"))
  if(onto == "GO:BP")
    load(system.file("", "ICsBPhumanall.rda", package = "CoExpNets"),.GlobalEnv)
  else if(onto == "GO:MF")
    load(system.file("", "ICsMFhumanall.rda", package = "CoExpNets"),.GlobalEnv)
  else
    load(system.file("", "ICsCChumanall.rda", package = "CoExpNets"),.GlobalEnv)
}

getICReport = function(gprof){
  ontologies = c("GO:BP","GO:CC","GO:MF")
  ics = NULL
  for(onto in ontologies){
    loadICsGOSim(onto)
    terms = gprof$term_id[gprof$source == onto]
    ics[[onto]] = IC[terms]
    names(ics[[onto]]) = terms
  }
  return(ics)
}




#' phenoExam - It returns the results of a phenotype enrichment analysis applied on each module of the selected network.
#'
#' @param net.file the pathway of the network selected for the functional enrichment analysis
#' @return it returns the data frame containing the results of the phenotype enrichment analysis, where each row represents
#' one annotation associated with a specific module
#' @export
#' @examples

phenoExam <- function(net.file) {

  require(PhenoExam)

  if(typeof(net.file) == "character")
    net <- readRDS(net.file)
  else
    net = net.file

  modules = unique(net$moduleColors)

  databases <- getdbnames()

  enrichResults <- as.data.frame(matrix(numeric(), nrow = 0, ncol = 9))

  for (module in modules) {
    myModuleGenes <- names(net$moduleColors[net$moduleColors==module])
    enrichModule <- PhenoExam::PhenoEnrichGenes(genes= myModuleGenes, database=databases, url=F)
    enrichModule <- enrichModule$alldata
    enrichModule <- enrichModule[enrichModule$adjust_pvalue<0.05, ]
    enrichModule <- as.data.frame(cbind(enrichModule, rep(module, nrow(enrichModule))))

    enrichResults <- rbind(enrichResults, enrichModule)
  }

  colnames(enrichResults)[10] <- "module"

  # Removing not interesting terms
  index <- which(enrichResults$term_name %in% c("No diseases associated", "No CRB phenotype", "No HPO phenotype"))

  if(length(index)>0) {
    enrichResults <- enrichResults[-index, ]
  }

  return(enrichResults)

}




#' myGenAnnotationCellType - It returns the results of the cell type markers enrichment for each module of a selected network using a Fisher's exact test.
#' It is based on genAnnotationCellType function from CoExpNets.
#'
#' @param tissue A label to use to refer to the results in files and downstream results
#' @param which.one if we would like to use a CoExpNets network already created, which one we would like to use
#' @param markerspath the path where the markers files are located
#' @param net.in the path or the variable name that contains the selected network
#' @param legend
#' @param doheatmap
#' @param notHuman
#' @param plot.file
#' @param threshold
#' @param return.processed
#' @param getMarkerSize
#' @param getOnlyOverlap
#' @param getMyOwnTest
#' @param applyFDR
#' @param display.cats
#' @return it returns a data frame with one column per module, one row per cell type markers list and each cell containing
#' the result of a Fisher's exact test.
#' @export
#' @examples

myGenAnnotationCellType = function(tissue="None",
                                   which.one="new",
                                   markerspath=system.file("ctall", "",
                                                           package = "CoExpNets"),
                                   net.in=NULL,
                                   legend=NULL,
                                   doheatmap=T,
                                   notHuman=F,
                                   plot.file=NULL,
                                   threshold=20,
                                   return.processed=F,
                                   getMarkerSize=F,
                                   getOnlyOverlap=F,
                                   getMyOwnTest=F,
                                   applyFDR=F,
                                   display.cats=NULL){ #This last flag FALSE when you want the raw p-values

  # require(gplots)

  cat("Entering genAnnotationCellType with",which.one,"\n")

  if(which.one == "new"){
    if(typeof(net.in) == "character")
      net = readRDS(net.in)
    else
      net = net.in
  }else{
    net = getNetworkFromTissue(tissue=tissue,which.one=which.one)
  }
  modules = unique(net$moduleColors)
  if(notHuman)
    names(net$moduleColors) = toupper(names(net$moduleColors))
  #So the 1st heatmap
  #will have a column for each cell.types element and a row for
  #each module and we will show -log10(p-values) in a scale

  files = list.files(path=markerspath,full.names = T)
  files = files[grep(pattern=".txt$",files)]
  if(getMarkerSize){
    sizes = lapply(files,function(x){
      nrow(read.delim(x,header=T))
    })
    return(cbind(markerset=files,genes=sizes))
  }


  markernames = gsub(".txt","",basename(files))


  ctypes = markernames
  ctypedata = matrix(ncol=length(modules),nrow=length(files))
  ctypedata[,] = 1
  rownames(ctypedata) = ctypes
  colnames(ctypedata) = modules

  if(getOnlyOverlap){
    myfunction = Vectorize(function(m,f){
      submarkers = read.delim(f,stringsAsFactors=F,header=T)[,1]
      #cat("Module is",m,"file is",f,"\n")
      #print(names(net$moduleColors)[net$moduleColors == m])
      #print(submarkers)
      return(sum(submarkers %in% names(net$moduleColors)[net$moduleColors == m]))
    })

    overlap = outer(modules,files,myfunction)
    colnames(overlap) = gsub(".txt","",basename(files))
    rownames(overlap) = modules
    if(applyFDF){
      overlap = p.adjust(overlap,method="fdr")
    }
    return(overlap)

  }

  if(getMyOwnTest){
    myfunction = Vectorize(function(m,f){
      submarkers = read.delim(f,stringsAsFactors=F,header=T)[,1]
      ov = sum(submarkers %in% names(net$moduleColors)[net$moduleColors == m])
      nm = sum(net$moduleColors == m)

      #cat("Module is",m,"file is",f,"\n")
      #print(names(net$moduleColors)[net$moduleColors == m])
      #print(submarkers)
      return(testGeneSet(n.module=nm,total.specific=length(submarkers),
                         total.net=length(net$moduleColors),n.module.and.specific=ov)$p.value)
    })

    overlap = outer(modules,files,myfunction)
    colnames(overlap) = gsub(".txt","",basename(files))
    rownames(overlap) = modules
    return(overlap)

  }



  all.gene.names = names(net$moduleColors)
  all.gene.names = CoExpNets::fromAny2GeneName(names(net$moduleColors))

  tmp.enr.f = paste0(markerspath,"enrichment_tmp.csv")
  if(file.exists(tmp.enr.f))
    file.remove(tmp.enr.f)
  ukbec.en = WGCNA::userListEnrichment(all.gene.names,net$moduleColors,files,
                                       nameOut=tmp.enr.f)
  if(file.exists(tmp.enr.f)){
    enrichment = read.csv(tmp.enr.f,
                          stringsAsFactors=F)

    for(i in 1:nrow(enrichment)){
      module = enrichment$InputCategories[i]
      for(j in 1:length(files)){
        if(grepl(ctypes[j],enrichment$UserDefinedCategories[i]))
          break
      }
      category = ctypes[j]
      ctypedata[category,module] = enrichment$CorrectedPvalues[i]
    }
  }
  rownames(ctypedata) = gsub("cell_type_TableS1.csv.","",rownames(ctypedata))
  uctypedata = ctypedata
  ctypedata = ctypedata[,apply(ctypedata,2,function(x){ any(x < 1)}),drop=FALSE]
  ctypedata = -log10(ctypedata)
  ctypedata[is.infinite(ctypedata)] = max(ctypedata[!is.infinite(ctypedata)])

  if(threshold > 0)
    ctypedata[ctypedata > threshold] = threshold
  #Order by cell type

  ctypedata = ctypedata[order(rownames(ctypedata)),,drop=FALSE]

  ctypedata = ctypedata[apply(ctypedata,1,function(x){ any(x > 2)}),]


  if(doheatmap){
    ctypedata = as.data.frame(ctypedata)
    if((ncol(ctypedata) >= 2) && nrow(ctypedata) >=2){
      ctypedata = as.matrix(ctypedata)
      if(is.null(legend))
        legend = tissue
      if(!is.null(plot.file))
        pdf(plot.file,width=15,height=8)
        gplots::heatmap.2(ctypedata,
                  trace="none",
                  col=heat.colors(100)[100:1],
                  cexCol=1.0,
                  cexRow=1.0,
                  Rowv=F,
                  Colv=T,
                  main=paste0("Cell type enrichment for ",legend),
                  key=F,srtCol=45,dendrogram="none",
                  margins=c(6,20))
      y= ctypedata
      if(!is.null(plot.file))
        dev.off()
    }
  }

  if(return.processed)
    return(ctypedata)
  return(uctypedata)
}




#' myGetDownstreamNetwork - It creates and annotates a single-cell gene co-expression network
#' It is based on getDownstreamNetwork function from CoExpNets.
#'
#' @param tissue A label to use to refer to the results in files and downstream results
#' @param n.iterations Number of iterations in the k-means algorithm to refine the
#' clustering process
#' @param expr.data A data frame with cell type cluster expression data obtained from svaCorrection function (genes at cols and cells at rows)
#' @param beta The smoothing parameter of the WGCNA pipeline. If -1, then the method
#' will suggest one automatically by looking at the R2 between gene connectivity and
#' Scale Free Topology features
#' @param job.path This method will generate a number of files. It needs a folder to write
#' results
#' @param min.cluster.size Minimum number of genes for a group to be considered as cluster
#' for the tree cutting algorithm to convert from a dendogram to a cluster.
#' @param net.type Whether a signed ("signed") or unsigned ("unsigned") network type will be created
#' @param debug Set this to true if you want a quick run of the method to test how it works with
#' a small amount of your genes
#' @param excludeGrey If WGCNA detects grey genes, set it to TRUE if you want them removed
#' from the network before applying k-means
#' @param fullAnnotation Set to TRUE if you want to annotate by cell type and GO, REACTOME and KEGG
#' @param fastKMeans Set to TRUE for a faster k-means use
#'
#' @return A file name that can be used to access your network
#' @return it returns a data frame with one column per module, one row per cell type markers list and each cell containing
#' the result of a Fisher's exacta test.
#' @export
#' @examples

myGetDownstreamNetwork = function(tissue="mytissue",
                                  n.iterations=50,	#Number of iterations for k-means, 50 recommended
                                  min.exchanged.genes=20,
                                  expr.data, 	#We expect a file name pointing to a dataframe (RDS format) with
                                  #genes in columns and samples in rows. Each gene name appears
                                  #in the column name. Better to use gene symbols as names
                                  beta=-1,	#If -1 the algorithm will seek for the best beta
                                  job.path="~/tmp/",	#Where to store all results
                                  min.cluster.size=30,		#Minimum number of genes to form a cluster
                                  net.type="signed",			#Leave it like that (see WGCNA docs)
                                  debug=F,
                                  blockTOM=F,
                                  save.tom=F,
                                  save.plots=F,
                                  excludeGrey=FALSE,
                                  fullAnnotation=T,
                                  silent=T){

  final.net=NULL
  distance.type="cor"
  centroid.type="pca"
  cor.type="pearson"

  if(debug){
    if(typeof(expr.data) == "character")
      expr.data = readRDS(expr.data)
    expr.data = expr.data[,1:1500]
    n.iterations=5
  }

  validgenes = unlist(apply(expr.data,2,function(x){
    return(var(x) != 0)
  }))
  if(sum(validgenes) < ncol(expr.data))
    cat("There are genes with 0 variance:",
        paste0(colnames(expr.data)[!validgenes],collapse=","),"\n")
  discGenes = colnames(expr.data)[!validgenes]
  expr.data = expr.data[,validgenes]

  net.and.tom = getAndPlotNetworkLong(expr.data=expr.data,
                                      beta=beta,
                                      tissue.name=tissue,
                                      min.cluster.size=min.cluster.size,
                                      save.plots=save.plots,
                                      excludeGrey=excludeGrey,
                                      additional.prefix=job.path,
                                      return.tom=T,
                                      cor.type=cor.type,
                                      silent=silent)

  if(is.null(net.and.tom))
    return(net.and.tom)

  if(is.null(final.net) & !is.null(job.path))
    final.net = paste0(job.path,"/","net",tissue,".",
                       net.and.tom$net$beta,".it.",n.iterations,".rds")
  if(is.null(job.path))
    final.net = NULL

  outnet = applyKMeans(tissue=tissue,
                       n.iterations=n.iterations,
                       net.file=net.and.tom$net,
                       expr.data=expr.data,
                       excludeGrey=excludeGrey,
                       min.exchanged.genes = min.exchanged.genes,
                       silent=silent)


  if(save.tom & !is.null(final.net)){
    if(blockTOM)
      saveTOM(tom=net.and.tom$tom,
              clusters=outnet$moduleColors,
              filepref=paste0(final.net,".tom."))
    else
      saveRDS(net.and.tom$tom,paste0(final.net,".tom.rds"))
  }

  outnet$beta = net.and.tom$net$beta
  outnet$file = final.net
  outnet$adjacency = net.and.tom$net$adjacency
  names(outnet$moduleColors) = colnames(expr.data)
  outnet$discGenes = discGenes

  # names(outnet$moduleColors) = unlist(lapply(str_split(names(outnet$moduleColors),
  #                                                     "\\."),
  #                                           function(x){return(x[[1]])}))

  if(!is.null(final.net)){
    saveRDS(outnet,final.net)
    if(save.tom){
      if(blockTOM)
        outnet$tom = paste0(final.net,".tom.")
      else
        outnet$tom = paste0(final.net,".tom.rds")
    }
    if(save.plots){
      cat("Generating mod sizes for",final.net,"\n")
      pdf(paste0(final.net,".mod_size.pdf"))
      plotModSizes(which.one="new",tissue=final.net)
      dev.off()
      pdf(paste0(final.net,".Eigengenes_clustering.pdf"))
      plotEGClustering(which.one="new",tissue=final.net)
      dev.off()
    }
    if(fullAnnotation){
      go = myGetGProfilerOnNet(net.file=final.net,
                               out.file=paste0(final.net,"_gprof.csv"))

      write.csv(myGenAnnotationCellType(net.in=final.net,
                                        return.processed = F,
                                        doheatmap=save.plots,
                                        tissue=tissue,
                                        plot.file=paste0(final.net, "_celltype.pdf")),
                paste0(final.net,"_celltype.csv"))

      write.table(CoExpNets::getMM(final.net,expr.data,genes=NULL),paste0(final.net, "_getMM.csv"),quote=F,row.names=F, sep="\t")
      write.csv(phenoExam(final.net), paste0(final.net, "_PEG.csv"))
    }


    return(final.net)
  }
  outnet
}




#' getNet - It applies the whole pipeline to create and annotate a cell-type-specific sginel-cell gene co-expression network

#' @param exprData a normalized cell type cluster gene expression matrix
#' @param covs a cell type cluster metadata table
#' @param batchCov the variable name that represents the batch if exists
#' @param idCov the variable name that represents the individuals ID
#' @param bioCovsToCorrect the biological or technical covariates that we are not interested in and we want to remove their effect on our data
#' @param name the name of the cluster
#' @param filteringGenes if filteringGenes=TRUE, we want to apply gene filtering based on a minimum expression value or cutoff on a percentage of cells
#' @param filteringCutoff if filteringGenes=TRUE, filteringCutoff represents the minimum expression value
#' @param filteringProportion if filteringGenes=TRUE, filteringProportion represents the percentage of cells
#' @param plots if plots=TRUE, plots will be saved
#' @param path the path where we want that files are saved
#' @return it returns a new folder containing theses files: the residuals, the network, the functional enrichment analysis results,
#' the phenotype enrichment analysis results, the cell type markers enrichment and module membership.
#' @export
#' @examples

getNet <- function(exprData,
                   covs,
                   batchCov=NULL,
                   idCov,
                   bioCovsToCorrect=NULL,
                   name=deparse(substitute(exprData)),
                   filteringGenes=F,
                   filteringCutoff=NULL,
                   filteringProportion=NULL,
                   plots=F,
                   path=getwd()) {

  # Check arguments
  if (missing(exprData)) { stop("Please provide an expression matrix\n")}
  if (missing(covs)) { stop("Please provide a a table of covariates\n")}
  if (missing(idCov)) { stop("Please provide the name of the variable that represents the id of the samples\n")}
  if (missing(bioCovsToCorrect)) { stop("Please provide the name of the biological variables for which you want to correct the expression data\n")}

  # Load libraries
  require(stringr)
  require(gplots)
  require(fs)
  require(dplyr)
  require(preprocessCore)
  require(WGCNA)
  require(limma)
  require(biomaRt)
  require(CoExpNets)
  require(sva)
  require(swamp)
  require(gprofiler2)
  require(readr)

  name = gsub("_exprData", "", name)

  # Create directories
  if(!(dir.exists(paste0(path, "/", name)))) {
    cat("Creating", paste0(path, "/", name), "directory\n")
    dir.create(paste0(path, "/", name))
    dir.create(paste0(path, "/", name, "/Net/"))
    dir.create(paste0(path, "/", name, "/Residuals/"))

    if (plots==T) {
      dir.create(paste0(path, "/", name, "/Plots/"))
      dir.create(paste0(path, "/", name, "/Plots/Plots_MDS_Batch/"))
      dir.create(paste0(path, "/", name, "/Plots/Plots_PCA_Batch/"))
      dir.create(paste0(path, "/", name, "/Plots/Plots_SVAs/"))
      dir.create(paste0(path, "/", name, "/Plots/Plots_Net/"))
    }
  }

  # Preparing data
  preparedData <- prepareData(exprData, covs, batchCov, idCov, filteringGenes, filteringCutoff, filteringProportion)
  preparedExprData <- preparedData[["normalized_expression"]]
  preparedCovs <- preparedData[["covariates"]]

  # MDS plots
  if (plots==T && !is.null(batchCov)) {
    plotMDS_Batch(preparedExprData, preparedCovs, batchCov,
                  save = T,
                  outputDir = paste0(path, "/", name, "/Plots/Plots_MDS_Batch/"),
                  name=name)
  }

  if (!is.null(batchCov)) {
    # Combat correction
    combatExprData = combatCorrect(preparedExprData,
                                   preparedCovs,
                                   batchCov,
                                   bioCovsToCorrect,
                                   plots,
                                   name,
                                   paste0(path, "/", name, "/Plots/Plots_PCA_Batch/"))
  } else {
    combatExprData = preparedExprData
  }

  if (!is.null(bioCovsToCorrect)) {
    # SVA correction
    resids = svaCorrection(combatExprData = combatExprData,
                           preparedCovs = preparedCovs,
                           bioCovsToCorrect = bioCovsToCorrect,
                           plots = plots,
                           idCov = idCov,
                           batchCov = batchCov,
                           name = name,
                           save = T,
                           outputDir = paste0(path, "/", name, "/Plots/Plots_SVAs/"),
                           residsOutput = paste0(path, "/", name, "/Residuals/"))

  } else {

    # Without SVA
    resids <- t(combatExprData)
    saveRDS(resids, paste0(path, "/", name, "/Residuals/", name, "_resids.rds"))
  }

  # Network creation and full annotation
  net = myGetDownstreamNetwork(tissue=name,
                               n.iterations=50,
                               net.type = "signed",
                               debug=F,
                               expr.data=resids,
                               fullAnnotation=T,
                               save.plots=plots,
                               job.path=paste0(path, "/", name, "/Net"))

  cat("Done.\n")

  # Correlation between modules and covariates
  path=paste0(path, "/", name, "/")

  # Moving files to their corresponding directories
  if (plots==T) {

    for (files in list.files(path=path, pattern = ".pdf")) {
      fs::file_move(paste0(path, files), paste0(path, "Plots/Plots_Net/"))
    }

    if (!(is.null(net)) ) {

      for (files in list.files(path=paste0(path, "Net/"), pattern = ".pdf")) {
        fs::file_move(paste0(path, "Net/", files), paste0(path, "Plots/Plots_Net/"))
      }
    }
  }

}

