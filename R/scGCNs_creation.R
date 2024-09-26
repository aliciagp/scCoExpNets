#' pair - It combines single-cell metadata and corresponding individual covariates in the same data frame
#' It is an auxiliary function for pseudoCells function
#'
#' @param indexes it represents a range of numbers from 1 to the number of columns of the gene expression matrix
#' @param pair it represents the ID of the cells. To group the cells per individual, it is necessary that individual ID is included in cell ID
#' @return a data frame containing pairs of cells grouped per individual ID. The data frame contains three columns: cell1, cell2, sampleID
#' @export

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



#' getModules - It creates and annotates a single-cell gene co-expression network
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

getModules = function(tissue="mytissue",
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
                      fullAnnotation=F,
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

    return(list(final.net=final.net, name=paste0(job.path,"/","net",tissue,".",
                                            net.and.tom$net$beta,".it.",n.iterations,".rds")))
  }
  outnet
}




#' get_Ti_GCN - It creates and annotates a Ti GCN of a specific cell type cluster. It is an auxiliary function of get_all_GCNs function.

#' @param exprData a normalized cell type cluster gene expression matrix
#' @param metadata a cell type cluster metadata table
#' @param gcnName the name of the network
#' @param idColname the name of the column that contains donors' ID
#' @param batchColname the name of the column that contains batch effect information
#' @param covsToCorrect the biological or technical covariates names that we are not interested in and we want to remove their effect on our data
#' @param filteringGenes if filteringGenes=TRUE, we want to apply gene filtering based on a minimum expression value or cutoff on a percentage of cells
#' @param filteringCutoff if filteringGenes=TRUE, filteringCutoff represents the minimum expression value
#' @param filteringProportion if filteringGenes=TRUE, filteringProportion represents the percentage of cells
#' @param plots if plots=TRUE, plots will be saved
#' @param path the path where the files will be saved
#' @param moduleTraitCorr, if moduleTraitCorr=T, the asociation between module eigengenes and traits will be estimated
#' @param covsToCorr=NULL, if moduleTraitCorr=T, covsToCorr includes the traits to be evaluated
#' @param phenotypeEnrich=T, if phenotypeEnrich=T, phenotype enrichment analysis will be carried out for each module of the network
#' @param correction.method the correction method used for functional enrichment analysis (deafult method is gSCS)
#' @param sources the databases we want to use for functional enrichment analysis (default databases are GO, KEGG and REAC)
#' @param organism the specie of the organism
#' @param exclude.iea if exclude.iea=TRUE, we remove automatically annotations (not manually curated)
#' @return it returns a new folder containing theses files: the residuals, the network, the functional enrichment analysis results,
#' the phenotype enrichment analysis results, the cell type markers enrichment and module membership.
#' @export

get_Ti_GCN <- function(exprData,
                       metadata,
                       gcnName="GCN1",
                       idColname=NULL,
                       batchColname=NULL,
                       covsToCorrect=NULL,
                       filteringGenes=F,
                       filteringCutoff=NULL,
                       filteringProportion=NULL,
                       plots=F,
                       path=getwd(),
                       moduleTraitCorr=F,
                       covsToCorr=NULL,
                       phenotypeEnrich=T,
                       correction_method="gSCS",
                       sources=c("GO", "KEGG", "REAC"),
                       organism="hsapiens",
                       exclude_iea=T) {

  # Check arguments
  if (missing(exprData)) { stop("Please provide an expression matrix\n")}
  if (missing(metadata)) { stop("Please provide a a table of covariates\n")}
  if (missing(idColname)) { stop("Please provide the name of the variable that represents the id of the samples\n")}

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

  # Create directories
  if(!(dir.exists(paste0(path, "/", gcnName)))) {
    cat("Creating", paste0(path, "/", gcnName), "directory\n")
    dir.create(paste0(path, "/", gcnName))
    dir.create(paste0(path, "/", gcnName, "/Net/"))
    dir.create(paste0(path, "/", gcnName, "/Residuals/"))

    if (plots==T) {
      dir.create(paste0(path, "/", gcnName, "/Plots/"))
      dir.create(paste0(path, "/", gcnName, "/Plots/Plots_MDS_Batch/"))
      dir.create(paste0(path, "/", gcnName, "/Plots/Plots_PCA_Batch/"))
      dir.create(paste0(path, "/", gcnName, "/Plots/Plots_SVAs/"))
      dir.create(paste0(path, "/", gcnName, "/Plots/Plots_Net/"))
    }
  }

  if(is.character(exprData)) {
    exprData <- readRDS(exprData)
  }

  if(is.character(metadata)) {
    metadata <- readRDS(metadata)
  }

  # Preparing data
  preparedData <- prepareData(exprData=exprData,
                              covs=metadata,
                              batchCov=batchColname,
                              idCov=idColname,
                              filteringGenes=filteringGenes,
                              filteringCutoff=filteringCutoff,
                              filteringProportion=filteringProportion)

  preparedExprData <- preparedData[["normalized_expression"]]
  preparedCovs <- preparedData[["covariates"]]

  # MDS plots
  if (plots==T && !is.null(batchColname)) {
    plotMDS_Batch(preparedExprData, preparedCovs, batchColname,
                  save = T,
                  outputDir = paste0(path, "/", name, "/Plots/Plots_MDS_Batch/"),
                  name=name)
  }

  if (!is.null(batchColname)) {
    # Combat correction
    combatExprData = combatCorrect(preparedExprData,
                                   preparedCovs,
                                   batchColname,
                                   covsToCorrect,
                                   plots,
                                   gcnName,
                                   paste0(path, "/", name, "/Plots/Plots_PCA_Batch/"))
  } else {
    combatExprData = preparedExprData
  }

  if (!is.null(covsToCorrect)) {
    # SVA correction
    resids = svaCorrection(combatExprData = combatExprData,
                           preparedCovs = preparedCovs,
                           bioCovsToCorrect = covsToCorrect,
                           plots = plots,
                           idCov = idColname,
                           batchCov = batchColname,
                           name = gcnName,
                           save = T,
                           outputDir = paste0(path, "/", gcnName, "/Plots/Plots_SVAs/"),
                           residsOutput = paste0(path, "/", gcnName, "/Residuals/"))

  } else {

    # Without SVA
    resids <- t(combatExprData)
    saveRDS(resids, paste0(path, "/", gcnName, "/Residuals/", gcnName, "_resids.rds"))
  }

  rownames <- ave(rownames(resids), rownames(resids), FUN = function(i) paste0(i, '_', seq_along(i)))
  stopifnot(identical(gsub("_.*", "", rownames), rownames(resids)))
  rownames(resids) <- rownames

  # Network creation
  net = getModules(tissue=gcnName,
                   n.iterations=50,
                   net.type = "signed",
                   debug=F,
                   expr.data=resids,
                   fullAnnotation=T,
                   save.plots=plots,
                   job.path=paste0(path, "/", gcnName, "/Net"))

  final.net <- net$name

  # Modules annotation
  cat("Get functional enrichment\n")
  write.csv(getFunctionalEnrichment(net=final.net,
                                    correction_method=correction_method,
                                    sources=sources,
                                    organism=organism,
                                    exclude_iea=exclude_iea),
            paste0(final.net, "_functionalEnrich.csv"))

  cat("Get module membership for each gene\n")
  write.csv(CoExpNets::getMM(final.net,resids,genes=NULL),
            paste0(final.net, "_moduleMembership.csv"),quote=F,row.names=F, sep="\t")

  if(moduleTraitCorr==T & !is.null(covsToCorr)) {
    cat("Get module trait corr\n")
    write.csv(getTraitCorr(net=final.net, metadata=preparedCovs, covsToCorr=covsToCorr),
              paste0(final.net, "_traitCorr.csv"))
  }

  if(phenotypeEnrich==T) {
    cat("Get phenotype enrichment\n")
    write.csv(getPhenotypeTerms(net=final.net),
              paste0(final.net, "_phenotypeEnrich.csv"))
  }

  cat("Done.\n")

  path=paste0(path, "/", gcnName, "/")

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


#' get_all_GCNs - It creates and annotates a set of GCNs for a specific cell type cluster
#' @param exprData a normalized cell type cluster gene expression matrix
#' @param metadata a cell type cluster metadata table
#' @param gcnName the name of the network
#' @param idColname the name of the column that contains donors' ID
#' @param batchColname the name of the column that contains batch effect information
#' @param minNumCells the minimum number of cells to create a GCN
#' @param covsToCorrect the biological or technical covariates names that we are not interested in and we want to remove their effect on our data
#' @param filteringGenes if filteringGenes=TRUE, we want to apply gene filtering based on a minimum expression value or cutoff on a percentage of cells
#' @param filteringCutoff if filteringGenes=TRUE, filteringCutoff represents the minimum expression value
#' @param filteringProportion if filteringGenes=TRUE, filteringProportion represents the percentage of cells
#' @param plots if plots=TRUE, plots will be saved
#' @param path the path where the files will be saved
#' @param moduleTraitCorr, if moduleTraitCorr=T, the asociation between module eigengenes and traits will be estimated
#' @param covsToCorr=NULL, if moduleTraitCorr=T, covsToCorr includes the traits to be evaluated
#' @param phenotypeEnrich=T, if phenotypeEnrich=T, phenotype enrichment analysis will be carried out for each module of the network
#' @param correction.method the correction method used for functional enrichment analysis (deafult method is gSCS)
#' @param sources the databases we want to use for functional enrichment analysis (default databases are GO, KEGG and REAC)
#' @param organism the specie of the organism
#' @param exclude.iea if exclude.iea=TRUE, we remove automatically annotations (not manually curated)
#' @return it returns a new folder containing theses files: the residuals, the network, the functional enrichment analysis results,
#' the phenotype enrichment analysis results, the cell type markers enrichment and module membership.
#' @export
#'
get_all_GCNs <- function(exprData,
                         metadata,
                         gcnName="GCN1",
                         idColname=NULL,
                         batchColname=NULL,
                         minNumCells=200,
                         covsToCorrect=NULL,
                         filteringGenes=F,
                         filteringCutoff=NULL,
                         filteringProportion=NULL,
                         plots=F,
                         path=getwd(),
                         moduleTraitCorr=F,
                         covsToCorr=NULL,
                         phenotypeEnrich=T,
                         correction_method="gSCS",
                         sources=c("GO", "KEGG", "REAC"),
                         organism="hsapiens",
                         exclude_iea=T) {

  if(is.character(exprData)) {
    exprData <- readRDS(exprData)
  }

  if(is.character(metadata)) {
    metadata <- readRDS(metadata)
  }

  iter <- 0
  newExprData <- exprData
  newMetadata <- metadata

  while(ncol(newExprData) > minNumCells) {
    newName <- paste0(gcnName, "_iter", iter)

    if(iter>0) {
      myPseudoCells <- pseudoCells(newExprData, newMetadata)
      newExprData <- myPseudoCells[[1]]
      newMetaData <- myPseudoCells[[2]]
    }

    gcn <- get_Ti_GCN(exprData=newExprData,
                      metadata=newMetadata,
                      gcnName=newName,
                      idColname=idColname,
                      batchColname=batchColname,
                      covsToCorrect=covsToCorrect,
                      filteringGenes=filteringGenes,
                      filteringCutoff=filteringCutoff,
                      filteringProportion=filteringProportion,
                      plots=plots,
                      path=path,
                      moduleTraitCorr=moduleTraitCorr,
                      covsToCorr=covsToCorr,
                      phenotypeEnrich=phenotypeEnrich,
                      correction_method=correction_method,
                      sources=sources,
                      organism=organism,
                      exclude_iea=exclude_iea)

    iter <- iter+1
  }

}




#' getPreservedModules  - Check whether one network module is preserved in a tissue gene expression profile
#'
#' @param network The network you want to assess about preservation (full path and file name or a network object).
#' @param expr.data.files A list with the expression data used to create the network (1st element) and the
#' expression data onto which we want the network assessed of preservation (2nd element)
#' @param tissues A vector of strings, to name the tissues
#' @param permutations The number of random permutations used to construct the null hypothesis. 200 may be a
#' reasonable number
#' @param maxModuleSize All modules over this size in genes will be considered as having as much size as this limit
#' @param maxGoldModuleSize The size of the random module to be used as a reference of a totally random module values
#' @param randomSeed Seed for the random number generator
#'
#' @return A data frame with modules in the rows and statistics at the columns.
#' @export

getPreservedModules <- function (network, expr.data.files = NULL, tissues = c("snig","putm"),
                                 permutations = 200, maxModuleSize = 5000, maxGoldModuleSize = 400,
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
