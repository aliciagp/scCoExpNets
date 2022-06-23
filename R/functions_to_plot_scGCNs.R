#' getOVerlap - It returns the number of genes that remain together at the end after successive pseudo-cells iterations
#' using all possible combinations of modules from different iterations
#'
#' @param nets_dir the path where we can find all the networks belonging to the same cell type
#' @param results_dir the path where we can find the results of this function stored
#' @return a data frame containing the number of genes that remain together at the end after successive pseudo-cells iterations
#' using all possible combinations of modules from different iterations
#' @export
#' @examples

getOverlap <- function(nets_dir,
                       results_dir=getwd()) {

  # Check if the file has been already created
  nets <- list.files(nets_dir, full.names=T, recursive=T, pattern=".50.rds$")
  name <- str_split(str_split(nets[1], "/Net/")[[1]], "/")[[1]]
  name <- name[length(name)]

  file <- list.files(results.dir, full.names=T, pattern=name)

  if (length(index)>0) {
    dat <- readRDS(file)
  } else {

    # Step 1: create a list with the genes of each module
    modulesList <- list()
    iter <- 0

    for (net in nets) {
      net <- readRDS(net)
      modules <- unique(net$moduleColors)

      for (module in modules) {
        modulesList[[paste0("Iter",iter, "-", module)]]<- names(net$moduleColors)[net$moduleColors==module]
      }
      iter <- iter+1
    }



    # Step 2: create a list of all possible combinations
    names <- names(modulesList)
    namesList <- list()

    for (i in 0:(length(nets)-1)) {
      namesList[[i+1]] <- names[grep(paste0("Iter", i), names)]
    }

    df <- as.data.frame(expand.grid(namesList))
    colnames(df) <- paste0("Iter", 0:(length(nets)-1))


    # Step 3: check the genes overlap for each combination
    results <- data.frame()

    for (i in 1:nrow(df)) {
      comb <- as.vector(unlist(df[i, ]))
      poolOfGenes <- table(unlist(modulesList[which(names(modulesList) %in% comb)]))
      overlap <- length(poolOfGenes[poolOfGenes==length(nets)])

      if(overlap>0) {
        mydf <- as.data.frame(t(c(comb, overlap)))
        colnames(mydf) <- c(paste0("Iter", 0:(length(nets)-1)), "overlap")
        results <- rbind(results, mydf)
      }
    }


    # Step 4: create a data frame with two columns
    results <- as.data.frame(apply(results, 2, function(x) gsub("-", ": ", x)))
    path <- apply(results[, -ncol(results)], 1 , function(x) paste0(x, collapse="-"))
    dat <- as.data.frame(cbind(path, as.numeric(as.character(results$overlap))))
    dat$path <- gsub("Iter", "", dat$path)

    saveRDS(dat, paste0(results_dir, paste0(gsub("_", "", name), "_overlap.rds")))
  }
  return(dat)
}




#' getColors - It returns the colors to be used for the sunburst plot based on the correlation of each module with clinical diagnosis, age or sex
#' using all possible combinations of modules from different iterations
#'
#' @param nets_dir the path where we can find all the networks belonging to the same cell type
#' @param metadata the path or the data frame where we can find the metadata of all the cells belonging to this cell type
#' @param cellID the name of the variable that contains the ID of each cell
#' @param covariates the name of the biological covariates we are interested in
#' @param covariateSelected the name of the covariate selected to define the colors of the sunburst
#' @return a list containing both the colors assigned and the names of the modules
#' @export
#' @examples

getColors <- function(nets_dir,
                      metadata,
                      cellID="cell_id",
                      covariates=c("Clinical_Dx", "Donor_age_years", "Sex"),
                      covariateSelected="Clinical_Dx") {


  if(is.character(metadata)) {
    metadata <- readRDS(metadata)
  }

  mydf <- getModuleTraitCorr(nets_dir,
                             metadata,
                             cellID=cellID,
                             covariates=covariates)

  # Color sunburst based on modules correlation with selected covariate
  if(covariateSelected=="Clinical_Dx") {
    cov <- "dx_pval"
  } else if(covariateSelected=="Donor_age_years") {
    cov <- "age_pval"
  } else {
    cov <- "sex_pval"
  }

  mydf2 <- mydf[, c("iter", "module", cov)]
  modules <- paste0(mydf2$iter, ": ", mydf2$module)
  modules <- gsub("^iter", "Iter", modules)
  col <- rep("grey", length(modules))

  log <- -log10(as.numeric(as.character(mydf2[, cov])))

  # Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('darkgreen','red'))
  index <- which(log < -log10(0.05))

  # This adds a column of color values based on the y values
  cols <- rbPal(30)[as.numeric(cut(log, breaks = 30))]
  cols[index] <- "lightgrey"

  return(list(cols=cols, modules=modules))
}



#' getSunburst - It returns a sunburst plot where the inner ring represent the T0 modules (the size of each segment is proportional to the number of genes)
#'  and the rest of rings represent how the genes of each T0 module remains or not together after successive pseudo-cells iterations
#'
#' @param nets_dir the path where we can find all the networks belonging to the same cell type
#' @param metadata the path or the data frame where we can find the metadata of all the cells belonging to this cell type
#' @param cellID the name of the variable that contains the ID of each cell
#' @param covariates the name of the biological covariates we are interested in
#' @param covariateSelected the name of the covariate selected to define the colors of the sunburst
#' @param results_dir the path where we can find the results of this function stored
#' @return a list containing both the colors assigned and the names of the modules
#' @export
#' @examples

getSunburst <- function(nets_dir,
                        metadata,
                        cellID="cell_id",
                        covariates=c("Clinical_Dx", "Donor_age_years", "Sex"),
                        covariateSelected="Clinical_Dx",
                        results_dir=getwd()) {

  require(stringr)
  require(sunburstR)
  require(htmltools)
  require(d3r)

  dat <- getOverlap(nets_dir, results_dir)
  res <- getColors(nets_dir, metadata, cellID, covariates, covariateSelected)

  cols <- res[["cols"]]
  modules <-res[["modules"]]
  modules <- gsub("Iter", "", modules)
  sb3 <- sund2b(dat, width="100%", colors= list(range = c("lightblue", cols), domain = c("ROOT", modules)), showLabels=TRUE, rootLabel="ALL")

  # div(
  #   style="display: flex; align-items:center;",
  #   sb3
  # )

  return(sb3)

}



