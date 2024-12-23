
getIC <- function() {

  require(GOSim)
  setEvidenceLevel(evidences = "IC", organism="hsapiens")
  setOntology("CC", loadIC=FALSE) # important: setOntology assumes that the IC file already exists. To prevent an error message we need the second argument
  calcICs()
  rm(GOSimEnv)

  setEvidenceLevel(evidences = "IC", organism="hsapiens")
  setOntology("BP", loadIC=FALSE) # important: setOntology assumes that the IC file already exists. To prevent an error message we need the second argument
  calcICs()
  rm(GOSimEnv)

  setEvidenceLevel(evidences = "IC", organism="hsapiens")
  setOntology("MF", loadIC=FALSE) # important: setOntology assumes that the IC file already exists. To prevent an error message we need the second argument
  calcICs()
  rm(GOSimEnv)

}



getJaccardIndex <- function(x, y) {
  x <- str_split(x, ", ")[[1]]
  y <- str_split(y, ", ")[[1]]

  int <- length(intersect(x,y))
  j <- int/length(y)

  return(j)
}



addIC <- function(df) {

  # Add IC information to each subgraph
  load(file="C:/Users/alici/Downloads/paperSeb/aarmdv2/aaCBM/GO_subgraphs/ICsCChsapiensIC.rda")
  ic <- IC
  load(file="C:/Users/alici/Downloads/paperSeb/aarmdv2/aaCBM/GO_subgraphs/ICsBPhsapiensIC.rda")
  ic <- c(ic, IC)
  load(file="C:/Users/alici/Downloads/paperSeb/aarmdv2/aaCBM/GO_subgraphs/ICsMFhsapiensIC.rda")
  ic <- c(ic, IC)

  df$IC <- rep(NA, nrow(df))

  for(i in 1:nrow(df)) {
    ids <- str_split(as.vector(unlist(df[i, "IDs"])), ", ")[[1]]
    myic <- ic[match(ids, names(ic))]
    # stopifnot(identical(names(myic), ids))
    mymin <- round(min(myic, na.rm=T), 2)
    mymean <- round(mean(myic, na.rm=T), 2)
    mymax <- round(max(myic, na.rm=T), 2)
    myic <- paste0(mymean, " (", mymin, " - ", mymax, ")")
    df$IC[i] <- myic
  }

  return(df)
}



##############################################################################################################################################

getCombsT0vsTi <- function(df) {

  # Get all possible combinations between T0 subgraphs and the ones from the following iterations
  dt0 <- df[df$iter==0, ]
  dti <- df[df$iter!=0, ]
  combs <- expand.grid(paste0(dt0$iter, "_", dt0$IDs), paste0(dti$iter, "_", dti$IDs))

  # Get T0 and Ti subgraphs
  t0 <- gsub("_.*", "", combs$Var1)
  t0_graph <- gsub(".*_", "", combs$Var1)
  ti <- gsub("_.*", "", combs$Var2)
  ti_graph <- gsub(".*_", "", combs$Var2)

  # Add size and IC
  dt0_df <- dt0[match(t0_graph, dt0$IDs), ]
  stopifnot(identical(dt0_df$IDs, t0_graph))
  t0_size <- dt0_df$size
  t0_IC <- dt0_df$IC

  dti_df <- dti[match(ti_graph, dti$IDs), ]
  stopifnot(identical(dti_df$IDs, ti_graph))
  ti_size <- dti_df$size
  ti_IC <- dti_df$IC

  # Create the new table
  combs <- data.frame(t0, t0_graph, t0_size, t0_IC, ti, ti_graph, ti_size, ti_IC)
  combs <- combs[order(combs$t0_graph), ]
  combs$jc <- apply(combs, 1, function(x) getJaccardIndex(x["t0_graph"], x["ti_graph"]))

  return(combs)

}


getSubgraphsTypes <- function(combs, info="summary") {

  results <- data.frame()

  # Select a GO subgraph
  for (s0 in unique(combs$t0_graph)) {
    dfs0 <- combs[combs$t0_graph==s0, ]

    # Select the GO subgraph with the maximum jc for each iteration
    iters <- as.numeric(as.character(unique(combs$ti)))
    types <- c()

    for(iter in iters) {
      dfsi <- dfs0[dfs0$ti==iter, ]
      dfsim <- dfsi[which.max(dfsi$jc), ]
      dif <- as.numeric(gsub(" \\(.*", "", dfsim$ti_IC)) - as.numeric(gsub(" \\(.*", "", dfsim$t0_IC))
      if(is.na(dif)) {
        dif <- 0
      }
      type <- "non-maintained"

      # Classify GO subgraphs
      if(dfsim$jc>=0.7 & dif>=0.1) {
        type <- "specialized"
      } else if(dfsim$jc>=0.7 & dfsim$ti_size>dfsim$t0_size & dif<0.1) {
        type <- "expanded"
      } else if (dfsim$jc>=0.7 & dfsim$ti_size<=dfsim$t0_size & dif<0.1) {
        type <- "maintained"
      }

      if(info=="all") {
        results <- rbind(results, data.frame(dfsim, type=type))
      } else {
        types <- c(types, type)
      }
    }

    if(info=="summary") {
      results <- rbind(results, c(dfsim$t0_graph, dfsim$t0_size, dfsim$t0_IC, types))
    }

  }

  if(info=="summary") {
    colnames(results) <- c("t0_graph", "t0_size", "t0_IC", paste0("type_", 1:(ncol(results)-3)))
  }

  return(results)
}


getSpAndExpSubgraphs <- function(results) {

  results$first_type <- rep(NA, nrow(results))
  results$first_iter <- rep(NA, nrow(results))

  for (i in 1:nrow(results)) {
    types <- unlist(results[i, -1])
    index_sp <- min(match("specialized", types))
    index_ex <- min(match("expanded", types))

    if(is.na(index_sp)) { index_sp <- 10 }
    if(is.na(index_ex)) { index_ex <- 10 }

    if(index_sp<index_ex) {
      results$first_type[i] <- "specialized"
      results$first_iter[i] <- index_sp
    } else if (index_ex<index_sp) {
      results$first_type[i] <- "expanded"
      results$first_iter[i] <- index_ex
    }

  }

  # Get the number of specialized and expanded subgraphs
  iters <- colnames(results)[-c(1, ncol(results), ncol(results)-1)]
  results$first_iter <- factor(results$first_iter, levels=iters)
  results$first_type <- factor(results$first_type, levels=c("expanded", "specialized"))

  tab <- table(results$first_iter, results$first_type)

  return(tab)

}



compareT0vsTisubgraphs <- function(df) {

  df <- addIC(df)
  combs <- getCombsT0vsTi(df)
  results <- getSubgraphsTypes(combs)
  tab <- getSpAndExpSubgraphs(results)

  return(tab)
}




##################################################################################################################################################



getCombsTivsT0 <- function(df) {

  results <- data.frame()
  iters <- unique(df$iter)
  iters <- iters[iters!=0]

  for (iter in iters) {
    # Get Ti subgraphs
    current <- df[df$iter==iter, ]

    # Get subgraphs from previous iters
    previous <- df[df$iter<iter, ]

    # Get all possible combinations between Ti and subgraphs from previous iterations
    combs <- expand.grid(paste0(current$iter, "_", current$IDs), paste0(previous$iter, "_", previous$IDs))
    colnames(combs) <- c("current_graph", "previous_graph")

    # Split graph and iter info
    current_iter <- gsub("_.*", "", combs$current_graph)
    current_graph <- gsub(".*_", "", combs$current_graph)
    previous_iter <- gsub("_.*", "", combs$previous_graph)
    previous_graph <- gsub(".*_", "", combs$previous_graph)

    # Add IC and terms names
    current_df <- df[match(current_graph, df$IDs), ]
    stopifnot(identical(current_df$IDs, current_graph))
    current_IC <- current_df$IC
    current_size <- current_df$size

    previous_df <- df[match(previous_graph, df$IDs), ]
    stopifnot(identical(previous_df$IDs, previous_graph))
    previous_IC <- previous_df$IC
    previous_size <- previous_df$size

    # Create the final table
    combs <- data.frame(current_iter, current_size, current_IC, current_graph,
                        previous_iter, previous_size, previous_IC, previous_graph)

    # Get Jaccard index for each pair of subgraphs
    combs$jaccard <- apply(combs, 1, function(x) getJaccardIndex(x["current_graph"], x["previous_graph"]))
    results <- rbind(results, combs)
  }

  results <- results[order(results$current_iter, results$current_graph), ]

  return(results)

}


getNewSubgraphs <- function(combs) {

  # Detect new subgraphs
  results <- data.frame()
  iters <- unique(combs$current_iter)

  for(iter in iters) {
    # Compare Ti subgraphs with subgraphs from previous iterations
    res <- combs[combs$current_iter==iter, ]

    # Detect completely new Ti subgraphs not found in previous iterations
    for(g in unique(res$current_graph)) {
      g_df <- res[res$current_graph==g, ]
      selected <- g_df[which.max(g_df$jaccard), ]

      if (selected$jaccard==0) {
        type <- "new"
      } else {
        type <- "nr"
      }

      results <- rbind(results, cbind(selected, type))
    }
  }

  return(results)

}


compareTivsT0subgraphs <- function(df) {

  df <- addIC(df)
  combs <- getCombsTivsT0(df)
  tab <- getNewSubgraphs(combs)

  return(tab)
}




plotNewKnowledge <- function(df) {

  tab1 <- compareT0vsTisubgraphs(df)
  tab2 <- compareTivsT0subgraphs(df)
  data <- cbind(tab1, new=tab2[, "new"])

  mydata <- melt(data, id.vars="iter")
  colnames(mydata) <- c("iter", "type", "value")
  mydata$value <- as.numeric(as.character(mydata$value))

  g <- ggplot(data=mydata, aes(x=iter, y=value, group=type, color=type)) +
              geom_line()+
              geom_point() +
              theme_classic() +
              theme(text=element_text(size=16)) +
              xlab("Pseudo-cells iteration") +
              ylab("GO subgraphs with\n new knowledge") +
              labs(color="Type of knowledge") +
              scale_y_continuous(limits=c(0,max(mydata$value)+2), breaks=seq(0,max(mydata$value)+2,2)) +
              scale_x_continuous(limits=c(0,max(mydata$iter)), breaks=seq(0, max(mydata$iter), 1))

  return(list(plot=g, data=mydata))
}
