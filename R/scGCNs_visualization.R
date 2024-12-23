#' getModulesCompositionPerIter  - It returns the genes that made up each module in each iteration.
#' It is an auxiliary function to create a sunburst plot, sankey plot or enrichment plot.
#'
#' @param nets_dir the path where we can find all the networks belonging to the same cell type
#' @return a data frame containing the number of genes that remain together at the end after successive pseudo-cells iterations
#' using all possible combinations of modules from different iterations
#' @export

getModulesCompositionPerIter <- function(nets.dir) {

  # Load libraries
  require(stringr)

  # Check if the file has been already created
  nets <- list.files(nets.dir, full.names=T, recursive=T, pattern=".50.rds$")
  name <- gsub("\\..*", "", gsub("net", "", basename(nets[1])))

  # Get modules
  modulesList <- list()

  for (net in nets) {
    iter <- grep("iter", basename(net))
    if(length(iter)==0) {
      iter <- 0
    } else {
      iter <- as.numeric(gsub("\\..*", "", str_split(basename(net), "iter")[[1]][2]))
    }
    net <- readRDS(net)
    modules <- unique(net$moduleColors)

    for (module in modules) {
      modulesList[[paste0("Iter",iter, "-", module, "_")]]<- names(net$moduleColors)[net$moduleColors==module]
    }
  }

  modulesList <- unlist(modulesList)
  names(modulesList) <- gsub("_.*", "", names(modulesList))

  # Get gene names
  genes <- names(readRDS(nets[1])$moduleColors)
  results <- data.frame()

  for(g in genes) {
    path <- modulesList[grep(paste0("^", g, "$"), modulesList)]
    row <- c(g, names(path))
    results <- rbind(results, row)
  }

  # Step 4: create a data frame with two columns
  results[, -1] <- as.data.frame(apply(results[, -1], 2, function(x) gsub("-", ": ", x)))
  colnames(results) <- c("gene", paste0("Iter", 0:(length(nets)-1)))
  results$path <- apply(results[, -1], 1 , function(x) paste0(x, collapse="-"))
  results$path <- gsub("Iter", "", results$path)

  return(results)
}


#' getColors - It returns the colors to be used for the sunburst plot
#'
#' @param criteria the criteria used to color the modules: "T0_modules", "all_modules" or "covariate".
#' @param modulesGenes if criteria!="covariate", a data frame with one row per gene and one column per iteration is needed. Each
#' @param nets_dir if criteria=="covariate", the path with all the networks that belong to the same cell type and different iterations
#' @param metadata if criteria=="covariate", the path or the data frame with the metadata of all the cells that belong to a specific cell type
#' @param cellID if criteria=="covariate", the name of the variable that contains the ID of each cell in the metadata table
#' @param covariate if criteria=="covariate", the name of the biological covariate we are interested in
#' @return a list containing both the colors assigned and the names of the modules
#' @export

getColors <- function(criteria=c("T0_modules", "all_modules", "covariate"),
                      modulesGenes=NULL,
                      nets_dir=NULL,
                      metadata=NULL,
                      cellID="cell_id",
                      covariate=c("Clinical_Dx")) {

  # Load library
  require(RColorBrewer)

  # Check arguments
  criteria <- match.arg(criteria)

  # If you want to color the modules based on the correlation with a covariate, criteria=="covariate"
  if(criteria=="covariate") {
    if(is.character(metadata)) {
      metadata <- readRDS(metadata)
    }

    mydf <- getModuleTraitCorr(nets_dir,
                               metadata,
                               cellID=cellID,
                               covariates=covariate)

    if(covariate=="Clinical_Dx") {
      cov <- "dx_pval"
    } else if(covariate=="Donor_age_years") {
      cov <- "age_pval"
    } else {
      cov <- "sex_pval"
    }

    mydf2 <- mydf[, c("iter", "module", cov)]
    modules <- paste0(mydf2$iter, ": ", mydf2$module)
    modules <- gsub("^iter", "Iter", modules)
    col <- rep("grey", length(modules))
    log <- -log10(as.numeric(as.character(mydf2[, cov])))
    rbPal <- colorRampPalette(c('darkgreen','red'))
    index <- which(log < -log10(0.05))
    colors <- rbPal(30)[as.numeric(cut(log, breaks = 30))]
    colors[index] <- "lightgrey"
    names(colors) <- modules
  }

  # If you want to color all modules based on T0 module genes, criteria=="T0_modules"
  overlap <- modulesGenes
  dat <- as.data.frame(table(overlap$path))
  colnames(dat) <- c("path", "overlap")
  all_modules <- as.vector(unlist(apply(overlap[, -c(1, ncol(overlap))], 2, function(x) unique(x))))
  all_modules <- gsub("Iter", "", all_modules)

  if(criteria=="T0_modules") {
    colors <- rep("lightgrey", length(all_modules))
    names(colors) <- all_modules
    T0_modules <- gsub("Iter", "", unique(overlap[, "Iter0"]))

    for (m in T0_modules) {
      mydat <- dat[grep(m, dat$path), ]
      mydat <- mydat[order(mydat$overlap, decreasing=T), ]
      selected <- unlist(str_split(mydat$path[1], "-"))
      mycolor <- gsub(".*\\: ", "", m)
      colors[match(selected, names(colors))] <- mycolor
    }
  }

  # If you want to assign one random color for each module except for T0 modules, which are colored with the same color
  # as their names, criteria="all_modules"
  if(criteria=="all_modules") {
    colors <- gsub(".*\\: ", "", all_modules)

    set.seed(1235)
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    colors <- rep(NA, length(all_modules))
    colors[grep("0: ", all_modules)] <- gsub("0: ", "", all_modules[grep("0: ", all_modules)])
    colors[-grep("0: ", all_modules)] <- sample(col_vector, length(colors[-grep("0: ", all_modules)]))
    names(colors) <- all_modules
  }

  return(colors)
}



#' getSunburstPlot - It returns a sunburst plot where the inner ring represent the T0 modules (the size of each segment is proportional to the number
#' of genes) and the rest of rings represent how the genes of each T0 module remains or not together after successive pseudo-cells iterations
#'
#' @param criteria the criteria used to color the modules: "T0_modules", "all_modules" or "covariate".
#' @param nets_dir the path with all the networks that belong to the same cell type and different iterations.
#' @param metadata if criteria=="covariate", the path or the data frame with the metadata of all the cells that belong to a specific cell type
#' @param cellID if criteria=="covariate", the name of the variable that contains the ID of each cell in the metadata table
#' @param covariate if criteria=="covariate", the name of the biological covariate we are interested in
#' @return a list containing both the colors assigned and the names of the modules
#' @export

getSunburstPlot <- function(criteria=c("T0_modules", "all_modules", "covariate"),
                            nets_dir=NULL,
                            metadata=NULL,
                            cellID="cell_id",
                            covariate=c("Clinical_Dx"),
                            showLabels=F) {

  require(stringr)
  require(sunburstR)
  require(htmltools)
  require(d3r)

  mg <- getModulesCompositionPerIter(nets_dir)

  colors <- getColors(criteria=criteria,
                     modulesGenes=mg,
                     nets_dir=nets_dir,
                     metadata=metadata,
                     cellID=cellID,
                     covariate=covariate)

  cols <- as.vector(colors)
  modules <- names(colors)
  modules <- gsub("Iter", "", modules)
  dat <- as.data.frame(table(mg$path))
  colnames(dat) <- c("path", "overlap")
  sb3 <- sund2b(dat, width="150%", colors= list(range = c("lightblue", cols), domain = c("ROOT", modules)), showLabels=showLabels, rootLabel="ALL")

  div(
    style="display: flex; align-items:center;",
    sb3
  )

  return(list(sunburst=sb3, overlap=mg, colors=colors))

}




#' getSankeyPlot - It returns a sankey plot that shows how the genes of a T0 module of interest are distributed in other modules through pseudo-cells
#' iterations. The plot contains as many columns as pseudo-cells iterations and each column contains as many boxes as modules in that iteration.
#' @param module_name the name of the T0 module selected to get the sankey plot.
#' @param nets_dir if modulesGenes is not provided, the path where all the networks of the same cell type and different iterations are located.
#' @return a sankey plot
#' @export

getSankeyPlot <- function(module_name="pink",
                          nets_dir=NULL) {

  # Load libraries
  require(stringr)
  require(ggsankey)
  require(dplyr)
  require(ggplot2)

  # Get the table
  mg <- getModulesCompositionPerIter(nets_dir)

  # Step 0
  df <- as.data.frame(table(mg$path))
  colnames(df) <- c("path", "overlap")
  newdf <- as.data.frame(do.call("rbind", str_split(df$path, "-")))
  t <- newdf
  colnames(t) <- paste0("T", 0:(ncol(t)-1))
  t$overlap <- df$overlap
  t <- t[t$overlap>10, ]

  # Step 1
  t$color <- rep("white", nrow(t))
  t$color[t$T0==paste0("0: ", module_name)] <- module_name
  df <- t %>% make_long(paste0(colnames(t)[1:(ncol(t)-2)]))
  df$color <- rep("white", nrow(df))

  for (i in 1:nrow(df)) {
    iter <- as.vector(unlist(df[i, "x"]))
    module <- as.vector(unlist(df[i, "node"]))
    myt <- t[, match(iter, colnames(t))]
    color <- unique(t$color[grep(module, myt)])

    if(length(color)>1) {
      color <- setdiff(color, "white")
    }
    df$color[i] <- color
  }

  # Step 2
  o <- mg[mg$Iter0==paste0("Iter0: ", module_name), ]
  lens <- apply(o[, -c(1,ncol(o))], 2, function(x) table(x))
  lens <- unlist(lens)
  names(lens) <- gsub(".*\\.", "", names(lens))
  names(lens) <- gsub("Iter", "", names(lens))
  min <- round((5*nrow(o))/100, 0)
  lens <- lens[lens>min]

  dagg <- df%>%
    dplyr::group_by(node)%>%
    tally()

  for (i in 1:nrow(dagg)) {
    index <- match(dagg$node[i], names(lens))

    if(!is.na(index)) {
      dagg$n[i] <- as.vector(lens[index])
    }
  }

  TotalCount = nrow(o)

  dagg <- dagg%>%
    dplyr::group_by(node)%>%
    dplyr::mutate(pct = n/TotalCount)

  # Step 3
  df2 <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
  df2$node <- paste0(" ", df2$node)
  df2$next_node <- paste0(" ", df2$next_node)

  p <- ggplot(df2, aes(x = x
                       , next_x = next_x
                       , node = node
                       , next_node = next_node
                       , fill= ifelse(node %in% paste0(" ", names(lens)), module_name, NA)
                       , color="black"
                       , label = ifelse(node %in% paste0(" ", names(lens)), paste0(node, "\nn=", n, ' (',  round(pct* 100,1), '%)' ), NA))) +
    geom_sankey(flow.alpha = .5,
                node.color = "gray30",
                node.size=0.5) +
    geom_sankey_label(size = 3.2, color = "black", fill = "white", width=2, hjust = -0.1, alpha=0.65) +
    theme_sankey(base_size = 18) +
    labs(x=NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5, size=14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill='transparent')) +
    scale_fill_manual(values=c(module_name, "white"), na.value=NA) +
    scale_color_manual(values=c("transparent"))

  return(p)
}




#' getEnrichmentPlot - It returns an enrichment plot for a T0 module of interest and show how the enrichment changes through pseudo-cells iterations
#' for this same set of genes.
#' @param source the name of the source to be used for the enrichment analysis (GO, KEGG, REAC).
#' @param modulesGenes a data frame with as many rows as genes and as many columns as iterations. For each gene (row), you can see the module where
#' it belongs in each pseudo-cell iteration.
#' @param nets_dir if modulesGenes is not provided, the path where all the networks of the same cell type and different iterations are located.
#' @return an enrichment plot
#' @export

getEnrichmentPlot <- function(module_name="pink",
                              source="GO:BP",
                              modulesGenes=NULL,
                              nets_dir=NULL) {

  # Load library
  require(ComplexHeatmap)
  require(circlize)
  require(gprofiler2)

  # Get maintained genes through pseudo-cells iterations
  if(is.null(modulesGenes)) {
    mg <- getModulesGenes(nets_dir)
  } else {
    if(is.character(modulesGenes)) {
      modulesGenes <- readRDS(modulesGenes)
    }
    mg <- modulesGenes
  }

  o <- mg[mg$Iter0==paste0("Iter0: ", module_name), ]
  min <- round((nrow(o)*5)/100, 0)
  modules <- module_name
  all.genes <- list()
  index <- 2:(ncol(o)-1)

  for (i in index) {
    tab <- table(o[, i])
    selected <- names(tab)[tab>=min]

    for(s in selected) {
      iter <- gsub("\\:.*", "", s)
      myo <- o[, c("gene", iter)]
      all.genes[[s]] <- myo$gene[myo[, 2]==s]
    }
  }

  # Get enrichment for each module
  enrich <- gprofiler2::gost(all.genes,
                             correction_method="g_SCS",
                             sources=source,
                             organism = "hsapiens",
                             exclude_iea = T,
                             highlight = T)$result

  if(source!="KEGG") {
    enrichh <- enrich[enrich$highlighted==TRUE, ]
  } else {
    enrichh <- enrich
  }

  terms <- unique(enrichh$term_name)
  modules <- unique(enrich$query)
  results <- data.frame()

  for(term in terms) {
    myrow <- rep(0, length(modules))
    names(myrow) <- modules
    mydf <- enrichh[enrichh$term_name==term, ]
    index <- match(mydf$query, names(myrow))
    myrow[index] <- -log10(mydf$p_value)
    results <- rbind(results, myrow)
  }

  colnames(results) <- modules
  rownames(results) <- terms
  lens <- lengths(all.genes)
  lens <- lens[match(colnames(results), names(lens))]
  stopifnot(identical(names(lens), colnames(results)))

  results <- rbind(results, as.vector(lens))
  results <- as.matrix(results)
  rownames(results)[nrow(results)] <- "overlap"

  max <- round(max(results[-nrow(results), ]),0) + 5

  groups <- gsub("\\:.*", "", colnames(results))
  groups <- gsub("Iter", "T", groups)

  col_fun = colorRamp2(c(0, 10, 20), c("white", "blue", "red"))

  colors <- rep("white", ncol(results))
  names(colors) <- colnames(results)

  for (iter in 0:(ncol(mg)-2)) {
    modules <- colnames(results)[grep(iter, colnames(results))]
    selected <- names(sort(results[nrow(results), modules], decreasing=T))[1]

    if(is.null(selected)) {
      selected <- modules
    }
    colors[match(selected, names(colors))] <- module_name
  }

  colnames(results) <- rep(" ", ncol(results))

  h <- Heatmap(results[-nrow(results), ],
               col = col_fun,
               cluster_columns = F,
               column_split = groups,
               row_km=3,
               column_title_gp = gpar(fontsize = 9),
               name="-log10P",
               border = TRUE,
               column_names_rot = 45,
               # show_column_names = F,
               row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9),
               top_annotation = HeatmapAnnotation(`mantained genes` = anno_barplot(results[nrow(results), ],
                                                                                   extend =0.15,
                                                                                   gp=gpar(fill=colors)),
                                                  annotation_name_gp= gpar(fontsize = 9)),
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           legend_width = unit(3, "cm"),
                                           title_gp=gpar(fontsize=9),
                                           labels_gp=gpar(fontsize=7)))

  draw(h, heatmap_legend_side="bottom", annotation_legend_side="right")

  return(list(heatmap=h, enrichment=enrichh))

}


#' getRobustAndEvolvedModules - For each module of the T0 network, it returns the percentage of the genes maintained together at Tn.
#' @param nets.dir the path where the networks of the same cell type and different iterations are located.
#' @param path the directory where the statistics will be saved.
#' @param name the name of the cell type.
#' @return it returns a table with one row for each T0 module and some columns that show how the module composition has changed from T0 to Tn.
#' @export

getRobustAndEvolvedModules <- function(nets.dir, path=getwd(), name="CT1") {

  file <- list.files(paste0(path, name, "_stats.rds"))

  if(length(file)==0) {
    dirs <- list.dirs(nets.dir, recursive=F)
    stats <- data.frame()
    overlap <- getModulesCompositionPerIter(nets.dir=nets.dir)


    for(module in unique(overlap$Iter0)) {
      # Get T0 module genes distribution
      o <- overlap[overlap$Iter0==module, ]
      ms <- nrow(o)
      iters <- colnames(o)[3:(ncol(o)-1)]

      for(iter in iters) {
        # Get stats
        min <- round((5*ms)/100,0)
        tab <- table(o[, iter])
        tab_filt <- tab[tab>=min]

        # Get entropy
        genes <- as.vector(tab/ms)
        e <- round(- sum(genes * log2(genes))/length(genes),4)

        # Save results
        row <- data.frame(name=name, T0_module=module, module_size=ms, iter=gsub("Iter", "", iter), nModulesFilt=length(tab_filt),
                          remained_genes=max(tab_filt), remained_genes_per=round((max(tab_filt)*100)/ms, 1), nModules=length(tab), entropy=e)

        stats <- rbind(stats, row)
      }
    }

    saveRDS(stats, paste0(path, name, "_stats.rds"))

    stats$lost_genes_per <- 100-stats$remained_genes_per
    stats$type <- rep("partiallyEvolved", nrow(stats))
    stats$type [stats$remained_genes_per>=70] <- "maintained"
    stats$type [stats$remained_genes_per<20] <- "evolved"
    stats$type <- factor(stats$type , levels=c("maintained", "partiallyEvolved", "evolved"))

  } else {
    stats <- readRDS(file)
  }

  return(stats)

}


#' summariseModulesAnnotations - For each module of the T0 network, it returns the main annotations, including the size, the percentage of genes maintained together at Tn, the association with clinical covariates and the association with different list of markers.
#' @param nets.dir the path where the networks of the same cell type and different iterations are located.
#' @param path the directory where the statistics will be saved.
#' @param name the name of the cell type.
#' @return it returns a table with one row for each T0 module and some columns that shows tha main annotations of each module.
#' @export


summariseModulesAnnotations <- function(nets.dir=nets.dir, path=getwd(), name="CT1") {

  # Load libraries
  require(stringr)

  nets.dir2 <- list.dirs(nets.dir, recursive=F)
  net.dir <- nets.dir2[grep("iter0", nets.dir2)]
  net.files <- list.files(paste0(net.dir, "/Net/"), full.names=T)

  # Separate the network from the rest of the files
  net <- readRDS(net.files[grep(".50.rds$", net.files)])
  files <- net.files[-grep(".50.rds$", net.files)]
  files <- files[-grep("functionalEnrich", files)]
  files <- files[-grep("phenotypeEnrich", files)]
  files <- files[-grep("moduleMembership", files)]

  # Get name of the iteration
  iter <- str_split(files[1], "iter")[[1]]
  iter <- paste0("iter", gsub("\\..*", "", iter[length(iter)]))

  # Add the basic modules' information to the final table
  tab <- sort(table(net$moduleColors), decreasing=T)
  results <- data.frame(iter=rep(iter, length(tab)), modules=names(tab), size=as.vector(tab))

  # Add the percentage of genes maintained together at Tn
  stats <- getRobustAndEvolvedModules(nets.dir=nets.dir, path=path, name=name)
  stats2 <- stats[match(results$modules, gsub("Iter0: ", "", stats$T0_module)), ]
  stopifnot(identical(gsub("Iter0: ", "", stats2$T0_module), results$modules))

  results <- cbind(results, genesTogetherAtTn_per=stats2$remained_genes_per)

  # Add the covariates p-val to the final table
  traitCorr <- read.csv(files[grep("traitCorr", files)])[, -1]
  traitCorr <- traitCorr[match(results$modules, traitCorr$module), ]
  stopifnot(identical(traitCorr$module, results$modules))
  results <- cbind(results, traitCorr[grep("pval", colnames(traitCorr))])
  files <- files[-grep("traitCorr", files)]

  # For each fgsea file
  for (file in files) {

    # Get name of the markers list/source of annotation
    name <- str_split(file, "_")[[1]]
    name <- gsub(".csv", "", name[length(name)])

    # Read file
    file <- read.csv(file)[, -1]

    # Add this information to the final table
    all_values <- rep(NA, nrow(results))
    names(all_values) <- results$modules

    for(m in unique(file$module)) {
      df <- file[file$module==m, ]
      newvalue <- paste0(paste0(df$pathway," P<(", format(df$padj2, digits=3), ")"), collapse=", ")
      all_values[match(m, names(all_values))] <- newvalue
    }
    results <- cbind(results, all_values)
    colnames(results)[ncol(results)] <- paste0("markers_", name)
  }

  return(results)
}




#' getInteractiveBrowser - It creates an interactive browser with T0 modules annotations and highlights the most relevant modules based on user's criteria.
#' @param nets.dir the path where the networks of the same cell type and different iterations are located.
#' @param path the directory where the statistics will be saved.
#' @param name the name of the cell type.
#' @param criteria the criteria to highlight the most relevant modules.
#' @return it returns an interative table with one row for each T0 module and some columns that shows the main annotations of each module.
#' @export

getInteractiveBrowser <- function(nets.dir, path, name, criteria=NULL) {

  # Load libraries
  library(reactable)
  library(htmltools)
  library(fontawesome)
  library(shiny)
  library(dplyr)
  library(reactablefmtr)


  # Get modules annotations
  results <- summariseModulesAnnotations(nets.dir=nets.dir, path=path, name=name)

  if(is.null(criteria)) {
    results$criteria_maintained[results$genesTogetherAtTn_per>=70] <- "YES"
    results$criteria_maintained[results$genesTogetherAtTn_per<70] <- "NO"

    results$criteria_evolved[results$genesTogetherAtTn_per>=70] <- "NO"
    results$criteria_evolved[results$genesTogetherAtTn_per<70] <- "YES"
  }

  r <- reactable(
    results,
    # defaultSorted = `Mantained genes at Tn (%)`,
    defaultColDef = colDef(
      header = function(value) gsub(".", " ", value, fixed = TRUE),
      cell = function(value) format(value, nsmall = 1),
      align = "left",
      minWidth = 100
    ),
    columnGroups = list(
      colGroup(name = "genesTogether", columns = colnames(results)[3:4]),
      colGroup(name = "covariates", columns = colnames(results)[grep("pval", colnames(results))]),
      colGroup(name = "Markers", columns = colnames(results)[grep("markers", colnames(results))]),
      colGroup(name = "Criteria", columns = colnames(results)[grep("criteria", colnames(results))])
    ),
    style = list(fontSize = '12px'),
    fullWidth = FALSE,
    height = 800,
    resizable = TRUE,
    wrap = T,
    highlight = TRUE,
    filterable = TRUE,
    compact = TRUE,
    pagination = FALSE,
    onClick = "expand",
    columns = list(
      `Iteration`= colDef(minWidth=80),
      `Module`=colDef(minWidth=80),
      `Module size`=colDef(filterable = FALSE, minWidth = 100, align = "left", cell = function(value) {
        width <- paste0(value / max(results$`Module size`) * 100, "%")
        bar_chart(value, width = width)
      }),
      `Mantained genes at Tn (%)`= colDef(filterable=FALSE, minWidth = 120, align = "left", cell = function(value) {
        width <- paste0(value / max(results$`Mantained genes at Tn (%)`) * 100, "%")
        bar_chart(value, width = width, fill = "#fc5185", background = "#e1e1e1")
      }),
      `Criteria Maintained`= colDef(filterable=FALSE, minWidth = 90, cell = function(value) {
        if (value=="NO") "\u274c NO" else "\u2714\ufe0f YES"}),
      `Criteria Evolved`= colDef(filterable=FALSE, minWidth = 80, cell = function(value) {
        if (value=="NO") "\u274c NO" else "\u2714\ufe0f YES"})),

    elementId = "download-table",
    showSortIcon = FALSE,
    borderless = TRUE,
    class = "standings-table",
    # theme = spotify_theme(),
    rowStyle = JS("function(rowInfo) {
      if (rowInfo.values['Criteria Maintained'] =='YES') {
        return { background: 'lightblue' }
      }
    }"),
    # rowClass = JS("function(rowInfo) {
    #   if (rowInfo.values['Criteria Maintained'] =='YES') {
    #     return 'bold'
    #   }
    # }"),
    rowClass = JS("
      function(rowInfo, state) {
        const firstSorted = state.sorted[0]
        if (firstSorted && firstSorted.id === 'Iteration') {
          const nextRow = state.pageRows[rowInfo.viewIndex + 1]
          if (nextRow && rowInfo.values['group'] !== nextRow.group) {
            return 'group-last'
          }
        }
      }"
    ),
    # selection = "multiple"
  )


  tbl <- htmltools::browsable(
    tagList(
      tags$button(
        tagList(fontawesome::fa("download"), "Download as CSV"),
        onclick = "Reactable.downloadDataCSV('DNs-download-table', 'DNs_interesting_modules.csv')"
      ), r
    ))


  div(class = "standings",
      div(class = "title",
          h2("Modules chracterization through pseudo-cells iterations"),
          "scCoExpNets"
      ),
      tbl,
      # "Forecast from before 3rd group matches"
  )

}
