
# Install scCoExpNets

``` r
devtools::install_github("aliciagp/scCoExpNets")
```


# First steps with scCoExpNets R package

By calling a single function, scCoExpNets performs the following tasks for each cell type:

    - Create a set of single-cell gene co-expression networks (scGCNs): scCoExpNets creates a scGCN from the original expression matrix. Then, it reduces the sparsity of the gene expression by creating pseudo-cells while the number of pseudo-cells remains above a certain threshold. For each pseudo-cell matrix, the corresponding scGCN is created.
    
    - Estimate the importance of each gene in each module: for each network, scCoExpnets calculates the module membership of each gene in each module.
    
    - Get the functional terms of each module: scCoExpNets returns the functional enrichment of each module of each network using GO, KEGG and REAC datasets. 
    
    - Get the phenotype terms of each module: if phenotypeEnrich==T, scCoExpNets also returns the phenotype enrichment of each module of each network using all databases available at PhenoExam R package. 
    
    - Get the association between modules and donors' traits: if moduleTraitCorr==T, scCoExpNets estimates the correlation between the eigengenes of the module and the traits selected (included in covsToCorr) 


``` r
matrix_files <- list.files(exprData_path, full.names=T)
metadata_files <- list.files(metadata_path, full.names=T)
cell_types <- gsub("_exprData.rds", "", gsub(".*/", "", matrix_files))

for (c in cell_types) {
  exprData <- matrix_files[grep(c, matrix_files)]
  metadata <- readRDS(metadata_files[grep(c, metadata_files)])
  metadata$Sex <- factor(metadata$Sex, levels=c(0,1))
  metadata$Clinical_Dx <- factor(metadata$Clinical_Dx, levels=c(0,1))

  get_all_GCNs(exprData=exprData,
              metadata=metadata,
              gcnName=c,
              idColname="cell_id",
              batchColname=NULL,
              covsToCorrect=NULL,
              filteringGenes=F,
              filteringCutoff=NULL,
              filteringProportion=NULL,
              plots=F,
              moduleTraitCorr=T,
              covsToCorr=c("Clinical_Dx", "Donor_age_years", "Sex"),
              phenotypeEnrich=T,
              path=save_path)
}

```


# scCoExpNets to annotate the modules with specific sources

scCoExpNets also includes a function to evaluate module enrichment for user-provided gene lists, allowing for more specific characterization according to the study objective or the starting hypothesis.

``` r
# Note: save the fgsea results in the same folder where the network is stored
net <- paste0(networks_path, "netDA_neurons_6.17.it.50.rds")
exprData <- paste0(exprData_path, "DA_neurons_6_exprData.rds")
markersList <- "cellTypeMarkers"

write.csv(getMarkersListEnrich(net=net,
                               markersList=markersList,
                               exprData=exprData), paste0(results_path, "_cellTypeMarkers.csv"))
    
```


# scCoExpNets to detect modules preserved in an independent dataset

scCoExpNets also includes a function to check whether one network module of a specific cell type is preserved in the gene expression profile of an independent dataset of the same tissue and cell type.

``` r
results <- getPreservedModules(network=dns_net,
                               expr.data.files = list(dns_resids, dns_ind_resids),
                               tissues = c("dns_discovery", "dns_replication"),
                               permutations = 200,
                               maxModuleSize = 5000,
                               maxGoldModuleSize = 400,
                               randomSeed = 1)
```


# scCoExpNets provides new ways to visualize module composition through pseudo-cells iterations

scCoExpNets provides an interactive sunburst plot to represent an overview of all the T0 module gene composition changes through pseudo-cells iterations.


``` r
nets_dir <- paste0(networks_path, "DA_neurons_6/")
metadata <- paste0(metadata_path, "DA_neurons_6_metadata.rds")

getSunburstPlot(criteria=c("T0_modules"),
                nets_dir=nets_dir,
                metadata=metadata,
                cellID="cell_id",
                covariate=NULL,
                showLabels=F)
```


Additionally, scCoExpNets includes a sankey plot to represent the changes in a T0 module gene composition through pseudo-cells iterations.

``` r
nets_dir <- paste0(networks_path, "DA_neurons_6/")
metadata <- paste0(metadata_path, "DA_neurons_6_metadata.rds")

getSankeyPlot(module_name="pink",
              nets_dir=nets_dir)
```



# scCoExpNets provides an interactive browser to detect interesting modules

Finally, scCoExpNets summarise the annotations of the T0 modules in an interactive table. Additionally, it highlights the most relevant modules based on users' criteria. 

``` r
nets.dir <- paste0(networks_path, "DA_neurons_6/")
path <- paste0(networks_path, "DA_neurons_6/")
name <- "DA_neurons_6"

getInteractiveBrowser(nets.dir, path, name)
```
