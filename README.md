
# scCoExpNets

<!-- badges: start -->
<!-- badges: end -->

This package is based on the CoExpNets R package, but adapted to create and analyze cell-type-specific gene co-expression networks from single-cell/single-nucleus RNA-seq data. 

For each cell type, we combat the sparsity of the expression matrix by collapsing pairs of cells that come from the same individuals while the number of cells> 200, therefore, obtaining multiple expression matrices and multiple co-expression GCNs. In this way, we have different points of view for the interpretation of the same data set, so that we can see the insights that are not affected by sparsity but also those that only appear when we reduce sparsity. We provide the necessary tools to analyze this large number of co-expression networks, detecting the modules with the highest correlation with clinical diagnosis, but also the ones enriched with cell type markers or neurological disorders gene panels. Moreover, we detect the modules that are well-preserved through pseudo-cells iterations, that is, the ones non-affected by the sparsity, and vice versa, and represent the results with a sunburst.

We applied this new tool on a snRNA-seq dataset human post-mortem substantia nigra pars compacta tissue of 13 controls and 14 Parkinson’s disease (PD) cases (18 males and 9 females) with 30-99 years. For more details, please consult the link to the manuscript:

## Installation

Recommended that you install first CoExpNets package like so:

``` r
devtools::install_github('juanbot/CoExpNets')
```

You can install the development version of scCoExpNets like so:

``` r
devtools::install_github('aliciagp/scCoExpNets')
```


## Credits

The development of this package has been possible under the supervision of Juan A. Botía from the University of Murcia and Sebastian Guelfi and many collaborators from Verge Genomics start-up.



