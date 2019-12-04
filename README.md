# RadialGO
## A radial tree visualization for Gene Ontology enrichment analysis

The goal of this package is to provide a visualization of GO enrichment analysis
in a radial tree format.

## Installation

You can install the latest version using
``` r 
require("devtools")
install_github("brucosper/RadialGO", build_vignettes=TRUE)
library("radialgo")
```

To run the package use
```r
radialgo()
```
## Overview

The package takes a list of Gene Ontology IDs and their corresponding 
P-values, and takes the graph implied by the top nodes (i.e. gets their ancestors),
removing the ancestors which have a corresponding P-value smaller than the user-specified cutoff.
It then graphs these nodes, along with the edge relationships, and outputs it.
![](./inst/extdata/PEREIRA_B_A1.png)

## Contributions

The author of the package is Bruno Pereira. The functions available within this 
package include. 

``` r 
library("RadialGO")
lsf.str("package:RadialGO")
```

- generateGraph

The generateGraph function was authored by Bruno Pereira, as well as the internal helper functions, 
which build the list of nodes implied by the input, as well as the edges. It uses the DiagrammeR
package for generating graphical output.
The GO.db and the AnnotationDbi packages were used to retrieve information about
GO terms. The invert function in the searchable package was used for 
facilitating data retrieval from a named numerical vector. The shiny package was 
used to provide an interactive app for running the function.
