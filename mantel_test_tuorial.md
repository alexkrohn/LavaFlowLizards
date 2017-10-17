Using Mantel Tests to Test for IBD and IBE
================

This script will take a genetic distance matrix, a geographic distance matrix and an ecological distance matrix to run Mantel and partial Mantel tests to determine the effect of geographic distance, and, separately, ecological distance on the observed genetic distances between individuals. In this case, the ecological distance matrix is binary (0 for two individuals from the same environment, 1 for two individuals from different environments).

Import Genetic Distance Matrix
==============================

See the [IBE\_IBE\_plots\_tutorial](https://github.com/alexkrohn/LavaFlowLizards/blob/master/IBE_IBD_plots_tutorial.md) for more information on how to generate this matrix using ANGSD.

``` r
geneticmatrix <- as.matrix(read.table("~/Documents/Berkeley/IBDIBE/reDone_analyses/genetics/sceloporus/redone_justPA_individuals/ngsdist/ngsdist_matrix_k70u5", 
    skip = 2))

geneticmatrix <- geneticmatrix[, 2:13]  #Remove the name column
```

Create Geographic Distance Matrix from Lat Longs
------------------------------------------------

``` r
metadata <- read.table("~/Documents/Berkeley/IBDIBE/reDone_analyses/genetics/sceloporus/redone_justPA_individuals/ngsdist/justPAindividuals_metadata.txt", 
    sep = "\t")  # Import the data
metadata <- metadata[, c(1, 2, 4)]  # Keep only the relevant columns
names(metadata) <- c("lat", "long", "Sample")
```

Give the lat longs a projection and create the distance matrix

``` r
require(RgoogleMaps)
require(raster)
require(gdata)
```

``` r
locs <- subset(metadata, select = c("long", "lat"))
coordinates(locs) <- c("long", "lat")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(locs) <- crs.geo

geodist <- pointDistance(locs, longlat = TRUE)
geodist <- geodist/1000  # Convert the distances to km instead of m
```

Import the ecological distance matrix
-------------------------------------

Don't forget to check the dimensions of the matrices to make sure they're equal.

``` r
eco <- as.matrix(read.table("~/Documents/Berkeley/IBDIBE/reDone_analyses/genetics/sceloporus/redone_justPA_individuals/ngsdist/ecodist_justpaindividuals.txt"))

dim(geodist)
```

    ## [1] 12 12

``` r
dim(geneticmatrix)
```

    ## [1] 12 12

``` r
dim(eco)
```

    ## [1] 12 12

Run the Mantel tests
--------------------

First, the IBD question: does genetic distance correlate with geographic distance?

``` r
require(vegan)
```

``` r
mantel(ydis = geneticmatrix, xdis = geodist)
```

    ## 
    ## Mantel statistic based on Pearson's product-moment correlation 
    ## 
    ## Call:
    ## mantel(xdis = geodist, ydis = geneticmatrix) 
    ## 
    ## Mantel statistic r:  0.62 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.161 0.227 0.309 0.370 
    ## Permutation: free
    ## Number of permutations: 999

It does!

Next, the IBE question: does genetic distance correlate with ecological distance when you control for the effect of geographic distance?

``` r
mantel.partial(ydis = geneticmatrix, xdis = eco, zdis = geodist)
```

    ## 
    ## Partial Mantel statistic based on Pearson's product-moment correlation 
    ## 
    ## Call:
    ## mantel.partial(xdis = eco, ydis = geneticmatrix, zdis = geodist) 
    ## 
    ## Mantel statistic r: 0.187 
    ##       Significance: 0.11 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.189 0.249 0.282 0.340 
    ## Permutation: free
    ## Number of permutations: 999

It doesn't!
