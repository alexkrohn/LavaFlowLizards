Making IBD/IBE Plots from ANGSD Datasets
================

Generating the Genotypes from ANGSD
-----------------------------------

First, you should follow the numerous instructions online to get bam files and corresponding reference genome for all of your individuals. I also created a "keep file," a file of sites that were high quality enough to keep in the analysis. I used only those sites from the keep file using the -sites flag. To generate the correct input for ngsDist I used this command:

``` bash
angsd -ref Reference.fasta -sites Sites_passed_filters.keep -only_proper_pairs 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -doPost 1 -doGeno 8 -gl 1 -bam 12inds_bamlist -out ngsdist_input 
```

Create the Genetic Distance Matrix
----------------------------------

Next, following the great tutorial by Matt Fumagalli [here](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), I used the -doGeno 8 output to create a genetic distance matrix.

First, save the number of sites that were included in your analysis.

``` bash
NSITES=`zcat ngsdist_input.mafs.gz | tail -n+2 | wc -l`
```

Next, I created a list of labels for each individual. These must be in the same order as the individuals in your bamlist.

In my bamlist (12inds\_bamlist), my individual names were separated by two underscores. This made making a label file very easy.

``` r
bamlist <- read.table("12inds_bamlist", sep = "_")
head(bamlist)  #the individual names are in the 7th column
```

    ##                                                V1                      V2
    ## 1 /global/home/kebi/data/Alex/IBDIBE/startingover 8Oct2015/correctCleanSE
    ## 2 /global/home/kebi/data/Alex/IBDIBE/startingover 8Oct2015/correctCleanSE
    ## 3 /global/home/kebi/data/Alex/IBDIBE/startingover 8Oct2015/correctCleanSE
    ## 4 /global/home/kebi/data/Alex/IBDIBE/startingover 8Oct2015/correctCleanSE
    ## 5 /global/home/kebi/data/Alex/IBDIBE/startingover 8Oct2015/correctCleanSE
    ## 6 /global/home/kebi/data/Alex/IBDIBE/startingover 8Oct2015/correctCleanSE
    ##                         V3              V4             V5
    ## 1 output/Sceloporus/justPA inds/individual contigs/PopGen
    ## 2 output/Sceloporus/justPA inds/individual contigs/PopGen
    ## 3 output/Sceloporus/justPA inds/individual contigs/PopGen
    ## 4 output/Sceloporus/justPA inds/individual contigs/PopGen
    ## 5 output/Sceloporus/justPA inds/individual contigs/PopGen
    ## 6 output/Sceloporus/justPA inds/individual contigs/PopGen
    ##                              V6    V7         V8
    ## 1 reference/alignment/EBRARK001 BD153 sorted.bam
    ## 2 reference/alignment/EBRARK001 BD154 sorted.bam
    ## 3 reference/alignment/EBRARK001 BD155 sorted.bam
    ## 4 reference/alignment/EBRARK001 BD158 sorted.bam
    ## 5 reference/alignment/EBRARK001 BD159 sorted.bam
    ## 6 reference/alignment/EBRARK001 CW100 sorted.bam

``` r
individuals <- bamlist[, 7]
write.table(individuals, "individual_labels", quote = F, row.names = F, col.names = F)
```

Finally, use ngsDist to actually make the matrix.

``` bash
ngsDist -geno ngsdist_input.geno.gz -n_ind 12 -n_sites $NSITES -labels individual_labels -o 12inds_ngs_dist_matrix
```

Make the IBD/IBE Plots
----------------------

In the next part, we will import the genetic distance matrix, create a matrix of geographic distances, and import a matrix of ecological distances.

First, import the ngsDist matrix, and remove the label column and header rows. We'll also change the upper triangle and diagnal row to be NA, just to ensure we're not plotting everything twice.

``` r
library(gdata)
geneticmatrix <- as.matrix(read.table("ngsdist_matrix_k70u5", skip = 2))
geneticmatrix <- geneticmatrix[, 2:13]

upperTriangle(geneticmatrix) <- NA
diag(geneticmatrix) <- NA
```

Next, import the geographic data. I have lats and longs for each of my lizards, so I'll start with that.

``` r
library(raster)

metadata <- read.table("justPAindividuals_metadata.txt", sep = "\t")
metadata <- metadata[, c(1, 2, 4)]  #only keep the lat, long and sample columns
names(metadata) <- c("lat", "long", "Sample")  #name the columns

# Give the lats and longs a projection.
locs <- subset(metadata, select = c("long", "lat"))
coordinates(locs) <- c("long", "lat")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(locs) <- crs.geo

# Create a geographic distance matrix from the lats and longs
geodist <- pointDistance(locs, longlat = TRUE)

# Convert the values into km instead of m, and, again, replace the diagnal with
# NA
geodist <- geodist/1000
diag(geodist) <- NA
```

Finally, import the ecodist matrix. For my data, these are just 0's and 1's. 0 if the comparison is between two individuals in the same environment (either lava-lava comparison, or non-lava to non-lava comparison), and 1 if the two individuals are from a different environment (lava to non-lava comparisons).

``` r
ecodist <- as.matrix(read.table("ecodist_justpaindividuals.txt"))
```

Now we've got the basic elements of our plot.

We will eventually color our pairwise comparisons by whether the two individuals are in the same environment or not. For the published paper, I did this by making the same environment comparisons hollow circles, and the different environments filled circles. At the end I'll show you how easy it is to make a color plot

``` r
ecodist[which(ecodist == 0)] <- 19
ecodist[which(ecodist == 1)] <- 1
upperTriangle(ecodist) <- NA
diag(ecodist) <- NA
```

A quick check to make sure all of our matricies are the same size:

``` r
dim(ecodist)
```

    ## [1] 12 12

``` r
dim(geodist)
```

    ## [1] 12 12

``` r
dim(geneticmatrix)
```

    ## [1] 12 12

The first attempt plot is very easy:

``` r
plot(x = geodist, y = geneticmatrix, pch = ecodist, main = "Sceloporus cowlesi", 
    xlab = "Pairwise distance (km)", ylab = "Genetic Distance", cex.lab = 1.5, cex.axis = 1.5, 
    cex = 2, ylim = c(0.2, 0.35))
legend(x = "topright", col = "black", pch = c(19, 1), legend = c("Same environment", 
    "Different environment"))
```

![](IBE_IBD_plots_tutorial_files/figure-markdown_github/unnamed-chunk-11-1.png)

We can add the "IBD" line to this plot, but doing a simple linear regression of genetic distance and pairwise distance. Notice that this is not the most appropriate way to calculate IBD because the data may violate assumptions of linear regressions. In the paper we used Mantel tests to determine significant patterns of IBD.

``` r
linear.model <- lm(lowerTriangle(geneticmatrix) ~ lowerTriangle(geodist))

plot(x = geodist, y = geneticmatrix, pch = ecodist, main = "Sceloporus cowlesi", 
    xlab = "Pairwise distance (km)", ylab = "Genetic Distance", cex.lab = 1.5, cex.axis = 1.5, 
    cex = 2, ylim = c(0.2, 0.35))
legend(x = "topright", col = "black", pch = c(19, 1), legend = c("Same environment", 
    "Different environment"))
abline(linear.model, lwd = 2)
```

![](IBE_IBD_plots_tutorial_files/figure-markdown_github/unnamed-chunk-12-1.png)

Changing this plot to color is as simple as modifying the eco matrix.

``` r
ecodist[which(ecodist == 19)] <- "dodgerblue3"
ecodist[which(ecodist == 1)] <- "tomato"

plot(x = geodist, y = geneticmatrix, col = ecodist, pch = 19, main = "Sceloporus cowlesi", 
    xlab = "Pairwise distance (km)", ylab = "Genetic Distance", cex.lab = 1.5, cex.axis = 1.5, 
    cex = 2, ylim = c(0.17, 0.3))
legend(x = "topright", col = c("dodgerblue3", "tomato"), pch = 19, legend = c("Same environment", 
    "Different environment"))
abline(linear.model, lwd = 2)
```

![](IBE_IBD_plots_tutorial_files/figure-markdown_github/unnamed-chunk-13-1.png)
