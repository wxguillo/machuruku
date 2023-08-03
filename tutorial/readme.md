# Machuruku: the Tutorial 2.0
This is the updated tutorial coinciding with the release of Machuruku 2.0. I originally developed Machuruku ~2019 with [Jason Brown](https://www.jasonleebrown.org/), my MS advisor (2017-2020). After the release of the original version, we published [a paper](https://academic.oup.com/sysbio/article/70/5/1033/6171196) in Systematic Biology in 2021 describing the method. Machuruku was my first major R project and as such the original version bore some of the hallmarks of a novice R coder, along with some other quirks that continued to bother me as I moved on into a doctoral program. A couple years later in 2023, with much more coding experience under my belt, I decided to revisit Machuruku which quickly ballooned into rewriting the entire thing, both in service to my own projects and the community at large. The scope of the changes merits a nomenclatural leap to Machuruku version 2.0.

The purpose of Machuruku remains unchanged: the package is for modeling, reconstructing, and visualizing niches in past and present climates, while accounting for evolution. Machuruku reconstructs climatic niches along a time-calibrated phylogeny and can project these niches into paleoclimatic data to visualize the potential geographic origins of lineages. The advances provided by Machuruku 2.0 are mostly practical. For example, many of the for-loops I was previously reliant on have been replaced with apply statements, which offer significant boosts in computation speed. Also boosting speed is Machuruku's new integration with the [Terra](https://github.com/rspatial/terra) R package, a faster upgrade to its predecessor Raster. The workings of the various functions provided in the package have been streamlined to be more user-friendly and flexible. There is also a potentially significant conceptual change to the way Machuruku calculates and visualizes uncertainty in ancestral niche reconstructions, which I will discuss in greater detail later. 

In this new tutorial, I will first provide a simple quick-start guide to using Machuruku 2.0 at its most basic. I then provide a deeper step-by-step dive into the each function and the possibilities they offer for reconstructing and visualizing ancestral niches. 
## Contents
* [Introduction](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#introduction)
  * [Installing Machuruku](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#installing-machuruku)
  * [Downloading and exploring tutorial data](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#downloading-and-exploring-tutorial-data)
    * [Loading tree](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#loading-tree)
    * [Loading occurrence data](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#loading-occurrence-data)
    * [Loading climate data](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#loading-climate-data)
* [Quick-start Guide](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#quick-start-guide)
  * [Estimating tip response curves](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#estimating-tip-response-curves)
  * [Estimating ancestral niches at a time-slice](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#estimating-ancestral-niches-at-a-time-slice)
  * [Projecting ancestral models into paleoclimatic data](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#projecting-ancestral-models-into-paleoclimatic-data)
* [Detailed Guide](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#detailed-guide)
  * [Estimating tip response curves](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#estimating-tip-response-curves-1)
    * [Accounting for spatial autocorrelation by rarefying occurrence data](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#accounting-for-spatial-autocorrelation-by-rarefying-occurrence-data)
    * [Visualizing climate response curves](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#visualizing-climate-response-curves)
  * [Estimating ancestral niches](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#estimating-ancestral-niches)
    * [Estimating niches at a time-slice](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#estimating-niches-at-a-time-slice)
    * [Estimating niches for all nodes](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#estimating-niches-for-all-nodes)
    * [Accounting for uncertainty in ancestral character estimation](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#accounting-for-uncertainty-in-ancestral-character-estimation)
    * [Accounting for uncertainty in divergence time estimation given one tree with error bars](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#accounting-for-uncertainty-in-divergence-time-estimation-given-one-tree-with-error-bars)
    * [Accounting for uncertainty in divergence time estimation given a Bayesian posterior distribution of trees](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#accounting-for-uncertainty-in-divergence-time-estimation-given-a-bayesian-posterior-distribution-of-trees)
    * [Saving and loading ace output](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#saving-and-loading-ace-output)
  * [Projecting ancestral models into paleoclimatic data](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#projecting-ancestral-models-into-paleoclimatic-data-1)
    * [Clipping response curve tails to produce cleaner models](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#clipping-response-curve-tails-to-produce-cleaner-models)
    * [Creating binary models](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#creating-binary-models)
    * [Incorporating uncertainty from ancestral character estimation](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#incorporating-uncertainty-from-ancestral-character-estimation)
    * [Limiting model area with inverse-distance weighting](https://github.com/wxguillo/machuruku/blob/main/tutorial/readme.md#limiting-model-area-with-inverse-distance-weighting)
## Introduction
### Installing Machuruku
To install Machuruku, simply use the `install_github()` function from `devtools`:
```
install.packages("devtools")
devtools::install_github("wxguillo/machuruku")
```
Load up Machuruku by running:
```
library(machuruku)
```
You can test if it worked by typing `machu` into the console (in RStudio) and seeing if all of the functions appear in autofill.
### Downloading and exploring tutorial data
I opted to provide the raw tutorial data in this repository rather than within the R package itself so that I can demonstrate how to load it. For now you will need 3 of the files: 
* `basslerigroup.treefile` - A nexus-formatted treefile containing the three members of the *Ameerega bassleri* group, a small clade of Amazonian dendrobatid poison frogs. This tree was time-calibrated in [BEAST 2](https://www.beast2.org/) and contains 95% HPD intervals for each node height (divergence time). It is a small subset of the UCE phylogeny of *Ameerega* presented in [Guillory et al. 2020](https://www.sciencedirect.com/science/article/pii/S1055790319304609), restricted to the *bassleri* group (*Ameerega bassleri, pepperi,* and *yoshina*) plus *A. silverstonei*, an in-genus outgroup taxon with a similar distribution.
* `basslerigroup.csv` - Occurrence data for each species in our phylogeny, in decimal-degree format. The first column is the species ID (which must match that in the tree), the second is longitude (x), and the third is latitude (y).
* `climate.zip` - A zipped file containing the climate data we will be using. If you unzip the file, you will get a folder called `climate/`, which in turn contains four subfolders: `current/`, `mis19/`, `mpwp/`, and `m2/`. These four subfolders correspond to different time periods (present, Marine Isotope Stage 19 (Pleistocene; 0.787 Ma), Mid-Pliocene Warming Period (3.205 Ma), and Marine Isotope Stage M2 (Pliocene; 3.3 Ma). Each subfolder contains 14 raster climate layers from [Paleoclim](http://www.paleoclim.org/), which I've cropped to an area surrounding the *bassleri* group's distribution and reduced in resolution to save space.  

Download and place these three files in a folder anywhere on your machine, and unzip `climate.zip`. Load up R and set your working directory to this location with the `setwd()` function.

#### Loading tree
There are multiple ways to load trees in R, such as `ape::read.tree`, but here let's use `ape::read.nexus` to load in our time-calibrated phylogeny. Install Ape if you haven't already using `install.packages("ape")`, then run: 
```
# load tree
library(ape)
tree <- read.nexus("basslerigroup.treefile")
```
Note that you shouldn't need to load in the ape package yourself, since it's exported in the namespace of machuruku already.

Simply use `plot()` to visualize the tree:
```
plot(tree)
```
![plot() tree](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/plot%20tree.png?raw=true)

#### Loading occurrence data
Loading occurrence data is simple. Use the `read.csv` or `read.delim` function:
```
# load occurrence data
occ <- read.csv("basslerigroup.csv")
```
We can visualize the occurrence data in a minute, after we've loaded the climate data as well.

#### Loading climate data
In Machuruku 2.0, I have mostly migrated the code to use the Terra function rather than its slower predecessor Raster. Terra is so much faster because, rather than trying to store the contents of a raster, which can get very large, on your machine's RAM, it simply reads them straight from the disk when need be. It also contains some new functions that can handle raster data more elegantly than the old Raster functions could. Here let's use `terra::rast`, the equivalent of the old `raster::raster`, to load our four climate datasets. Again, install Terra if you haven't already, then run:
```
# load climate data
library(terra)
current <- rast(list.files("climate/current", full.names=T))
mis19 <- rast(list.files("climate/mis19", full.names=T))
mpwp <- rast(list.files("climate/mpwp", full.names=T))
m2 <- rast(list.files("climate/m2", full.names=T))
```
Let's visualize both the climate and the occurrence data. First plot the current-climate Bio1 raster, then the occurrence data, and finally a legend.
```
# plot climate
plot(current$bio_1)
# add occurrence data
taxa <- unique(occ$species)
cols <- c("cyan", "yellow", "red", "orange")
for (i in taxa) points(subset(occ, species==i)$long, subset(occ, species==i)$lat, 
                       col = cols[which(taxa==i)], 
                       pch=19, cex=0.75)
legend(x = "bottomleft", legend = taxa, col = cols, pch=19, bty="n")
```
![bio1 w/ points](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/bio1%20w%20pts.png?raw=true)

Now we can see that the species' ranges do not really overlap, but are very close to each other. They evolved in a place with lots of topographic heterogeneity (the foothills of the Andes), which probably promoted their diversification. We might be able to see some niche evolution in this group since they occupy subtly different habitats.
## Quick-start Guide
Machuruku is composed of three principal steps. In this section I'll demonstrate them without any extraneous bells-and-whistles. Let's set a phylogenetic niche modeling goal for ourselves: To see what the ancestors of the *bassleri* group and *A. silverstonei* were doing in the Late Pliocene, approximately 3.3 million years ago (the putative age of our M2 climate layers). What's important to note here is that it's highly unlikely that any of our nodes will exactly line up with the age of our paleoclimate data, e.g., none of the nodes will be 3.3 million years old. To accommodate this, Machuruku can interpolate ancestral niches along the branches subtending the desired time-slice. This can be kind of tricky to visualize, so Machuruku has a built-in function, `machu.treeplot()`, that can do it for us.
```
# visualize tree with timeslice at 3.3 Ma
machu.treeplot(tree, timeslice=3.3)
```
![treeplot 3.3](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/treeplot.png?raw=true)

The time-slice at 3.3 Ma recovers *two* taxa, not the four at the tips, so we'll be making two ancestral niche models.
### Estimating tip response curves
The first step in Machuruku is the `machu.1.tip.resp()` function. This function takes occurrence data and present-day climate data, constructs a Bioclim niche model for each taxon, and then characterizes the response of each taxon to each climate variable as a skew-normal distribution. 
```
# estimate tip response curves
resp <- machu.1.tip.resp(occ, current)
``` 
### Estimating ancestral niches at a time-slice
The second step is using `machu.2.ace()` to estimate ancestral niches at our desired time-slice. This step takes the response table and uses the `ace()` function from [Ape](http://ape-package.ird.fr/) to estimate the niche parameters for whatever taxa exist at our 3.3 Ma timeslice.
```
# estimate ancestral niches at timeslice
ace <- machu.2.ace(resp, tree, timeslice=3.3, unc=T)
```
The `unc=T` argument is used to incorporate ancestral character estimation uncertainty into the final models.

### Projecting ancestral models into paleoclimatic data
The final major step is to use `machu.3.anc.niche()` to project the ancestral niche models into our 3.3 Ma paleoclimate data. This step loops through each taxon in `ace` and constructs a Bioclim model for it in the provided paleoclimate data.
```
# project ancestral niches into paleoclimate data
mod <- machu.3.anc.niche(ace, m2)
```
This produces two ancestral niche models, one for each branch-taxon present at 3.3 Ma, and projected into our 3.3 Ma paleoclimate data. We can visualize the results with a final function, `machu.plotmap()`.
```
# visualize ancestral niches
machu.plotmap(mod, plot="together", plot.asp=20/9, axes=F, to.scale=T)
```
![2 models](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/2%20models.png?raw=true)

We can see that the ancestor of the *bassleri* group may have existed further to the south than the present distribution of the group. On the other hand, the ancestor of *silverstonei* seems constrained to its present-day distribution.
## Detailed guide 
### 1. Estimating tip response curves
Given the occurrence data and present-day climate rasters, the first major step in Machuruku, `machu.1.tip.resp`, constructs niche models for each present-day taxon in the dataset, then characterizes the response of each taxon to each climate variable as a skew-normal distribution. First, the `dismo::bioclim` function is used to extract the climate values of each variable at each occurrence point for each species. Then, the `sn::selm` function fits a skew-normal distribution to each variable to describe the relationship between occurrence and climate (`fGarch::snormFit` was used before, but I found some issues with it; now the sn package is used throughout Machuruku). The skew-normal distribution is itself characterized by three components: mean, standard deviation, and skew (the upper and lower 95% confidence limits are no longer utilized). This is a modified version of the [Bioclim](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12144) niche modeling algorithm (widely regarded as the first such algorithm), which uses uniform distributions for each variable instead of skew-normal. At its most basic, the function looks like this:
```
### 1. estimate tip response curves
resp <- machu.1.tip.resp(occ, current, verbose=T)
```
With the `verbose=T` option turned on, it will print the following output to the screen:
```
[1] "'clim' is a SpatRaster, attempting to convert to RasterStack for compatibility with dismo::bioclim..."
[1] "Processed taxon bassleri"
[1] "Processed taxon pepperi"
[1] "Processed taxon silverstonei"
[1] "Processed taxon yoshina"
```
On my machine this process took about 10 seconds to complete, the vast majority of it being the conversion of the SpatRasters to a RasterStack. `machu.1.anc.niche` will accept SpatRaster format, but will automatically convert it to the older RasterStack format for compatibility with `dismo::bioclim`. 

To view the results for the first two climate variables, run `resp[,1:6]`:
```
             bio_1_mean bio_1_stdev  bio_1_skew bio_10_mean bio_10_stdev bio_10_skew
bassleri       238.2260    18.64478   -1.987488    243.8598     18.48886   -1.893566
pepperi        253.8157    19.63661 -183.446148    260.4739     20.25066 -183.446148
silverstonei   183.9266    26.75841  183.446148    189.7586     26.46845  183.446148
yoshina        249.0784    14.71177   -3.913331    254.3711     14.71131   -3.594554
```
This "response table" shows the mean, standard deviation, and skew for each climate variable, for each taxon. We can visualize these distributions using the function `machu.respplot`:
```
# visualize response curves for two variables
machu.respplot(resp, clim=c("bio_1", "bio_10"), plot="t")
```
![respplot1](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/respplot1.png?raw=true)

The `machu.respplot` function (response-plot) can plot response curves in various combinations. Here I specify the first two climate variables in `resp` by name, as well as to plot the curves "together" in one window. The curves are nearly identical because the bio1 (mean annual temp) and bio10 (mean temp of warmest quarter) bioclimatic variables are very similar. To view the response curves for *every* climate variable, we can simply omit the `clim` argument:
```
# visualize response curves for all variables
machu.respplot(resp, plot="t", legend.cex=0.6)
```
![respplot2](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/respplot2.png?raw=true)

Another feature of `machu.1.anc.niche`, new to Machuruku 2.0, is the visualization of present-day niche models with `dismo::predict`. Activate this feature with the `plot` argument:
```
# visualize niche models
machu.1.tip.resp(occ, current, plot="t")
```
![currentmodels1](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/currentmodels1.png?raw=true)

This is similar to `machu.respplot` where `plot="t"` (equivalent to `plot="together"`) will show every map in one window. You can also add the occurrence points to confirm that the estimated niche model lines up with the known distribution of the species:
```
dev.off()
machu.1.tip.resp(occ, current, plot="s", plot.points=T)
```
![currentmodels2](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/currentmodels2.png?raw=true)

Here `dev.off()` is used to reset the plot window. The `plot="s"` (equivalent to `plot="separately"`) option shows every map in its own window; here I'm only showing the first one for *A. bassleri*. 

Finally, the function can also output these niche models instead of the response table with the `output.bioclim` option. This can be useful if you'd like to plot, save, or analyze the present-day niche models for each species on your own. 
```
mod <- machu.1.tip.resp(occ, current, output.bioclim=T)
par(mfrow=c(2,2), mar=c(0,0,2,0))
lapply(1:4, function(x) plot(mod[[x]], axes=F, legend=F, box=F, main=names(mod)[x]))
```
![currentmodels3](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/currentmodels3.png?raw=true)

The output `mod` is a list where each element is a RasterLayer; the `lapply` function loops through each one and plots it after setting up a 2x2 plotting window with `par`.
#### Accounting for spatial autocorrelation by rarefying occurrence data
Occurrence data are often spatially autocorrelated due to systematic biases in the way we collect specimens, i.e. they tend to be clustered in easily accessible locations. Spatial autocorrelation can induce downstream impacts on niche modeling, necessitating the use of algorithms to "rarefy", or thin, occurrence data beforehand. In Machuruku, the function `machu.occ.rarefy` can do this for you. 
```
# rarefy occurrence data
occ.r <- machu.occ.rarefy(occ, rarefy.dist=10, plot=T)
```
The `rarefy.dist` option specifies the amount of rarefication to perform; in my example, I've set it to 10 km (the default unit; decimal degrees are also available with `rarefy.units="dd"`) such that all points must be at least 10 km from any other point. Adjust this value to fit the spatial scale of your data. Also, the `plot=T` option (new for Machuruku 2.0) visualizes the extent of rarefication performed:
![rarefied](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/rarefied.png?raw=true)

The function will tell you that we have gone from 73 to 33 points after rarefication. However, there is a problem: `occ.r` contains fewer than 10 occurrence data points for all species except *bassleri* now. The `selm` function that `machu.1.tip.resp` uses to fit a skew-normal distribution to the climate responses requires at least 10 points for each taxon. The best solution is simply to find more occurrence data ([GBIF](https://www.gbif.org/) is a great resource), but `machu.1.tip.resp` will automatically generate random occurrence points up to n=10 for each taxon that requires them. The random points are selected from within a minimum convex polygon surrounding the occurrence points for each species. This is a decent solution for species with a single contiguous range; less so for those with strongly disjunct distributions. 
If we run `machu.1.tip.resp` again, with `occ.r`, the function will generate random points. 
```
resp <- machu.1.tip.resp(occ.r, current, plot="t", plot.points=T, verbose=T)
```
It will warn us that it is doing so:
```
[1] "The following taxa have fewer than 10 occurrence points: pepperi, silverstonei, yoshina"
[1] "Warning: adding random points up to n=10 for each of these species, contained within MCP defined by points."
[1] "'clim' is a SpatRaster, attempting to convert to RasterStack for compatibility with dismo::bioclim..."
[1] "Plotting all plots together. n2mfrow chose 2 rows and 2 columns based on an aspect ratio of 16/9"
[1] "Processed taxon bassleri"
[1] "Processed taxon pepperi"
[1] "Processed taxon silverstonei"
[1] "Processed taxon yoshina"
```
When the function plots the niche models and occurrence points, the randomly generated points are highlighted in red:
![currentmodels4](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/currentmodels4.png?raw=true)

Additionally, the saved output `resp` is a list containing both the normal output (in this case, a response table) as well as the modified occurrence data table with a new column `"rand"` specifying which points were randomly generated. 
#### Selecting the most important climate variables
Another important task to do prior to creating the present-day niche models with `machu.1.top.env` is to reduce the climate dataset. Multiple correlated climate rasters being provided to the niche modeling algorithm can induce downstream biases just as spatial autocorrelation in the occurrence data can. Machuruku contains the function `machu.top.env` for identifying the most important climate variables for each taxon by constructing generalized boosted regression niche models. There are three algorithms for selecting the best variables: "contribution", "nvars", and "estimate". 
The "contribution" method only retains those variables that contribute more than x% of relative importance to each species. This is the default algorithm. The default minimum relative importance, specified by `contrib.greater`, is 5%. 
```
# select most important climate variables
# contribution method
micv.contrib <- machu.top.env(occ.r, current, method="contrib", contrib.greater=5, verbose=T)
```
This prints progress to the screen; the function takes a bit of time to run. It requires occurrence data and present-day climate data to construct the niche models with. For occurrence data, I've provided the spatially rarefied output from before, so as not to introduce bias into the niche models the function constructs. For climate data, I've provided the same SpatRaster `current`, which the function has to convert to a RasterStack for compatibility. Once the function is finished, running `micv.contrib' displays which variables were retained:
```
[1] "bio_1"  "bio_11" "bio_12" "bio_13" "bio_14" "bio_15" "bio_16" "bio_18" "bio_19" "bio_4"  "bio_8"  "bio_9"
```
This is 12 of the 14 starting variables, which is probably too many. If we want a specific number of the "best" variables, we can use the "nvars" algorithm. 
```
# nvars method
micv.nvars <- machu.top.env(occ.r, current, method="nvars", nvars.save=7, verbose=T)
```
The default value of `nvars.save`, which specifies the number of best variables to retain, is 5; here I've set it to 7. Running `micv.nvars` displays the seven variables retained:
```
[1] "bio_4"  "bio_15" "bio_18" "bio_14" "bio_16" "bio_13" "bio_12"
```
The final algorithm is "estimate", which chooses the number of best variables by systematically removing variables until average change in the model exceeds the original standard error of deviance explained. This algorithm is more computationally intensive than the others and takes longer to finish. 
```
micv.est <- machu.top.env(occ.r, current, method="estimate", verbose=T)
```
On my machine, it took a minute or two to finish. Running `micv.est` displays the variables retained:
```
[1] "bio_1"  "bio_11" "bio_12" "bio_13" "bio_14" "bio_15" "bio_16" "bio_17" "bio_18" "bio_19" "bio_4"  "bio_8"  "bio_9"
```
In this case it retained 13 of the 14 original variables, which also isn't very helpful. Here I personally found the "nvars" method most useful, so I'm going to reduce the climate variable set to just those, and use this reduced set for the rest of the tutorial.
```
# reduce climate datasets to only most important variables
current.reduced <- current[[micv.nvars]]
mis19.reduced <- mis19[[micv.nvars]]
mpwp.reduced <- mpwp[[micv.nvars]]
m2.reduced <- m2[[micv.nvars]]
```
It's important to make sure that all time periods have the same set of climate variables. Next I'll re-run `machu.1.tip.resp` using the newly reduced climate dataset, as well as the spatially rarefied occurrence data from before.
```
# re-run tip response table with reduced climate dataset and rarefied occurrence data
resp <- machu.1.tip.resp(occ.r, current.reduced, plot="t", plot.points=T)
```
![currentmodels5](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/currentmodels5.png?raw=true)

These new models are actually somewhat broader than the previous ones (this is especially visible with *yoshina*), as we've removed some of the noise from the climate dataset. The models will be slightly different anyway since we re-sampled new random occurrence points to get to 10 for each species. Now we're ready to move to the ancestral character estimation step.
### 2. Estimating ancestral niches
#### Visualizing timeslices
The second major step in Machuruku is using `machu.2.ace` to estimate ancestral niches for a set of taxa along the phylogenetic tree. The primary inputs are the tree and the response table from `machu.1.tip.resp`. The function uses the `ape::ace` to estimate the climate response parameters (mean, stdev, skew) for each climate variable at each node with Brownian motion. There are two possible outputs: the reconstructed values for every node in the tree (along with the tips, which remain unchanged from the `resp` input), or reconstructed values at one or more timeslices. If a timeslice is specified, `machu.2.ace` will identify whatever branches existed at that point in time, and interpolate the climate response values from the nodes or tips at either end of the branch to that point. To visualize this, use the function `machu.treeplot` and specify the dates of our three paleoclimate datasets for the `timeslice` parameter.
```
# visualize time-slices
dev.off()
options(warn=-1)
machu.treeplot(tree, timeslice = c(0.787,3.205,3.3))
```
Use `dev.off()` to clear the plot window and `par` parameters, then use `options(warn=-1)` to turn off warnings since this function tends to generate ones that don't really matter. For the `timeslice` option, I've specified a `c` vector with the putative ages (in Ma) of each of our paleoclimate datasets. Here is the result:

![treeplot1](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/treeplot1.png?raw=true)

This looks pretty good, but there are a lot of graphical options in `machu.treeplot` that can be tweaked to make the figure look better. 
```
machu.treeplot(tree, timeslice = c(0.787,3.205,3.3), x.l.lim = 13, x.u.lim=-3, nodelabsize = 0.35, col = "skyblue", timelaboffset = -0.2)
```
Here I manually changed the vertical size of the window to make the tree more compact, and adjusted the horizontal scaling with the `x.l.lim` and `x.u.lim` parameters to better fill the space. This necessitated adjusting the size of the circles at each node with `nodelabsize`, and adjusting the position of the time labels with `timelaboffset`. Finally, I changed the color of the timeslice lines to "skyblue" with `col`, just for fun.

![treeplot2](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/treeplot2.png?raw=true)

The `machu.treeplot` function was previously based in the ggtree R package, but this caused some issues with the code and often required a separate installation, so for Machuruku 2.0 the function is now based mostly in the [phytools](http://blog.phytools.org/) package. The function will automatically select a scale and tick mark scheme based on the order of magnitude of the age of the tree, and in general produces cleaner results. 

There are three timeslices visible: two near 3 Ma, which are the m2 and mpwp paleoclimate datasets, and one near 0.7 Ma, which is the mis19 dataset. The two older timeslices "slice through" two branches: the one from node 1 to *silverstonei*, and the one from node 1 to node 2. The more recent mis19 timeslice "slices through" four branches: node1-*silverstonei*, node2-*pepperi*, node3-*yoshina*, and node3-*bassleri*. As such, when projecting ancestral niche models into the m2 or mpwp dataset, there should only be two ancestral niche models, one for each of these ancestral "branch-taxa", and when projecting into the mis19 dataset, there should be four. The `machu.2.ace` function takes care of this.
#### Estimating ancestral niches at timeslices
To run `machu.2.ace` with one or more timeslices, specify the `timeslice` argument:
```
# estimate ancestral niches at timeslices
resp <- resp[[2]]
ace.ts <- machu.2.ace(resp, tree, timeslice=c(0.787,3.205,3.3))
```
The first line of code is simply reassigning `resp` to solely the response table, to minimize confusion. The `machu.2.ace` function takes two primary inputs: the response table (`resp`) and a time-calibrated phylogeny (`tree`). The tree must be ultrametric (i.e., all tips at the same height), and the names at each tip must match the taxa in the response table exactly (the function will make these checks for you). The `timeslice` argument is identical to that from `machu.treeplot`. 

The ability to specify multiple timeslices is a significant improvement for Machuruku 2.0. In the previous version only one time-slice at a time was possible, necessitating calling the function in a for loop or other repeating construct if one wanted to use multiple timeslices. This was enabled by re-coding the function as a series of nested `lapply` statements; as such the output of `machu.2.ace` is now formatted as a list, where each element corresponds to one timeslice. Running `ace.ts` shows this; each element is named after a timeslice. Running `names(ace.ts)` shows:
```
[1] "timeslice_0.787" "timeslice_3.205" "timeslice_3.3"
```
Let's examine this output more closely. Running `ace.ts[[1]][,1:8]` shows the first eight columns of the first element, "`timeslice_0.787`" (mis19):
```
                   branch_start branch_end bio_4_mean bio_4_stdev bio_4_skew bio_15_mean bio_15_stdev bio_15_skew
Node3-bassleri                7          1   488.8174    61.24433   45.76368    24.84982     3.640581    24.50394
Node3-yoshina                 7          4   464.7665    45.17915   159.7225    23.44164     7.757321    141.1891
Node2-pepperi                 6          2   561.1887    17.35599   165.7438    29.50367     2.949623    15.94719
Node1-silverstonei            5          3   602.4069    53.03126   8.241115    40.01165     5.556379   -174.4965
```
Each row in this table corresponds to a "branch-taxon" (the same ones visualized by the third (rightmost) line in the `machu.treeplot` images above). The first two columns give the indices of the nodes or tips at either end of the branch-taxon's branch, and subsequent columns represent climate response parameters for two climate variables, Bio4 (temperature seasonality) and Bio15 (precipitation seasonality). These values were linearly interpolated from the values at the nodes or tips at either end of the branch; for a timeslice exactly halfway between two nodes, the interpolated parameter value would be halfway between the nodes' values. 

Moving on to the second timeslice, "`timeslice_3.205`" (mpwp), running `ace.ts[[2]][,1:8]` shows that there are only two branch-taxa in this timeslice. We already knew this from visualizing the tree and timeslices with `machu.treeplot` before.
```
                   branch_start branch_end bio_4_mean bio_4_stdev bio_4_skew bio_15_mean bio_15_stdev bio_15_skew
Node1-Node2                   5          6   514.5594    40.36691   119.6053    27.01833     4.706292    43.41313
Node1-silverstonei            5          3   591.3218     51.4332   22.29359    38.37208      5.44911   -146.9995
```
#### Estimating ancestral niches at each node
The alternative mode of `machu.2.ace` is returning ancestral niche parameters for each node, instead of those branch-taxa at one or more timeslices. This is simply accomplished by omitting the `timeslice` argument:
```
# estimate ancestral niches at each node
ace.n <- machu.2.ace(resp, tree)
```
The output of this function is a list with a single element called `"tips_and_nodes"`. This name is important in triggering downstream effects in the `machu.3.anc.niche` function. Running `ace.n[[1]][,1:7]` allows us to examine the output as before:
```
                 times bio_4_mean bio_4_stdev bio_4_skew bio_15_mean bio_15_stdev bio_15_skew
bassleri      0.000000   483.4028   69.894895   8.911683    24.38023     2.883884    4.735996
pepperi       0.000000   581.1957    8.074401 183.446148    30.65780     2.239828    2.407736
silverstonei  0.000000   606.0149   53.551384   3.667379    40.54528     5.591293 -183.446148
yoshina       0.000000   446.5674   45.290139 183.446148    22.22352     9.188911  183.446148
Node1        11.577178   552.9404   45.900022  70.949409    32.69519     5.077696  -51.792643
Node2         2.710415   512.2920   40.040041 122.479628    26.68297     4.684351   49.037393
Node3         2.267554   499.0037   44.970337 115.091966    25.73324     5.064125   61.692618
```
The `"tips_and_nodes"` table contains 7 rows: 4 for the tips and 3 for the nodes. The climate response values for the tips are identical to those from the response table; they are included in the `machu.2.ace` output for the sake of comparison with the nodes. The first column `times` shows the divergence time for each node; this can be useful when deciding which paleoclimate dataset to project a certain node's ancestral niche model into. 

We can visualize the evolution of climatic responses along the tree with the `machu.respplot` function, which works with the output of `machu.2.ace` as well as `machu.1.tip.resp`. 
```
# visualize niche evolution
machu.respplot(ace.n[[1]], clim="bio_12", fill=T)
```
![respplot3](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/respplot3.png?raw=true)

Here I've just visualized the Bio12 variable, annual precipitation, which is most intuitive to look at. I've also specified `fill=T` (a new option for Machuruku 2.0) to fill in the space beneath each response curve and make them easier to see. In this case we can see that *silverstonei* prefers climates with higher precipitation, while the actual members of the *bassleri* group are clustered at the lower end of the spectrum. The interesting part is how the nodal taxa's climate responses are intermediate between these two extremes.
#### Characterizing uncertainty in ancestral niche estimation
The values offered by ancestral character estimation algorithms like `ace` are subject to considerable uncertainty, which increases further back in time. Relying on a single point estimate of a given value may be too precise for some to stomach, especially for something as nebulous as the standard deviation of a species' response to precipitation. To that end, Machuruku 2.0 includes a new method of characterizing this uncertainty to produce more conservative estimates of the ancestral niche. If the user specifies `unc=T`, `machu.2.ace` retains the upper and lower 95% confidence intervals for each parameter, provided by `ace`. When passed to `machu.3.anc.niche` later on, these "uncertainty values" produce broader, less specific ancestral niche models. 
```
# include uncertainty
ace.ts.u <- machu.2.ace(resp, tree, timeslice=c(0.787,3.205,3.3), unc=T)
```
Other than including `unc=T`, this command is the same as two sections ago where we took 3 timeslices of the tree. Running `ace.ts.u[[1]][,1:8]` shows how `unc=T` changes the `machu.2.ace` output:
```
                   branch_start branch_end bio_4_mean bio_4_mean_lCI bio_4_mean_uCI bio_4_stdev bio_4_stdev_lCI bio_4_stdev_uCI bio_4_skew bio_4_skew_lCI bio_4_skew_uCI
Node3-bassleri                7          1   488.8174       457.3668        520.268    61.24433        47.95568        74.53297   45.76368       14.31306       77.21431
Node3-yoshina                 7          4   464.7665       433.3159       496.2171    45.17915         31.8905        58.46779   159.7225       128.2718       191.1731
Node2-pepperi                 6          2   561.1887       532.6862       589.6913    17.35599        6.146849        29.39901   165.7438       137.2413       194.2464
Node1-silverstonei            5          3   602.4069       586.6109        618.203    53.03126        53.46503        59.70547   8.241115      -7.554937       24.03717
```
Columns 3 through 8 are the climate response values for a single climate variable, Bio4. Each original parameter (mean, stdev, skew) now has an additional two columns directly following it, "`lCI`" (lower confidence interval) and "`uCI`" (upper confidence interval). The original parameter values represent the medians. When this output is passed to `machu.3.anc.niche` later on, the function detects the presence of these uncertainty values and automatically incorporates them into the analysis. Rather than use a single skew-normal distribution for a given climate variable, the function takes *every combination* of median, lCI, and uCI values for each parameter, constructs a skew-normal distribution describing the relationship of suitability and climate for each combination, converts these to direct estimates of suitability for each pixel in the niche model, then sums and rescales all niche models to produce the result. There are 3<sup>3</sup>=27 combinations for each climate variable. We can visualize these with `machu.respplot`:
```
# visualize uncertainty
machu.respplot(ace.ts.u[[1]], clim="bio_12")
```
![respplot4](https://github.com/wxguillo/machuruku/blob/machuruku-2.0/tutorial/images/respplot4.png?raw=true)

When uncertainty values are passed to this function, it automatically plots the 27 skew-normal distributions formed by the different parameter combinations. There isn't a ton of variability in this case, so the default `fill=F` actually makes it easier to see. Node1-*silverstonei* is very constrained, while Node3-*bassleri* and Node3-*yoshina* have more variation. Thick dotted lines are drawn to show the median distribution for each taxon that would be used were `unc=F`. 
#### Characterizing uncertainty in divergence times

Let's confirm that our output contains the ancestral taxa that we'd expect. Remember that Mis19 corresponds to the third time-slice pictured above, so we'd expect four ancestral taxa. The `ace.mis19` output looks like this:
```
$bio_1_mean
                name branch_start branch_end    value
1     Node7-bassleri            7          1 229.7862
2      Node7-yoshina            7          4 233.7089
3      Node6-pepperi            6          2 240.4937
4 Node5-silverstonei            5          3 192.6381

$bio_1_stdev
                name branch_start branch_end    value
1     Node7-bassleri            7          1 16.15572
2      Node7-yoshina            7          4 13.93988
3      Node6-pepperi            6          2 14.06111
4 Node5-silverstonei            5          3 13.88672
...
```
(I've only shown the first two elements of the list, to save space.) We can see that the output is formatted as a list, where each element corresponds to one of the climatic response parameters (e.g., a column from the output of `machu.1.tip.resp()`). Each element is a table, where the value for that parameter is shown for each taxon. Each row corresponds to a taxon; we see right away that each table has four rows, so Machuruku got the right number of taxa at this time period. 
> Also shown are the nodes/tips that lie on either side of (subtend) the branch that that taxon represents (`branch_start` is the older one and `branch_end` is the newer one). The numbers refer to node and tip IDs that are encoded into the tree when it is imported by the `read.nexus()` function. The tips are given IDs first, as 1-through-4 (or however many taxa you have). The nodes are then added; in this case, there are 3 nodes, so the nodes are given IDs 5-through-7. Referring to these numbers in conjunction with a tree visualization from `machu.treeplot()` can help you get your bearings in the tree.

Let's check one of the Pliocene outputs, because for these (the first two time-slices in the figure above), we'd expect two taxa each.
```
$bio_1_mean
                name branch_start branch_end    value
1        Node5-Node6            5          6 232.5243
2 Node5-silverstonei            5          3 197.6712

$bio_1_stdev
                name branch_start branch_end    value
1        Node5-Node6            5          6 14.62475
2 Node5-silverstonei            5          3 13.97985
...
```
These are the first two elements of `ace.mpwp`. We can see that there are now two taxa retained, one along the branch leading from node 5 to node 6 (the ancestor of the *bassleri* group), and one leading from node 5 to *A. silverstonei.*
#### Estimating niches for all nodes
Let's say you're not interested in a particular time-slice, but you want to compare the niches at each node in a phylogeny with the niches of the tips (descendant taxa) in the same climate-space. This can be easily accomplished simply by running `machu.2.ace` without specifying a time-slice.
```
ace.all <- machu.2.ace(resp, bassleritree)
```
The output of this, `ace.all`, will provide us with niche models for the nodes, as well as retaining those of the tips (to make it easy to compare them later):
```
$bio_1_mean
          name      time NA    value
1     bassleri  0.000000 NA 228.1429
2      pepperi  0.000000 NA 243.3333
3 silverstonei  0.000000 NA 191.0000
4      yoshina  0.000000 NA 234.1507
5        Node5 11.577178 NA 215.0977
6        Node6  2.710415 NA 233.5538
7        Node7  2.267554 NA 232.8777

$bio_1_stdev
          name      time NA    value
1     bassleri  0.000000 NA 16.87901
2      pepperi  0.000000 NA 13.82269
3 silverstonei  0.000000 NA 13.85641
4      yoshina  0.000000 NA 13.48532
5        Node5 11.577178 NA 14.30230
6        Node6  2.710415 NA 14.64380
7        Node7  2.267554 NA 14.79502
...
```
These are the first two elements of `ace.all`. There are seven taxa, the first four of which are the extant taxa (bearing no difference from their values in `resp`), and the last three are each of the nodes; Node5 is the root node, Node6 is the immediate common ancestor of the *bassleri* group, and Node7 is the common ancestor of *A. bassleri* and *A. yoshina*. 
> Rather than tell you about branch starts and ends, this output shows you the age of each taxon. The third `NA` column is a placeholder (four columns are required in this output).

We can visualize how these taxa have evolved in climate space with respect to each climate variable by using `machu.respplot()`:
```
machu.respplot(ace.all, comb = T, legend.pos = "topright")
```
![visualization of ancestors](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/ancestor%20vis.png?raw=true)

There's a lot to see in this plot, partially occluded by the legend. But right away we can see how Node5, the putative common ancestor of the *bassleri* group and *A. silverstonei*, occupied an intermediate climate space between these two descendants. (Of course, including all ~32 *Ameerega* species in our analyses would change these results quite a bit!)
#### Accounting for uncertainty in ancestral character estimation
As any good phylogeneticist will tell you (and a bunch of classic papers), using ancestral character estimation (=ancestral state reconstruction) can be fraught because of the uncertainty inherent in this technique, and our assumptions that evolution is essentially happening randomly (Brownian motion). Absent any other information, these seem like fine null assumptions, so there is not much we can do. But if you wish, Machuruku can take some of the uncertainty inherent in ancestral character estimation into account when constructing niche models in the final step. This is actually a recommended approach for our program, as it'll lead to less stringent and more interpretable models. You can activate this feature by using the `n.unc` argument and specifying a number:
```
ace.M2.unc <- machu.2.ace(resp, bassleritree, T = 3.3, n.unc = 4)
```
Note that here we're taking a time-slice corresponding for our M2 paleoclimate data. What the `n.unc = 4` specification does is take `n.unc` evenly-spaced samples of the 95% confidence intervals surrounding the median reconstructed value. If you look at `ace.M2.unc`, you can see that there are a bunch of new columns, each corresponding to one of these samples:
```
$bio_1_mean
                name branch_start branch_end    value     unc1     unc2     unc3     unc4
1        Node5-Node6            5          6 232.3266 220.9530 228.5354 236.1177 243.7001
2 Node5-silverstonei            5          3 197.8689 190.8343 195.5240 200.2138 204.9035

$bio_1_stdev
                name branch_start branch_end    value     unc1     unc2     unc3     unc4
1        Node5-Node6            5          6 14.62110 12.92806 14.05675 15.18544 16.31413
2 Node5-silverstonei            5          3 13.98351 12.93633 13.63445 14.33256 15.03068
...
```
There is still the `value` column, with the median value, and after that, four (`n.unc`) columns of samples from the 95% confidence interval distribution. The first of these columns (`unc1`) is essentially the lower 95% confidence limit, and the last (`unc4`) is the upper limit. `unc2` and `unc4` will be evenly spaced between the limits. Imagine that each of these columns represents an independent niche model, at a certain part of the distribution of possible models. Basically, in step 3, you'd be constructing 4 (in this case) independent niche models, then summing and rescaling them. This accounts for uncertainty much better than simply using the median.
#### Accounting for uncertainty in divergence time estimation given one tree with error bars
Another source of uncertainty in Machuruku comes from our time-calibrated tree itself. In a phylogeny, there are going to be two main areas where uncertainty occurs: the topology, and the branch lengths. We can't do a lot about the topology here - if you aren't certain about the topology of your tree, just run Machuruku on two or more trees with slightly different topologies. For uncertainty in branch lengths, which for a time-calibrated tree corresponds to divergence times, we can do something, though.

Many (if not all) programs used for divergence time estimation will provide each node with error bars that usually describe the 95% probability distribution of ages for that node. Our *bassleri* group tree was calibrated using [BEAST 2](https://www.beast2.org/), a popular Bayesian phylogenetics program. It comes with error bars called 95% Highest Posterior Density (HPD) intervals. (Maximum likelihood programs like RelTime will give you 95% confidence intervals.) However, when we load our tree with the `read.nexus()` function from `ape`, we lose that error bar information when the tree is converted to a Phylo object. So, we need to use another function: `read.beast()` from the [`treeio`](https://bioconductor.org/packages/release/bioc/html/treeio.html) package.
```
# install.packages("treeio")
# library(treeio)
# you may have to install and load treeio separately
# load tree as treedata object
bassleribeast <- read.beast("basslerigroup.treefile")
```
The tree is now stored as a Treedata object, which contains all of the information regarding node heights (and more) from the original file. 
> To extract this information for your own purposes, use the construct `bassleribeast@data$height` (for example) or convert the tree to a "Tibble" object using the `as_tibble()` function.

> Currently I have only ensured that `machu.tree.unc()` only works with output from BEAST 2. Because of the way RelTime stores its node height data, it will not work with this function absent some significant changes.

Now let's run the `machu.tree.unc()` function, which will take our tree and produce two additional trees using the 95% HPD intervals at each node: one tree where the node heights correspond to the upper 95% HPD limit, and another where they correspond to the lower. Then, let's use `machu.treeplot()` to visualize all three trees at once:
```
trees <- machu.tree.unc(bassleribeast)
machu.treeplot(trees, upperX = 3, timelabeloffset = 1.25, nodecirclecolor = "darkseagreen2")
```
![3 trees](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/3trees.png?raw=true)

Three trees have been plotted on the same time scale. The top tree is our "lower confidence interval" (LCI) tree, the middle tree is our median starting tree, and the bottom tree is our "upper confidence interval" (UCI) tree. The node heights of the LCI and UCI trees have been pushed backwards or forwards to match the corresponding confidence limits of the median tree. The three trees are stored together in a list (for us, called `trees`) that contains them as its elements: `LCI tree`, `input tree`, and `UCI tree`.

So what's the point of this? We can use a single time-slice to illustrate why `machu.tree.unc()` might be useful.
```
machu.treeplot(trees, upperX = 3, timelabeloffset = 1.25, nodecirclecolor = "darkseagreen2", timeslice = 2.49)
```
![3 trees w/ time-slice](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/3trees%20w%20timeslice.png?raw=true)

If you are time-slicing the phylogeny in a certain place, you may recover different branch-taxa depending on the heights of the nodes. In the LCI tree (top), we are getting two taxa: Node5-Node6 and Node5-silverstonei. In the median input tree, we are getting three: Node6-Node7, Node6-pepperi, and Node5-silverstonei. And in the UCI tree (bottom), we get four: Node7-yoshina, Node7-bassleri, Node6-pepperi, and Node5-silverstonei. These will dramatically change your ancestral niche modeling results, because you're literally making niche models for a different number of taxa. The idea is to construct these additional trees, and then run Machuruku for them three separate times, and compare the results. That would look something like this:
```
ace.LCI <- machu.2.ace(resp, as.phylo(trees$`LCI tree`), T = 2.49)
ace.input <- machu.2.ace(resp, as.phylo(trees$`input tree`), T = 2.49)
ace.UCI <- machu.2.ace(resp, as.phylo(trees$`UCI tree`), T = 2.49)
```
Note that we need to use the `as.phylo` function to convert the treedata objects stored in the `trees` list to get them to work properly. If you look at each of the three objects we made, you'll see that `ace.LCI` has two taxa, `ace.input` has three, and `ace.UCI` has four - just as we predicted.
> You may be wondering why we characterize divergence time uncertainty this way and not by utilizing the whole posterior distribution of trees, or a "tree cloud", where nearly every possible tree within the 95% HPD intervals is represented. The short answer is that this is perhaps a good idea, but we don't think it particularly useful because of our emphasis on visualization. We can't really think of a way to visualize a "cloud" of reconstructed niches that is easy to interpret, especially since for a particular time-slice you'd have all manner of taxa present in different trees.

#### Accounting for uncertainty in divergence time estimation given a Bayesian posterior distribution of trees
Very astute readers may notice that the above method for accounting for divergence time uncertainty isn't perfect, because the two additional trees constructed by the `machu.tree.unc` function are *chimeric,* meaning their node-height information comes from multiple trees that were ultimately summarized by TreeAnnotator when we turned the Bayesian posterior distribution of trees into a single summary tree (that we then imported into R with `read.beast`). Because of this, the additional trees are not ultrametric, e.g. the different paths from root to tip are not the same length (though they are close); we solve this problem using the `force.ultrametric` function from [`phytools`](http://blog.phytools.org/). 

If this bugs you, there is yet another function for accounting for uncertainty in divergence time estimation: `machu.trees.unc`. It's called *trees* as opposed to *tree* because this time the input is an entire posterior distribution of trees, as is commonly output from most Bayesian phylogenetics programs (for this reason, `machu.trees.unc` should work with output from programs other than BEAST, as well). Because these posterior distributions can contain thousands or even millions of trees, using them in R can be next to impossible. For this reason, `machu.trees.unc` works *without* loading the trees into R; you just call the file location instead:
```
trees <- machu.trees.unc("posterior.trees")
```
Similarly to `machu.tree.unc`, this function will output a list of three trees, only this time the trees are actually sampled directly from the posterior, contained in the example file `posterior.trees`. This file is called simply as a text string and contains a very small "posterior distribution" of 19 trees from our [Guillory et al. 2020](https://www.sciencedirect.com/science/article/pii/S1055790319304609) *Ameerega* paper, with 35 taxa in each tree (1 per *Ameerega* species). Using `machu.treeplot` again, we can visualize the trees like so:
```
machu.treeplot(trees, upperX = 7, tiplabelsize = 2.5,
               timelabeloffset = 1.01, timelabelsize = 2,
               nodecirclecolor = "darkseagreen2", nodecirclesize = 5, nodelabelsize = 3)
```
![3 trees from machu.trees.unc](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/machu%20trees.png?raw=true)

The function works by loading each tree in the file one at a time in a series of three passes. The first pass simply counts how many trees are present, and then calculates the number to be skipped at a desired burn-in level (by default this is set to 10%, rounded up; in our case, with 19 trees, the first 2 trees were skipped). 
> Burn-in refers to the phenomenon where generally a chain of posterior samples takes a while to stabilize around a certain value, so when summarizing the posterior (as in the program TreeAnnotator, etc.), the first X% of samples are generally ignored. In `machu.trees.unc`, the exact burn-in percentage ignored can be specified with the `burnin` parameter, by default set to 0.1 (=10%).

The second pass then calculates the distribution of tree-heights in the posterior (i.e., it calculates tree-height for every tree), and then calculates the 95% HPD interval and mean of this tree-height distribution. The third and final pass identifies those trees with tree-heights closest to the upper and lower 95% HPD limits and average, and puts them into the output list (`trees` in our case).
> Much like the burn-in percentage, the HPD confidence limit percentage can also be manipulated using the `conf` parameter, which is by default set to 0.95 (=95%). To identify the most extreme trees (min and max heights), you would simply set `conf = 1`. However, these are likely to be total outliers and not very useful for your project.

#### Saving and loading ace output
We also provide a function that allows you to save your `machu.2.ace()` output to an external file and reload it for later. If you use the `savename` parameter in the `machu.2.ace()` function:
```
# save machu.2.ace() output
ace.all <- machu.2.ace(resp, bassleritree, savename = "ace_all")
```
Machuruku will save the output as a comma-separated value file called `ace_all.csv`. This can be opened in Excel or whatever program you like to use, if you want. If you load it back into R with `read.csv`, it looks like this:
```
    X         name      time NA.        value       climvar
1   1     bassleri  0.000000  NA 2.281429e+02    bio_1_mean
2   2      pepperi  0.000000  NA 2.433333e+02    bio_1_mean
3   3 silverstonei  0.000000  NA 1.910000e+02    bio_1_mean
4   4      yoshina  0.000000  NA 2.341138e+02    bio_1_mean
5   5        Node5 11.577178  NA 2.150910e+02    bio_1_mean
6   6        Node6  2.710415  NA 2.335421e+02    bio_1_mean
7   7        Node7  2.267554  NA 2.328640e+02    bio_1_mean
8   8     bassleri  0.000000  NA 1.687901e+01   bio_1_stdev
9   9      pepperi  0.000000  NA 1.382269e+01   bio_1_stdev
10 10 silverstonei  0.000000  NA 1.385641e+01   bio_1_stdev
11 11      yoshina  0.000000  NA 1.351322e+01   bio_1_stdev
12 12        Node5 11.577178  NA 1.430706e+01   bio_1_stdev
13 13        Node6  2.710415  NA 1.465221e+01   bio_1_stdev
14 14        Node7  2.267554  NA 1.480499e+01   bio_1_stdev
...
```
(I've only printed the first 14 lines.) As you can see, the structure of the .csv is quite different than the `ace.all` object, which was a list where each element was a table with the values of each separate climate response parameter for each taxon. Here, it's been reformatted so that the `climvar` column takes the place of separating each parameter into its own element of a list. The issue is that, for this to actually be useful, we need to be able to load the output back into Machuruku. We can use the function `machu.ace.load()` here:
```
# load ace_all.csv
ace.all <- machu.ace.load("ace_all.csv")
```
If you look at `ace.all` here, you'll see it's identical to the original and formatted correctly as a list.
### Projecting ancestral models into paleoclimatic data
The final major step in Machuruku is to use `machu.3.anc.niche()` to project the ancestral niche models estimated and parsed by `machu.2.ace()` into paleoclimate data. This step loops through each taxon in the `machu.2.ace()` output and constructs a Bioclim model in the provided paleoclimate data. Here I'm projecting our `ace.M2` models, consisting of two taxa interpolated to 3.3 Ma, into the M2 climate data:
```
mod <- machu.3.anc.niche(ace.M2, ClimM2, verbose = T)
```
This step is the most time-consuming part of Machuruku, though it shouldn't take too long in this instance since we only have two taxa and two climate variables. As you add taxa and climate variables, it'll begin to take longer.

If we look at `mod`, we see that it is a list where each element is a Raster layer corresponding to a niche model for each taxon:
```
$`Node5-Node6`
class      : RasterLayer 
dimensions : 328, 253, 82984  (nrow, ncol, ncell)
resolution : 0.04166667, 0.04166667  (x, y)
extent     : -81.49866, -70.95699, -14.0478, -0.3811359  (xmin, xmax, ymin, ymax)
crs        : +proj=longlat +datum=WGS84 +no_defs 
source     : memory
names      : layer 
values     : 0, 0.9999929  (min, max)

$`Node5-silverstonei`
class      : RasterLayer 
dimensions : 328, 253, 82984  (nrow, ncol, ncell)
resolution : 0.04166667, 0.04166667  (x, y)
extent     : -81.49866, -70.95699, -14.0478, -0.3811359  (xmin, xmax, ymin, ymax)
crs        : +proj=longlat +datum=WGS84 +no_defs 
source     : memory
names      : layer 
values     : 0, 0.9939741  (min, max)
```
We can use the accessory function `machu.plotmap()` to actually show what these models look like, which is the whole point of this process anyway:
```
par(mfrow=c(1,2))
machu.plotmap(mod)
```
![2 models w rarefied occ data](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/2%20mods%20rare.png?raw=true)

Note that these are nearly identical to the models I constructed in the Quick-start Guide above; the only differences are due to the fact that I used spatially rarefied occurrence data here, whereas in the Quick-start Guide I did not.
> You can specify different color ramps using the `col` argument in this function and specifying a number (1:6). 1 = base, 2 = plasma, 3 = viridis, 4 = l17, 5 = r3, 6 = white-to-black. Default = 1.

In the case where we ran `machu.2.ace()` without specifying a time-slice, getting `ace.all`, we can project both the extant taxa and the three nodes into the same climate layers to see how they directly compare in identical climate space:
```
# build present-day models for all taxa including ancestors
mod.all <- machu.3.anc.niche(ace.all, ClimCur, verbose = T)
par(mfrow=c(2,4))
machu.plotmap(mod.all)
plot(bassleritree)
```
![7 models in current climate w/ tree](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/7%20mod.png?raw=true)

I've used a plot of the tree to fill up the eighth panel. Remember, this is all in modern climate, so the three nodal taxa would not necessarily have lived in these areas, but perhaps they would today if they were still extant. (Of course, these are terrible models since we're only using two climate variables in this tutorial; the only even somewhat accurate one is the one for *A. silverstonei*, though of course the ancestral models are speculative). 

If you're only interested in a few taxa at a time, you can use the `taxa` argument to specify which ones:
```
# build present-day models for only extant taxa
mod.extant <- machu.3.anc.niche(ace.all, ClimCur, verbose = T, taxa = 1:4)
par(mfrow=c(1,4))
machu.plotmap(mod.extant)
```
![4 models of extant spp only](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/4%20mod%20extant.png?raw=true)

#### Clipping response curve tails to produce cleaner models
The `machu.3.anc.niche()` function includes a variety of additional arguments that allow you to further explore your models. One is the `clip.Q` argument, which "trims" or "clips" the tails of the response curves by setting the portions beyond the reconstructed `lowerQ` and `upperQ` parameters equal to zero. By default, the `clip.Q` argument is actually turned on, because it tends to produce cleaner models. However, to fully characterize the tails of these distributions, we can turn `clip.Q` off:
```
# turn clip.Q off
mod.clip.Q <- machu.3.anc.niche(ace.M2, ClimM2, verbose = T, clip.Q = F)
par(mfrow=c(2,2))
machu.plotmap(mod)
machu.plotmap(mod.clip.Q)
```
![4 models clip Q off bottom](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/4%20mod%20clipQ.png?raw=true)

The top two models are from the `mod` object we previously made, where `clip.Q = TRUE`. The bottom two models are from `mod.clip.Q`, where `clip.Q = FALSE`. You'll notice that it captures quite a bit more area, albeit of lower suitability.
#### Creating binary models
Another option is the `resp.curv` argument, which tells Machuruku that you want to use response curve information to characterize the niche. If you set this argument to FALSE, Machuruku instead sets all values outside of the bounds defined by `lowerQ` and `upperQ` for each climate response curve equal to 0, and all values inside of them to 1 - effectively creating a uniform distribution defined by `lowerQ` and `upperQ`. This is closer in nature to the original way Bioclim worked. In cases where you want a binary suitability surface, it may be useful.
```
# use resp.curv = F
mod.resp.curv <- machu.3.anc.niche(ace.M2, ClimM2, verbose = T, resp.curv = F)
par(mfrow=c(2,2))
machu.plotmap(mod)
machu.plotmap(mod.resp.curv)
```
![4 models resp curv on](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/4%20mod%20respcurv.png?raw=true)

Once again, the top two models are from `mod`, with `resp.curv = TRUE` (default). The bottom two models are with `resp.curv = FALSE`. As you can see, they capture a similar area, but all of it has the same (maximum) suitability. This is more like saying the taxon is equally likely to inhabit any of these pixels.
> Note that `resp.curv` and `clip.Q` can never be turned on at the same time because they're redundant; `resp.curv` already clips the tails since it uses `upperQ` and `lowerQ` as the bounds for the uniform distribution. `clip.Q` simply does not set the intervening values equal to one like `resp.curv`.

#### Incorporating uncertainty from ancestral character estimation
You might remember the `n.unc` argument from `machu.2.ace()`, that we used to sample the distribution of possible ancestral character estimation values for each climate response parameter, and then sort them into additional niche models. Here is where we incorporate that information. Setting `n.unc` to TRUE (it's FALSE by default) will essentially construct separate models for all of the samples you took (so, five in our case), then sum them, and rescale them so that the maximum value is still 1. 
```
# use calc.unc = T
mod.calc.unc <- machu.3.anc.niche(ace.M2.unc, ClimM2, verbose = T, calc.unc = T)
par(mfrow=c(2,2))
machu.plotmap(mod)
machu.plotmap(mod.calc.unc)
```
![4 models calc unc on](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/4%20mod%20calcunc.png?raw=true)

Again, the top two models are from `mod`, with `calc.unc = FALSE`, and the bottom two models are from `mod.calc.unc`, with `calc.unc = TRUE` (also notice that we used `ace.M2.unc` as our input, not `ace.M2`). The models incorporating uncertainty are actually a bit narrower in terms of the highly suitable habitat, but should characterize the total variability in suitability better than the other models.

#### Limiting model area with inverse-distance weighting
Often, niche models will recover areas obviously outside the known range of a species as suitable. While this can be useful, depending on the question it may be superfluous or misleading. To that end, we added a function called `machu.geo.idw` that uses inverse-distance weighting to restrict niche models to the areas surrounding the taxon's distribution. Essentially, this is a way for accounting for migration limitations by removing from consideration suitable areas that are nonetheless separated from the taxon's actual distribution by geographic barriers or other factors. An obvious hurdle for this function is that *we don't know* what the actual distributions of ancestral taxa were; figuring that out is the whole point of Machuruku in the first place. To (sort of) get around this, the user just specifies an extant taxon distribution (as occurrence data) to represent the ancestral taxon. This is obviously an imperfect system, but it's the best that can be done absent voluminous fossil data. The function is used as such:
```
clip <- machu.geo.idw(mod[[2]], occ, taxa = "silverstonei", buffer.dist = 100, kernel.size = 2, MCP.percent = 50)
```
> Note that here we're using `occ` rather than `occ.rarefied` because the latter only contains 3 occurrence points for *A. silverstonei* and the function needs at least 5 to work.

In its current implementation, only one taxon can be run through `machu.geo.idw` at once. This is because the user will likely want to specify different occurrence data and buffer distances, etc., for each taxon. In this example, we are using `mod[[2]]` as our taxon to clip, which represents the ancestral taxon "node5-silverstonei" projected into M2 climate data. We specify our `occ.rarefied` occurrence data, and further specify *A. silverstonei* as the extant taxon whose data we wish to use. 
> If you want to specify a *combination* of extant taxa to use, you can do that too. For instance, if you're clipping `mod[[1]]`, the "node5-node6" taxon (ancestor to the *bassleri* group, you can specify `taxa=c("pepperi","bassleri","yoshina")` to include the entire distribution of the *bassleri* group as a proxy for the distribution of its MRCA. 

What this function does is basically construct a minimum convex polygon (MCP) around the input points, whose size is determined by the `buffer.dist` parameter, which is by default set to 300 kilometers. It then creates a second buffer around the MCP that smooths the model's probability from 1 (adjacent to the MCP) to 0 at the edge of this second buffer. The second smooth buffer's size is determined by the `kernel.size` parameter (default = 1), which is multiplied against the `buffer.dist` parameter to calculate the size; i.e., if `buffer.dist = 100` and `kernel.size = 2`, the second buffer will be 200 km wide. The `mcp.percent` (default = 100) determines the percentage of outlier occurrence points to be removed before construction of the original MCP. 

To visualize what `machu.geo.idw` does, do the following:
```
par(mfrow=c(1,2))
plot(mod[[2]])
plot(clip)
```
![normal model vs. linear-distance weighted model](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/idw.png?raw=true)

As you can see, the clipped model (right) is now restricted to the area surrounding the extant *A. silverstonei* distribution. The left model is the original node5-silverstonei model in M2 climate. 
