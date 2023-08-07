# Machuruku: the Tutorial 2.0
This is the updated tutorial coinciding with the release of Machuruku 2.0. I originally developed Machuruku ~2019 with [Jason Brown](https://www.jasonleebrown.org/), my MS advisor (2017-2020). After the release of the original version, we published [a paper](https://academic.oup.com/sysbio/article/70/5/1033/6171196) in Systematic Biology in 2021 describing the method. Machuruku was my first major R project and as such the original version bore some of the hallmarks of a novice R coder, along with some other quirks that continued to bother me as I moved on into a doctoral program. A couple years later in 2023, with much more coding experience under my belt, I decided to revisit Machuruku which quickly ballooned into rewriting the entire thing, both in service to my own projects and the community at large. The scope of the changes merits a nomenclatural leap to Machuruku version 2.0.

The purpose of Machuruku remains unchanged: the package is for modeling, reconstructing, and visualizing niches in past and present climates, while accounting for evolution. Machuruku reconstructs climatic niches along a time-calibrated phylogeny and can project these niches into paleoclimatic data to visualize the potential geographic origins of lineages. The advances provided by Machuruku 2.0 are mostly practical. For example, many of the for-loops I was previously reliant on have been replaced with apply statements, which offer significant boosts in computation speed. Also boosting speed is Machuruku's new integration with the [Terra](https://github.com/rspatial/terra) R package, a faster upgrade to its predecessor Raster. The workings of the various functions provided in the package have been streamlined to be more user-friendly and flexible. There is also a potentially significant conceptual change to the way Machuruku calculates and visualizes uncertainty in ancestral niche reconstructions, which I will discuss in greater detail later. 

In this new tutorial, I will first provide a simple quick-start guide to using Machuruku 2.0 at its most basic. I then provide a deeper step-by-step dive into the each function and the possibilities they offer for reconstructing and visualizing ancestral niches. 
## Contents
* [Introduction](https://github.com/wxguillo/machuruku/tree/main/tutorial#introduction)
  * [Installing Machuruku](https://github.com/wxguillo/machuruku/tree/main/tutorial#installing-machuruku)
  * [Downloading and exploring tutorial data](https://github.com/wxguillo/machuruku/tree/main/tutorial#downloading-and-exploring-tutorial-data)
    * [Loading tree](https://github.com/wxguillo/machuruku/tree/main/tutorial#loading-tree)
    * [Loading occurrence data](https://github.com/wxguillo/machuruku/tree/main/tutorial#loading-occurrence-data)
    * [Loading climate data](https://github.com/wxguillo/machuruku/tree/main/tutorial#loading-climate-data)
* [Quick-start Guide](https://github.com/wxguillo/machuruku/tree/main/tutorial#quick-start-guide)
  * [Estimating tip response curves](https://github.com/wxguillo/machuruku/tree/main/tutorial#estimating-tip-response-curves)
  * [Estimating ancestral niches at a time-slice](https://github.com/wxguillo/machuruku/tree/main/tutorial#estimating-ancestral-niches-at-a-time-slice)
  * [Projecting ancestral models into paleoclimatic data](https://github.com/wxguillo/machuruku/tree/main/tutorial#projecting-ancestral-models-into-paleoclimatic-data)
* [Detailed Guide](https://github.com/wxguillo/machuruku/tree/main/tutorial#detailed-guide)
  * [1. Estimating tip response curves](https://github.com/wxguillo/machuruku/tree/main/tutorial#1-estimating-tip-response-curves)
    * [Accounting for spatial autocorrelation by rarefying occurrence data](https://github.com/wxguillo/machuruku/tree/main/tutorial#accounting-for-spatial-autocorrelation-by-rarefying-occurrence-data)
    * [Selecting the most important climate variables](https://github.com/wxguillo/machuruku/tree/main/tutorial#selecting-the-most-important-climate-variables)
  * [2. Estimating ancestral niches](https://github.com/wxguillo/machuruku/tree/main/tutorial#2-estimating-ancestral-niches)
    * [Visualizing timeslices](https://github.com/wxguillo/machuruku/tree/main/tutorial#visualizing-timeslices)
    * [Estimating ancestral niches at a time-slice](https://github.com/wxguillo/machuruku/tree/main/tutorial#estimating-ancestral-niches-at-timeslices)
    * [Estimating ancestral niches at each node](https://github.com/wxguillo/machuruku/tree/main/tutorial#estimating-ancestral-niches-at-each-node)
    * [Characterizing uncertainty in ancestral niche estimation](https://github.com/wxguillo/machuruku/tree/main/tutorial#characterizing-uncertainty-in-ancestral-niche-estimation)
    * [Characterizing uncertainty in divergence times from a single tree](https://github.com/wxguillo/machuruku/tree/main/tutorial#characterizing-uncertainty-in-divergence-times-from-a-single-tree)
    * [Characterizing uncertainty in divergence times from a posterior distribution of trees](https://github.com/wxguillo/machuruku/tree/main/tutorial#characterizing-uncertainty-in-divergence-times-from-a-posterior-distribution-of-trees)
    * [Saving and loading ace output](https://github.com/wxguillo/machuruku/tree/main/tutorial#saving-and-loading-ace-output)
  * [3. Projecting ancestral niche models into paleoclimate](https://github.com/wxguillo/machuruku/tree/main/tutorial#3-projecting-ancestral-niche-models-into-paleoclimate)
    * [Projecting multiple timeslices into multiple paleoclimates](https://github.com/wxguillo/machuruku/tree/main/tutorial#projecting-multiple-timeslices-into-multiple-paleoclimates)
    * [Projecting into paleoclimate while accounting for uncertainty](https://github.com/wxguillo/machuruku/tree/main/tutorial#projecting-into-paleoclimate-while-accounting-for-uncertainty)
    * [Clipping niche models](https://github.com/wxguillo/machuruku/tree/main/tutorial#clipping-niche-models)
    * [Making binary niche models](https://github.com/wxguillo/machuruku/tree/main/tutorial#making-binary-niche-models)
    * [Projecting a single scenario into multiple paleoclimates](https://github.com/wxguillo/machuruku/tree/main/tutorial#projecting-a-single-scenario-into-multiple-paleoclimates)
    * [Projecting multiple scenarios into a single paleoclimate](https://github.com/wxguillo/machuruku/tree/main/tutorial#projecting-multiple-scenarios-into-a-single-paleoclimate)
    * [Saving niche models to a folder](https://github.com/wxguillo/machuruku/tree/main/tutorial#saving-niche-models-to-a-folder)
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
* `climate.zip` - A zipped file containing the climate data we will be using. If you unzip the file, you will get a folder called `climate/`, which in turn contains five subfolders: `current/`, `mis19/`, `mpwp/`, and `m2/`. The subfolders correspond to different time periods (present, Marine Isotope Stage 19 (Pleistocene; 0.787 Ma), Mid-Pliocene Warming Period (3.205 Ma), and Marine Isotope Stage M2 (Pliocene; 3.3 Ma). Each contains 14 raster climate layers from [Paleoclim](http://www.paleoclim.org/), which I've cropped to an area surrounding the *bassleri* group's distribution and reduced in resolution to save space. 

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
![simpletree](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/simpletree.png?raw=true)

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
legend(-81, -11, legend = taxa, col = cols, pch=19, bty="n")
```
![bio1 w/ points](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/bio1%20w%20pts.png?raw=true)

Now we can see that the species' ranges do not really overlap, but are very close to each other. They evolved in a place with lots of topographic heterogeneity (the foothills of the Andes), which probably promoted their diversification. We might be able to see some niche evolution in this group since they occupy subtly different habitats.
## Quick-start Guide
Machuruku is composed of three principal steps. In this section I'll demonstrate them without any extraneous bells-and-whistles. Let's set a phylogenetic niche modeling goal for ourselves: To see what the ancestors of the *bassleri* group and *A. silverstonei* were doing in the Late Pliocene, approximately 3.3 million years ago (the putative age of our M2 climate layers). What's important to note here is that it's highly unlikely that any of our nodes will exactly line up with the age of our paleoclimate data, e.g., none of the nodes will be 3.3 million years old. To accommodate this, Machuruku can interpolate ancestral niches along the branches subtending the desired time-slice. This can be kind of tricky to visualize, so Machuruku has a built-in function, `machu.treeplot()`, that can do it for us.
```
# visualize tree with timeslice at 3.3 Ma
machu.treeplot(tree, timeslice=3.3)
```
![treeplot 3.3](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/treeplot.png?raw=true)

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
![2 models](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/2%20models.png?raw=true)

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
![respplot1](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/respplot1.png?raw=true)

The `machu.respplot` function (response-plot) can plot response curves in various combinations. Here I specify the first two climate variables in `resp` by name, as well as to plot the curves "together" in one window. The curves are nearly identical because the bio1 (mean annual temp) and bio10 (mean temp of warmest quarter) bioclimatic variables are very similar. To view the response curves for *every* climate variable, we can simply omit the `clim` argument:
```
# visualize response curves for all variables
machu.respplot(resp, plot="t", legend.cex=0.6)
```
![respplot2](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/respplot2.png?raw=true)

Another feature of `machu.1.anc.niche`, new to Machuruku 2.0, is the visualization of present-day niche models with `dismo::predict`. Activate this feature with the `plot` argument:
```
# visualize niche models
machu.1.tip.resp(occ, current, plot="t")
```
![currentmodels1](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/currentmodels1.png?raw=true)

This is similar to `machu.respplot` where `plot="t"` (equivalent to `plot="together"`) will show every map in one window. You can also add the occurrence points to confirm that the estimated niche model lines up with the known distribution of the species:
```
dev.off()
machu.1.tip.resp(occ, current, plot="s", plot.points=T)
```
![currentmodels2](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/currentmodels2.png?raw=true)

Here `dev.off()` is used to reset the plot window. The `plot="s"` (equivalent to `plot="separately"`) option shows every map in its own window; here I'm only showing the first one for *A. bassleri*. 

Finally, the function can also output these niche models instead of the response table with the `output.bioclim` option. This can be useful if you'd like to plot, save, or analyze the present-day niche models for each species on your own. 
```
mod <- machu.1.tip.resp(occ, current, output.bioclim=T)
par(mfrow=c(2,2), mar=c(0,0,2,0))
lapply(1:4, function(x) plot(mod[[x]], axes=F, legend=F, box=F, main=names(mod)[x]))
```
![currentmodels3](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/currentmodels3.png?raw=true)

The output `mod` is a list where each element is a RasterLayer; the `lapply` function loops through each one and plots it after setting up a 2x2 plotting window with `par`.
#### Accounting for spatial autocorrelation by rarefying occurrence data
Occurrence data are often spatially autocorrelated due to systematic biases in the way we collect specimens, i.e. they tend to be clustered in easily accessible locations. Spatial autocorrelation can induce downstream impacts on niche modeling, necessitating the use of algorithms to "rarefy", or thin, occurrence data beforehand. In Machuruku, the function `machu.occ.rarefy` can do this for you. 
```
# rarefy occurrence data
dev.off()
occ.r <- machu.occ.rarefy(occ, rarefy.dist=10, plot=T)
```
The `rarefy.dist` option specifies the amount of rarefication to perform; in my example, I've set it to 10 km (the default unit; decimal degrees are also available with `rarefy.units="dd"`) such that all points must be at least 10 km from any other point. Adjust this value to fit the spatial scale of your data. Also, the `plot=T` option (new for Machuruku 2.0) visualizes the extent of rarefication performed:
![rarefied](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/rarefied.png?raw=true)

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
![currentmodels4](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/currentmodels4.png?raw=true)

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
![currentmodels5](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/currentmodels5.png?raw=true)

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

![treeplot1](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/treeplot1.png?raw=true)

This looks pretty good, but there are a lot of graphical options in `machu.treeplot` that can be tweaked to make the figure look better. 
```
machu.treeplot(tree, timeslice = c(0.787,3.205,3.3), x.l.lim = 13, x.u.lim=-3, nodelabsize = 0.35, col = "skyblue", timelaboffset = -0.2)
```
Here I manually changed the vertical size of the window to make the tree more compact, and adjusted the horizontal scaling with the `x.l.lim` and `x.u.lim` parameters to better fill the space. This necessitated adjusting the size of the circles at each node with `nodelabsize`, and adjusting the position of the time labels with `timelaboffset`. Finally, I changed the color of the timeslice lines to "skyblue" with `col`, just for fun.

![treeplot2](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/treeplot2.png?raw=true)

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
![respplot3](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/respplot3.png?raw=true)

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
![respplot4](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/respplot4.png?raw=true)

When uncertainty values are passed to this function, it automatically plots the 27 skew-normal distributions formed by the different parameter combinations. There isn't a ton of variability in this case, so the default `fill=F` actually makes it easier to see. Node1-*silverstonei* is very constrained, while Node3-*bassleri* and Node3-*yoshina* have more variation. Thick dotted lines are drawn to show the median distribution for each taxon that would be used were `unc=F`. 
#### Characterizing uncertainty in divergence times from a single tree
Another source of uncertainty in the Machuruku pipeline comes from the time-calibrated phylogeny itself. Just like estimates of a character value at each node, the timing of the node itself is uncertain, more so for older nodes. Dealing with this in Machuruku is, perhaps, "not my job", but I wrote a function to do it anyway. `machu.tree.unc` takes a time-calibrated phylogeny, or a Bayesian posterior distribution of phylogenies, and produces two additional phylogenies, one in general older than the input, one in general younger, that characterize the uncertainty in divergence times and can be used as separate inputs for `machu.2.ace`. 

The `machu.tree.unc` function has two modes depending on the input. If the input is a single tree, it produces the two additional trees from the divergence time uncertainty contained in the tree file (i.e., the error bars). If the input is a posterior distribution of trees, it produces the two additional trees from the corresponding distribution of tree-heights; more detail on this mode later. For the former, single-tree mode, a different type of input is necessary than what we've been using, which is the "phylo" format from Ape. Phylo format does not record divergence time uncertainty data, even if the read-in input file contains it. Instead, we have to use the "treedata" format from the package [Treeio](https://bioconductor.org/packages/release/bioc/html/treeio.html), which does record this information in "tibble" format (part of the Tidyverse world).
```
# characterize divergence time uncertainty
# reload tree as 'treedata' object
library(treeio)
beast <- read.beast("basslerigroup.treefile")
```
The `treeio::read.beast` function creates a treedata object. Running `beast` will display the tibble containing the divergence time uncertainty data; they are stored in the "`height_0.95_HPD`" column, which we can view in compact form by running `do.call(rbind, as_tibble(beast)$height_0.95_HPD)`:
```
         [,1]         [,2]
[1,] 0.000000 7.105427e-15
[2,] 0.000000 7.105427e-15
[3,] 0.000000 7.105427e-15
[4,] 0.000000 7.105427e-15
[5,] 7.112935 1.620734e+01
[6,] 1.312191 4.255582e+00
[7,] 1.030775 3.665719e+00
```
Each row in this makeshift table represents the upper and lower 95% HPD (highest posterior density; analogous to a confidence interval) limits for each tip and node in the tree. The first four rows represent the tips; they are essentially zero since it doesn't make sense for present-day taxa to have any divergence time uncertainty. The last three rows represent the nodes; the left column is the lower 95% HPD limit, and the right column is the upper 95% HPD limit. `machu.tree.unc` uses these values to construct the additional trees, essentially setting the node heights of the trees equal to these 95% HPD limits.

Note that this function is only tested and designed for outputs from a few programs: [BEAST](https://www.beast2.org/) and BEAST-adjacent software such as [SNAPP](https://www.beast2.org/snapp/) (specifically the consensus output produced by TreeAnnotator), and the RelTime method implemented in [MEGA](https://www.megasoftware.net/). Other programs may produce incompatible outputs; let me know and I can try to include them.
```
# run tree-uncertainty utility from a single input tree
beast.trees <- machu.tree.unc(beast)
```
Running `names(beast.trees)` shows the names of the three trees output by `machu.tree.unc`:
```
[1] "lCItree"   "inputtree" "uCItree"
```
The first tree, `"lCItree"` (lower confidence interval tree) has node heights based on the first column of the above table; the second tree, `"inputtree"` is the same as the input tree `beast` and represents the median node heights; the third tree `"uCItree"` (upper confidence interval tree) has node heights based on the second column of the above table. Running `beast.trees` itself will show the structure of the output:
```
$lCItree

Phylogenetic tree with 4 tips and 3 internal nodes.

Tip labels:
  bassleri, pepperi, silverstonei, yoshina

Rooted; includes branch lengths.

$inputtree

Phylogenetic tree with 4 tips and 3 internal nodes.

Tip labels:
  bassleri, pepperi, silverstonei, yoshina

Rooted; includes branch lengths.

$uCItree

Phylogenetic tree with 4 tips and 3 internal nodes.

Tip labels:
  bassleri, pepperi, silverstonei, yoshina

Rooted; includes branch lengths.
```
The output is a list where each element is one of the three trees discussed above. We can visualize these trees all at once using `machu.treeplot`:
```
# visualize trees
machu.treeplot(beast.trees, timeslice=c(0.787,3.205,3.3), nodelabsize=0.5, timelaboffset=-0.15, col="skyblue")
```
![3trees](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/3trees.png?raw=true)

This figure illustrates why one may want to use the different trees output from `machu.tree.unc` as input to `machu.2.ace`. When the trees have different node heights, the same timeslice can cut through different branches. In this figure, the two older timeslices (corresponding to the m2 and mpwp paleoclimate datasets) slice through two branch-taxa (Node1-Node2 and Node1-*silverstonei*) in the first two trees (the lCI and input trees). However, in the third tree (the uCI tree), they slice through four branch-taxa. If we re-run `machu.2.ace` with these new trees and take a single timeslice at 3.3 Ma, we see that the output `ace.uCItree` contains data for two branch-taxa, and the others contain data for four.
```
# re-run machu.2.ace with the new trees for comparison
ace.lCItree <- machu.2.ace(resp, beast.trees$lCItree, timeslice=3.3)
ace.inputtree <- machu.2.ace(resp, beast.trees$inputtree, timeslice=3.3)
ace.uCItree <- machu.2.ace(resp, beast.trees$uCItree, timeslice=3.3)
# count number of taxa in each of the new outputs
c(nrow(ace.lCItree[[1]]), nrow(ace.inputtree[[1]]), nrow(ace.uCItree[[1]]))
```
The output of the last line is `2 2 4`, which is what we expect.
#### Characterizing uncertainty in divergence times from a posterior distribution of trees
The other mode of `machu.tree.unc` is activated when the input is a posterior distribution of trees, i.e. the raw output of a Bayesian phylogenetics program like [BEAST](beast2.org), rather than a single consensus tree from TreeAnnotator as used previously. In the previous version of Machuruku, this mode was contained in the separate function `machu.trees.unc`; for Machuruku 2.0 I've combined them into a single function. Some users may want to use this mode instead of the previous one because it has greater flexibility and may characterize the divergence time uncertainty more accurately. Unlikely the previous mode, this mode actually returns trees that were directly computed by the original Bayesian software.

Because a posterior can contain thousands or even millions of trees, loading it all into R's memory would be a bad idea. Instead, simply specify by name a NEXUS file that contains the posterior distribution. You can download the file "posterior.trees" from the tutorial repository and examine the format to make sure the function will work with your data (the function is untested with other formats and almost certainly won't work). Essentially, each separate tree in the distribution is written on a separate line near the end of the file; the top of the file contains the list of taxa and a translation table giving each taxon name a numeric ID that appears in the trees to save space. The example file "posterior.trees" is a heavily truncated version of the output from our [2020 paper](https://www.sciencedirect.com/science/article/pii/S1055790319304609) on *Ameerega* phylogenetics, containing all of the taxa in the study but just 19 trees for the sake of computational speed. Place it in your tutorial folder and run the following:
```
# run tree-uncertainty utility from a posterior
post.trees <- machu.tree.unc("posterior.trees")
```
As `verbose=T` by default, the function will provide some progress updates as it works:
```
[1] "Multiple trees in file 'posterior.trees' detected. Burnin: 0.1; Confidence level: 0.95"
[1] "FIRST PASS: Calculating trees to remove at 10% burnin."
[1] "Found 19 trees, skipping the first 2 under 10% burnin."
[1] "SECOND PASS: Calculating 95% HPD of tree heights at 10% burnin."
[1] "95% HPD min: 182.44; Median: 282.65; 95% HPD max: 348.88; from 17 trees under 10% burnin."
[1] "THIRD PASS: Finding trees with closest heights to 95% HPD limits and median at 10% burnin."
```
This also describes how the algorithm actually works. The function conducts three passes through the dataset, loading each tree directly from the file one at a time to save memory. The first pass simply counts the number of trees in the dataset and calculates how many to remove based on the `burnin` parameter. Removing some trees from the beginning of the posterior is a common practice because the MCMC algorithms that Bayesian phylogenetics software use take a while to converge (i.e. "burn in") so the first x% of the distribution may be inconsistent with the rest; by default the `burnin` parameter is set to 0.1 (i.e., discarding the first 10% of the posterior). 

The function's second pass skips the first 10% of trees (when `burnin=0.1`) and calculates the tree height for each one (i.e., the distance from the root to the tip) to produce a distribution of tree heights. It then calculates the 95% HPD limits of this distribution; in our example, the lower 95% HPD limit was 182.44 and the upper one was 348.88. The confidence level can be set with the `conf` parameter, which by default is 0.95. Setting it lower will return trees closer to the median tree; on the other hand setting `conf=1` will identify the most extreme trees.

In the third pass, the function then identifies which trees have heights closest to the HPD limits and median, and returns these in a list similar to the other mode of `machu.tree.unc`. We can visualize these trees again with `machu.treeplot`:
```
# visualize trees
machu.treeplot(post.trees, timelabs=F, nodelabs=F)
```
![3bigtrees](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/3bigtrees.png?raw=true)

These trees are much larger than the ones we've been using (they contain all *Ameerega* species, not just the *bassleri* group), so I've set `timelabs=F` and `nodelabs=F` to reduce clutter. Similarly to before, the top tree is the lCI tree, the middle one is the median tree, and the bottom tree is the uCI tree.
#### Saving and loading ace output
Because the output of `machu.2.ace` in Machuruku 2.0 is now formatted differently, as a list, the output when saved to a CSV file is also different. You can save output from `machu.2.ace` simply by specifying a filename with the `csv.name` parameter (this was called `savename` before Machuruku 2.0):
```
# save ace output
ace.ts <- machu.2.ace(resp, tree, timeslice=c(0.787,3.205,3.3), csv.name="ace_ts.csv")
```
A CSV file called "ace_ts.csv" has been saved to your working directory (`getwd()`). We can examine its structure by viewing the first 7 columns with `read.csv("ace_ts.csv")[,1:7]`:
```
         scenario              taxon branch_start branch_end bio_4_mean bio_4_stdev bio_4_skew
1 timeslice_0.787     Node3-bassleri            7          1   488.8174    61.24433  45.763683
2 timeslice_0.787      Node3-yoshina            7          4   464.7665    45.17915 159.722456
3 timeslice_0.787      Node2-pepperi            6          2   561.1887    17.35599 165.743824
4 timeslice_0.787 Node1-silverstonei            5          3   602.4069    53.03126   8.241115
5 timeslice_3.205        Node1-Node2            5          6   514.5594    40.36691 119.605293
6 timeslice_3.205 Node1-silverstonei            5          3   591.3218    51.43320  22.293585
7   timeslice_3.3        Node1-Node2            5          6   514.9949    40.42969 119.053189
8   timeslice_3.3 Node1-silverstonei            5          3   590.8863    51.37041  22.845688
```
Each row in this table represents a single taxon. The first column lists the timeslice of each taxon (called a "scenario" for generality), and the second column lists the taxon name. The remaining columns contain the same data we saw before with the `machu.2.ace` output. Saving these files may be important for reproducibility and other custom analyses.

Loading the CSV file into R and trying to put it into the `machu.3.anc.niche` for downstream analysis will not work because it's formatted differently than the `machu.2.ace` output. To load a saved `machu.2.ace` output, use the function `machu.ace.load`:
```
# load ace output
ace.ts <- machu.ace.load("ace_ts.csv")
```
Examining `ace.ts` will show that the loaded file is now in its proper format and ready for downstream analysis.
### 3. Projecting ancestral niche models into paleoclimate
The final major step in Machuruku is to use `machu.3.anc.niche` to project the ancestral niche models estimated by `machu.2.ace` into paleoclimate data. For each taxon in each timeslice, the function calculates suitability in the timeslice's associated paleoclimate data by reconstructing a skew-normal distribution for each paleoclimate variable with the `sn::dsn` function. Various parameters and input schemes modify this basic blueprint in various ways. The major advance for Machuruku 2.0 is in allowing the input of multiple timeslices and paleoclimate layers at once, and in the new method of characterizing uncertainty in the ancestral niche.
#### Projecting multiple timeslices into multiple paleoclimates
Because it has been rewritten as a series of nested `lapply` statements, `machu.3.anc.niche` can now handle multiple timeslices and multiple paleoclimates as input. The main inputs are the output of `machu.2.ace` (ancestral climate response values) and paleoclimate raster layers. The paleoclimate layers can be provided in various formats, including both SpatRaster and RasterStack. When using paleoclimate from multiple time periods as input, it is best to combine them into a list with the `list` function. Ensure the datasets are listed in the same order as their corresponding timeslices in the `machu.2.ace` output, and then provide it as input to `machu.3.anc.niche`:
```
### 3. project ancestral niche models into paleoclimate
# multiple timeslices into multiple paleoclimates
clim <- list(mis19.reduced, mpwp.reduced, m2.reduced)
mod.ts <- machu.3.anc.niche(ace.ts, clim, verbose=T)
```
I've made sure to provide the reduced climate datasets, because they contain the same climate variables that were used in creating `ace.ts`. On my machine, the function takes only about a second to finish. Because `verbose=T`, it prints the following output:
```
[1] "Did not detect uncertainty samples in 'ace' ('unc' in machu.2.ace())."
[1] "Associating scenarios ('ace') with rastersets ('clim') (numbers are indices, rows are associations):"
     clim ace
[1,]    1   1
[2,]    2   2
[3,]    3   3
[1] "Processed scenario 1"
[1] "Processed scenario 2"
[1] "Processed scenario 3"
```
The table in the center of this output is the most important part: it is a sanity check describing the association of each paleoclimate dataset ("`clim`") with each timeslice ("`ace`"). Each row represents one such association, with each number being an index of the corresponding input (i.e., "`1`" in the "`ace`" column is the first element of `ace.ts`). With this table, we know that the first timeslice (`ace.ts[[1]]`) is being projected into the first paleoclimate dataset (`clim[[1]]`), and so on.

Running `mod.ts` shows the structure of the `machu.3.anc.niche` output. It is a list of lists, where each element corresponds to one timeslice, itself a list of ancestral niche models, one per taxon. Running `names(mod.ts)` shows `[1] "timeslice_0.787" "timeslice_3.205" "timeslice_3.3"`, while examining the names of the first element of `mod.ts` with `names(mod.ts[[1]])` yields `[1] "Node3-bassleri"     "Node3-yoshina"      "Node2-pepperi"      "Node1-silverstonei"`. From our knowledge of `ace.ts`, we'd expect `mod.ts` to have four niche models for the first timeslice (mis19), and two for the others (mpwp and m2). Running `sapply(mod.ts, length)` yields what we expect:
```
timeslice_0.787 timeslice_3.205   timeslice_3.3 
              4               2               2
```
We can visualize these niche models using the `machu.plotmap` function. 
```
machu.plotmap(mod.ts, plot="t", axes=F, title.cex=0.85, to.scale=T, plot.asp=20/9)
```
This works similarly to the other plotting functions in Machuruku 2.0 where the plots can be shown `"together"` or `"separately"`. I've set `axes=F` to reduce visual clutter, `title.cex=0.85` so that the titles fit in the window, and `plot.asp=20/9` to plot the maps in a compact horizontally-oriented format. The `to.scale` parameter identifies the maximum suitability value across all niche models and plots all of them on the same scale so that they can be directly compared. Here is the result:
![ancniche1](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/ancniche1.png?raw=true)

It turns out that with the settings and inputs we used, there is actually no suitable habitat for our ancestral taxa except in the oldest time period! Because the taxa are still extant, we know this can't be the case, so we can try being more conservative with our inputs.
#### Projecting into paleoclimate while accounting for uncertainty
Earlier, we accounted for uncertainty in the reconstruction of the ancestral niche parameters by setting `unc=T` in `machu.2.ace`. When passing a `machu.2.ace` output where `unc=T` to `machu.3.anc.niche`, the function automatically detects the presence of uncertainty and includes it in the projections into paleoclimate. This is done by taking all 3<sup>3</sup>=27 possible combinations of niche parameters for each climate variable, reconstructing a skew-normal distribution for each one, and converting to suitability across the climate raster to essentially produce 27 different niche models. These are then summed and rescaled to produce a single model that is generally broader than one not accounting for uncertainty. To run this analysis, we will use the `ace.ts.u` result from before as input:
```
# projecting while accounting for uncertainty
mod.ts.u <- machu.3.anc.niche(ace.ts.u, clim, verbose=T)
```
The function acknowledges the presence of uncertainty information when printing progress updates to screen:
```
[1] "Detected uncertainty samples in 'ace' ('unc' in machu.2.ace())."
[1] "Associating scenarios ('ace') with rastersets ('clim') (numbers are indices, rows are associations):"
     clim ace
[1,]    1   1
[2,]    2   2
[3,]    3   3
[1] "Processed scenario 1"
[1] "Processed scenario 2"
[1] "Processed scenario 3"
```
In this case the function took significantly longer to finish than without uncertainty because it is constructing 27x as many niche models. Let's view the results again:
```
machu.plotmap(mod.ts.u, plot="t", axes=F, title.cex=0.85, to.scale=T, plot.asp=20/9)
```
![ancniche2](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/ancniche2.png?raw=true)

This shows some improvement, as we now have some suitable area in the second timeslice (maps #5 and 6). However, the first timeslice still shows no suitable habitat. This is what happens when you choose a very range-limited taxon as the example for your tutorial. 

#### Clipping niche models
By default, the `machu.3.anc.niche` function clips the tails of the skew-normal distributions describing the relationships between each climate variable and suitability (i.e., setting suitability (y) below the lower confidence limit (x) to zero, and setting suitability above the upper confidence limit to the upper confidence limit). This produces "cleaner" models. This behavior is controlled by the `clip.Q` parameter, which is `TRUE` by default. The confidence limits to clip by are in turn set by the `clip.amt` parameter, which is set to 0.95 by default. 

In the previous version of Machuruku, the confidence limits to clip by were themselves climate response variables parsed in `machu.1.tip.resp` and passed through `machu.2.ace`. To simplify matters, these variables are omitted from Machuruku 2.0 and the confidence limits are calculated directly within `machu.3.anc.niche` by randomly sampling from the corresponding skew-normal distribution with `sn::rsn` and calculating the confidence limits from the sample with `quantile`. The number of samples taken by `rsn` is set with the parameter `clip.samples`, which is set by default to 10,000; decreasing this number may improve the speed of the function, while sacrificing accuracy (for most applications 10,000 will be fine).

In the following code block, I demonstrate the functionality of `clip.Q` by projecting ancestral niches into the m2 timeslice, first with no clipping, then clipping at the default level (95%), then performing very stringent clipping. 
```
# clipping niche models
mod.ts.u.nc <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], clip.Q=F) # no clipping
mod.ts.u.c95 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], clip.Q=T) # default clipping (95%)
mod.ts.u.c50 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], clip.Q=T, clip.amt=0.5) # stringent clipping (50%)
```
When subsetting the `machu.2.ace` output `ace.ts.u`, single brackets must be used so that the input is still in "list" format. Now let's combine these results into a single list with `c` and visualize them with `machu.plotmap`:
```
# combine and visualize
clip.demo <- c(mod.ts.u.nc, mod.ts.u.c95, mod.ts.u.c50)
names(clip.demo) <- c("no clipping", "clipping at 95%", "clipping at 50%")
machu.plotmap(clip.demo, plot="t", axes=F, to.scale=T, plot.asp=1)
```
I reset the names of each timeslice/scenario to better reflect what I'm trying to show. In `machu.plotmap`, I specified `plot.asp=1`, which controls the aspect ratio of the output window that the function uses to determine how many rows and columns to use. I'm after a 3x2 grid.
![clipdemo](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/clipdemo.png?raw=true)

Here we see that the models are broadest at the top, and much narrower at the bottom where I've clipped the tails of the skew-normal distributions halfway. This can be a useful way of making models more or less conservative a posteriori.
#### Making binary niche models
In some applications one may want to produce a binary (0,1) niche model rather than a continuous one. This is accommodated in `machu.3.anc.niche` via the `resp.curv` parameter, which specifies whether to use response curves in constructing the niche model. By default, `resp.curv=TRUE`. When `FALSE`, the function will instead produce binary models. With this method, the values of pixels with suitabilities within or greater than confidence limits determined by `clip.amt` are converted to 1; the values of pixels with suitabilites lower than the lower confidence limit are converted to 0. Below, we produce binary models at the default confidence level (95%) and at a more stringent one (75%), then display the results.
```
# producing binary models
mod.ts.u.b95 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], resp.curv=F) # default clipping (95%)
mod.ts.u.b75 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], resp.curv=F, clip.amt=0.75) # stringent clipping (75%)
# combine and visualize
binary.demo <- c(mod.ts.u.b95, mod.ts.u.b75)
names(binary.demo) <- c("clipping at 95%", "clipping at 75%")
machu.plotmap(binary.demo, plot="t", title.cex=0.8, axes=F, plot.asp=1)
```
![binarydemo](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/binarydemo.png?raw=true)

The models are now in binary format, where a pixel either is suitable or it isn't. The bottom row with more stringent clipping produces less suitable area.
#### Projecting a single scenario into multiple paleoclimates
Because it is written as a series of nested apply statements, `machu.3.anc.niche` can handle any number of timeslices or paleoclimate datasets. The question is how it associates a timeslice to a paleoclimate dataset. In the previous sections, each of three timeslices was associated with one of three paleoclimates. But in some cases, a user may want to project one timeslice into multiple paleoclimates, perhaps for comparison. In fact, this is analogous to the usual practice of projecting present-day niche models into past or future climates, which does not account for niche evolution (Machuruku was developed to overcome this issue). In any case, all this takes is specifying a single timeslice and multiple paleoclimates:
```
# one timeslice into multiple paleoclimates
clim <- list(mis19.reduced, mpwp.reduced, m2.reduced)
mod.multipc <- machu.3.anc.niche(ace.n, clim, taxa=1:4, clip.Q=F, verbose=T)
```
Here I'm using the same three climate datasets as before, and `ace.n`, which you may recall doesn't correspond to any timeslice but instead contains the reconstructed climate response variables for the three nodes in our phylogeny, as well as the original climate response values of the four taxa at the tips of the phylogeny. The `ace.n` object is a list with a single element called `"tips_and_nodes"`. In this case, I'm only interested in projecting the four modern taxa (*bassleri, yoshina, pepperi,* and *silverstonei*) into the three paleoclimate datasets, so I specified `taxa=1:4`, which tells the function only to run the first four taxa in each scenario (taxa can also be specified by name). Since `verbose=T`, the function prints the following progress report to screen:
```
[1] "Retained the following scenarios and taxa:"
[1] "tips_and_nodes: bassleri, pepperi, silverstonei, yoshina"
[1] "Did not detect uncertainty samples in 'ace' ('unc' in machu.2.ace())."
[1] "Associating scenarios ('ace') with rastersets ('clim') (numbers are indices, rows are associations):"
     clim ace
[1,]    1   1
[2,]    2   1
[3,]    3   1
[1] "Processed scenario 1"
[1] "Processed scenario 2"
[1] "Processed scenario 3"
```
When specific taxa are specified, `machu.3.anc.niche` prints a sanity check describing which scenarios (=timeslices) and taxa were retained in the analysis. This is because not all taxa may be present in all scenarios when timeslices are taken. If a certain timeslice contains none of the specified taxa, it will be dropped from the analysis. In this case, the second line of the progress report tells us that for our one `"tips_and_nodes"` scenario, only the extant taxa have been retained, as intended.

Now look at the association table in the middle of the report. Contrary to the previous case, where each column went `1 2 3`, "`ace`" now reads `1 1 1`. This is showing that the one and only scenario in `ace.n` (`"tips_and_nodes"`) is being projected into each of the three paleoclimates. 

Let's visualize the results:
```
# visualize
machu.plotmap(mod.multipc, col=2, plot="t", axes=F, to.scale=T)
```
![multipleclimates](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/multipleclimates.png?raw=true)

This shows that none of the extant taxa have any suitable habitat in mis19 ("rasterSet1"), and have varyingly small amounts in mpwp and m2. But again, these are based on the extant taxa, not their ancestors, which doesn't account for evolution in their physioclimatic tolerances.

Note that the maps are drawn in a new color palette, which I specified with `col=2` in the above `machu.plotmap` command. This function comes with 6 premade color palettes, inspired by the [Viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html) R package for colorblind-friendly graphing. Palette `2` is the "plasma" ramp from Viridis, which can also be specified by name. You can also specify your own color palette with the `colorRampPalette` function. 
#### Projecting multiple scenarios into a single paleoclimate
One can also do the opposite of the above and project multiple timeslices into a single paleoclimate dataset. This doesn't strike me as particularly useful except for comparing taxa from multiple time periods in the same climatic context. 
```
# multiple timeslices into a single paleoclimate
mod.multi.ts <- machu.3.anc.niche(ace.ts.u, m2.reduced, clip.Q=F, verbose=T)
```
Here I specify as input `ace.ts.u`, which contains all three timeslices from before with uncertainty data, plus the m2 paleoclimate dataset by itself. This prints the following progress report to screen:
```
[1] "Detected uncertainty samples in 'ace' ('unc' in machu.2.ace())."
[1] "Warning: solo SpatRaster detected, with no raster.sets specified; treating as a single timeslice."
[1] "Associating scenarios ('ace') with rastersets ('clim') (numbers are indices, rows are associations):"
     clim ace
[1,]    1   1
[2,]    1   2
[3,]    1   3
[1] "Processed scenario 1"
[1] "Processed scenario 2"
[1] "Processed scenario 3"
```
The association table is now reversed from the previous section: the "`clim`" column reads `1 1 1`, showing that the one and only paleoclimate dataset has been associated with all three timeslices. Now let's visualize the results (using `col=3`, which is the "viridis" color palette):
```
# visualize
machu.plotmap(mod.multi.ts, col=3, plot="t", axes=F, to.scale=T)
```
![multipletimeslices](https://github.com/wxguillo/machuruku/blob/main/tutorial/images/multipletimeslices.png?raw=true)

Recall that there are 8 taxa total in our 3 timeslices, 4 in the mis19 timeslice, and 2 each in the mpwp and m2 timeslices. The titles of these plots refer to these three timeslices, but it's important to note that they are *only projected into the m2 paleoclimate data* (i.e. "timeslice_3.3). The scenario component of each title refers to the taxon's timeslice of origin, rather than the paleoclimate scenario it is projected into.

`machu.3.anc.niche` can also accommodate the odd case where there are differing amounts of timeslices and paleoclimates, all greater than 1. In this situation, the function will do a 1:1 association of timeslices and paleoclimates *up to* the length of the shorter input. In other words, when 2 timeslices and 3 paleoclimates are provided, timeslice 1 and paleoclimate 1 will be associated, timeslice 2 and paleoclimate 2 will be associated, and paleoclimate 3 will be dropped from the analysis.
```
# differing numbers of timeslices and paleoclimates (all >1)
machu.3.anc.niche(ace.ts[1:2], clim, verbose=T) %>% invisible
```
Here I've specified only the first two timeslices from `ace.ts`, and recall that `clim` contains three paleoclimate datasets. The `%>% invisible` part is using a dplyr pipe (`%>%`) to send the output of the function to `invisible` so that it isn't printed to screen at the end. Since `verbose=T`, the function still prints the following progress report:
```
[1] "Did not detect uncertainty samples in 'ace' ('unc' in machu.2.ace())."
[1] "Associating scenarios ('ace') with rastersets ('clim') (numbers are indices, rows are associations):"
     clim ace
[1,]    1   1
[2,]    2   2
[1] "Processed scenario 1"
[1] "Processed scenario 2"
```
Only the first two paleoclimate datasets were used since there were fewer timeslices than paleoclimates.
#### Saving niche models to a folder
The final feature of `machu.3.anc.niche` to cover is saving outputs. Of course, it isn't very hard to save rasters to a file in whatever format you like with the `writeRaster` function, but Machuruku can save you the trouble with the `output.folder` parameter.
```
# save outputs to folder
machu.3.anc.niche(ace.ts, clim, output.folder=getwd(), verbose=T) %>% invisible
```
I've set `output.folder` to the current working directory. Checking for the files with `list.files(pattern=".tif")` shows their names:
```
[1] "timeslice_0.787_Node1-silverstonei.tif" "timeslice_0.787_Node2-pepperi.tif"      "timeslice_0.787_Node3-bassleri.tif"    
[4] "timeslice_0.787_Node3-yoshina.tif"      "timeslice_3.205_Node1-Node2.tif"        "timeslice_3.205_Node1-silverstonei.tif"
[7] "timeslice_3.3_Node1-Node2.tif"          "timeslice_3.3_Node1-silverstonei.tif"
```
The output file format is always TIFF. The name format is the "scenario_taxon.tif". In certain cases, such as when subsetting your `ace` input, this format can get screwed up. You can use `machu.plotmap` to "preview" the titles in that sense because it works similarly. If your titles are getting screwed up, you might as well create your own function to write the rasters to file. 












