library(dismo)
library(rgdal)
require(raster)
library(phytools)
library(ape)
library(SuppDists)
library(fGarch)
library(schoolmath)
library(ggtree)
library(cowplot)
library(ggplot2)
library(treeio)
library(tidytree)
library(tcltk)
library(spatstat)
library(stringr)
library(RColorBrewer)

######## load the data

# set working directory
setwd("YOUR/DIRECTORY/tutorial")

# current climate data
bio1 <- raster("climate/current/bio_1.tif")
bio12 <- raster("climate/current/bio_12.tif")

# early pleistocene data
bio1mis19 <- raster("climate/mis19/bio_1.tif")
bio12mis19 <- raster("climate/mis19/bio_12.tif")

# mid pliocene warm period
bio1mpwp <- raster("climate/mpwp/bio_1.tif")
bio12mpwp <- raster("climate/mpwp/bio_12.tif")

# pliocene data
bio1m2 <- raster("climate/m2/bio_1.tif")
bio12m2 <- raster("climate/m2/bio_12.tif")

# stack climate data
ClimCur <- stack(bio1, bio12)
ClimMis19 <- stack(bio1mis19, bio12mis19)
ClimMpwp <- stack(bio1mpwp, bio12mpwp)
ClimM2 <- stack(bio1m2, bio12m2)

# load occurrence data
occ	<- read.delim("basslerigroup.csv", h=T, sep=",")

# load tree
bassleritree <- read.nexus("basslerigroup.treefile")

######### exploring the data
# visualize raster layer
plot(bio1)

# visualize occurrence data w/ legend
taxa <- c(as.character(unique(occ$species)))
for (i in taxa){
  points(subset(occ, species==i)$long_DD, subset(occ, species==i)$lat_DD, col = which(taxa==i)+3, pch=16, cex=0.75)
}
legend("topright", legend = c(as.character(unique(occ$species))), col = 4:8, pch=16)

# visualize phylogenetic tree
machu.treeplot(bassleritree, upperX = 2, timelabeloffset = 1, timeslice = 3.3)

######### quick start
resp <- machu.1.tip.resp(occ, ClimCur, verbose = T)
ace <- machu.2.ace(resp, bassleritree, T=3.3)
mod <- machu.3.anc.niche(ace, ClimM2, verbose = T)
par(mfrow=c(1,2))
machu.plotmap(mod)

######### detailed guide
# calculate tip responses with rarefied occurrence data
occ.rarefied <- machu.occ.rarefy(in.pts = occ, rarefy.dist = 10)
resp <- machu.1.tip.resp(occ.rarefied, ClimCur, verbose = T)
# plot rarefied points
plot(bio1)
points(occ$long_DD, occ$lat_DD, pch = 16, cex = 0.75)
taxa <- c(as.character(unique(occ.rarefied$species)))
for (i in taxa){
  points(subset(occ.rarefied, species==i)$long_DD, subset(occ.rarefied, species==i)$lat_DD, col = which(taxa==i)+3, pch=16, cex=0.75)
}
legend("topright", legend = c(as.character(unique(occ.rarefied$species))), col = 4:8, pch=16)

# plot rarefied points
plot(bio1, axes = F, xlim = c(-78, -74), ylim = c(-10,-5))
points(occ$long_DD, occ$lat_DD, pch = 16, cex = 0.75)
for (i in c(as.character(unique(occ.rarefied$species)))){
  points(subset(occ.rarefied, species==i)$long_DD, subset(occ.rarefied, species==i)$lat_DD, col = which(b==i)+3, pch=16, cex=0.75)
}
legend("topright", legend = c(as.character(unique(occ$species))), col = 4:8, pch=16)

# visualize response curves
machu.respplot(resp)
machu.respplot(resp, taxa = "bassleri", linewidth = 2, linecolor = "coral3")
machu.respplot(resp, taxa = c(1,4), clim = "bio_1", linewidth = 2, linecolor = "coral3")
machu.respplot(resp, taxa = c(1,4), clim = "bio_1", linewidth = 2, comb = T)
machu.respplot(resp, linewidth = 2, comb = T, legend.pos = "topright")
# visualize time-slices
machu.treeplot(bassleritree, upperX = 2, timelabeloffset = 1, timeslice = c(0,0.787,3.205,3.3))
# run one time-slice for each time period
ace.mis19 <- machu.2.ace(resp, bassleritree, T=0.787)
ace.mpwp <- machu.2.ace(resp, bassleritree, T=3.205)
ace.M2 <- machu.2.ace(resp, bassleritree, T = 3.3)
# reconstruct niches for all nodes
ace.all <- machu.2.ace(resp, bassleritree)
# visualize climate evolution
machu.respplot(ace.all, comb = T, legend.pos = "topright")
# reconstruct w/ uncertainty
ace.M2.unc <- machu.2.ace(resp, bassleritree, T = 3.3, n.unc = 4)
# build new trees with 95 HPD
# load tree as treedata object
bassleribeast <- read.beast("basslerigroup.treefile")
trees <- machu.tree.unc(bassleribeast)
machu.treeplot(trees, upperX = 3, timelabeloffset = 1.25, nodelabelcolor = "darkseagreen2", timeslice = 2.49)
# run machu.2.ace with each tree separately
ace.LCI <- machu.2.ace(resp, as.phylo(trees$`LCI tree`), T = 2.49)
ace.input <- machu.2.ace(resp, as.phylo(trees$`input tree`), T = 2.49)
ace.UCI <- machu.2.ace(resp, as.phylo(trees$`UCI tree`), T = 2.49)
# save machu.2.ace() output
ace.all <- machu.2.ace(resp, bassleritree, savename = "ace_all")
# load ace_all.csv
ace.all <- machu.ace.load("ace_all.csv")
# build models with machu.3.anc.niche()
mod <- machu.3.anc.niche(ace.M2, ClimM2, verbose = T)
par(mfrow=c(1,2))
machu.plotmap(mod)
# build present-day models for all taxa including ancestors
mod.all <- machu.3.anc.niche(ace.all, ClimCur, verbose = T)
par(mfrow=c(2,4))
machu.plotmap(mod.all)
plot(bassleritree)
# build present-day models for only extant taxa
mod.extant <- machu.3.anc.niche(ace.all, ClimCur, verbose = T, taxa = 1:4)
par(mfrow=c(1,4))
machu.plotmap(mod.extant)
# turn clip.Q off
mod.clip.Q <- machu.3.anc.niche(ace.M2, ClimM2, verbose = T, clip.Q = F)
par(mfrow=c(2,2))
machu.plotmap(mod)
machu.plotmap(mod.clip.Q)
# use resp.curv = F
mod.resp.curv <- machu.3.anc.niche(ace.M2, ClimM2, verbose = T, resp.curv = F)
par(mfrow=c(2,2))
machu.plotmap(mod)
machu.plotmap(mod.resp.curv)
# use calc.unc = T
mod.calc.unc <- machu.3.anc.niche(ace.M2.unc, ClimM2, verbose = T, calc.unc = T)
par(mfrow=c(2,2))
machu.plotmap(mod)
machu.plotmap(mod.calc.unc)
