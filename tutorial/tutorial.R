library(machuruku)

######## load the data
# set working directory
setwd("YOUR/DIRECTORY/tutorial")

# load tree
library(ape)
tree <- read.nexus("basslerigroup.treefile")
plot(tree)

# load occurrence data
occ <- read.csv("basslerigroup.csv")

# load climate data
library(terra)
current <- rast(list.files("climate/current", full.names=T))
mis19 <- rast(list.files("climate/mis19", full.names=T))
mpwp <- rast(list.files("climate/mpwp", full.names=T))
m2 <- rast(list.files("climate/m2", full.names=T))

# plot climate
plot(current$bio_1)
# add occurrence data
taxa <- unique(occ$species)
cols <- c("cyan", "yellow", "red", "orange")
for (i in taxa) points(subset(occ, species==i)$long, subset(occ, species==i)$lat, 
                       col = cols[which(taxa==i)], 
                       pch=19, cex=0.75)
legend(x = "bottomleft", legend = taxa, col = cols, pch=19, bty="n")

######### quick start
# visualize tree with timeslice at 3.3 Ma
machu.treeplot(tree, timeslice=3.3)
# estimate tip response curves
resp <- machu.1.tip.resp(occ, current)
# estimate ancestral niches at timeslice
ace <- machu.2.ace(resp, tree, timeslice=3.3, unc=T)
# project ancestral niches into paleoclimate data
mod <- machu.3.anc.niche(ace, m2)
# visualize ancestral niches
machu.plotmap(mod, plot="together", plot.asp=20/9, axes=F, to.scale=T)

######### detailed guide
### 1. estimate tip response curves
resp <- machu.1.tip.resp(occ, current, verbose=T)
resp[,1:6]
# visualize response curves for two variables
machu.respplot(resp, clim=c("bio_1", "bio_10"), plot="t")
# visualize response curves for all variables
machu.respplot(resp, plot="t", legend.cex=0.6)
# visualize niche models
machu.1.tip.resp(occ, current, plot="t")
dev.off()
machu.1.tip.resp(occ, current, plot="s", plot.points=T)
# output niche models instead of response table
mod <- machu.1.tip.resp(occ, current, output.bioclim=T)
par(mfrow=c(2,2), mar=c(0,0,2,0))
lapply(1:4, function(x) plot(mod[[x]], axes=F, legend=F, box=F, main=names(mod)[x]))

# rarefy occurrence data
occ.r <- machu.occ.rarefy(occ, rarefy.dist=10, plot=T)
resp <- machu.1.tip.resp(occ.r, current, plot="t", plot.points=T, verbose=T)

# select most important climate variables
# contribution method
micv.contrib <- machu.top.env(occ.r, current, method="contrib", contrib.greater=5, verbose=T)
micv.contrib
# nvars method
micv.nvars <- machu.top.env(occ.r, current, method="nvars", nvars.save=7, verbose=T)
micv.nvars
# estimate method
micv.est <- machu.top.env(occ.r, current, method="estimate", verbose=T)
micv.est
# reduce climate datasets to only most important variables
current.reduced <- current[[micv.nvars]]
mis19.reduced <- mis19[[micv.nvars]]
mpwp.reduced <- mpwp[[micv.nvars]]
m2.reduced <- m2[[micv.nvars]]
# re-run tip response table with reduced climate dataset and rarefied occurrence data
resp <- machu.1.tip.resp(occ.r, current.reduced, plot="t", plot.points=T)

### 2. estimate ancestral niches
# visualize time-slices
dev.off()
options(warn=-1)
machu.treeplot(tree, timeslice = c(0.787,3.205,3.3))
machu.treeplot(tree, timeslice = c(0.787,3.205,3.3), x.l.lim = 13, x.u.lim=-3, nodelabsize = 0.35, col = "skyblue", timelaboffset = -0.2)

# estimate ancestral niches at timeslices
resp <- resp[[2]]
ace.ts <- machu.2.ace(resp, tree, timeslice=c(0.787,3.205,3.3))
ace.ts
names(ace.ts)
ace.ts[[1]][,1:8]
ace.ts[[2]][,1:8]

# estimate ancestral niches at each node
ace.n <- machu.2.ace(resp, tree)
names(ace.n)
ace.n[[1]][,1:7]
# visualize niche evolution
machu.respplot(ace.n[[1]], clim="bio_12", fill=T)

# include uncertainty
ace.ts.u <- machu.2.ace(resp, tree, timeslice=c(0.787,3.205,3.3), unc=T)
ace.ts.u[[1]][,1:11]
# visualize uncertainty
machu.respplot(ace.ts.u[[1]], clim="bio_12")

# characterize divergence time uncertainty
# reload tree as 'treedata' object
library(treeio)
beast <- read.beast("basslerigroup.treefile")
beast
do.call(rbind, as_tibble(beast)$height_0.95_HPD)
# run tree-uncertainty utility from a single input tree
beast.trees <- machu.tree.unc(beast)
names(beast.trees)
beast.trees
# visualize trees
machu.treeplot(beast.trees, timeslice=c(0.787,3.205,3.3), nodelabsize=0.5, timelaboffset=-0.15, col="skyblue")
# re-run machu.2.ace with the new trees for comparison
ace.lCItree <- machu.2.ace(resp, beast.trees$lCItree, timeslice=3.3)
ace.inputtree <- machu.2.ace(resp, beast.trees$inputtree, timeslice=3.3)
ace.uCItree <- machu.2.ace(resp, beast.trees$uCItree, timeslice=3.3)
# count number of taxa in each of the new outputs
c(nrow(ace.lCItree[[1]]), nrow(ace.inputtree[[1]]), nrow(ace.uCItree[[1]]))

# run tree-uncertainty utility from a posterior
post.trees <- machu.tree.unc("posterior.trees")
# visualize trees
machu.treeplot(post.trees, timelabs=F, nodelabs=F)

# save ace output
ace.ts <- machu.2.ace(resp, tree, timeslice=c(0.787,3.205,3.3), csv.name="ace_ts.csv")
read.csv("ace_ts.csv")[,1:7]
# load ace output
ace.ts <- machu.ace.load("ace_ts.csv")

### 3. project ancestral niche models into paleoclimate
# multiple timeslices into multiple paleoclimates
clim <- list(mis19.reduced, mpwp.reduced, m2.reduced)
mod.ts <- machu.3.anc.niche(ace.ts, clim, verbose=T)
mod.ts
names(mod.ts)
names(mod.ts[[1]])
sapply(mod.ts, length)
machu.plotmap(mod.ts, plot="t", axes=F, title.cex=0.85, to.scale=T, plot.asp=20/9)

# projecting while accounting for uncertainty
mod.ts.u <- machu.3.anc.niche(ace.ts.u, clim, verbose=T)
machu.plotmap(mod.ts.u, plot="t", axes=F, title.cex=0.85, to.scale=T, plot.asp=20/9)

# clipping niche models
mod.ts.u.nc <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], clip.Q=F) # no clipping
mod.ts.u.c95 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], clip.Q=T) # default clipping (95%)
mod.ts.u.c50 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], clip.Q=T, clip.amt=0.5) # stringent clipping (50%)
# combine and visualize
clip.demo <- c(mod.ts.u.nc, mod.ts.u.c95, mod.ts.u.c50)
names(clip.demo) <- c("no clipping", "clipping at 95%", "clipping at 50%")
machu.plotmap(clip.demo, plot="t", axes=F, to.scale=T, plot.asp=1)

# producing binary models
mod.ts.u.b95 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], resp.curv=F) # default clipping (95%)
mod.ts.u.b75 <- machu.3.anc.niche(ace.ts.u[3], clim[[3]], resp.curv=F, clip.amt=0.75) # stringent clipping (75%)
# combine and visualize
binary.demo <- c(mod.ts.u.b95, mod.ts.u.b75)
names(binary.demo) <- c("clipping at 95%", "clipping at 75%")
machu.plotmap(binary.demo, plot="t", title.cex=0.8, axes=F, plot.asp=1)

# one timeslice into multiple paleoclimates
clim <- list(mis19.reduced, mpwp.reduced, m2.reduced)
mod.multipc <- machu.3.anc.niche(ace.n, clim, taxa=1:4, clip.Q=F, verbose=T)
# visualize
machu.plotmap(mod.multipc, col=2, plot="t", axes=F, to.scale=T)

# multiple timeslices into a single paleoclimate
mod.multi.ts <- machu.3.anc.niche(ace.ts.u, m2.reduced, clip.Q=F, verbose=T)
# visualize
machu.plotmap(mod.multi.ts, col=3, plot="t", axes=F, to.scale=T)

# differing numbers of timeslices and paleoclimates (all >1)
machu.3.anc.niche(ace.ts[1:2], clim, verbose=T) %>% invisible

# save outputs to folder
machu.3.anc.niche(ace.ts, clim, output.folder=getwd(), verbose=T) %>% invisible
list.files(pattern=".tif")

