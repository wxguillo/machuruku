##############################################################
#                Machuruku v0.1.5
# Includes modeling with uncertainty in ancestral state
# reconstructions and in divergence time estimations
#
# machu.2.ace() now also has the option to return models
# for each node rather than each branch at a time-slice
#
# added time-slice plotting to machu.treeplot()
#
# machu.3.anc.niche() now has a verbose argument and can
# construct models for both types of output from machu.2.ace()
#
#                Machuruku v0.1.6
# machu.tree.unc now accounts for cases in which a node's
# error bars extend past those of its parent or descendant
# node, which would create a paradoxical topology
#
# machu.tree.unc now prints an error message if the input
# tree is not ultrametric
#
# added machu.occ.rarefy (from humboldt) to rarefy spatial
# occurrence data and account for spatial autocorrelation
#
# added improved verbosity to machu.3.anc.niche
#
# added the "taxa" argument to machu.3.anc.niche so that
# you can run it on select taxa
#
# added a "verbose" option to machu.1.tip.resp
#
# added a "savename" option to machu.2.ace so that you can
# save output as a .csv from this step
#
# added machu.ace.load() function to load saved output from
# machu.2.ace()
#
#                Machuruku v0.1.7.1
# formatted for roxygen2, with extraneous commands removed
# (e.g., everything not part of a function)
#
# added machu.respplot() function
#
# fixed a bug where using the "savename" parameter in
# machu.2.ace() didn't save a .csv unless a time-slice
# was specified
#
#                Machuruku v0.1.7.2
# changed n.unc to work with any whole number, and re-charac-
# terized it to be about "n samples" rather than "n+1 values
# delimiting n bins" which was way too confusing
#
# changed machu.plotmap() to use custom color palettes from
# Humboldt
#
# added machu.top.env() to identify the most important climate
# variables for the taxon set
#
#                Machuruku v0.1.7.3
# first major reformatting as a package and commit to Github
#
#                Machuruku v0.1.8
# added the machu.trees.unc() function to characterize node
# height uncertainty from a posterior dist of trees
#
# rewrote machu.treeplot() to be more efficient and flexibly
# accept input from both machu.tree.unc() and machu.trees.unc()
#
# added the machu.geo.idw() function to clip extraneous areas
# from niche models via inverse-distance weighting
#
#                Machuruku v0.1.8.1
# changed importation of "spatstat" to "spatstat.geom" as the
# 'nndist' function was moved to the latter package.

##############################################################

##############################################################
#######################    STEP 1    #########################
#   get response curves for ancestral state reconstruction   #
##############################################################

#' Calculate tip climate response curves
#'
#' For each taxon, construct a BIOCLIM species distribution model
#' and estimate a response curve for each climate variable.
#'
#' @param occ occurrence data for each species, formatted as a dataframe. Format should be species, x/long, and y/lat in that order.
#' @param ClimCur present-day climate data, formatted as a RasterStack
#' @param verbose if TRUE, print progress (by-taxon) to the screen
#'
#' @return table consisting of the response of each species to the climate
#' data. Each response is represented by a skew-normal distribution.
#'
#' @examples
#' response.table <- machu.1.tip.resp(occ, ClimCur)
#' @import dismo
#' @import fGarch
#' @export
machu.1.tip.resp <- function(occ, ClimCur, verbose = FALSE){

  # check if occurrence data is a dataframe
  if (is.data.frame(occ) == TRUE){
    # check if ClimCur is a raster stack
    if (class(ClimCur)[1] == "RasterStack"){

      # get list of taxa
      species.in <- apply(occ[1], 2, list)
      # reduce to unique species names
      species.in <-rapply(species.in,function(x)(unique(x)))

      # initialize output list
      output.env <- list()
      # sample climate for each species' localities
      for (i in species.in){
        sp.temp.in <- subset(occ, species == i)
        sp.temp.xy <- sp.temp.in[,2:3]
        sp_name    <- i
        Z<-bioclim(ClimCur, sp.temp.xy)
        output.env[[(i)]] <- Z@presence
      }
      # initialize some vectors that'll come in handy later
      parm_names <- c("mean", "stdev", "skew", "lowerQ", "upperQ")
      response.list <- list()

      # for each species:
      for (j in species.in){

        sp_name   <- j
        in.pres   <- output.env[[j]]
        n.bios    <- ncol(in.pres)
        name.bios <- colnames(in.pres)

        # initialize vector that will later be each row of the table
        row.vector <- c()

        # for each climate variable for each species:
        for (i in 1:n.bios){

          # mean, sd, xi
          suppressWarnings(normFit <- snormFit(in.pres[,i]))
          mean <- normFit$par[1]
          sd   <- normFit$par[2]
          xi   <- normFit$par[3]

          # if skew estimate is poor, mainly due to low sample size, fix skew to 1 and measure tradional mean/sd
          if(normFit$convergence==1){
            mean <- mean(in.pres[,i])
            sd   <- sd(in.pres[,i])
            xi   <- 1
          }
          # upper and lower quantiles
          quants <- quantile(in.pres[,i], probs = c(0.025, 0.975), names = FALSE)
          lQ <- quants[1]
          uQ <- quants[2]

          parms <- c(mean, sd, xi, lQ, uQ)

          # name the elements of the vector according to the name of the climate variable and the response variable
          # (e.g., "bio_1_q1")
          for (h in 1:5){
            names(parms)[h] <- paste0(name.bios[i], "_", parm_names[h])
          }
          # build a giant vector for the whole species, that will eventually be turned into a row for the big table later
          row.vector <- append(row.vector, parms)
        }
        # add each species' row.vector to a list to be stored for later
        response.list[[j]] <- row.vector
        if (verbose == TRUE){ print(paste0("Processed taxon ", j)) }
      }

      # initialize the final table of response curve values
      response.table <- matrix(nrow = 1, ncol = (5*(n.bios))) # the multiplier preceding n.bios corresponds
      # to the number of variables estimated.
      # iterate through each species' response values and add them by row to the table
      for (i in 1:(length(response.list))){
        response.table <- rbind(response.table, response.list[[i]])
      }

      # remove the first row of the table (which is just NAs)
      response.table <- response.table[2:nrow(response.table),]

      # add row (species) names to the table
      rownames(response.table) <- species.in

      return(response.table)
    } else {
      print("climate layers not a RasterStack")
    }
  } else {
    print("occ is not a data frame.")
  }
}
#' Select top environmental variables
#'
#' Perform a boosted regression tree analysis to identify the most important climate variables for your taxon set. This function is a modified version of humboldt.top.env() from humboldt. It runs generalized boosted regression models (a machine learning ENM algorithm) to select top parameters for inclusion your analyses. This is important because you want the models to reflect variables that are relevant to the species' distribution. Alternatively, you can run Maxent outside of R and manually curate the variables you include (also recommended).
#'
#' @param occ occurrence data for all taxa. Identical to input for machu.1.tip.resp(). Dataframe, with columns in the order of species, x/long, y/lat.
#' @param clim climate data for all taxa. Identical to input for machu.1.tip.resp(). A rasterstack of corresponding climate variables.
#' @param learning.rt value from 0.01 to 0.001 for building the ENMs, start with 0.01 and if prompted, change to 0.001. the default value is 0.01
#' @param steps numbers of trees to add at each cycle for modelling each taxon. Start with 50 and if you run into problems gradually decrease, stopping at 1. The default value is 50
#' @param method this determines how important environmental variables are selected.There are three options: "estimate", "contrib", "nvars". If method="estimate", the boosted regression tree algorithm will choose the number of variables to include by systematically removing variables until average change in the model exceeds the original standard error of deviance explained. This is the most computationally intensive method. If method="contrib", variables above a relative influence value will be kept. See associated parameter 'contrib.greater'. If method="nvars", a fixed number of user specified variables will be kept. See associated parameter 'nvars.save'. The kept variables are selected by their relative influence. The 'nvars.save'-highest contributing variables for each taxon are retained and pooled, then ranked, and the 'nvars.save'-highest contributing variables for the whole pool are finally retained.
#' @param nvars.save if method="nvars",this variable is required. It is the number of the top variables to save. The kept variables are selected by their relative influence in predicting the species distribution, selecting for the highest contributing variables. Often the total variables retained is lower due to identical variables select among both species. The default value is 5.  This value will be ignored if method="estimate" or "contrib"
#' @param contrib.greater if method="contrib", this variable is required. The kept variables are selected for their relative influence in predicting the species' distribution. Here, users select variables equal to or above an input model contribution value. The default value for this method is 5 (= variables with 5 percent or higher contribution to model of either species are kept). This value will be ignored if method="estimate" or "nvars"
#' @param pa.ratio ratio of pseudoabsences to occurrence points, typically this is 4. The default value is 4. There have to be at least 50 total points (occ+pseudoabsences) for the model to work; if the sum does not total to 50, the difference is taken as the number of pseudoabsences, rather than the value of occ*pa.ratio.
#' @param verbose print progress to the screen if TRUE, do not if FALSE
#' @return prints the important variables to the screen. It is suggested to then combine them into a new rasterstack.
#' @import dismo
#' @import sp
#' @import raster
#' @import gbm
#' @examples
#' # identify the top 6 climate variables across all taxa
#' machu.top.env(occ, clim, method = "nvars", nvars.save = 6)
#' # identify all climate variables with a contribution greater than 10%
#' machu.top.env(occ, clim, method = "contrib", contrib.greater = 10)
#' @export

machu.top.env <- function(occ, clim, learning.rt = 0.01, pa.ratio = 4, steps = 50, method = "contrib", nvars.save = 5, contrib.greater = 5, verbose = T) # method=='estimate' method=='contrib' method=='nvars'

{
  if (method== "ESTIMATE"){method = "estimate"}
  if (method== "Estimate"){method = "estimate"}
  if (method== "CONTRIB"){method = "contrib"}
  if (method== "Contrib"){method = "contrib"}
  if (method== "contribution"){method = "contrib"}
  if (method== "Nvars"){method = "nvars"}
  if (method== "NVARS"){method = "nvars"}

  # this is necessary bc the gbm.step function works by setting "silent" rather than "verbose", which are opposites
  if (verbose == T){verbose <- F} else if (verbose == F){verbose <- T}

  #################
  # specify occ and clim
  #################
  # determine indices (i.e., which columns, corresponding to no. climate variables) of predictor vars
  e.var <- 4:(3+nlayers(clim))

  # identify which species has the most occurrence data
  n.occ.max <- 0
  for (i in unique(occ$species)){
    n.sp <- nrow(subset(occ, species==i))
    if (n.sp > n.occ.max){
      n.occ.max <- n.sp
    }
  }

  # randomly sample n.occ.max*30 points from the raster climate data : these will be the pool from
  # which pseudoabsences will be drawn
  env.points <- rasterToPoints(clim[[1]], fun=NULL, spatial=FALSE)[,1:2]

  if (verbose == F) { print("Extracting climate values from rasters. This could take a minute")}
  rast.val <- data.frame(env.points, extract(clim, env.points))
  if (verbose == F) { print("Finished extracting climate values")}

  ID <- rep(0, nrow(rast.val))
  rast.val <- cbind(ID, rast.val)
  pa.pool <- rast.val[sample(nrow(rast.val), n.occ.max*30), ]

  # initialize output list
  imp.vars <- list()
  # loop over each taxon
  for (i in 1:length(unique(occ$species))){

    # get species name
    sp.name <- as.character(unique(occ$species)[i])

    # get number of points for that taxon
    n.occ <- nrow(subset(occ, species==sp.name))
    # randomly sample n.occ*4 pseudoabsence (pa) points from pa.pool
    # there have to be at least 50 points for the gbm model to work
    if (n.occ+(n.occ*pa.ratio) < 50){
      n.pa <- 50 - n.occ
    } else {
      n.pa <- n.occ*pa.ratio
    }
    sp.pa <- pa.pool[sample(nrow(pa.pool), n.pa), ]

    # get occurrence points for that species
    sp.occ <- subset(occ, species==sp.name)[,2:3]
    colnames(sp.occ) <- c("x", "y")
    ID <- rep(1, nrow(sp.occ))

    # extract climate values for actual occurrence points and put it all in a table
    sp.data <- data.frame(ID, sp.occ, extract(clim, sp.occ))

    # combine pseudoabsence and species occurrence tables
    datamodel <- rbind(sp.pa, sp.data)

    if (method == "nvars" || method == "contrib") {
      modelOut <- gbm.step(data = datamodel, gbm.x = e.var, gbm.y = 1, family = "bernoulli",
                            tree.complexity = 5, learning.rate = learning.rt, step.size = steps, bag.fraction = 0.5,
                            plot.main = FALSE, plot.folds = FALSE, silent = verbose)
    }
    if (method == "estimate") {
      modelOut1 <<- gbm.step(data = datamodel, gbm.x = e.var, gbm.y = 1, family = "bernoulli",
                              tree.complexity = 5, learning.rate = learning.rt, step.size = steps, bag.fraction = 0.5,
                              plot.main = FALSE, plot.folds = FALSE, silent = verbose)
      mod.simp <<- gbm.simplify(modelOut1)
      min.y <- min(c(0, mod.simp$deviance.summary$mean))
      n.vars1 <- match(min.y, c(0, mod.simp$deviance.summary$mean)) - 1
      vars1.use <- mod.simp$pred.list[[n.vars1]]
      names(env1)[vars1.use]
      modelOut <- gbm.step(data = datamodel, gbm.x = mod.simp$pred.list[[n.vars1]], gbm.y = 1,
                            family = "bernoulli", tree.complexity = 5, learning.rate = learning.rt, step.size = steps,
                            bag.fraction = 0.5, plot.main = FALSE, plot.folds = FALSE, silent = verbose)
    }
    if (method == "nvars") {
      varInc <- as.vector(modelOut$contributions[1:nvars.save, ])
    }
    if (method == "estimate") {
      varInc <- as.vector(modelOut$contributions[, 1])
    }
    if (method == "contrib") {
      varInc <- as.vector(modelOut$contributions$var[modelOut$contributions$rel.inf > contrib.greater])
    }
    imp.vars[[i]] <- varInc
  }
  # sort the combined variable importances and take the top nvars.save across all taxa
  if (method == "nvars") {
    h <- as.data.frame(bind_rows(imp.vars))
    k <- as.vector(unique(h$var))
    for (i in k){
      j<-h[h$var==i,]
      if (match(i, k) == 1){
        l <- j[j$rel.inf==max(j$rel.inf),]
      } else {
        l <- rbind(l, j[j$rel.inf==max(j$rel.inf),])
      }
    }
    l <- l[order(-l$rel.inf),]
    imp.var.out <- as.vector(l[1:nvars.save,][,1])
  }
  if (method == "contrib" || method == "estimate"){
    imp.var.out <- sort(unique(unlist(imp.vars)))
  }
  print("Important variables:")
  print(imp.var.out)
}

#' Rarefy occurrence points
#'
#' Rarefy spatial occurrence points to reduce spatial autocorrelation and model bias. This is a modified version of humboldt.occ.rarefy() from humboldt.
#'
#' @param in.pts input data frame
#' @param colxy columns corresponding to longitude, then latitude
#' @param rarefy.dist distance to rarefy points (values need to be in km (recommended) or decimal degrees).  See associated parameter rarefy.units.
#' @param rarefy.units the units of rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees
#' @param verbose if verbose=T, texts boxes displaying progress will be displayed
#' @return A script to systematically select localities within a specified area at specified spatial resolution.  The outcome is always the same and is not random.  This reduces sampling biases in downstream analyses- you should do it! Output is a reduced dataset with less spatial autocorrelation.
#' @importFrom spatstat.geom nndist
#' @importFrom tcltk setTkProgressBar
#' @examples
#' ##remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.
#' occ <- machu.occ.rarefy(in.pts = occ, colxy = 2:3, rarefy.dist = 50, rarefy.units = "km")
#' @export

# use removeCollinearity() from virtualspecies to pick out co-correlated variables

machu.occ.rarefy <- function(in.pts, colxy = 2:3, rarefy.dist = 0, rarefy.units = "km", verbose = T) {
  switch(Sys.info()[['sysname']],
         Windows= {userOS=1},
         Linux  = {userOS=2},
         Darwin = {userOS=2})

  if (rarefy.units == "KM"){rarefy.units = "km"}
  if (rarefy.units == "Km"){rarefy.units = "km"}
  if (rarefy.units == "DD"){rarefy.units = "dd"}
  if (rarefy.units == "Dd"){rarefy.units = "dd"}

  if (rarefy.units == "km") {
    min.dist <- rarefy.dist * 1000  #values in km
    sp2p1 <- SpatialPoints(in.pts[, colxy], CRS("+proj=longlat +datum=WGS84"))
    sp2p2 <- spTransform(sp2p1, CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    xy <- data.frame(x4r = coordinates(sp2p2)[, 1], y4r = coordinates(sp2p2)[, 2])
    xy <- data.frame(cbind(xy,in.pts))#new
  }

  if (rarefy.units == "dd") {
    yy <- colxy[2]
    maxLat <- max(in.pts[, yy])
    minLat <- min(in.pts[, yy])
    # estimage dd to km for study area
    rare.dd <- (mean(c(maxLat, minLat)))
    adjKm <- (-0.0139 * (rare.dd * rare.dd)) + (0.0898 * rare.dd) + 111.1
    min.dist <- adjKm * rarefy.dist * 1000  #values in km
    print(paste("Value used for rarefying:", round((min.dist/1000), 2), "km. Remember that the length of a decimal degrees changes latitudinally due to the convergence of the lines of longitude at the poles. The value used here is the average distance of decimal-degrees within your study area. Alternatively, simply input distance as km value and change rarefy.units='km'"))
    sp2p1 <- SpatialPoints(in.pts[, colxy], CRS("+proj=longlat +datum=WGS84"))
    sp2p2 <- spTransform(sp2p1, CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    xy <- data.frame(x4r = coordinates(sp2p2)[, 1], y4r = coordinates(sp2p2)[, 2])
  }

  nPts<-nrow(xy)
  if (verbose == T & userOS==1){pb <- winProgressBar(title = "Initializing",min = 0,max =nPts, width = 300)}
  if (verbose == T & userOS==2){pb <- tkProgressBar(title = "Initializing", label = "", min = 0, max = nPts, initial = nPts, width = 300)}
  # setup env- sepaerate from distance
  spName <- in.pts[, 1][1]
  new.data <- NULL

  to <- 1
  del.min.dist <- xy

  repeat {
    nn1 <- nndist(del.min.dist[, "x4r"], del.min.dist[, "y4r"])  # calculate distance nearest neighbour
    if (sum(nn1 < min.dist) == 0) {
      break
    }

    # iteratively removing points starting with the one having the minimal distance to the
    # nearest neighbour
    nn2 <- nndist(del.min.dist[, "x4r"], del.min.dist[, "y4r"], k = 2)

    del1 <- nn1 == min(nn1)
    del2 <- nn2 == min(nn2[del1])
    delk <- del1 & del2
    if (sum(del2) > 1) {
      for (k in 3:6) {
        nn <- nndist(del.min.dist[, "x4r"], del.min.dist[, "y4r"], k = k)
        delk <- delk & nn == min(nn[delk])
        if (verbose == T & userOS==1){setWinProgressBar(pb, length(delk), title=paste(" Rarefying:",length(delk),"remaining localities (of", nPts,"input)"))}
        if (verbose == T & userOS==2){setTkProgressBar(pb, length(delk), title=paste(length(delk),"remaining pts (of", nPts))}
        if (sum(nn[delk] == min(nn[delk])) > 1) {
          break
        }
      }
    }
    # from the two points which are the nearest neighbours of the whole set, remove the one
    # closest to the second neighbour
    del.min.dist <- del.min.dist[-(which(delk)[1]), ]
  }
  if (verbose == T){close(pb)}
  new.data <- rbind(new.data, del.min.dist)
  nc<-(ncol(new.data))
  col.orig<-c(3:nc)
  data.out<-new.data[,col.orig]
  print(paste("Starting points =", nrow(xy), ", Final rarefied points =", nrow(new.data)))
  return(data.out)
}

###############################################################
#######################    STEP 2    ##########################
#       ancestral state reconstuction of response curves      #
###############################################################

#'Account for uncertainty in divergence time estimation from a single tree
#'
#'This function generates two additional trees from a single input tree with uncertainty
#'(e.g. 95% confidence intervals) in divergence time estimations. This function is written to
#'work with output from the BEAST 2 program TreeAnnotator and may not work with other formats.
#'One tree will be built using the lower confidence interval values for each node (i.e. a
#'more recent tree), and the other tree will be built using the upper confidence interval
#'values for each node (i.e. an older tree). The output of the function is a list of three
#'treedata objects, each of which can be separately used as tree input for Machuruku.
#'As the values of 95% confidence intervals do not generally make an ultrametric tree when
#'used as node heights, the phytools function force.ultrametric() is used to force the
#'two generated trees to be ultrametric using the "nnls" function, which may result in
#'some shifting of node heights from the 95% CI values.
#'
#'@param tree input treedata object, imported with the read.beast() function from treeio.
#'Only tested with output from the BEAST 2 utility TreeAnnotator.
#'
#'@return returns a list of three treedata objects better characterizing the uncertainty
#'of the input tree's node heights.
#'
#'@examples
#'trees <- machu.tree.unc(beast)
#'@import ape
#'@import tidytree
#'@import phytools
#'@export
machu.tree.unc <- function(tree){
  if (is.ultrametric(as.phylo(tree)) == TRUE){
    ##### create upper CI tree#####
    # convert the tree to a "tibble"
    UCItreetib <- as_tibble(tree)
    r <- nrow(UCItreetib)
    t <- Ntip(tree)
    # for every row in the tibble, copy the upper CI value for height, length, and rate
    # into the true value slot. this makes the upper CI tree.
    for (x in 1:r){
      nodeUCI <- UCItreetib$height_0.95_HPD[[x]][2]
      pnode <- UCItreetib$parent[[x]]
      pnodeUCI <- UCItreetib$height_0.95_HPD[[pnode]][2]
      # if the upper CI for the node in question is higher in value than that of the
      # parent node, copy the parental upper CI into the node in question's height
      # slot instead.
      if (nodeUCI > pnodeUCI && x != pnode){
        print(paste0("node ", x, " upper CI (", nodeUCI, ") overlaps parent node upper CI (", pnodeUCI, ")"))
        print(paste0("resetting node ", x, " upper CI to ", pnodeUCI))
        UCItreetib$height[x] <- pnodeUCI
      } else {
        # else, just copy the upper CI value for the node's height into the true value
        # slot. This is what should usually happen.
        UCItreetib$height[x] <- UCItreetib$height_0.95_HPD[[x]][2]
      }
    }
    # for every branch in the tree, rescale it so that the tree remains ultrametric.
    for (z in 1:r){
      parNode <- UCItreetib$parent[z]
      UCItreetib$branch.length[z] <- UCItreetib$height[parNode] - UCItreetib$height[z]
    }
    # The tree may not be ultrametric due to the confidence intervals not corresponding to the node heights
    # of an ultrametric tree. We use the phytools function force.ultrametric(), with the nnls method from Phangorn,
    # to force the tree to be ultrametric.
    if (is.ultrametric(as.phylo(as.treedata(UCItreetib))) == TRUE){
      UCItree <- as.treedata(UCItreetib)
    } else {
      UCItree <- as.treedata(force.ultrametric(as.phylo(as.treedata(UCItreetib)), method="nnls"))
      print("UCI tree was rescaled to be ultrametric")
    }

    ##### create lower CI tree#####
    # convert the tree to a "tibble"
    LCItreetib <- as_tibble(tree)
    r <- nrow(LCItreetib)
    t <- Ntip(tree)
    # for every row in the tibble, copy the lower CI value for height, length, and rate
    # into the true value slot. this makes the lower CI tree.
    for (x in 1:r){
      nodeLCI <- LCItreetib$height_0.95_HPD[[x]][1]
      pnode <- LCItreetib$parent[[x]]
      pnodeLCI <- LCItreetib$height_0.95_HPD[[pnode]][1]
      # if the lower CI for the node in question is higher in value than that of the
      # parent node, copy the node's lower CI into the parent node's height
      # slot instead. Also copy the node's lower CI into its own height slot.
      if (nodeLCI > pnodeLCI && x != pnode){
        print(paste0("parent node ", pnode, " lower CI (", pnodeLCI, ") overlaps child node ", x, " lower CI (", nodeLCI, ")"))
        print(paste0("resetting node ", pnode, " lower CI to ", nodeLCI))
        LCItreetib$height[pnode] <- nodeLCI
        LCItreetib$height[x] <- LCItreetib$height_0.95_HPD[[x]][1]
      } else {
        # else, just copy the lower CI value for the node's height into the true value
        # slot. This is what should usually happen.
        LCItreetib$height[x] <- LCItreetib$height_0.95_HPD[[x]][1]
      }
    }
    # for every branch in the tree, rescale it so that the tree remains ultrametric.
    for (z in 1:r){
      parNode <- LCItreetib$parent[z]
      LCItreetib$branch.length[z] <- LCItreetib$height[parNode] - LCItreetib$height[z]
    }
    # The tree is likely not ultrametric due to the confidence intervals not corresponding to the node heights
    # of an ultrametric tree. We use the phytools function force.ultrametric(), with the nnls method from Phangorn,
    # to force the tree to be ultrametric.
    if (is.ultrametric(as.phylo(as.treedata(LCItreetib))) == TRUE){
      LCItree <- as.treedata(LCItreetib)
    } else {
      LCItree <- as.treedata(force.ultrametric(as.phylo(as.treedata(LCItreetib)), method="nnls"))
      print("LCI tree was rescaled to be ultrametric")
    }

    treelist <- list(LCItree, tree, UCItree)
    names(treelist) <- c("LCI tree", "input tree", "UCI tree")
    return(treelist)
  } else {
    print("Input tree is not ultrametric.")
  }
}

#'Account for uncertainty in divergence time estimation from a posterior distribution of trees
#'
#'This function creates a list of three trees that characterize the uncertainty in divergence
#'times from a Bayesian posterior distribution of trees. It is written to work with a NEXUS file
#'containing multiple trees and a translation table of taxon names, specifically the output of
#'BEAST 2 but the format is generic enough to likely work with other program outputs as well.
#'The function is written to use a file without loading it into R and overloading memory, as
#'posterior treefiles can be very large in size.
#'The function conducts three passes through the posterior:
#'1) Calculates the overall number of trees and the amount skipped as burn-in
#'2) Calculates the distribution of tree-heights and attendant 95% HPD interval
#'3) Finds the trees with tree-heights closest to the upper 95% HPD, lower 95% HPD, and average
#'The output of the function is a list of three phylo objects, each of which can be separately
#'used as tree input for Machuruku.
#'
#'@param file input file, provided as a filepath string. File must be NEXUS format with multiple
#'trees and a translation table of taxon names (e.g., 1 Homo sapiens). Only tested with BEAST 2 output.
#'@param burnin a numeric value between 0 and 1, denoting the percentage of trees to skip. MCMC chains do
#'not usually converge right away, so it is common practice to skip the first 10% or so of samples as "burn-in".
#'The default value is 10%.
#'@param conf a numeric value between 0 and 1, denoting the desired tree-height HPD interval to sample trees from. Values
#'closer to 1 will result in more extreme trees. The default value is 95%.
#'@param inc the increment value at which the function will shout its progress at you, if verbose=TRUE. Default
#'is 10,000.
#'@param verbose shout progress, turned on by default
#'
#'@return returns a list of three phylo objects characterizing the uncertainty in divergence time of the posterior
#'distribution of trees and enabling the user to use the trees as separate inputs for Machuruku.
#'
#'@examples
#'trees <- machu.trees.unc(file)
#'@import ape
#'@import stringr
#'@import TeachingDemos
#'@export
machu.trees.unc <- function(file, burnin = 0.1, conf = 0.95, inc = 10000, verbose = TRUE) {
  ##########################
  ##### ERROR CHECKING #####
  ##########################
  # Check if burnin is a percentage between 0 and 100
  if (burnin > 1 | burnin < 0){ stop("Specified burnin is not between 0 and 1.")}
  # Check if conf is a percentage between 0 and 100
  if (conf > 1 | conf < 0){ stop("Specified confidence level is not between 0 and 1.")}
  #########################################################################
  ##### FIRST PASS: Calculating the number of trees and burnin amount #####
  #########################################################################
  if (verbose == TRUE){ print(paste0("FIRST PASS: Calculating trees to remove at ", burnin*100, "% burnin")) }
  con <- file(file, "r")
  # initiate tree counter and translation table
  treenumber <- 1
  tlntab <- rep(NA, Ntip(read.nexus(file))[[1]])
  tlnnumber <- 1
  # read in file line by line
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    # get tip translation table
    if (grepl("^\\s+\\d+\\s.+", line) == TRUE) {
      tlntab[tlnnumber] <- line
      tlnnumber <- tlnnumber+1
    }
    # if the line is actually a tree, do this:
    if (grepl("^tree", line) == TRUE) {
      # shout to screen
      if (treenumber %% inc == 0 & verbose == TRUE) { print(paste0("First pass (burnin): Tree ", treenumber)) }
      treenumber <- treenumber+1
    }
  }
  # trim commas from tip names (except for last one)
  for (i in 1:(length(tlntab)-1)){ tlntab[i] <- str_sub(tlntab[i], end=-2) }
  # create translation table
  tlntab <- strsplit(tlntab,"\\s+")
  tlntab <- do.call(rbind, tlntab)[,2:3]
  # calculate burnin
  burn.no <- round(burnin*treenumber)
  # report burnin and treenumber
  if (verbose == TRUE){ print(paste0("Found ", treenumber, " trees, skipping the first ", burn.no, " under ", burnin*100, "% burnin.")) }
  close(con)
  #####################################################################
  ##### SECOND PASS: Calculating HPD distribution of tree heights #####
  #####################################################################
  if (verbose == TRUE){ print(paste0("SECOND PASS: Calculating ", conf*100, "% HPD of tree heights at ", burnin*100, "% burnin")) }
  con <- file(file, "r")
  # re-initiate tree counter and a vector of tree heights
  treenumber <- 1
  treeheights <- c()
  # read in file line by line
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    # if the line is actually a tree, do this:
    if (grepl("^tree", line) == TRUE) {
      # if treenumber is still below the burnin threshold, skip the tree
      if (treenumber <= burn.no){
        # shout to screen
        if (treenumber %% inc == 0 & verbose == TRUE){ print(paste0("Second pass (HPD): Tree ", treenumber, " (skipped at ", burnin*100, "% burnin)")) }
        treenumber <- treenumber+1
        # if treenumber is above burnin threshold, calculate tree height and add to vector
      } else if (treenumber > burn.no){
        height <- sum(read.tree(text=line)$edge.length)
        treeheights <- c(treeheights, height)
        # shout to screen
        if (treenumber %% inc == 0 & verbose == TRUE){ print(paste0("Second pass (HPD): Tree ", treenumber)) }
        treenumber <- treenumber+1
      }
    }
  }
  # Calculate HPDs and average
  HPD.min <- emp.hpd(treeheights)[1]
  HPD.max <- emp.hpd(treeheights)[2]
  HPD.avg <- mean(treeheights)
  # report HPDs and average
  if (verbose == TRUE){ print(paste0(conf*100, "% HPD min: ", round(HPD.min,digits=2), "; Average: ", round(HPD.avg,digits=2), "; ",
                                     conf*100, "% HPD max: ", round(HPD.max,digits=2), " from ", length(treeheights), " trees under ",
                                     burnin*100, "% burnin.")) }
  close(con)
  ##########################################################
  ##### THIRD PASS: Find closest trees to HPDs and avg #####
  ##########################################################
  if (verbose == TRUE){ print(paste0("THIRD PASS: Finding trees with closest heights to ", conf*100, "% HPD limits and average at ", burnin*100, "% burnin")) }
  con <- file(file, "r")
  # re-re-initiate tree counter
  treenumber <- 1
  # read in file line by line
  #mintree <- "1"; avgtree <- "2"; maxtree <- "3"
  # remove(lowest.from.avg)
  # remove(lowest.from.max)
  # remove(lowest.from.min)
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    # if the line is actually a tree, do this:
    if (grepl("^tree", line) == TRUE) {
      # if treenumber is still below the burnin threshold, skip the tree
      if (treenumber <= burn.no){
        # shout to screen
        if (treenumber %% inc == 0 & verbose == TRUE){ print(paste0("Third pass (trees): Tree ", treenumber, " (skipped at ", burnin*100, "% burnin)")) }
        treenumber <- treenumber+1
      }
      # if tree number = burn.no+1, set all values to this tree
      # this gives us a place to start from
      if (treenumber == burn.no+1){
        height <- sum(read.tree(text=line)$edge.length)
        mintree <- line
        avgtree <- line
        maxtree <- line
        lowest.from.min <- abs(height-HPD.min)
        lowest.from.avg <- abs(height-HPD.avg)
        lowest.from.max <- abs(height-HPD.max)
        # shout to screen
        if (treenumber %% inc == 0 & verbose == TRUE){ print(paste0("Third pass (trees): Tree ", treenumber)) }
        treenumber <- treenumber+1
      }
      # if tree number > burn.no+1, test whether treeheight distance is lower than the current lowest for all values
      # if so, set new values and trees
      if (treenumber > burn.no+1) {
        height <- sum(read.tree(text=line)$edge.length)
        # HPD.min
        if (abs(height-HPD.min) < lowest.from.min){
          lowest.from.min <- abs(height-HPD.min)
          mintree <- line }
        # HPD.avg
        if (abs(height-HPD.avg) < lowest.from.avg){
          lowest.from.avg <- abs(height-HPD.avg)
          avgtree <- line }
        # HPD.max
        if (abs(height-HPD.max) < lowest.from.max){
          lowest.from.max <- abs(height-HPD.max)
          maxtree <- line }
        # shout to screen
        if (treenumber %% inc == 0 & verbose == TRUE){ print(paste0("Third pass (trees): Tree ", treenumber)) }
        treenumber <- treenumber+1
      }
    }
  }
  # Put trees into a list
  treelist <- list(read.tree(text=mintree), read.tree(text=avgtree), read.tree(text=maxtree))
  names(treelist) <- c("min.tree", "avg.tree", "max.tree")
  # need to restore the original tip labels, since nexus represents each with a number
  # sort the translation table made in pass 1 by the order of the numbered tips in the trees
  sorted.tip.labs <- tlntab[match(treelist[[1]]$tip.label, tlntab[,1]),][,2]
  for (i in 1:length(treelist)){
    treelist[[i]]$tip.label <- sorted.tip.labs
  }
  return(treelist)
  close(con)
}

#'Plot trees to visualize each step in the Machuruku process
#'
#'Plot trees with node labels and divergence dates on a time scale. Accepts trees as a
#'single "phylo" or "treedata" object, or a list of "treedata" objects (output from the
#'machu.tree.unc() or machu.trees.unc() functions). Visualizing a single tree allows for the
#'assessment of the proper value for T encompassing the desired number of lineages, and
#'shows node numbers so that the name of a given lineage (which for an internal branch may
#'be something like node5-node6) can be easily attributed to a particular branch by the user.
#'Visualizing a list of trees on the same timescale visualizes the different branch-taxa recoverable
#'by a single time-slice for trees with multiple heights.
#'
#'@param tree tree input. Can be a single tree in treedata or phylo format, or a list of trees in either
#'format. Output of the plotting function will change depending on the type of input.
#'@param upperX set the upper (right-side) limit of the timescale for the treeplot. Must be a
#'positive number. The real purpose of this value is to scale the treeplot so that the tip labels
#'are visible. At a value of 0 the tip labels will likely be cut off by the right border of the plot.
#'@param timelabeloffset set the horizontal offset for the divergence time labels. At different scales
#'the tree plot will place the divergence time labels at varying distances from their corresponding
#'nodes. This value shifts them horizontally so that the user can make an aesthetically pleasing image
#'at their preferred scale.
#'@param nodecirclecolor set the color for the node circles. Set to "cadetblue1" by default.
#'@param nodecirclesize set the size for the node circles. Set to 7 by default.
#'@param nodelabelsize set the size for the node label text (node numbers). Set to 5 by default.
#'@param tiplabelsize set the size for the tip label text. Set to 5 by default.
#'@param timelabelsize set the size for the divergence time label text. Set to 5 by default.
#'@param timeslice display a vertical line through the phylogeny at the desired time-slice for a
#'machu.2.ace analysis. By default this feature is turned off.
#'@param timeslicecolor set the color for the time-slice line. Set to "coral1" by default.
#'@param timeslicelinetype set the line-type for the time-slice line. Set to "solid" by default. Other
#'options include "dashed" and "dotted".
#'@param timelabeldigits set the number of node label (divergence time) decimal points displayed. Default is 2.
#'
#'@return Returns a tree plot with different properties depending on whether the input is a single tree
#'or a list of trees from the functions machu.tree.unc() and machu.trees.unc().
#'
#'@examples
#'# plot three trees derived from machu.tree.unc (a list of three treedata objects)
#'machu.treeplot(trees, upperX=3, timelabeloffset=1.25, nodecirclecolor="orangered1")
#'# plot three trees derived from machu.tree.unc with a timeslice at 2.5 mya
#'machu.treeplot(trees, upperX=3, timelabeloffset=1.25, nodecirclecolor="orangered1", timeslice = 2.5)
#'# plot a single treedata object with no time-slice
#'machu.treeplot(tree, upperX=3, nodecirclecolor = "lightblue", timelabeloffset = 1)
#'# plot a single treedata object with time-slice at 3 mya
#'machu.treeplot(tree, upperX=3, nodecirclecolor = "lightblue", timelabeloffset = 1, timeslice=3, timeslicecolor = "blue", timeslicelinetype = "dashed")
#'@import ape
#'@import ggtree
#'@import ggplot2
#'@import gridExtra
#'@import cowplot
#'@export
machu.treeplot <- function(tree, upperX=20,
                           nodecirclecolor="cadetblue1", nodecirclesize=7,
                           nodelabelsize=5, tiplabelsize=5,
                           timelabeloffset=3, timelabelsize=5, timelabeldigits=2,
                           timeslice=NULL, timeslicecolor="coral1", timeslicelinetype="solid"){
  # if we're plotting a single tree:
  if (length(tree) == 1 || class(tree) == "phylo" ) {

    ### get branching times for the single input tree
    nt <- Ntip(as.phylo(tree))
    nn <- Nnode(as.phylo(tree))
    singletreetimes <- c()
    singletreetimes[1:nt] <- NA
    singletreetimes[(nt+1):(nt+nn)] <- round(branching.times(as.phylo(tree)), digits=timelabeldigits)

    ### plot single tree
    xlimit <- ceiling(max(branching.times(as.phylo(tree))))
    ### plot tree without time-slice line if time-slice argument is not passed
    if (is.null(timeslice) == TRUE){
      revts(ggtree(tree)) +
        geom_nodepoint(color=nodecirclecolor, size=nodecirclesize) +
        geom_text2(aes(subset=!isTip, label=node), size=nodelabelsize) +
        geom_tiplab(size=tiplabelsize) +
        theme_tree2() +
        geom_text2(aes(label=singletreetimes), nudge_x = timelabeloffset, size=timelabelsize) +
        xlim(-xlimit,upperX)
      ### plot tree with time-slice line if time-slice argument is passed
    } else {
      revts(ggtree(tree)) +
        geom_nodepoint(color=nodecirclecolor, size=nodecirclesize) +
        geom_text2(aes(subset=!isTip, label=node), size=nodelabelsize) +
        geom_tiplab(size=tiplabelsize) +
        theme_tree2() +
        geom_text2(aes(label=singletreetimes), nudge_x = timelabeloffset, size=timelabelsize) +
        xlim(-xlimit,upperX) +
        geom_vline(xintercept = -timeslice, color = timeslicecolor, linetype = timeslicelinetype)
    }
    ### if the user has more than one tree (usually 3), plot them all on the same x-axis
  } else if (class(tree) == "list" && length(tree) > 1){

    # get the max tree height value across all trees to set the highest X for the plot
    phylos <- lapply(tree, as.phylo)
    xlimit <- ceiling(max(do.call(rbind, lapply(phylos, branching.times))))

    # if a time-slice is not to be plotted:
    if (is.null(timeslice) == TRUE){
      z<-lapply(phylos,
                function(x){
                  nt <-  Ntip(x)
                  nn <- Nnode(x)
                  times <- rep(NA, nt)
                  times[(nt+1):(nt+nn)] <- round(branching.times(x), digits=timelabeldigits)
                  revts(ggtree(x) +
                          geom_nodepoint(color=nodecirclecolor, size=nodecirclesize) +
                          geom_text2(aes(subset=!isTip, label=node), size=nodelabelsize) +
                          geom_tiplab(size=tiplabelsize) +
                          theme_tree2()) +
                    geom_nodelab(aes(label=times), nudge_x = timelabeloffset, size=timelabelsize) +
                    xlim(-xlimit,upperX)
                }
      )
      # if a time-slice is to be plotted:
    } else if (is.null(timeslice) == FALSE){
      z<-lapply(phylos,
                function(x){
                  nt <-  Ntip(x)
                  nn <- Nnode(x)
                  times <- rep(NA, nt)
                  times[(nt+1):(nt+nn)] <- round(branching.times(x), digits=timelabeldigits)
                  revts(ggtree(x) +
                          geom_nodepoint(color=nodecirclecolor, size=nodecirclesize) +
                          geom_text2(aes(subset=!isTip, label=node), size=nodelabelsize) +
                          geom_tiplab(size=tiplabelsize) +
                          theme_tree2()) +
                    geom_nodelab(aes(label=times), nudge_x = timelabeloffset, size=timelabelsize) +
                    xlim(-xlimit,upperX) +
                    geom_vline(xintercept = -timeslice, color = timeslicecolor, linetype = timeslicelinetype)
                }
      )
    }
    # plot
    gridExtra::grid.arrange(grobs=z, nrow=length(z))
  }
}

#'Visualize response curves of taxa to climate variables
#'
#'Visualize species response curves in a variety of different ways. A "response table" (the output
#'of machu.1.tip.resp()) or an "ace list" (the output of machu.2.ace()) can both be provided. Any
#'combination of taxa can be displayed either in separate plots or on the same plot (using matplot()),
#'and any combination of climate variables used in constructing the models can also be displayed. A
#'few arguments for adjusting line width, color, etc. are provided.
#'
#'@param resp input. This can either be output from machu.1.tip.resp() or from machu.2.ace().
#'@param taxa specify which taxa you wish to plot. This can either be done with text (e.g., "species1", or
#'c("species1,"species2")) or integers (e.g., 1 or 1:3, or c(1,4,6:10)). Note that with integers you must
#'refer to the "resp" input to know which integer corresponds to which taxon. When no taxa are specified,
#'all taxa are displayed by default.
#'@param clim specify which climate variables you wish to plot. This can only be done via text, and the
#'provided text must match the climate variable name exactly, e.g. "bio_1" or c("bio_1", "bio_10"). When
#'no climate variables are specified, all variables are displayed by default.
#'@param comb specify whether to combine taxon plots for the same climate variable. This option will
#'plot only as many charts as there are climate variables, and all taxa will be displayed in each chart,
#'allowing for more direct comparison. Different taxa are colored using the "Set2" color palette from
#'RColorBrewer. Color interpolation is used in cases of more than 8 taxa (the number of colors in the
#'palette). Default is set to FALSE.
#'@param linetype specify the desired line type for the plots. Default is set to 1, a solid line.
#'@param linewidth specify the desired line width for the plots. Default is set to 1.
#'@param linecolor specify the desired line color for the plots. Default is set to black. Note that this
#'only applies in the case where comb = FALSE. When comb = TRUE, the "Set2" color palette from RColorBrewer
#'is used.
#'@param legend.pos specify the desired legend position. You may wish to change this depending on the shape
#'of the curves, so as not to occlude them. The default setting is "topleft".
#'@param par.override override the default panel settings by providing a vector (e.g., c(2,3)) of the
#'desired configuration of rows and columns for the panels. By default, when comb = FALSE, each column
#'corresponds to a climate variable, and each row corresponds to a taxon; when comb = TRUE, all plots
#'are provided in a single row. For large numbers of climate variables, you may wish to specify multiple
#'rows of panels with this parameter.
#'
#'@return Returns a plot or plots showing the response of the specified taxa to the specified climate variable(s).
#'
#'@examples
#'# Plot two separate plots, the responses of both taxa to the climate variable bio_1
#'machu.respplot(resp, taxa = c("species1","species2"), clim = "bio_1")
#'# Plot the same two taxa, but on the same plot
#'machu.respplot(resp, taxa = 1:2, clim = "bio_1", comb = TRUE)
#'# Plot the first five taxa (except for taxon 2), in all climate variables, in separate plots
#'machu.respplot(resp, taxa = c(1,3:5))
#'# Plot all taxa in four specified climate variables, and override paneling to print in a more compact format (2x2)
#'machu.respplot(resp, clim = c("bio_1","bio_2","bio_3","bio_4), par.override = c(2,2))
#'@import stringr
#'@import fGarch
#'@import RColorBrewer
#'@export
machu.respplot <- function(resp, taxa = NULL, clim = NULL, comb = FALSE, linetype = 1, linewidth = 1, linecolor = "black", legend.pos = "topleft", par.override = NULL){
  ##############################
  # if the output of machu.1.tip.resp is provided
  if (typeof(resp) == "double"){
    # use all taxa by default
    if (is.null(taxa) == TRUE){
      taxa <- 1:nrow(resp)
    }
    # use all climate variables by default
    if (is.null(clim) == TRUE){
      q <- colnames(resp)
      q <- str_remove(q, "_mean"); q <- str_remove(q, "_stdev"); q <- str_remove(q, "_skew"); q <- str_remove(q, "_lowerQ"); q <- str_remove(q, "_upperQ");
      clim <- unique(q)
    }
    # if comb is FALSE
    if (comb == FALSE){
      # set number of panels
      if (is.null(par.override) == TRUE){
        par(mfrow = c(length(taxa), length(clim)))
      } else if (is.null(par.override) == FALSE){
        par(mfrow = par.override)
      }
      # for every taxon specified
      for (i in taxa){
        # and for every climate variable specified
        for (j in clim){
          # get climate response values for taxon and climate variable
          mu <- resp[,paste0(j, "_mean")][i]
          sd <- resp[,paste0(j, "_stdev")][i]
          xi <- resp[,paste0(j, "_skew")][i]
          # make the plot
          climate.value <- seq(round(mu-sd*3), round(mu+sd*3))
          density <- dsnorm(climate.value, mean = mu, sd = sd, xi = xi)
          plot(climate.value, density, type="l", lty=linetype, lwd=linewidth, col=linecolor)
          if (typeof(i) == "character"){
            title(paste(i, j))
          } else {title(paste(rownames(resp)[i], j))}
        }
      }
      # if comb is TRUE
    } else if (comb == TRUE){
      # set number of panels
      if (is.null(par.override) == TRUE){
        par(mfrow = c(1, length(clim)))
      } else if (is.null(par.override) == FALSE){
        par(mfrow = par.override)
      }
      # for every climate variable specified
      for (i in clim){
        # create the x-bounds of the plot
        clim.mu <- resp[,paste0(i, "_mean")][taxa]
        clim.sd <- resp[,paste0(i, "_stdev")][taxa]
        lobound <- min(clim.mu-clim.sd*3)
        upbound <- max(clim.mu+clim.sd*3)
        climate.value <- seq(round(lobound), round(upbound))
        vec <- c()
        # for every taxon
        for (j in taxa){
          # get climate response values for taxon and climate variable
          mu <- resp[,paste0(i, "_mean")][j]
          sd <- resp[,paste0(i, "_stdev")][j]
          xi <- resp[,paste0(i, "_skew")][j]
          vec <- append(vec, dsnorm(climate.value, mean = mu, sd = sd, xi = xi))
        }
        # create data matrix
        mat <- matrix(vec, nrow = length(climate.value), ncol = length(taxa))
        # create color ramp
        # if there are more taxa than colors in the Set2 palette (8), use a colorRampPalette interpolation
        if (length(taxa) > 8){
          cols <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(length(taxa))
        } else {
          cols <- brewer.pal(n = length(taxa), name = "Set2")
        }
        # get names for the legend
        if (typeof(taxa[1]) == "character"){
          leg <- taxa
        } else {
          leg <- rownames(resp)[taxa]
        }
        # make the plot
        matplot(climate.value, mat, type = "l", col = cols, lty = linetype, lwd = linewidth, ylab = "density", xlab = "climate.value")
        title(i)
        legend(legend.pos, legend = leg, col = cols, pch=19)
      }
    }
  }
  ###############################
  # if the output of machu.2.ace is provided
  if (typeof(resp) == "list"){
    # use all taxa by default
    if (is.null(taxa) == TRUE){
      taxa <- 1:nrow(resp[[1]])
    }
    # use all climate variables by default
    if (is.null(clim) == TRUE){
      q <- names(resp)
      q <- str_remove(q, "_mean"); q <- str_remove(q, "_stdev"); q <- str_remove(q, "_skew"); q <- str_remove(q, "_lowerQ"); q <- str_remove(q, "_upperQ");
      clim <- unique(q)
    }
    # if comb is FALSE
    if (comb == FALSE){
      # set number of panels
      if (is.null(par.override) == TRUE){
        par(mfrow = c(length(taxa), length(clim)))
      } else if (is.null(par.override) == FALSE){
        par(mfrow = par.override)
      }
      # if the taxon set was provided as a character (not integers), convert to integers:
      if (typeof(taxa) == "character"){taxa <- match(taxa, as.character(resp[[1]]$name))}
      # for every taxon specified
      for (i in taxa){
        # and for every climate variable specified
        for (j in clim){
          # get climate response values for taxon and climate variable
          mu <- resp[[paste0(j, "_mean")]]$value[i]
          sd <- resp[[paste0(j, "_stdev")]]$value[i]
          xi <- resp[[paste0(j, "_skew")]]$value[i]
          # make the plot
          climate.value <- seq(round(mu-sd*3), round(mu+sd*3))
          density <- dsnorm(climate.value, mean = mu, sd = sd, xi = xi)
          plot(climate.value, density, type="l", lty=linetype, lwd=linewidth, col=linecolor)
          title(paste(as.character(resp[[1]]$name[i]), j))
        }
      }
      # if comb is TRUE
    } else if (comb == TRUE){
      # set number of panels
      if (is.null(par.override) == TRUE){
        par(mfrow = c(1, length(clim)))
      } else if (is.null(par.override) == FALSE){
        par(mfrow = par.override)
      }
      # if the taxon set was provided as a character (not integers), convert to integers:
      if (typeof(taxa) == "character"){taxa <- match(taxa, as.character(resp[[1]]$name))}
      # for every climate variable specified
      for (i in clim){
        # create the x-bounds of the plot
        clim.mu <- resp[[paste0(i, "_mean")]]$value[taxa]
        clim.sd <- resp[[paste0(i, "_stdev")]]$value[taxa]
        lobound <- min(clim.mu-clim.sd*3)
        upbound <- max(clim.mu+clim.sd*3)
        climate.value <- seq(round(lobound), round(upbound))
        vec <- c()
        # for every taxon
        for (j in taxa){
          # get climate response values for taxon and climate variable
          mu <- resp[[paste0(i, "_mean")]]$value[j]
          sd <- resp[[paste0(i, "_stdev")]]$value[j]
          xi <- resp[[paste0(i, "_skew")]]$value[j]
          vec <- append(vec, dsnorm(climate.value, mean = mu, sd = sd, xi = xi))
        }
        # create data matrix
        mat <- matrix(vec, nrow = length(climate.value), ncol = length(taxa))
        # create color ramp
        # if there are more taxa than colors in the Set2 palette (8), use a colorRampPalette interpolation
        if (length(taxa) > 8){
          cols <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(length(taxa))
        } else {
          cols <- brewer.pal(n = length(taxa), name = "Set2")
        }
        # get names for the legend
        leg <- as.character(resp[[1]]$name)[taxa]
        # make the plot
        matplot(climate.value, mat, type = "l", col = cols, lty = linetype, lwd = linewidth, ylab = "density", xlab = "climate.value")
        title(i)
        legend(legend.pos, legend = leg, col = cols, pch=19)
      }
    }
  }
  par(mfrow=c(1,1))
}

#'Perform ancestral character estimation of response curves
#'
#'Use ace() to generate climate response curves at each node of a
#'time-calibrated phylogeny, and extract the values along the branches
#'at a particular time if so desired. By default, the output is a set
#'of climate response curves for each node (and each tip). As each node
#'occurs at different times (making paleoclimatic projections tricky),
#'the user can also provide a time-slice value T. The function will
#'find each branch present at time T and interpolate response curves
#'along them to that point, then return those values. The function can
#'also estimate uncertainty with the n.unc parameter.
#'
#'@param tip.resp a table containing the character data for each tip.
#'Corresponds to the output of machu.1.tip.resp().
#'@param tree an ultrametric phylogeny, formatted as a tree object. Tip
#'labels must be identical to those in the tip.resp table. Note that
#'treedata objects do not work (use as.phylo() to convert).
#'@param T time from which to extract reconstructed values along any
#'branches present at T. Must be a value between 0 and the root age of
#'the tree. Should correspond to the estimated date of paleoclimatic
#'layers for best results. If T is not provided, the function will output
#'climate response curves for each node (as well as each tip).
#'@param n.unc a whole, positive number that specifies the number of evenly-spaced
#'samples of the 95% confidence interval to take for characterizing uncertainty in
#'ancestral character estimations. n.unc sampled values will be tacked on to the end
#'of the output, and in machu.3.ace(), specifying calc.unc = T will construct separate
#'models for each sample and average them to get the final model.
#'@param savename a string that, when provided, saves the output as a .csv file. It is not
#'necessary to add the file extension to the string. Use machu.ace.load() to load one of these
#'saved files.
#'
#'@return a list where each element corresponds to a response curve parameter for a
#'certain climate variable, and contains the values of that parameter for each taxon
#'in question. If T is not provided (default), the taxa returned are the tips and nodes
#'of the provided tree. If T is provided, the taxa returned are the branches present
#'at time T. If n.unc is provided, n.unc additional models will also be generated
#'based on the 95% confidence intervals of each value. These additional models will
#'be used in step 3 to incorporate uncertainty in ancestral character estimation
#'into the process by generating an ensemble model.
#'
#'@examples
#'# get model for each node, without uncertainty
#'ace.out <- machu.2.ace(response.table, bassleritree)
#'# get model for each lineage at T=2.5, without uncertainty
#'ace.out <- machu.2.ace(response.table, bassleritree, T=2.5)
#'# get model for each node, with 4 uncertainty samples
#'ace.out <- machu.2.ace(response.table, bassleritree, n.unc=4)
#'# get model for each lineage at T=2.5, with 4 uncertainty samples
#'ace.out <- machu.2.ace(response.table, bassleritree, T=2.5, n.unc=4)
#'@import ape
#'@export
machu.2.ace <- function(tip.resp, tree, T = NULL, n.unc=NULL, savename=NULL){

  # check if tree taxa match tip.response table
  if (identical(tree$tip.label, rownames(tip.resp[match(tree$tip.label, rownames(tip.resp)),])) == TRUE){
    # check if tree is ultrametric
    if (is.ultrametric(tree) == TRUE){

      ###################
      ### 1: Get data ###
      ###################

      # Sort the trait data file so that the order of rows matches the order of the tips
      # in the tree. This should prevent any misattribution errors down the line, because
      # I tested the program on a test dataset sorted in such a manner.
      response.table.sorted <- as.data.frame(tip.resp[match(tree$tip.label, rownames(tip.resp)),])
      traitnames <- c(colnames(response.table.sorted))

      # This creates a list of named vectors, where each vector is a column from the trait data file, e.g. a trait.
      tiptraits <- apply(response.table.sorted, 2, list)

      # This runs ace on each of the traits and puts the output for each trait into a list,
      # so that each element can be easily recalled later.
      acelist <- list()
      for (i in 1:length(tiptraits)) acelist[[i]] <- list(ace(tiptraits[[i]][[1]], tree, type = "continuous", method = "ML"))
      names(acelist) <- c(names(tiptraits))

      #############################
      ### 2: Create data matrix ###
      #############################

      ntips <- Ntip(tree)      # Get number of tips on the tree
      nnodes <- Nnode(tree)    # Get number of nodes in the tree

      treenames <- tree$tip.label  # improve readability with label vector
      for (i in (ntips+1):(ntips+nnodes)) treenames[i] <- paste0("Node", i)  # adds node number to names vector

      times <- c()                 # make empty vector (will be filled with branching times)
      times[1:ntips] <- 0          # all tips are assigned a branching time of 0
      times[(ntips+1):(ntips+nnodes)] <- branching.times(tree)   # gives branching time of nodes

      # This makes a new list of traits which is a modified version of tiptraits. Besides
      # containing trait values for tips, it appends the node values (retrieved from acelist)
      # for each trait as well.
      alltraits <- list()
      for (i in 1:length(tiptraits)) {
        alltraits[[i]] <- append(tiptraits[[i]][[1]], acelist[[i]][[1]]$ace)
      }
      names(alltraits) <- traitnames

      #calculate table with uncertainty (n.unc is provided)
      if (is.null(n.unc) == FALSE) {

        # initialize table and names vector
        dfCI <- matrix(nrow=(ntips+nnodes), ncol=2)
        dfCI.names <- c()

        # this loop generates both the table of confidence intervals for each trait, derived from
        # acelist[[i]][[1]]$CI95, as well as the list of names
        for (i in 1:length(acelist)){
          # generate the table, where each column is the uncertainty (lower CI95 or upper CI95)
          # at each node. Tips are included so the table can be melded with the big trait table later.
          # However, tips are given a value of zero.
          lCI <- c(tiptraits[[i]][[1]], acelist[[i]][[1]]$CI95[,1])
          uCI <- c(tiptraits[[i]][[1]], acelist[[i]][[1]]$CI95[,2])
          dfCI <- data.frame(dfCI, cbind(lCI, uCI))

          # generate the names for the table in vector form.
          tmp.names <- c(paste0(traitnames[i], "_ace_lCI"), paste0(traitnames[i], "_ace_uCI"))
          dfCI.names <- append(dfCI.names, tmp.names)
        }
        dfCI <- dfCI[,3:ncol(dfCI)]  # truncate the table because the first two columns will be NAs
        colnames(dfCI) <- dfCI.names # name the columns using the name vector from the loop

        # make table with uncertainty
        # initialize table and a temp version of the dfCI table, as well as name vector(s)
        table <- data.frame(treenames,times)
        dfCI.tmp <- dfCI
        table.names <- c()
        dfCI.names.tmp <- dfCI.names

        # this loop makes the table with uncertainty. The structure is three columns per trait:
        # the first column is the trait value at each node and tip, the second column is the lower CI95
        # at each node (with 0 at each tip), and the third column is the upper CI95 at each node.
        for (i in 1:length(alltraits)) {

          # this loop only takes the first two columns of the dfCI.tmp table. Also starts the name vector.
          table <- data.frame(table, alltraits[[i]], dfCI.tmp[,1:2])
          names.tmp <- c(names(alltraits[i]), dfCI.names.tmp[1:2])
          names.tmp
          table.names <- append(table.names, names.tmp)
          table.names

          # after taking the first two columns (corresponding to the lower and upper CI95s of the
          # current trait, "trim" these columns from the dfCI.tmp table, so that for every loop iteration,
          # the first two columns correspond to the trait denoted by the loop's index.
          # The if statement is there to prevent error messages.
          if (ncol(dfCI.tmp) > 2) {
            dfCI.tmp <- dfCI.tmp[,3:ncol(dfCI.tmp)]
          }
          if (length(dfCI.names.tmp) > 2) {
            dfCI.names.tmp <- dfCI.names.tmp[3:length(dfCI.names.tmp)]
          }
        }
        # Name the columns in the table
        colnames(table)[3:ncol(table)] <- table.names

      } else {
        # calculate table without uncertainty (n.unc is not provided)
        # This starts the data matrix and adds in the data contained in the alltraits list.
        table <- data.frame(treenames,times)
        for (i in 1:length(alltraits)) {
          table <- data.frame(table, alltraits[[i]])
        }
        # Name the columns in the table, each corresponding to a trait
        colnames(table)[3:ncol(table)] <- traitnames
      }

      # if T is not provided, don't take time-slice:
      if (is.null(T) == TRUE) {

        #############################################
        ### 3a: Return trait values for each node ###
        #############################################

        # This is the default option, which just returns a list of dataframes, where each dataframe corresponds
        # to a trait. Each trait's dataframe contains the values of that trait for all tips and all ancestors.
        # Later this can be given to machu.3.anc.niche to reconstruct node or tip models at different times.

        # initialize list
        ace.output <- list()

        # for every trait column (starting at column 3) in "table":
        for (j in 3:(ncol(table))){
          # grab relevant info and trait values and put them into a data frame for each trait.
          ace.output[[(j-2)]] <- data.frame(table$treenames, table$times, c(rep(NA, nrow(table))), table[,j])
          # change column names
          colnames(ace.output[[(j-2)]]) <- c("name", "time", "NA", "value")
        }

        # incorporate uncertainty if n.unc is provided
        if (is.null(n.unc) == FALSE) {

          # applies trait names to output list (including CIs for now)
          names(ace.output) <- table.names

          # transforms output to usable format for next step
          s<-seq(from=1, to=length(ace.output), by=3)
          ace.modified <- list()

          # for every trait (excluding lCI and uCI values)
          for (l in 1:length(s)){

            # copy the trait data (excluding lCI/uCI) to the modified ace list
            ace.modified[l] <- ace.output[s[l]]

            # for every taxon in ace.output
            for (j in 1:nrow(ace.output[[1]])){
              lCI <- ace.output[[(s[l]+1)]][j,4] # lower confidence value
              uCI <- ace.output[[(s[l]+2)]][j,4] # upper confidence value
              diff <- abs(uCI-lCI) # difference between confidence values
              inc <- diff/(n.unc-1)    # amount to increment by

              # for every bin, the number of which is assigned by the user
              for (i in 0:(n.unc-1)){
                ace.modified[[l]][j,(i+5)] <- lCI+(i*inc) # add an incrementing value to the corresponding trait table
                colnames(ace.modified[[l]])[(i+5)] <- paste0("unc", i+1)
              }
            }
          }
          ace.output <- ace.modified
          names(ace.output) <- traitnames

        } else {
          # rename dataframes to correspond to traits
          names(ace.output) <- traitnames
        }
        if (is.null(savename) == FALSE){
          for (i in 1:length(ace.output)){
            climvar <- rep(names(ace.output[i]), nrow(ace.output[[i]]))
            dfred <- data.frame(ace.output[[i]], climvar)
            if (i == 1){
              savethis <- dfred
            } else {
              savethis <- rbind(savethis, dfred)
            }
          }
          # naming the file
          write.csv(savethis, file=paste0(savename, ".csv"))
        }
        return(ace.output)

      # if T is provided, take time-slice:
      } else {

        ######################################################
        ### 3b: Calculate trait values at a desired time T ###
        ######################################################

        edge <- tree$edge        # Save edge matrix from ape

        # check if T is greater than branch height
        if (T <= max(branching.times(tree))) {

          # Below, the outer loop loops over each trait (corresponding to a column in "table").
          # The inner loop, for each row in the edge matrix (where each row consists of a node and a node/tip it is
          # connected to), tests whether T is between the branching time values for those nodes/tips. Essentially,
          # it is making sure to only return an interpolated ancestral trait value for branches that existed at time
          # T. The loop calculates values of v, which is the branch distance from T to the node/tip on either side
          # of it. It then plugs v and x (trait values) into an equation to calculate a, the interpolated value of
          # the ancestral state at that point in time.
          # Eventually, the output is given as a list of matrices, where the first two rows correspond to the origins
          # and ends of each branch that crosses time T, and the third row is the interpolated trait value at T.

          # if n.unc is provided, the confidence interval values are also interpolated at T. However, the
          # calculation assumes that the trait at the tip is known without error, as the confidence interval
          # values for each tip are set to the tip value itself. So any interpolation that occurs between a
          # node and a tip is calculated between the node confidence interval value and the tip trait value.

          # initialize list
          ace.output <- list()

          for (j in 3:(ncol(table))) {

            tempv_br_name  <- c()  # initialize four data vectors
            tempv_br_start <- c()  # with which we will
            tempv_br_end   <- c()  # later create a data frame
            tempv_a        <- c()  # containing our data

            for (i in 1:nrow(edge)) {

              if (table[edge[i,1],2] > T && T > table[edge[i,2],2]) {

                xi <- table[edge[i,1],j]               # trait value for tip/node i
                xj <- table[edge[i,2],j]               # trait value for tip/node j
                vi <- abs(T-table[edge[i,1],2])        # distance from tip/node i to T
                vj <- abs(T-table[edge[i,2],2])        # distance from tip/node j to T
                a <- ((xi/vi)+(xj/vj))/((1/vi)+(1/vj)) # calculates a, the value of the trait along each branch that exists at T

                # the following 4 lines append, to the four temp vectors we initilialized at the start of the loop, the respective
                # values of each (branch name, branch start, branch end (these from the edge matrix), and the value of a at T).
                tempv_br_name  <- append(tempv_br_name, paste0(table[edge[i,1],1], "-", table[edge[i,2],1]))
                tempv_br_start <- append(tempv_br_start, edge[i,1])
                tempv_br_end   <- append(tempv_br_end, edge[i,2])
                tempv_a        <- append(tempv_a, a)
              }
            }
            ace.output[[(j-2)]] <- data.frame(tempv_br_name, tempv_br_start, tempv_br_end, tempv_a)        # create output data frame
            colnames(ace.output[[(j-2)]]) <- c("name", "branch_start", "branch_end", "value")  # apply column names to each frame
          }
          # apply uncertainty if n.unc is provided
          if (is.null(n.unc) == FALSE) {

            # applies trait names to output list (including CIs for now)
            names(ace.output) <- table.names

            # transforms output to usable format for next step
            s<-seq(from=1, to=length(ace.output), by=3)
            ace.modified <- list()

            # for every trait (excluding lCI and uCI values)
            for (l in 1:length(s)){
              ace.modified[l] <- ace.output[s[l]] # copy the trait data (excluding lCI/uCI) to the modified ace list

              # for every ancestor at T
              for (j in 1:nrow(ace.output[[1]])){
                lCI <- ace.output[[(s[l]+1)]][j,4] # lower confidence value
                uCI <- ace.output[[(s[l]+2)]][j,4] # upper confidence value
                diff <- abs(uCI-lCI) # difference between confidence values
                inc <- diff/(n.unc-1)    # amount to increment by

                # for every bin, the number of which is assigned by the user
                for (i in 0:(n.unc-1)){
                  ace.modified[[l]][j,(i+5)] <- lCI+(i*inc) # add an incrementing value to the corresponding trait table
                  colnames(ace.modified[[l]])[(i+5)] <- paste0("unc", i+1)
                }
              }
            }
            ace.output <- ace.modified
            names(ace.output) <- traitnames
          } else {
            names(ace.output) <- traitnames
          }
          if (is.null(savename) == FALSE){
            for (i in 1:length(ace.output)){
              climvar <- rep(names(ace.output[i]), nrow(ace.output[[i]]))
              dfred <- data.frame(ace.output[[i]], climvar)
              if (i == 1){
                savethis <- dfred
              } else {
                savethis <- rbind(savethis, dfred)
              }
            }
            # naming the file
            write.csv(savethis, file=paste0(savename, ".csv"))
          }
          return(ace.output)
        } else {
          print("T is greater than the total tree height.")
        }
      }
    } else {
      print("tree is not ultrametric.")
    }
  } else {
    print("tree names and response table names do not match. Make sure they match exactly.")
  }
}

#'Load saved ancestral character estimations from .csv file
#'
#'Using the savename parameter in machu.2.ace(), you can save output
#'as a .csv. This function allows you to load this back into the special
#'list format which is necessary for machu.3.anc.niche().
#'
#'@param file filename (or filepath) of the .csv saved by machu.2.ace().
#'
#'@return correctly-formatted output of machu.2.ace().
#'
#'@examples
#'load.ace <- machu.ace.load("saved_ace.csv")
#'
#'@export
machu.ace.load <- function(file){
  # loading ace list from file and converting back to proper format
  aceload <- read.csv(file=file)
  climvars <- unique(aceload$climvar)
  acelist <- list()
  for (j in 1:length(climvars)){
    acelist[[j]] <- filter(aceload, climvar == climvars[j])[,2:5]
  }
  names(acelist) <- as.vector(climvars)
  return(acelist)
}

###############################################################
#######################    STEP 3   ###########################
#    convert ancestral states to spatial predictions          #
###############################################################

#'Convert ancestral response values to niche models
#'
#'Create a suitability map for each ancestor at time T, then convert
#'into a Bioclim model for each ancestor.
#'
#'@param ace a list of tables with the reconstructed response values
#'for each ancestor at T. Corresponds to the output of machu.2.ace().
#'@param clim paleoclimatic data, formatted as a RasterStack. Should
#'be as close to time T as possible for best results.
#'@param taxa optionally select specific taxa from ace to create models
#'for. This is useful if you wish to visualize only one specific taxon or
#'node at a certain time. The form is a row index or vector of indices. For
#'example, taxa = 1 will only create a model for the first taxon, and
#'taxa = 1:3 will create models for the first three. taxa = c(1,3:5) will
#'create models for the first taxon and for taxa 3 to 5.
#'@param resp.curv if TRUE, create ancestral models using response curves.
#'If false, assume response curve is a uniform distribution between the LCI
#'and UCI response values. Default = TRUE. Setting resp.curv = FALSE will
#'also automatically set clip.Q = FALSE.
#'@param clip.Q if TRUE, clip the tails of each response curve at .025 and
#'0.975. May produce cleaner models. Default = TRUE.
#'@param calc.unc if TRUE, incorporate uncertainty bins from using n.unc
#'in machu.2.ace(). Visualizes uncertainty in the model. Default = FALSE.
#'@param verbose if TRUE, print to screen the current taxon and climate
#'layer being processed. Default = FALSE.
#'
#'@return a list consisting of one raster per ancestor at T, showing the
#'habitat suitability for each pixel on the landscape.
#'
#'@examples
#'#run w/ default settings (resp.curv & clip.Q = TRUE, do not account for uncertainty)
#'output.models <- machu.3.anc.niche(ace.output, clim)
#'#assume response curves are uniform distributions between min and max
#'output.models <- machu.3.anc.niche(ace.output, clim, resp.curv=F)
#'#Do not clip response curve tails
#'#Note that clip.Q cannot be true if resp.curv is false.
#'output.models <- machu.3.anc.niche(ace.output, clim, clip.Q=F)
#'#incorporate uncertainty (make sure to specify a value for n.unc in machu.2.ace() beforehand)
#'output.models <- machu.3.anc.niche(ace.out.N, ClimMis19, calc.unc = T)
#'@import fGarch
#'@import raster
#'@export
machu.3.anc.niche <- function(ace, clim, taxa = NULL, resp.curv = TRUE, clip.Q = TRUE, calc.unc = FALSE, verbose = FALSE) {
  #check if the paleoclimate data is a raster stack
  if (class(clim)[1] == "RasterStack"){

    # initialize output
    output.bioclim <- list()

    # input rasterstack (climate layers), extract names, and make sequence corresponding to number of layers
    clim.names <- names(clim)
    n.in <- c(1:length(clim.names))
    traitnames <- names(ace)

    # if specific taxa are selected, modify ace to only include those taxa
    if (is.null(taxa) == FALSE){
      ace <- lapply(ace, function(x) return(x[taxa,]))
    }

    # for each ancestral taxon at T (corresponding to a row in each entry of the ace list):
    for (j in 1:nrow(ace[[1]])){

      # initialize list
      output.tmp <- list()

      # and for each climate layer, do:
      for (i in n.in){

        # get the correct names for each variable
        t.va <- c((((i-1)*5)+1):(i*5))     # first get the columns of the trait table corresponding to the current climate layer
        name.rc.parms <- traitnames[t.va]  # use these to extract names from the trait table
        out.bio <- clim.names[i]           # store output climate layer name

        # incorporate uncertainty into analysis
        if (calc.unc == TRUE) {

          # convert the raster corresponding to the current climate layer to a matrix
          in.raster.p <- paste0("clim$",out.bio)
          raster <- as.matrix(eval(parse(text=in.raster.p)))
          raster[is.na(raster)] <- 0
          raster <- c(raster)

          # initialize zero vector with the same number of components as raster
          projVal <- rep(0, length(raster))

          # for each uncertainty "bin" (specified in machu.2.ace by n.unc and given by each column after column 4 in ace.output):
          for (h in 5:ncol(ace[[1]])){
            # mean
            in.m.v <- paste0("ace$", name.rc.parms[1], "[", j, ",", h, "]")
            m.val <- eval(parse(text = in.m.v))
            # stdev
            in.sd.v <- paste0("ace$", name.rc.parms[2], "[", j, ",", h, "]")
            sd.val <- eval(parse(text = in.sd.v))
            # skew
            in.xi.v <- paste0("ace$", name.rc.parms[3], "[", j, ",", h, "]")
            xi.val <- eval(parse(text = in.xi.v))
            # lower quantile
            in.qL.v <- paste0("ace$", name.rc.parms[4], "[", j, ",", h, "]")
            qL.val <- eval(parse(text = in.qL.v))
            # upper quantile
            in.qU.v <- paste0("ace$", name.rc.parms[5], "[", j, ",", h, "]")
            qU.val <- eval(parse(text = in.qU.v))

            # project raster into suitability (is there a better explanation for this step?)
            # if resp.curve = T, reconstruct normal distribution from suitability values using m, sd, and xi
            if(resp.curv == TRUE){
              projVal.tmp <- dsnorm(raster, mean=m.val, sd=sd.val, xi=xi.val)}
            # if resp.curve = F, clip values below/above lower/upper quantiles and set all values inside quantiles to 1 (uniform dist.)
            if(resp.curv == FALSE){
              projVal.tmp <- raster
              projVal.tmp[raster <= qL.val] <- 0
              projVal.tmp[raster >= qU.val] <- 0
              projVal.tmp[raster <  qU.val & raster >qL.val] <- 1
              clip.Q = FALSE}

            # add projVal.tmp to the zero vector projVal generated before this loop. Basically this loop
            # generates a suitability vector for every "uncertainty bin" and sums them all. After the loop
            # the total suitability vector projVal will be rescaled to 1.
            projVal <- projVal + projVal.tmp
          }
          # if resp.curve = T was selected, and you wish to clip the suitability values by quantiles, do so
          # if resp.curve = F was selected, quantiles are automatically clipped in the previous step
          if(clip.Q == TRUE){projVal[raster <= qL.val] <- 0
          projVal[raster >= qU.val] <- 0}

          # rescale suitability values from 0-1 so that they can be compared
          max <- max(projVal)
          min <- min(projVal)
          pred_value <- 1/(max-min)*(projVal-min)

          #plots of values- not necssary later on
          ##plot(raster,pred_value, main=q.val)

          # turn list back into matrix
          xy <- matrix(pred_value,eval(parse(text=in.raster.p)) @nrows, eval(parse(text=in.raster.p)) @ncols)

          # turn the matrix into a raster
          rast <- raster(xy)
          extent(rast) <- eval(parse(text=in.raster.p)) @extent
          projection(rast) <- eval(parse(text=in.raster.p)) @crs

          #add each raster to the temp output list (this corresponds to one ancestor)
          output.tmp[[(i)]] <- rast

          # when each layer for each ancestor is processed, print the ancestor and layer to the screen
          if (verbose == TRUE){
            print(paste(ace[[1]]$name[j], out.bio, "--- with uncertainty"))
          }

        } else {

          # extract the reconstructed value for each parameter from the ace list
          # mean
          in.m.v  <- paste0("ace$",name.rc.parms[1],"$value")
          m.val   <- eval(parse(text = in.m.v))[j]
          # stdev
          in.sd.v <- paste0("ace$",name.rc.parms[2],"$value")
          sd.val  <- eval(parse(text = in.sd.v))[j]
          # skew
          in.xi.v <- paste0("ace$",name.rc.parms[3],"$value")
          xi.val  <- eval(parse(text = in.xi.v))[j]
          # lower quantile
          in.qL.v <- paste0("ace$",name.rc.parms[4],"$value")
          qL.val  <- eval(parse(text = in.qL.v))[j]
          # upper quantile
          in.qU.v <- paste0("ace$",name.rc.parms[5],"$value")
          qU.val  <- eval(parse(text = in.qU.v))[j]

          # convert the raster corresponding to the current climate layer to a matrix
          in.raster.p <- paste0("clim$",out.bio)
          raster <- as.matrix(eval(parse(text=in.raster.p)))
          raster[is.na(raster)] <- 0
          raster <- c(raster)

          # project raster into suitability (is there a better explanation for this step?)
          # if resp.curve = T, reconstruct normal distribution from suitability values using m, sd, and xi
          if(resp.curv == TRUE){
            projVal <- dsnorm(raster, mean=m.val, sd=sd.val, xi=xi.val)
          }
          # if resp.curve = F, clip values below/above lower/upper quantiles and set all values inside quantiles to 1 (uniform dist.)
          if(resp.curv==FALSE){
            projVal <- raster
            projVal[raster <= qL.val] <- 0
            projVal[raster >= qU.val] <- 0
            projVal[raster < qU.val & raster >qL.val] <- 1
            clip.Q = FALSE
          }

          # if resp.curve = T was selected, and you wish to clip the suitability values by quantiles, do so
          # if resp.curve = F was selected, quantiles are automatically clipped in the previous step
          if(clip.Q == TRUE){
            projVal[raster <= qL.val] <- 0
            projVal[raster >= qU.val] <- 0
          }

          # rescale suitability values from 0-1 so that they can be compared
          max <- max(projVal)
          min <- min(projVal)
          pred_value <- 1/(max-min)*(projVal-min)

          #plots of values- not necssary later on
          ##plot(raster,pred_value, main=q.val)

          # turn list back into matrix
          xy <- matrix(pred_value,eval(parse(text=in.raster.p)) @nrows, eval(parse(text=in.raster.p)) @ncols)

          # turn the matrix into a raster
          rast <- raster(xy)
          extent(rast) <- eval(parse(text=in.raster.p)) @extent
          projection(rast) <- eval(parse(text=in.raster.p)) @crs

          #add each raster to the temp output list (this corresponds to one ancestor)
          output.tmp[[(i)]] <- rast

          # when each layer for each taxon is processed, print the taxon and layer to the screen
          if (verbose == TRUE){
            print(paste(out.bio, "processed for", ace[[1]]$name[j]))
          }
        }
      }
      # add bioclim names to corresponding layers
      names(output.tmp) <- clim.names

      # "brick" each climate layer for each species
      rastbrick.tmp <- brick(output.tmp)

      # add each raster brick to the main output list
      output.bioclim[[(j)]] <- rastbrick.tmp
    }
    #add ancestor names to output.bioclim
    names(output.bioclim) <- as.vector(ace[[1]]$name)

    ###############################################################
    #######################    STEP 4    ##########################
    #     create BIOCLIM model from all spatial predictions       #
    ###############################################################

    # initialize model list
    output.models <- list()

    # for each raster brick (corresponding to an ancestor at T) in output.bioclim:
    for (i in 1:length(output.bioclim)){

      # find limiting values from all rasters
      anc.niche <- calc(output.bioclim[[i]], function(x) min(x, na.rm = TRUE))

      # create a classification scheme for masking oceans out
      rcl.mask <- as.matrix(c(0,max(raster),1))

      # create actual raster of backgroud using classication scheme
      background.mask <- reclassify((eval(parse(text=in.raster.p))), rcl=rcl.mask)

      # use mask on final ancestral distribution
      anc.niche.f <- background.mask*anc.niche  # find limiting values from all rasters

      # store the model for each ancestor in output.model list and add names
      output.models[[(i)]] <- anc.niche.f

      if (verbose == TRUE){
        print(paste("Finished processing", ace[[1]]$name[i]))
      }
    }
    # add names corresponding to ancestor to each model in output.model
    names(output.models) <- as.vector(names(output.bioclim))
  } else {
    print("clim is not a RasterStack.")
  }
  return(output.models)
}

#'Map model output from machu.3.anc.niche()
#'
#'Plot one or more maps using model output from machu.3.anc.niche().
#'User can specify color ramps from the viridis package.
#'
#'@param model output from machu.3.anc.niche(). A list of rasters representing
#'the models to display. Using the full output displays each model successively.
#'User can also specify a subset of models using single-bracket [x] notation to
#'only plot certain models or plot one model at a time.
#'@param col color ramp to use, a numeric value from 1 to 6. 1 = base, 2 =
#'plasma, 3 = viridis, 4 = l17, 5 = r3, 6 = white-to-black. Default = 1.
#'@param title if TRUE, print the taxon name above each corresponding model.
#'Default = TRUE.
#'
#'@return a series of maps corresponding to each output model specified.
#'
#'@examples
#'#plot all models with default color ramp and titles:
#'machu.plotmap(output.models)
#'#plot only the first model and use the "plasma" color ramp:
#'machu.plotmap(output.models[1], color=2)
#'#plot only the second and third models and do not display titles:
#'machu.plotmap(output.models[2:3], title=F)
#'#display three models at the same time
#'par(mfrow=c(1,3))
#'machu.plotmap(output.models)
#'#this function is designed for ease of use. For more control over plotting,
#'#use the image() function. Example below:
#'image(output.models[[1]], col=viridis(256), asp=1, axes=FALSE, xaxs="i", xaxt='n', yaxt='n', ann=F)
#'@import raster
#'@export
machu.plotmap <- function(model, col=1, title=TRUE) {
  # pick color
  #ORDER:1.base, 2.plasma, 3.viridis, 4.l17, 5.r3, 6. b&w
  colr.o <- list()
  if (col == 1) {
    #	base - pure
    colr.o <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  }
  if (col == 2) {
    #plasma - pure
    colr.o<- c("#000000", "#000C7D", "#000D7E", "#000D80", "#000E82", "#000E84", "#000E86", "#000F87", "#000F89", "#00108B", "#00108C", "#00118E", "#001190", "#001191", "#001293", "#001294", "#001296", "#001397", "#001399", "#00139A", "#00149B", "#00149D", "#00149E", "#00149F", "#0015A0", "#0015A1", "#0015A2", "#0015A3", "#0015A4", "#0015A5", "#0016A6", "#0016A7", "#0016A7", "#0016A8", "#0016A9", "#0016A9", "#0016A9", "#0A16AA", "#1516AA", "#1D15AA", "#2315AA", "#2915AA", "#2F15A9", "#3414A9", "#3914A8", "#3E13A7", "#4313A6", "#4712A5", "#4C12A4", "#5011A3", "#5311A2", "#5710A1", "#5A0FA0", "#5E0F9F", "#610E9F", "#640E9E", "#670D9D", "#6A0D9C", "#6C0C9B", "#6F0B9A", "#720B99", "#740A99", "#770A98", "#790997", "#7C0897", "#7E0896", "#800795", "#820795", "#850694", "#870693", "#890693", "#8B0592", "#8D0592", "#8F0491", "#910491", "#930490", "#950390", "#97038F", "#99038F", "#9B028E", "#9D028E", "#9F028E", "#A1018D", "#A3018D", "#A5018C", "#A7018C", "#A9008B", "#AB008B", "#AC008A", "#AE008A", "#B00089", "#B20089", "#B40088", "#B60088", "#B80087", "#B90087", "#BB0087", "#BD0086", "#BF0086", "#C10085", "#C30085", "#C40084", "#C60084", "#C80083", "#CA0083", "#CC0082", "#CE0082", "#CF0081", "#D10081", "#D30080", "#D50080", "#D6007F", "#D8007F", "#DA007E", "#DB007E", "#DD007D", "#DE007C", "#E0017C", "#E2027B", "#E3047B", "#E5067A", "#E6087A", "#E80B79", "#E90D78", "#EA1078", "#EC1277", "#ED1477", "#EE1676", "#F01875", "#F11A75", "#F21C74", "#F41E73", "#F52073", "#F62272", "#F72471", "#F82671", "#F92870", "#FB2A6F", "#FC2C6F", "#FD2E6E", "#FE306D", "#FF326C", "#FF346C", "#FF366B", "#FF386A", "#FF3A6A", "#FF3D69", "#FF3F68", "#FF4167", "#FF4366", "#FF4566", "#FF4765", "#FF4964", "#FF4B63", "#FF4D62", "#FF5062", "#FF5261", "#FF5460", "#FF565F", "#FF585E", "#FF5A5D", "#FF5D5C", "#FF5F5B", "#FF615B", "#FF635A", "#FF6559", "#FF6758", "#FF6A57", "#FF6C56", "#FF6E55", "#FF7054", "#FF7253", "#FF7452", "#FF7651", "#FF7850", "#FF7A4E", "#FF7C4D", "#FF7E4C", "#FF7F4B", "#FF814A", "#FF8349", "#FF8548", "#FF8747", "#FF8845", "#FF8A44", "#FF8C43", "#FF8E42", "#FF8F40", "#FF913F", "#FF933E", "#FF953C", "#FF963B", "#FF9839", "#FF9A38", "#FF9B36", "#FF9D35", "#FF9F33", "#FFA032", "#FFA230", "#FFA32F", "#FFA52E", "#FFA62C", "#FFA82B", "#FFA92A", "#FFAB29", "#FFAC28", "#FFAE27", "#FFAF26", "#FFB126", "#FFB225", "#FFB424", "#FFB523", "#FFB623", "#FFB822", "#FFB922", "#FFBB21", "#FFBC20", "#FFBD20", "#FFBF1F", "#FFC01F", "#FFC21F", "#FFC31E", "#FFC41E", "#FFC61E", "#FFC71D", "#FFC81D", "#FFCA1D", "#FFCB1D", "#FFCC1D", "#FFCE1D", "#FFCF1C", "#FFD01C", "#FFD21C", "#FFD31C", "#FFD41C", "#FFD61C", "#FFD71D", "#FFD81D", "#FFDA1D", "#FFDB1D", "#FFDC1D", "#FFDE1D", "#FFDF1E", "#FFE01E", "#FFE21E", "#FFE31E", "#FFE41F", "#FFE61F", "#FFE71F", "#FFE820", "#FFE920", "#FFEB21", "#FFEC21", "#FFED22", "#FFEF22", "#FFF123")
  }
  if (col == 3) {
    #viridis - pure
    colr.o <- colorRampPalette(c("#440154", "#440154ff", "#440558ff", "#450a5cff", "#450e60ff", "#451465ff", "#461969ff", "#461d6dff", "#462372ff", "#472775ff", "#472c7aff", "#46307cff", "#45337dff", "#433880ff", "#423c81ff", "#404184ff", "#3f4686ff", "#3d4a88ff", "#3c4f8aff", "#3b518bff", "#39558bff", "#37598cff", "#365c8cff", "#34608cff", "#33638dff", "#31678dff", "#2f6b8dff", "#2d6e8eff", "#2c718eff", "#2b748eff", "#29788eff", "#287c8eff", "#277f8eff", "#25848dff", "#24878dff", "#238b8dff", "#218f8dff", "#21918dff", "#22958bff", "#23988aff", "#239b89ff", "#249f87ff", "#25a186ff", "#25a584ff", "#26a883ff", "#27ab82ff", "#29ae80ff", "#2eb17dff", "#35b479ff", "#3cb875ff", "#42bb72ff", "#49be6eff", "#4ec16bff", "#55c467ff", "#5cc863ff", "#61c960ff", "#6bcc5aff", "#72ce55ff", "#7cd04fff", "#85d349ff", "#8dd544ff", "#97d73eff", "#9ed93aff", "#a8db34ff", "#b0dd31ff", "#b8de30ff", "#c3df2eff", "#cbe02dff", "#d6e22bff", "#e1e329ff", "#eae428ff", "#f5e626ff", "#fde725ff")) (256)
  }
  if (col == 4) {
    #l17 white and pure
    colr.o<- colorRampPalette(c("#FFFFFF", "#FFFEFC", "#FEFDF9", "#FEFDF7", "#FDFCF4", "#FDFBF1", "#FCFAEE", "#FCFAEB", "#FBF9E9", "#FAF7E3", "#FAF7E0", "#F9F6DD", "#F8F5DB", "#F8F4D8", "#F7F3D5", "#F7F3D2", "#F6F2CF", "#F6F1CD", "#F4F0C7", "#F4EFC4", "#F3EEC2", "#F2EDBF", "#F2ECBD", "#F2EBBB", "#F2EABA", "#F2E9B8", "#F3E8B6", "#F3E6B3", "#F3E5B1", "#F3E4AF", "#F3E3AE", "#F3E2AC", "#F3E1AA", "#F3E0A9", "#F2DFA7", "#F2DEA5", "#F2DCA2", "#F2DBA0", "#F2DA9F", "#F2D99D", "#F2D79C", "#F2D69A", "#F2D598", "#F2D497", "#F2D396", "#F2D193", "#F3D092", "#F3CF91", "#F3CD90", "#F3CC8E", "#F3CB8D", "#F3CA8C", "#F3C98B", "#F3C688", "#F4C587", "#F4C486", "#F4C385", "#F4C284", "#F4C183", "#F4BF81", "#F4BE80", "#F4BD7F", "#F4BB7D", "#F4BA7C", "#F4B87B", "#F4B77A", "#F4B67A", "#F5B579", "#F5B378", "#F5B277", "#F5B177", "#F5AE75", "#F5AD75", "#F5AC74", "#F5AB73", "#F5AA73", "#F5A872", "#F5A771", "#F5A670", "#F5A470", "#F5A26E", "#F5A16E", "#F59F6D", "#F59E6C", "#F59D6C", "#F59C6C", "#F59A6B", "#F5996B", "#F5986B", "#F5956A", "#F5946A", "#F5936A", "#F5916A", "#F5906A", "#F48F69", "#F48D69", "#F48C69", "#F48A69", "#F48868", "#F48768", "#F48668", "#F48468", "#F38367", "#F38167", "#F38067", "#F37F67", "#F27C67", "#F27B67", "#F27A68", "#F27868", "#F17768", "#F17668", "#F17468", "#F07368", "#F07269", "#EF6F69", "#EF6E69", "#EF6C69", "#EE6B6A", "#EE6A6A", "#EE686A", "#ED676A", "#ED666A", "#ED646A", "#EC616B", "#EC606B", "#EB5F6B", "#EA5D6C", "#EA5C6C", "#E95B6D", "#E85A6D", "#E8586E", "#E7576E", "#E6556F", "#E5536F", "#E55270", "#E45170", "#E34F71", "#E34E71", "#E24D72", "#E14B72", "#E04973", "#E04773", "#DF4674", "#DE4474", "#DE4375", "#DD4175", "#DC4075", "#DB3F76", "#DA3E77", "#D83C78", "#D73B79", "#D63A79", "#D5397A", "#D4377A", "#D3367B", "#D2357C", "#D1347C", "#D0337D", "#CE307E", "#CD2F7F", "#CC2E7F", "#CB2D80", "#CA2B80", "#C92A81", "#C82982", "#C72882", "#C52683", "#C32584", "#C12485", "#C02385", "#BE2386", "#BD2287", "#BC2287", "#BA2188", "#B92089", "#B72089", "#B41F8B", "#B31E8B", "#B11D8C", "#B01D8C", "#AE1C8D", "#AD1C8E", "#AB1B8E", "#AA1B8F", "#A71990", "#A51991", "#A31892", "#A21892", "#A01893", "#9E1993", "#9C1994", "#9A1A94", "#981A95", "#941B96", "#921B97", "#901B97", "#8E1C98", "#8C1C99", "#891C99", "#871D9A", "#851D9A", "#831D9B", "#7F1E9C", "#7C1E9D", "#7A1E9D", "#781E9E", "#751F9E", "#731F9F", "#701F9F", "#6E20A0", "#6B21A0", "#6522A1", "#6223A1", "#5F23A2", "#5C24A2", "#5924A2", "#5525A3", "#5225A3", "#4E26A3", "#4B26A4", "#4327A4", "#3E27A5", "#3A28A5", "#3428A6", "#2F29A6", "#2929A6", "#2129A7", "#182AA7", "#002AA8"))(256)
  }
  if (col == 5) {
    # r3 - pure
    colr.o<-colorRampPalette(c("#085CF8", "#0F5FF4", "#1361F1", "#1763ED", "#1965E9", "#1B67E5", "#1C6AE1", "#1D6CDE", "#1D6EDA", "#1D72D2", "#1D74CE", "#1C75CB", "#1B77C7", "#1979C3", "#187BBF", "#167DBB", "#147EB8", "#1380B4", "#1283AC", "#1385A8", "#1486A4", "#1788A0", "#1A899C", "#1D8A98", "#208C93", "#238D8F", "#278E8B", "#2D9082", "#2F927E", "#329379", "#349475", "#359570", "#37966C", "#389767", "#399862", "#3A9A5E", "#3B9C54", "#3B9D50", "#3C9E4B", "#3C9F46", "#3CA042", "#3DA13D", "#3DA239", "#3EA335", "#3FA431", "#42A62B", "#44A728", "#47A826", "#49A824", "#4CA922", "#4EAA21", "#51AA20", "#54AB20", "#5AAC1F", "#5CAD1E", "#5FAE1E", "#62AE1E", "#65AF1E", "#67AF1D", "#6AB01D", "#6CB11D", "#6FB11D", "#74B21D", "#77B31C", "#79B41C", "#7BB41C", "#7EB51C", "#80B51B", "#83B61B", "#85B61B", "#87B71B", "#8CB81A", "#8EB91A", "#91B91A", "#93BA1A", "#95BA19", "#97BB19", "#9ABB19", "#9CBC19", "#9EBD18", "#A3BE18", "#A5BE18", "#A7BF17", "#A9BF17", "#ABC017", "#AEC016", "#B0C116", "#B2C116", "#B4C215", "#B8C315", "#BBC314", "#BDC414", "#BFC514", "#C1C513", "#C3C613", "#C5C613", "#C7C712", "#CCC812", "#CEC811", "#D0C911", "#D2C910", "#D4CA10", "#D6CA0F", "#D8CB0F", "#DBCB0E", "#DDCB0E", "#E1CC0E", "#E3CD0E", "#E5CD0E", "#E7CD0F", "#E9CE10", "#EBCE12", "#EDCE14", "#EFCD16", "#F1CD19", "#F4CC1F", "#F6CB23", "#F7CA26", "#F8C929", "#F9C82D", "#F9C730", "#FAC633", "#FAC536", "#FBC339", "#FCC03F", "#FCBF41", "#FCBE44", "#FCBC46", "#FCBB48", "#FDB94B", "#FDB84D", "#FDB64F", "#FDB551", "#FDB255", "#FEB057", "#FEAF59", "#FEAD5B", "#FEAC5D", "#FEAA5F", "#FEA961", "#FEA763", "#FEA466", "#FEA368", "#FEA16A", "#FEA06C", "#FF9E6D", "#FF9D6F", "#FF9B71", "#FF9A72", "#FF9874", "#FF9577", "#FF9379", "#FF927A", "#FF907C", "#FF8F7D", "#FF8D7F", "#FE8B80", "#FE8A82", "#FE8884", "#FE8587", "#FE8388", "#FE818A", "#FE808B", "#FE7E8C", "#FE7C8E", "#FE7A8F", "#FD7991", "#FD7792", "#FD7395", "#FD7197", "#FD7098", "#FD6E99", "#FC6C9B", "#FC6A9C", "#FC689E", "#FC669F", "#FC64A0", "#FB60A3", "#FB5EA4", "#FB5CA5", "#FA5AA6", "#FA58A7", "#FA56A8", "#FA54A8", "#F952A8", "#F94EA7", "#F84CA6", "#F84AA5", "#F849A3", "#F747A1", "#F7469E", "#F7449C", "#F64399", "#F64296", "#F53F8F", "#F53E8C", "#F43D88", "#F43C85", "#F33A82", "#F3397E", "#F2387B", "#F13778", "#F13674", "#F0336D", "#EF326A", "#EE3167", "#EE2F63", "#ED2E60", "#EC2D5D", "#EB2C59", "#EB2A56", "#EA2953", "#E8264C", "#E82549", "#E72445", "#E62242", "#E5213F", "#E4203B", "#E31E38", "#E21D35", "#E21B31", "#E0182A", "#DF1626", "#DE1423", "#DD131F", "#DC111B", "#DB0F17", "#DA0C12", "#D90A0C", "#D70500"))(256)
  }
  #white to black
  if (col == 6) {
    colr.o <- colorRampPalette(c("#FFFFFF", "#404040", "#000000"))(256)
  }
  # successively plot each model
  for (i in 1:length(model)) {
    image(model[[i]], col=colr.o, asp=1, axes=FALSE, xaxs="i", xaxt='n', yaxt='n', ann=F)

  # if the user specifies they want titles, do this:
    if (title == TRUE){
      title(main=names(model)[i])
    }
  }
}
#'Clip models via inverse-distance weighting
#'
#'Niche modeling may recover areas that are distant and non-contiguous to a species' core
#'range as suitable habitat. This can be either a feature or a bug. In certain cases one
#'may wish to exclude these "extraneous" suitable areas from an analysis because they could
#'not realistically be migrated to and are thus superfluous. This function uses inverse-distance
#'weighting to clip extraneous suitable areas a certain distance from a core range, as
#'delineated by occurrence data. Because occurrence data for extinct lineages does not exist
#'(barring fossils, which do not exist for many taxa), extant taxon occurrence data is used by
#'proxy. We suggest using machu.treeplot() to discern the extant descendants of the lineage
#'in question, and then restricting the occurrence data to those taxa using the taxa argument
#'of this function. Currently, the function only takes one niche model as input at a time.
#'
#'@param in.sdm niche model, a single raster object. Corresponds to one element
#'of the list output of machu.3.anc.niche().
#'@param in.pts occurrence points, a dataframe. Default input format is taxon, long, lat,
#'though this can be finagled with other arguments.
#'@param buffer.dist distance (in km) by which to buffer a minimum convex polygon (MCP)
#'around the occurrence points. Any model values inside this MCP will remain unchanged.
#'@param kernel.size a multiplier that defines the broadness of the smoothing edge
#'surrounding the MCP created by buffer.dist. A kernel.size of 2 creates a smoothing
#'edge twice the width of the MCP buffer.dist. Thus, higher values create broader models.
#'Minimum value of 1.
#'@param MCP.percent defines proportion of outlier occurrence points to be excluded from
#'the MCP creation.
#'@param taxa.col the column in the occurrence data with taxon names. Set to 1 by default.
#'@param long.col the column in the occurrence data with longitudes (x). Set to 2 by default.
#'@param lat.col the column in the occurrence data with latitudes (y). Set to 3 by default.
#'@param taxa a string or string vector declaring which taxa are to be used in the creation of the
#'MCP. Must match exactly. If the argument is not called, the function uses all occurrence
#'points by default.
#'
#'@return a clipped niche model, formatted as a raster layer.
#'
#'@examples
#'# after running through the "quick-start" Machuruku tutorial:
#'clip1<-machu.geo.idw(mod[[1]], occ, taxa=c("pepperi","bassleri","yoshina"), buffer.dist = 100, kernel.size = 2, MCP.percent = 50)
#'plot(mod[[1]]);plot(clip1)
#'# repeat for the silverstonei ancestor:
#'clip2<-machu.geo.idw(mod[[2]], occ, taxa="silverstonei", buffer.dist = 100, kernel.size = 2, MCP.percent = 50)
#'plot(mod[[2]]);plot(clip2)
#'@import rgdal
#'@import raster
#'@import adehabitatHR
#'@export
machu.geo.idw <- function(in.sdm, in.pts, buffer.dist=300, kernel.size=2, MCP.percent=100,
                          taxa.col=1, long.col=2, lat.col=3, taxa=NULL) {

  # restrict occ data to a particular species (or several)
  if (is.null(taxa)==FALSE){
    in.pts <- in.pts[in.pts[,taxa.col] %in% taxa,]
  }

  ## required estimate km from dd
  if (kernel.size<1){ kernel.size=1 }
  #print(kernel.size)
  ext<-extent(in.sdm)
  latU<-ext[3]
  latL<-ext[4]

  # estimage dd to km for study area
  Env.dd.y <- mean(c(latU,latL))

  ### Jason's dirty dd to km for lat
  AvgKm <- abs(((-0.0139*(Env.dd.y^2))+(0.0898*Env.dd.y)+111.1))
  buffer.sp <- buffer.dist/AvgKm #values in km

  ## convert to shapefile
  OcSp <- SpatialPoints(in.pts[,c(long.col,lat.col)], CRS("+proj=longlat +datum=WGS84"))

  ## create MCP
  OcSp2 <- mcp(OcSp, percent = MCP.percent) # note this could be adjusted

  ## project shapefile
  OcSpT <- spTransform(OcSp2, CRS("+proj=longlat +datum=WGS84"))

  ## buffer MCP
  OcSpT2 <- raster::buffer(OcSpT, width = buffer.sp, dissolve = T)

  #### code for inverse distance weighting of model
  ## buffer 2 MCP - used to speed up distance function
  OcSpT2B <- raster::buffer(OcSpT2, width = (buffer.sp*2.1*kernel.size), dissolve = T)

  ## convert shapefile polygons of buffer to raster for use as a masks
  rr <- rasterize(OcSpT2, in.sdm,update=FALSE)
  tt <- rasterize(OcSpT2B,in.sdm,update=FALSE, background=0)

  #subtrack mask from eachother, replace nulls, add 1
  cc<-tt-rr
  cc[cc==1]<-NA
  cc=cc+1

  ## create distance to buffer function
  dist<-distance(cc)

  ## buffer 3 - distance to scale values outside of first buffer
  OcSpT3 <- raster::buffer(OcSpT2, width = (buffer.sp*kernel.size), dissolve = T)

  #convert final buffer to raster
  raster_mask1 <- rasterize(OcSpT3,dist,update=FALSE)

  #convent mask to distance calculation
  raster_mask2 <- raster_mask1 * dist
  raster_mask2_rescale <- 1- (raster_mask2/ (maxValue(raster_mask2)))

  #final mask
  ##this would be the step where you would loop through a bunch of SDMs per species - each 'in.sdm' could be a SDM from another time period
  final_model <-in.sdm* raster_mask2_rescale

  #output
  return(final_model)
}

