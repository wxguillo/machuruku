#                Machuruku v2.0
# updated machu.top.env() to allow for pseudoabsence sampling of
# low resolution rasters when the n.occ*30 or even n.occ exceeds
# the number of pixels. Also streamlined the code some
#
# added a plot argument to machu.occ.rarefy() that visualizes
# rarefication of points
#
# rewrote machu.treeplot(), now primarily based in phytools
# with more streamlined arguments and better optimized visuals
#
# added the ability to plot present-day niche models to
# machu.1.tip.resp() and streamlined the code. Changed the
# skew-normal dist estimation function from fGarch::snormFit
# to sn::selm, as the "xi" skew param offered by snormFit
# doesn't behave like a skew param should (values <1 basically
# create an upside-down mirror of the dist, rather than make the
# skew go the other way). selm produces very similar results but
# with a standard skew param for which a BM evolution model is
# more justifiable. This function can also now output present-day
# niche models with the output.bioclim argument. When there are
# fewer than 10 points for a given species, the function will now
# randomly generate up to 10 points within an MCP for that species
#
# machu.2.ace() output format is streamlined and now only returns
# uCI and lCI values when unc=T rather than considering multiple
# uncertainty "bins". This didn't make much sense since low means
# would be paired with low stdevs for no real reason, etc.
#
# rewrote machu.ace.load() for compatibility with new machu.2.ace()
# output file format.
#
# big changes to machu.3.anc.niche():
# 1. works with terra
# 2. clip.Q=T and resp.curv=F now work without needing reconstructed lCI and uCI values
# 3. works with multiple timeslices
# 4. uses sn::dsn to calculate suitability (when resp.curv=F)
# 5. organically detects whether uncertainty is included in the calculations
# 6. can output raster files
#
# subsumed machu.trees.unc() into machu.tree.unc(), which now detects whether it should
# create uCI/lCI trees from a file containing multiple trees or a single tree with
# confidence intervals automatically. Can also handle output from MEGA/RelTime as well.
#
# rewrote machu.respplot and machu.plotmap for greater flexibility and compatibility
#
# removed machu.geo.idw because conceptually using present-day distributions to restrict
# past niche models doesn't make a lot of sense, and I wasn't happy with the results.

#' Calculate tip climate response curves
#'
#' For each taxon, construct a Bioclim species distribution model and estimate a response curve for each climate variable.
#'
#' @param occ Occurrence data for each species, formatted as a dataframe. Each species must have at least 10 points.
#' @param clim A single RasterStack of present-day climate data. A SpatRaster object is also acceptable.
#' @param sp.col Integer specifying which column of the input occurrence data corresponds to species ID. Default = 1.
#' @param col.xy vector specifying long (x) and lat (y) of occurrence data. Default = 2:3.
#' @param plot Whether and how to make plots of present-day Bioclim models for each taxon. "separate": plot each model in its own window. "together": plot all models in the same window, arranged dynamically at a given aspect ratio. "n": default (no plot).
#' @param plot.points If TRUE, plot each taxon's occurrence data on top of its corresponding niche model. Default = F.
#' @param plot.asp Aspect ratio used to calculate the best arrangement of plots when plot="together". Default = 16/9.
#' @param output.bioclim If TRUE, output a Bioclim niche model raster for each taxon instead of climate response curves.
#' @param verbose If TRUE, print progress updates to the screen.
#'
#' @details
#' This function uses the function dismo::bioclim() to construct present-day Bioclim niche models. It is only compatible with the older Raster package, so a SpatRaster object (from the newer Terra package) will automatically be converted before being passed to the rest of the function.
#'
#' Unfortunately the constraints of the 'sn::selm()' function disallow any taxa having fewer than 10 occurrence points. To that end, this function contains a utility to randomly sample occurrence points within the minimum convex polygon comprised of the occurrence data for species, up to n=10. When the plotting functionality is activated (i.e. plot="s" or plot="t"), these random points are drawn in red. In this case, the output of the function will be a list that contains the normal output (response table or niche models) as well as the occurrence data table with the newly added random points. Obviously, it is better to have at least 10 real occurrence points; however for rare species, or range-limited species after spatial rarefication, that may be difficult or impossible.
#'
#' @return Table consisting of the response of each species to the climate data. Each response is represented by a skew-normal distribution. Alternatively, the function can output the actual Bioclim niche models.
#'
#' @examples
#' ## acceptable 'clim' formats
#' # Single RasterStack (raster) (preferred)
#' clim <- stack(list.files(rasterfolder, pattern="T0_", full.names=T))
#' # Single SpatRaster (terra)
#' clim <- c(rast(list.files(rasterfolder, pattern="T0_", full.names=T)))
#'
#' # basic
#' response.table <- machu.1.tip.resp(occ, clim)
#' # plot all plots separately
#' response.table <- machu.1.tip.resp(occ, clim, plot="s")
#' # plot all plots together, with points, at an aspect ratio of 16:10 (you have a 1600p screen)
#' response.table <- machu.1.tip.resp(occ, clim, plot="t", plot.points=T, plot.asp=16/10)
#' @import dismo
#' @import sn
#' @import raster
#' @import sf
#' @import dplyr
#' @export
machu.1.tip.resp <- function(occ, clim, sp.col=1, col.xy=2:3, plot="n", plot.points=F, plot.asp=16/9, output.bioclim=F, verbose=F){

  # get species' names
  sp <- unique(occ[,sp.col])

  # check whether each species has at least 10 occurrence points
  occ.counts <- sapply(sp, function(x) nrow(subset(occ, occ[sp.col] == x)))
  # if so, add points up to n=10 within the bounding box comprised of the species' points, and print warning messages since this isn't exactly best practice
  if (any(occ.counts < 10)){
    add.pts <- names(occ.counts[occ.counts < 10])
    print(paste("The following taxa have fewer than 10 occurrence points:", paste(add.pts, collapse=", ")))
    print(paste("Warning: adding random points up to n=10 for each of these species, contained within MCP defined by points."))
    # create a list where each element is a table for a single species (only species with <10 points)
    t <- lapply(add.pts, function(sp) occ[occ[,sp.col]==sp,]); names(t) <- add.pts
    # randomly select 10 minus the number of points for each taxon within a minimum convex polygon
    t <- lapply(t, function(x) x %>% st_as_sf(coords=col.xy) %>% summarize(geometry=st_union(geometry)) %>% st_convex_hull() %>% st_sample(10-nrow(x)))
    # convert to dataframes
    t <- lapply(t, function(x) unlist(x) %>% matrix(ncol=2, byrow=T) %>% as.data.frame())
    # add species name column and new "rand" column specifying which points are randomly generated
    t <- lapply(names(t), function(x) cbind(rep(x, nrow(t[[x]])), t[[x]], rep(T, nrow(t[[x]]))))
    # combine into new table
    t <- setNames(do.call(rbind, t), c(names(occ), "rand"))
    # combine with original table
    t <- rbind(setNames(cbind(occ, rep(F, nrow(occ))), c(names(occ), "rand")), t)
    # sort table
    occ <- t[order(t[,sp.col]),]
  }

  if (class(clim)=="SpatRaster"){
    if (verbose==T) print("'clim' is a SpatRaster, attempting to convert to RasterStack for compatibility with dismo::bioclim...")
    clim <- stack(clim)
  }
  if (class(clim)=="list" | class(clim)=="RasterBrick") stop("Provide a single RasterStack or SpatRaster for 'clim'.")

  # sample each species' occurrence data for each climate variable
  enms <- lapply(as.list(sp), function(x){
    sp.pts <- subset(occ, occ[,sp.col] == x)[,col.xy]
    bioclim(clim, sp.pts)
  })
  #########
  # PLOTS #
  #########
  # separate plots:
  if (plot=="separate" | plot=="separately" | plot=="sep" | plot=="s"){
    if (verbose==T) print(paste("Plotting each plot separately. Scroll back and forth to see each plot."))
    for (i in 1:length(enms)) {
      plot(dismo::predict(clim, enms[[i]]), main=sp[i])
      if (plot.points==T){
        p <- occ[occ[,sp.col]==sp[i],]
        points(p[,col.xy])
        if (any(colnames(occ)=="rand")) points(p[p$rand==T, col.xy], col="red")
      }
    }
  }
  if (plot=="together" | plot=="tog" | plot=="t"){
    arr <- n2mfrow(length(enms), asp=plot.asp) # calculate best plot arrangement
    par(mfrow=arr, mar=c(3.1, 2.1, 1.1, 1.1))
    if (verbose==T) print(paste("Plotting all plots together. n2mfrow chose", arr[1], "rows and", arr[2], "columns based on an aspect ratio of", MASS::fractions(plot.asp)))
    for (i in 1:length(enms)){
      plot(dismo::predict(clim, enms[[i]]), main=sp[i])
      if (plot.points==T){
        p <- occ[occ[,sp.col]==sp[i],]
        points(p[,col.xy])
        if (any(colnames(occ)=="rand")) points(p[p$rand==T, col.xy], col="red")
      }
    }
  }
  if (output.bioclim==T){
    if (verbose==T) print("Outputting Bioclim models instead of response curves.")
    ob <- lapply(enms, function(x) predict(clim, x))
    names(ob) <- sp
    # if 'rand' was added to 'occ', output 'occ' as well
    if (any(colnames(occ)=="rand")) return(list(occ, ob)) else return(ob)
  }
  #########
  # end plot section
  #########

  # initialize some objects that'll come in handy later
  response.list <- list()

  # prep enms list for next part
  enms <- lapply(enms, function(x) x@presence)
  names(enms) <- sp

  # for each species:
  for (j in sp){

    in.pres <- enms[[j]]

    # initialize vector that will later be each row of the table
    row.vector <- c()

    # for each climate variable for each species:
    for (i in colnames(in.pres)){

      # mean, sd, xi
      normFit <- extractSECdistr(selm(eval(parse(text=i)) ~ 1, data=in.pres))
      mean <- normFit@dp["xi"]
      sd   <- normFit@dp["omega"]
      skew <- normFit@dp["alpha"]
      params <- c(mean, sd, skew)

      # name the elements of the vector according to the name of the climate variable and the response variable
      names(params) <- paste0(i, "_", c("mean", "stdev", "skew"))

      # build a giant vector for the whole species, that will eventually be turned into a row for the big table later
      row.vector <- append(row.vector, params)
    }
    # add each species' row.vector to a list to be stored for later
    response.list[[j]] <- row.vector
    if (verbose == TRUE) print(paste0("Processed taxon ", j))
  }
  # combine list elements into a single output table
  response.table <- do.call(rbind, response.list)
  # if 'rand' was added to 'occ', output 'occ' as well
  if (any(names(occ)=="rand")) return(list(occ, response.table)) else return(response.table)
}

#' Select top environmental variables
#'
#' Perform a boosted regression tree analysis to identify the most important climate variables for your taxon set.
#'
#' @param occ Occurrence data for all taxa. Identical to input for machu.1.tip.resp(). Dataframe, with columns in the order of species, x/long, y/lat.
#' @param clim Climate data for all taxa. Identical to input for machu.1.tip.resp(). A RasterStack of corresponding climate variables. SpatRaster (from Terra) is also acceptable
#' @param sp.col Specify which column of the input occurrence data corresponds to species ID. Default = 1.
#' @param col.xy vector specifying long (x) and lat (y) of occurrence data. Default = 2:3.
#' @param learning.rt Value from 0.001 to 0.01 for building the ENMs, start with 0.01 and if prompted, change to 0.001. Default = 0.01.
#' @param steps Numbers of trees to add at each cycle for modelling each taxon. Start with 50 and if you run into problems gradually decrease, stopping at 1. Default = 50.
#' @param method This determines how important environmental variables are selected.There are three options: "estimate", "contrib", "nvars". If method="estimate", the boosted regression tree algorithm will choose the number of variables to include by systematically removing variables until average change in the model exceeds the original standard error of deviance explained. This is the most computationally intensive method. If method="contrib", variables above a relative influence value will be kept. See associated parameter 'contrib.greater'. If method="nvars", a fixed number of user specified variables will be kept. See associated parameter 'nvars.save'. The kept variables are selected by their relative influence. The 'nvars.save'-highest contributing variables for each taxon are retained and pooled, then ranked, and the 'nvars.save'-highest contributing variables for the whole pool are finally retained.
#' @param nvars.save If method="nvars",this variable is required. It is the number of the top variables to save. The kept variables are selected by their relative influence in predicting the species distribution, selecting for the highest contributing variables. Often the total variables retained is lower due to identical variables select among both species. The default value is 5.  This value will be ignored if method="estimate" or "contrib".
#' @param contrib.greater If method="contrib", this variable is required. The kept variables are selected for their relative influence in predicting the species' distribution. Here, users select variables equal to or above an input model contribution value. The default value for this method is 5 (= variables with 5 percent or higher contribution to model of either species are kept). This value will be ignored if method="estimate" or "nvars".
#' @param pa.ratio Ratio of pseudoabsences to occurrence points, typically this is 4. The default value is 4. There have to be at least 50 total points (occ+pseudoabsences) for the model to work; if the sum does not total to 50, the difference is taken as the number of pseudoabsences, rather than the value of occ*pa.ratio.
#' @param verbose Tf TRUE, print progress to the screen. Default = F.
#'
#' @details This function is a modified version of humboldt.top.env() from the package Humboldt. It runs generalized boosted regression models (a machine learning ENM algorithm) to select top parameters for inclusion your analyses. This is important because you want the models to reflect variables that are relevant to the species' distribution. Alternatively, you can run Maxent outside of R and manually curate the variables you include (also recommended).
#'
#' @return Prints the important climate variables to the screen. You can then combine them into a new RasterStack or SpatRaster object.
#'
#' @examples
#' ## acceptable 'clim' formats
#' # Single RasterStack (raster) (preferred)
#' clim <- stack(list.files(rasterfolder, pattern="T0_", full.names=T))
#' # Single SpatRaster (terra)
#' clim <- c(rast(list.files(rasterfolder, pattern="T0_", full.names=T)))
#'
#' # identify the top 6 climate variables across all taxa
#' machu.top.env(occ, clim, method = "nvars", nvars.save = 6)
#' # identify all climate variables with a contribution greater than 10%
#' machu.top.env(occ, clim, method = "contrib", contrib.greater = 10)
#' @import dismo
#' @import sp
#' @import raster
#' @import gbm
#' @import dplyr
#' @export
machu.top.env <- function(occ, clim, sp.col=1, col.xy=2:3, learning.rt=0.01, steps=50, method="contrib", nvars.save=5, contrib.greater=5, pa.ratio=4, verbose=F) # method=='estimate' method=='contrib' method=='nvars'
{
  # clim format check
  if (class(clim)=="SpatRaster"){
    if (verbose==T) print("'clim' is a SpatRaster, attempting to convert to RasterStack...")
    clim <- stack(clim)
  }
  if (class(clim)=="list" | class(clim)=="RasterBrick") stop("Provide a single RasterStack or SpatRaster for 'clim'.")

  # input adjustments
  if (method=="ESTIMATE" | method=="Estimate") method <- "estimate"
  if (method=="CONTRIB" | method=="Contrib" | method=="contribution") method <- "contrib"
  if (method=="Nvars" | method=="NVARS") method <- "nvars"
  # this is necessary bc the gbm.step function works by setting "silent" rather than "verbose", which are opposites
  if (verbose == T) verbose <- F else if (verbose == F) verbose <- T

  # determine indices (i.e., which columns, corresponding to no. climate variables) of predictor vars
  e.var <- 4:(3+nlayers(clim))

  # identify max number of occurrence points across all species
  n.occ.max <- max(table(occ$sp))
  # randomly sample n.occ.max*30 points from the raster climate data: these will be the pool from which pseudoabsences will be drawn
  # however if the data is low resolution there will likely be too few cells to sample w/o replacement
  # check if this is the case, and if so, change n.occ.max to n.cells
  n.cells <- nrow(rasterToPoints(clim[[1]], fun=NULL, spatial=FALSE)) # doesn't count ocean cells
  if (n.occ.max*30 > n.cells) n.occ.max <- n.cells
  if (verbose == F) { print("Extracting climate values from rasters. This could take a minute.")}
  pa.pool <- sampleRandom(clim, n.occ.max*30, xy=T)
  if (verbose == F) { print("Finished extracting climate values.")}
  ID <- rep(0, nrow(pa.pool))
  pa.pool <- cbind(ID, pa.pool)

  # initialize output list
  imp.vars <- list()
  # loop over each taxon
  for (i in 1:length(unique(occ$sp))){

    # get species name
    sp.name <- unique(occ$sp)[i]

    # get number of points for that taxon
    n.occ <- nrow(occ[occ[,sp.col]==sp.name,])

    # randomly sample n.occ*4 pseudoabsence (pa) points from pa.pool (4 is changeable with pa.ratio)
    # there have to be at least 50 points for the gbm model to work
    # if n.occ*pa.ratio > n.cells, set n.pa = n.cells
    # if n.occ > n.cells, set n.pa
    # perform various other checks and corrections
    n.pa <- n.occ*pa.ratio
    if (n.occ+(n.occ*pa.ratio) < 50) n.pa <- 50-n.occ
    if (n.occ*pa.ratio > n.cells) n.pa <- n.cells

    # sample pseudoabsences
    sp.pa <- pa.pool[sample(nrow(pa.pool), n.pa),]

    # get occurrence points for that species
    sp.occ <- occ[occ[,sp.col]==sp.name, col.xy]
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
      #names(env1)[vars1.use]
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
  return(imp.var.out)
}

#' Rarefy occurrence points
#'
#' Rarefy spatial occurrence points to reduce spatial autocorrelation and model bias. This is a modified version of humboldt.occ.rarefy() from humboldt.
#'
#' @param in.pts Input dataframe.
#' @param colxy Columns corresponding to longitude, then latitude. Default = 2:3.
#' @param rarefy.dist Distance to rarefy points (values need to be in km (recommended) or decimal degrees).  See associated parameter rarefy.units. Default = 0.
#' @param rarefy.units The units of the rarefy.dist parameter, either "km" for kilometers or "dd" for decimal degrees. Default = "km".
#' @param plot Creates an optional plot visualizing the points removed and kept. Default = F.
#' @param verbose If verbose=T, text boxes displaying progress will be displayed. Default = T.
#'
#' @details A script to systematically select localities within a specified area at specified spatial resolution.  The outcome is always the same and is not random.  This reduces sampling biases in downstream analyses- you should do it! Output is a reduced dataset with less spatial autocorrelation.
#'
#' @return A dataframe with rarefied occurrence data.
#'
#' @importFrom spatstat.geom nndist
#' @importFrom tcltk setTkProgressBar
#' @examples
#' ##remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.
#' occ <- machu.occ.rarefy(in.pts = occ, colxy = 2:3, rarefy.dist = 50, rarefy.units = "km")
#' @export
machu.occ.rarefy <- function(in.pts, colxy = 2:3, rarefy.dist = 0, rarefy.units = "km", plot=F, verbose = T) {
  switch(Sys.info()[['sysname']],
         Windows= {userOS=1},
         Linux  = {userOS=2},
         Darwin = {userOS=2})
  # input adjustments
  if (rarefy.units == "KM" | rarefy.units == "Km") rarefy.units <- "km"
  if (rarefy.units == "DD" | rarefy.units == "Dd") rarefy.units <- "dd"

  if (rarefy.units == "km") {
    min.dist <- rarefy.dist * 1000  #values in km
    sp2p1 <- SpatialPoints(in.pts[,colxy], CRS("+proj=longlat +datum=WGS84"))
    sp2p2 <- spTransform(sp2p1, CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    xy <- data.frame(x4r = coordinates(sp2p2)[,1], y4r = coordinates(sp2p2)[,2])
    xy <- data.frame(cbind(xy, in.pts))
  }

  if (rarefy.units == "dd") {
    yy <- colxy[2]
    maxLat <- max(in.pts[, yy])
    minLat <- min(in.pts[, yy])
    # estimate dd to km for study area
    rare.dd <- (mean(c(maxLat, minLat)))
    adjKm <- (-0.0139 * (rare.dd * rare.dd)) + (0.0898 * rare.dd) + 111.1
    min.dist <- adjKm * rarefy.dist * 1000  #values in km
    print(paste("Value used for rarefying:", round((min.dist/1000), 2), "km. Remember that the length of a decimal degree changes latitudinally due to the convergence of the lines of longitude at the poles. The value used here is the average distance of decimal-degrees within your study area. Alternatively, simply input distance as km value and change rarefy.units='km'"))
    sp2p1 <- SpatialPoints(in.pts[, colxy], CRS("+proj=longlat +datum=WGS84"))
    sp2p2 <- spTransform(sp2p1, CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    xy <- data.frame(x4r = coordinates(sp2p2)[, 1], y4r = coordinates(sp2p2)[, 2])
  }

  nPts <- nrow(xy)
  if (verbose==T & userOS==1) pb <- winProgressBar(title = "Initializing", min = 0, max =nPts, width = 300)
  if (verbose==T & userOS==2) pb <- tkProgressBar(title = "Initializing", label = "", min = 0, max = nPts, initial = nPts, width = 300)

  # setup env - separate from distance
  spName <- in.pts[,1][1]
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
  if (verbose == T) close(pb)
  new.data <- rbind(new.data, del.min.dist)
  nc <- (ncol(new.data))
  col.orig <- c(3:nc)
  data.out <- new.data[,col.orig]
  print(paste0("Starting points = ", nrow(xy), ", Final rarefied points = ", nrow(new.data)))

  # optional plot
  if (plot==T){
    plot(xy[,colxy+2], main="Before vs. After Spatial Rarefication", pch=19)
    points(new.data[,colxy+2], pch=19, col="red")
    legend("bottomright", c("Removed", "Kept"), pch=19, col=c("black", "red"),
           inset=c(0,-0.275), xpd=T, horiz=T, bty="n")
  }

  return(data.out)
}

#'Account for uncertainty in divergence time estimation from a single tree
#'
#'Create a list of three trees that represent uncertainty in divergence times for downstream analysis
#'
#'@param tree Primary input. Can be a single treedata object (loaded via treeio::read.beast) or a string indicating a nexus-style treefile containing multiple trees (i.e. a Bayesian posterior distribution).
#'@param burnin Percentage of trees to skip at the start of the multiple-trees file since MCMC algorithms take a while to converge. Must be between 0 and 1. Pertains only when the input is a Bayesian posterior of multiple trees. Default = 0.1.
#'@param conf Confidence level from which to calculate HPD limits for the tree heights in a multiple-trees file. Must be between 0 and 1. Pertains only when the input is a Bayesian posterior of multiple trees. Default = 0.1.
#'@param inc Increment at which to report progress, i.e. when every 'inc' trees are processed. Pertains only when the input is a Bayesian posterior of multiple trees, and verbose = TRUE. Default = 10000.
#'@param verbose Report progress and checks to screen. Default = TRUE.
#'
#'@details
#'This function generates a list of three trees that characterize the uncertainty of divergence time estimation. These trees can then be used in machu.2.ace() to explore scenarios in which varying numbers of taxa may have been present at a certain timeslice, given divergence time uncertainty.
#'
#'The function has two modes, depending on what input is provided. When a single tree is provided, the function outputs the input tree, plus two additional trees, constructed from the upper and lower confidence limits provided to characterize divergence time uncertainty by time calibration software such as BEAST. The 'phylo' object (from ape::read.tree) usually used to represent phylogenetic trees in R does not store this information, so the input must be in the 'treedata' format (from treeio::read.beast). Currently the function is compatible with output from BEAST and BEAST-adjacent software such as SNAPP, as well as output from MEGA calibrated with the RelTime method.
#'
#'The other mode characterizes uncertainty directly from a Bayesian posterior of trees. The function identifies the trees in the posterior with total tree heights closest to the median and the highest posterior density (HPD) limits of the distribution of tree heights. Because a tree posterior can consist of millions of trees, rather than overwhelm R by loading them all into an object (i.e. 'treedataList' from treeio), the function reads the trees directly from a file one at a time. Thus this mode is activated by specifying a filename; the file should be output from BEAST or BEAST-adjacent software, generally a .trees file in NEXUS format (open the file in a text editor to check, and examine the example file provided by the tutorial). The method consists of three "passes" through the posterior distribution of trees. In the first pass, the number of trees is simply counted, and the number of samples to skip as burn-in, given the burnin percentage specified by the user, is calculated. In the second pass, total tree heights are calculated for each tree, and the median and HPD limits at the confidence level specified are calculated for this distribution. In the third pass, the trees with total heights closest to the median and HPD limits are found, and sent to the output list.
#'
#'@return A list of three 'phylo' objects characterizing uncertainty in divergence times.
#'
#'@examples
#'# load a tree and store divergence time uncertainty info
#'tree <- treeio::read.beast("tree.treefile")
#'# characterize uncertainty
#'trees <- machu.tree.unc(tree)
#'# plot trees to visualize uncertainty
#'machu.treeplot(trees)
#'
#'# characterize uncertainty from a posterior of trees and treat the first 20% as burn-in
#'trees <- machu.tree.unc("posterior.trees", burnin=0.2)
#'
#'@import ape
#'@import treeio
#'@import TeachingDemos
#'@export
machu.tree.unc <- function(tree, burnin=0.1, conf=0.95, inc=10000, verbose=T){

  # format checks
  if (class(tree)=="treedataList") stop("'treedataList' object detected. For the multiple trees method, supply the filename rather than the object itself. This saves processing power.")
  if (class(tree)=="phylo") stop("'phylo' object detected. This format doesn't store uncertainty data. Use treeio::read.beast() to supply a 'treedata' object instead.")

  # single tree method
  if (class(tree)=="treedata"){
    # check for ultrametricity
    if (is.ultrametric(as.phylo(tree))==F) stop("Input tree must be ultrametric.")
    # check if this is a reltime tree
    if (length(grep("reltime", colnames(tree@data))) > 0){
      rt <- T
      if (verbose==T) print("Single RelTime tree detected.")
    } else {
      rt <- F
      if (verbose==T) print("Single BEAST-type tree detected.")
    }
    # lower CI tree
    LCItib <- as_tibble(tree)
    for (n in 1:nrow(LCItib)){
      # get upper CI for that node
      if (rt==F) nodeLCI <- LCItib$height_0.95_HPD[[n]][1] else if (rt==T) nodeLCI <- LCItib$divtime_0.95_CI[[n]][1]
      # get that node's parent node
      pnode <- LCItib$parent[n]
      # get upper CI for parent node
      if (rt==F) pnodeLCI <- LCItib$height_0.95_HPD[[pnode]][1] else if (rt==T) pnodeLCI <- LCItib$divtime_0.95_CI[[pnode]][1]
      # if the active row is a tip, set LCI to zero
      if (n <= Ntip(tree)) nodeLCI <- 0
      if (pnode <= Ntip(tree)) pnodeLCI <- 0
      # if any of these valuse is NA or NULL, set to zero
      if (is.null(nodeLCI)) nodeLCI <- 0; if (is.na(nodeLCI)) nodeLCI <- 0
      if (is.null(pnodeLCI)) pnodeLCI <- 0; if (is.na(pnodeLCI)) pnodeLCI <- 0
      # if the upper CI for the focal node is higher than that of its parent node, change its upper CI value to upper CI of parent node
      if (nodeLCI > pnodeLCI & n != pnode) nodeLCI <- pnodeLCI
      # set tree node heights to upper CI values
      LCItib$branch.length[n] <- pnodeLCI - nodeLCI
    }
    LCItree <- as.phylo(LCItib)
    if (verbose==T) print("Created LCI tree based on lower end of the 95% confidence intervals for divergence times.")

    # upper CI tree
    UCItib <- as_tibble(tree)
    for (n in 1:nrow(UCItib)){
      # get upper CI for that node
      if (rt==F) nodeUCI <- UCItib$height_0.95_HPD[[n]][2] else if (rt==T) nodeUCI <- UCItib$divtime_0.95_CI[[n]][2]
      # get that node's parent node
      pnode <- UCItib$parent[n]
      # get upper CI for parent node
      if (rt==F) pnodeUCI <- UCItib$height_0.95_HPD[[pnode]][2] else if (rt==T) pnodeUCI <- UCItib$divtime_0.95_CI[[pnode]][2]
      # if the active row is a tip, set UCI to zero
      if (n <= Ntip(tree)) nodeUCI <- 0
      if (pnode <= Ntip(tree)) pnodeUCI <- 0
      # if any of these valuse is NA or NULL, set to zero
      if (is.null(nodeUCI)) nodeUCI <- 0; if (is.na(nodeUCI)) nodeUCI <- 0
      if (is.null(pnodeUCI)) pnodeUCI <- 0; if (is.na(pnodeUCI)) pnodeUCI <- 0
      # if the upper CI for the focal node is higher than that of its parent node, change its upper CI value to upper CI of parent node
      if (nodeUCI > pnodeUCI & n != pnode) nodeUCI <- pnodeUCI
      # set tree node heights to upper CI values
      UCItib$branch.length[n] <- pnodeUCI - nodeUCI
    }
    UCItree <- as.phylo(UCItib)
    if (verbose==T) print("Created UCI tree based on upper end of the 95% confidence intervals for divergence times.")

    # create output list
    out <- list(LCItree, as.phylo(tree), UCItree)
    names(out) <- c("lCItree", "inputtree", "uCItree")
    return(out)

  # multiple trees method
  } else if (class(tree)=="character"){

    # announcements and checks
    if (verbose==T) print(paste0("Multiple trees in file '", tree, "' detected. Burnin: ", burnin, "; Confidence level: ", conf))
    if (burnin > 1 | burnin < 0) stop("'burnin' must be between 0 and 1.")
    if (conf > 1 | conf < 0 ) stop("'conf' must be between 0 and 1.")

    ##### FIRST PASS: Calculating the number of trees and burnin amount #####
    if (verbose==T) print(paste0("FIRST PASS: Calculating trees to remove at ", burnin*100, "% burnin."))
    con <- file(tree, "r")
    # initiate tree counter and translation table
    treenumber <- 1
    #tlntab <- rep(NA, Ntip(read.nexus(tree))[[1]])
    tlntab <- c()
    while (T) {
      line <- readLines(con, n=1)
      if (length(line)==0) break
      # get tip translation table
      if (grepl("^\\s+\\d+\\s.+", line)==T) tlntab <- c(tlntab, line)
      # if the line is actually a tree, do this:
      if (grepl("^tree", line)==T) {
        if (treenumber %% inc == 0 & verbose==T) print(paste0("First pass (burnin): Tree ", treenumber))
        treenumber <- treenumber+1
      }
    }
    # create and format translation table for later
    tlntab <- gsub(pattern=",$", replacement="", x=tlntab) # remove trailing commas
    tlntab <- strsplit(tlntab,"\\s+")
    tlntab <- do.call(rbind, tlntab)[,2:3]
    # calculate burnin
    burn.no <- ceiling(burnin*(treenumber-1))
    # report burnin and treenumber
    if (verbose==T) print(paste0("Found ", treenumber-1, " trees, skipping the first ", burn.no, " under ", burnin*100, "% burnin."))
    close(con)

    ##### SECOND PASS: Calculating HPD distribution of tree heights #####
    if (verbose==T) print(paste0("SECOND PASS: Calculating ", conf*100, "% HPD of tree heights at ", burnin*100, "% burnin."))
    con <- file(tree, "r")
    # re-initiate tree counter and a vector of tree heights
    treenumber <- 1
    treeheights <- c()
    # read in file line by line
    while (T) {
      line = readLines(con, n = 1)
      if (length(line)==0) break
      # if the line is actually a tree, do this:
      if (grepl("^tree", line)==T) {
        # if treenumber is still below the burnin threshold, skip the tree
        if (treenumber <= burn.no){
          # shout to screen
          if (treenumber %% inc == 0 & verbose==T){ print(paste0("Second pass (HPD): Tree ", treenumber, " (skipped at ", burnin*100, "% burnin)")) }
          treenumber <- treenumber+1
          # if treenumber is above burnin threshold, calculate tree height and add to vector
        } else if (treenumber > burn.no){
          height <- sum(read.tree(text=line)$edge.length)
          treeheights <- c(treeheights, height)
          # shout to screen
          if (treenumber %% inc == 0 & verbose==T) print(paste0("Second pass (HPD): Tree ", treenumber))
          treenumber <- treenumber+1
        }
      }
    }
    # Calculate HPDs and average
    HPD.min <- emp.hpd(treeheights, conf=conf)[1]
    HPD.max <- emp.hpd(treeheights, conf=conf)[2]
    HPD.med <- median(treeheights)
    # report HPDs and average
    if (verbose==T) print(paste0(conf*100, "% HPD min: ", round(HPD.min, digits=2), "; Median: ", round(HPD.med, digits=2), "; ",
                                 conf*100, "% HPD max: ", round(HPD.max, digits=2), "; from ", length(treeheights), " trees under ",
                                 burnin*100, "% burnin."))
    close(con)

    ##### THIRD PASS: Find closest trees to HPDs and median #####
    if (verbose==T) print(paste0("THIRD PASS: Finding trees with closest heights to ", conf*100, "% HPD limits and median at ", burnin*100, "% burnin."))
    con <- file(tree, "r")
    # re-re-initiate tree counter
    treenumber <- 1
    # read in file line by line
    while (T) {
      line = readLines(con, n=1)
      if (length(line)==0) break
      # if the line is actually a tree, do this:
      if (grepl("^tree", line)==T) {
        # if treenumber is still below the burnin threshold, skip the tree
        if (treenumber <= burn.no){
          if (treenumber %% inc == 0 & verbose == TRUE){ print(paste0("Third pass (trees): Tree ", treenumber, " (skipped at ", burnin*100, "% burnin)")) }
          treenumber <- treenumber+1
        # if tree number = burn.no+1, set all values to this tree; this gives us a place to start from
        } else if (treenumber == burn.no+1){
          height <- sum(read.tree(text=line)$edge.length)
          mintree <- line
          medtree <- line
          maxtree <- line
          lowest.from.min <- abs(height-HPD.min)
          lowest.from.med <- abs(height-HPD.med)
          lowest.from.max <- abs(height-HPD.max)
          # shout to screen
          if (treenumber %% inc == 0 & verbose==T) print(paste0("Third pass (trees): Tree ", treenumber))
          treenumber <- treenumber+1
        # if tree number > burn.no+1, test whether treeheight distance is lower than the current lowest for all values, and if so, set new values and trees
        } else if (treenumber > burn.no+1) {
          height <- sum(read.tree(text=line)$edge.length)
          # HPD.min
          if (abs(height-HPD.min) < lowest.from.min){
            lowest.from.min <- abs(height-HPD.min)
            mintree <- line }
          # HPD.med
          if (abs(height-HPD.med) < lowest.from.med){
            lowest.from.med <- abs(height-HPD.med)
            medtree <- line }
          # HPD.max
          if (abs(height-HPD.max) < lowest.from.max){
            lowest.from.max <- abs(height-HPD.max)
            maxtree <- line }
          # shout to screen
          if (treenumber %% inc == 0 & verbose==T) print(paste0("Third pass (trees): Tree ", treenumber))
          treenumber <- treenumber+1
        }
      }
    }
    # Put trees into a list
    treelist <- list(read.tree(text=mintree), read.tree(text=medtree), read.tree(text=maxtree))
    names(treelist) <- c(paste0("lower", conf*100, "HPDtree"), "mediantree", paste0("upper", conf*100, "HPDtree"))
    # Need to restore the original tip labels, since nexus represents each with a number
    # Sort the translation table made in pass 1 by the order of the numbered tips in the trees
    treelist <- lapply(treelist, function(x){
      sorted.tip.labs <- tlntab[match(x$tip.label, tlntab[,1]),2]
      x$tip.label <- sorted.tip.labs
      return(x)
    })
    return(treelist)
    close(con)
  }
}

#'Plot trees to visualize each step in the Machuruku process
#'
#'Visualize trees to help with interpreting Machuruku analysis and results. Will plot a single tree or multiple trees (i.e. the results of machu.tree.unc) on the same axis. Node numbers/IDs and divergence times can be visualized. Timeslices can be plotted to visualize which "branch-taxa" may be returned at a given time-slice.
#'
#'@param tree Tree input. Can be a single tree in phylo format, or a list of trees. Untested with other formats.
#'@param timeslice Numeric vector with the timeslices to plot, if any.
#'@param nodelabs Whether to plot labels identifying each node. Default = T.
#'@param nodelabsize This is given to the circle.exp arg of phytools::labelnodes(). Default = 0.5, which is optimized for simple trees.
#'@param col Color of the timeslices. Default = "red".
#'@param x.u.lim Upper limit of the x-axis (i.e. right side of plot). By default this is the order of magnitude of the tree height, which is usually enough to accommodate decently sized taxon names. If you have longer names you can manually adjust this value.
#'@param x.l.lim Lower limit of the x-axis (i.e. left side of plot). By default this is scaled by the order of magnitude of the tree height. If you have your own scale in mind you can override the default.
#'@param timelabs Whether to plot divergence time labels at each node. Default = T.
#'@param timelabsize Size of div time labels. Default = 0.8, which is optimized for simple trees.
#'@param timelaboffset Offset of div time labels from their respective nodes. Default = -0.35, which is optimized for simple trees. Making this number more negative will offset the labels further to the right.
#'
#'@details
#'This function tends to generate warnings that shouldn't affect the output. You can turn them off with options(warn=-1).
#'
#'@return Plots a tree or tree with node labels, divergence times, and time-slices.
#'
#'@examples
#'# plot tree by itself
#'machu.treeplot(tree)
#'# plot with three timeslices
#'machu.treeplot(tree, c(1e6, 2e6, 3e6))
#'@import phytools
#'@export
machu.treeplot <- function(tree, timeslice=NULL, nodelabs=T, nodelabsize=0.5, col="red", x.u.lim=NULL, x.l.lim=NULL, timelabs=T, timelabsize=0.8, timelaboffset=-0.35){
  # save user's previous par settings
  pars <- par()

  # plot single tree if a list of length=1 is given, or a simple phylo object
  if (length(tree) == 1 || class(tree) == "phylo" ){

    if (class(tree)=="list") tree <- tree[[1]] # if tree is a list with length=1, unlist it

    par(mar=c(3.1, 2.1, 1.1, 1.1)) # set plot margins

    th <- max(nodeHeights(tree)) # calculate tree height
    om <- 10^floor(log10(th)) # calculate order of magnitude of tree height
    xlims <- c(ceiling(th/om)*om, -om) # dynamically calculate a good x axis
    if (is.null(x.l.lim)==F) xlims[1] <- x.l.lim # override lower x-limit
    if (is.null(x.u.lim)==F) xlims[2] <- x.u.lim # override upper x-limit
    ylims <- c(0.5, 1.2*Ntip(tree)) # dynamically calculate a good y axis

    plot(NA, xlim=xlims, ylim=ylims, bty="n", axes=F, xlab="", ylab="") # set up plot
    axis(1, at=seq(0, xlims[1], by=om)) # draw dynamic timescale
    abline(v=timeslice, col=col) # draw timeslices, if provided
    plotTree(tree, direction="leftwards", xlim=xlims, ylim=ylims, mar=par()$mar, color="black", add=T)
    if (nodelabs==T) labelnodes(text=1:tree$Nnode, node=1:tree$Nnode+Ntip(tree), circle.exp=nodelabsize, interactive=F)
    if (timelabs==T) ape::nodelabels(text=ape::branching.times(tree), node=1:tree$Nnode+Ntip(tree),
                                     frame="none", adj=timelaboffset, cex=timelabsize)
  }
  # plot multiple trees if a list greater than length=1 is given
  if (class(tree) == "list" && length(tree) > 1){

    par(mar=c(0, 2.1, 0, 1.1), mfrow=c(length(tree), 1), xpd=NA) # set pars, including a column of 1 plot/tree

    th <- max(sapply(tree, function(x) max(nodeHeights(x)))) # get maximum node height from all trees
    om <- 10^floor(log10(th)) # calculate order of magnitude of tree height
    xlims <- c(ceiling(th/om)*om, -om) # dynamically calculate a good x axis
    if (is.null(x.l.lim)==F) xlims[1] <- x.l.lim # override lower x-limit
    if (is.null(x.u.lim)==F) xlims[2] <- x.u.lim # override upper x-limit

    # loop through each tree
    for (i in 1:length(tree)){
      # if on the last tree, add a bottom margin to give room for the timescale
      if (i==length(tree)) par(mar=c(3.1, 2.1, 0, 1.1))

      ylims <- c(0.5, 1.2*Ntip(tree[[i]])) # dynamically calculate a good y axis for each tree
      plotTree(tree[[i]], direction="leftwards", xlim=xlims, ylim=ylims, mar=par()$mar, color="black")
      if (nodelabs==T) labelnodes(text=1:tree[[i]]$Nnode, node=1:tree[[i]]$Nnode+Ntip(tree[[i]]),
                                  circle.exp=nodelabsize/2, interactive=F)
      if (timelabs==T) ape::nodelabels(text=ape::branching.times(tree[[i]]), node=1:tree[[1]]$Nnode+Ntip(tree[[i]]),
                                       frame="none", adj=timelaboffset/0.8, cex=timelabsize/0.9)
    }
    axis(1, at=seq(0, xlims[1], by=om)) # plot timescale
    abline(v=timeslice, col=col) # plot timeslices. goes across all plots because par(xpd=T)
  }
  # return user's previous par settings
  par(pars)
}

#'Visualize response curves of taxa to climate variables
#'
#'Visualize response curves for taxa and climate variables. The main input can be either a "response table" (machu.1.tip.resp() output), for visualizing present-day taxa, or "ace table" (subsetted machu.2.ace() output), for visualizing ancestral taxa. In the case when uncertainty for reconstructed parameters is retained by machu.2.ace() (when unc=T), the function will visualize this uncertainty by reconstructing every possible combination of skew-normal distribution parameters (n=27) plus the median distribution as a thick dashed line.
#'
#'@param x Main input, either output from machu.1.tip.resp() or from machu.2.ace(). In the case of the latter it must be a single subsetted list element, i.e. ace[[1]]. If uncertainty was included from the machu.2.ace analysis, the plots will visualize it.
#'@param taxa Taxa to be plotted, a character or numeric vector specifying the desired taxa names or indices. If NULL, all taxa will be plotted.
#'@param clim Climate variables to be plotted, a character vector. If NULL, all climate variables will be plotted.
#'@param lty Line type for the plots. Default = 1 (solid line).
#'@param lwd Line width for the plots. Default = 1.
#'@param col Colors for the plots. If NULL, defaults to the "Set2" palette from RColorBrewer.
#'@param fill Whether to fill the area beneath each response curve. Default = F.
#'@param plot How to arrange plots of multiple climate variables. When "separate", plot each in its own window. When "together", plot all in the same window, using n2mfrow to calculate the best arrangement given a specified aspect ratio ('plot.asp'). Default = "separate"
#'@param plot.asp Aspect ratio used to calculate the best arrangement of plots when plot="together". Default = 16/9.
#'@param legend.cex Size of the taxon legend. Default = 1.
#'
#'@return Displays a plot or plots showing the response of the specified taxa to the specified climate variable(s).
#'
#'@examples
#'# Plot the responses of every taxon to every climate variable, together in the same window
#'machu.respplot(resp, plot="t")
#'# Plot the response of only the first taxon to the first climate variable, filling in the area under the curve
#'dev.off()
#'machu.respplot(resp, taxa=1, clim="bio1", fill=T)
#'# Plot the response of the first three taxa to the first three climate variables, in separate windows
#'machu.respplot(resp, taxa=1:3, clim=c("bio1", "bio2", "bio3"), fill=T, plot="s")
#'
#'# Plot the response of ancestral taxa at 1 Ma, without uncertainty
#'ace <- machu.2.ace(resp, tree, timeslice=1e6)
#'machu.respplot(ace[[1]], fill=T, plot="t")
#'# Plot the response of ancestral taxa at 1 Ma, with uncertainty
#'dev.off()
#'ace <- machu.2.ace(resp, tree, timeslice=1e6, unc=T)
#'machu.respplot(ace[[1]], fill=T, plot="t")
#'
#'# Compare the responses of present-day to ancestral (nodal) taxa with respect to the first climate variable
#'dev.off()
#'ace <- machu.2.ace(resp, tree)
#'machu.respplot(ace[[1]], clim="bio1", fill=T)
#'@import mgsub
#'@import scales
#'@import RColorBrewer
#'@import sn
#'@export
machu.respplot <- function(x, taxa=NULL, clim=NULL, lty=1, lwd=1, col=NULL, fill=F, plot="separate", plot.asp=16/9, legend.cex=1){

  ########
  # 1. detect input type (resp or ace) and convert to a universal format (simple dataframe)
  ########
  # machu.1.tip.resp input
  if (class(x)[1]=="matrix"){
    x <- as.data.frame(x)
  # machu.2.ace input
  } else if (class(x)[1]=="data.frame"){
    # figure out how many accessory data columns to remove
    if (colnames(x)[1]=="branch_start") cr <- 1:2 else if (colnames(x)[1]=="times") cr <- 1 else stop("Format should be the output of machu.1.tip.resp() or machu.2.ace().")
    rn <- rownames(x)
    x <- x[,-cr]
    x <- as.data.frame(lapply(x, unlist))
    rownames(x) <- rn
  } else if (class(x)[1]=="list") stop("Subset input from machu.2.ace() with double brackets.") else stop("Format should be the output of machu.1.tip.resp() or machu.2.ace().")

  # prune input table based on which taxa and climates were selected
  if (is.null(taxa)==F){
    x <- x[taxa,]
    if (is.na(x[1,1])==T) stop("None of the provided 'taxa' matched any in the input.")
  }
  if (is.null(clim)==F){
    ci <- c(sapply(clim, function(z) grep(paste0(z, "_"), colnames(x)))) # get indices of matching climate layers for each input in 'clim'
    if (length(ci)==1) stop("None of the provided 'clim' matched any in the input.")
    if (length(ci)%%3) stop("Check 'clim' input for mismatches.")
    x <- x[,ci]
  }

  # turn the table into a list of 1 table per climate variable
  cvs <- unique(mgsub(colnames(x), pattern=c("_mean.*", "_stdev.*", "_skew.*"), replacement=rep("", 3))) # get climate variable names
  x <- lapply(cvs, function(z) x[,grep(paste0(z, "_"), colnames(x))])
  names(x) <- cvs

  # if 9 slots for each clim are detected (i.e., uncertainty is present), plot all possible combinations at 50% opacity plus median at 100% opacity
  if (ncol(x[[1]])==3) unc <- F else if (ncol(x[[1]])==9) unc <- T else if (ncol(x[[1]])!=3 & ncol(x[[1]])!=9) stop("Check 'clim' input for mismatches.")

  ########
  # 2. make plots
  ########
  # set up plotting space
  if (plot == "together" | plot == "tog" | plot == "t"){
    arr <- n2mfrow(length(x), asp=plot.asp) # calculate best plot arrangement
    par(mfrow=arr, mar=c(2.1, 2.1, 1.7, 1.1))
  } else if (plot == "separate" | plot == "separately" | plot == "sep" | plot == "s"){
    par(mar=c(2.1, 2.1, 1.7, 1.1))
  }

  # make one plot per climate variable 'var' (element in 'x')
  # loop through climate variables names
  for (var in names(x)){

    # no uncertainty
    if (unc==F){

      # make x-axis
      means  <- x[[var]][,paste0(var, "_mean")]
      stdevs <- x[[var]][,paste0(var,"_stdev")]
      lobound <- floor(min(means-stdevs*3)) # lower bound for x axis
      upbound <- ceiling(max(means+stdevs*3)) # upper bound for x axis
      climval <- seq(from=lobound, to=upbound, by=(upbound-lobound)/1000) # x axis
      # make y axis
      density <- apply(x[[var]], 1, function(z) dsn(climval, xi=z[1], omega=z[2], alpha=z[3]))
      # colors
      if (is.null(col)==T) colors <- brewer.pal(n=ncol(density), name="Set2") else if (is.null(col)==F) colors <- col
      # main plot
      matplot(climval, density, type="l", xlab=NA, ylab=NA, main=var, lty=lty, lwd=lwd, col=colors)
      # polygons if fill=T
      if (fill==T) for (i in 1:ncol(density)){
        polygon(x=c(min(climval), climval, max(climval)), y=c(0, density[,i], 0), col=alpha(colors[i], 0.5), border=NA)
        lines(x=climval, y=density[,i], lty=lty, lwd=lwd, col=colors[i])
      }
      # legend
      legend("topleft", rownames(x[[var]]), pch=19, pt.cex=1.5, bty="n", col=colors, cex=legend.cex)

    # yes uncertainty
    } else if (unc==T){

      # get combos
      comb <- list()
      for (i in 1:nrow(x[[var]])){
        means  <- as.vector(as.matrix(x[[var]][i ,grep("_mean",  colnames(x[[var]]))]))
        stdevs <- as.vector(as.matrix(x[[var]][i ,grep("_stdev", colnames(x[[var]]))]))
        skews  <- as.vector(as.matrix(x[[var]][i ,grep("_skew",  colnames(x[[var]]))]))
        comb[[i]] <- expand.grid(list(means=means, stdevs=stdevs, skews=skews))
      }
      names(comb) <- rownames(x[[var]])

      # make x axis
      bounds <- lapply(comb, function(co){
        lobound <- floor(min(co$means-co$stdevs*3))
        upbound <- ceiling(max(co$means+co$stdevs*3))
        return(c(lobound, upbound))
      })
      bounds <- do.call(c, bounds)
      climval <- seq(from=min(bounds), to=max(bounds), by=(max(bounds)-min(bounds))/1000) # x axis

      # make y axis
      density <- lapply(comb, function(co) apply(co, 1, function(z) dsn(climval, xi=z[1], omega=z[2], alpha=z[3])))

      # colors
      if (is.null(col)==T) colors <- brewer.pal(n=length(density), name="Set2") else if (is.null(col)==F) colors <- col

      # make plot
      matplot(climval, density[[1]], type="l", xlab=NA, ylab=NA, main=var, lty=lty, lwd=lwd, col=alpha(colors[1], 0.5))
      # additional taxa
      if (length(density) > 1) for (i in 2:length(density)){
        matlines(climval, density[[i]], lty=lty, lwd=lwd, col=alpha(colors[i], 0.25))
      }
      # polygons if fill=T
      if (fill==T) for (i in 1:length(density)){
        apply(density[[i]], 2, function(p) polygon(x=c(min(climval), climval, max(climval)), y=c(0, p, 0),
                                                   col=alpha(colors[i], 0.1), border=NA))
      }
      # draw median lines for emphasis
      for (i in 1:nrow(x[[var]])){
        lines(climval, dsn(climval,
                           xi=x[[var]][i, grep("_mean$", colnames(x[[var]]))],
                           omega=x[[var]][i, grep("_stdev$", colnames(x[[var]]))],
                           alpha=x[[var]][i, grep("_skew$", colnames(x[[var]]))]),
              lty=2, lwd=lwd*2, col=colors[i])
      }
      # legend
      legend("topleft", rownames(x[[var]]), pch=19, pt.cex=1.5, bty="n", col=colors, cex=legend.cex)
    }
  }
}

#'Perform ancestral character estimation of response curves
#'
#'Extract reconstructed response curve parameters at a given set of timeslices, or at each tip and node
#'
#'@param tip.resp Output of machu.1.tip.resp().
#'@param tree Phylogenetic tree, a phylo object. Must be ultrametric, and taxon names must match those of 'tip.resp' exactly.
#'@param timeslice A numeric vector giving the timeslices to extract reconstructed climate response values from. If NULL, provides values for the tips and nodes of the tree. All values must be less than the overall height of the tree.
#'@param unc If TRUE, output upper and lower confidence intervals for each climate response value. Default = F.
#'@param csv.name A string that specifies the output CSV file's location and name. Should include .csv extension. Will be saved to getwd() unless another path is specified.
#'@param ace.method A string that specifies to ape::ace() which 'method' to use. Options are "REML" and "ML". See ?ape::ace for more details. Default = "REML".
#'@param verbose Print progress updates and output file location/name to screen. Default = F.
#'
#'@details
#'This function uses ace() from the package "ape" to generate climate response curves at each node of a time-calibrated phylogeny, and extracts the values along the branches at a particular time if so desired. By default, the output is a set of climate response curve parameters for each node and each tip in the tree. However, each node occurs at different poitns in time, making paleoclimatic projections tricky. Thus the user can also provide one or more timeslices with the 'timeslice' parameter. The function will find each branch present at each timeslice and interpolate the response curve parameters along the branches to those points in time. The function can also record uncertainty in ancestral character estimation with the 'unc' parameter.
#'
#'You can use phytools::force.ultrametric() to fix small deviations from ultrametricity in your tree. If you are using an alternative tree format from ape's 'phylo' format, the function will automatically attempt to convert it to phylo. If it fails, the function will probably break down the line. Try to convert to 'phylo' yourself before using it in the function.
#'
#'You should use a dated tree with this function, otherwise your projections into paleoclimate data will not make much sense. If you're really desperate, you can try ape::chronos() to get a quick-and-dirty dated tree, but I doubt it would be acceptable for publication. Really you should be using the output of time-calibration software such as BEAST or RelTime.
#'
#'The 'timeslice' parameter uses linear interpolation to calculate values for each "branch-taxon". For example, timeslice=1e6 will detect which branches exist at 1 Ma, and perform a linear interpolation of climate response values using the values at the subtending nodes/tips. See Guillory & Brown (2021) for a more detailed explanation. Of course, make sure to specify timeslices in the relative scale your tree's node ages are coded in. For example, BEAST node ages are often in units of millions of years, so "1" would be equivalent to the "1e6" I use in my example.
#'
#'The 'unc' parameter behavior has been changed for machuruku 2.0. Now setting unc=T only extracts the upper and lower 95\% confidence intervals for each climate response parameter provided by ace(). These are automatically incorporated into the Bioclim calculations in machu.3.anc.niche(), producing more conservative niche models. When "tips and nodes" are returned (not timeslices), the tips do not have associated confidence intervals because they are held over from the 'tip.resp' input, not produced via ace().
#'
#'Use machu.ace.load() to reload outputs saved with 'csv.name'.
#'
#'@return A list of one element per scenario, either one or more timeslices or all tips and nodes, each a table of reconstructed climate response parameters.
#'
#'@examples
#'# return timeslices for 1 and 3 Ma
#'ace <- machu.2.ace(tip.resp, tree, timeslice=c(1e6, 3e6))
#'# same as above, with uncertainty
#'ace <- machu.2.ace(tip.resp, tree, timeslice=c(1e6, 3e6), unc=T)
#'# return only tips and nodes
#'ace <- machu.2.ace(tip.resp, tree)
#'# same as above, saving output to getwd()
#'ace <- machu.2.ace(tip.resp, tree, csv.name="ace.csv")
#'@import ape
#'@export
machu.2.ace <- function(tip.resp, tree, timeslice=NULL, unc=F, csv.name=NULL, ace.method="REML", verbose=F){

  # Check if tree is a phylo object
  if (class(tree)!="phylo"){
    if (verbose==T) print("Tree is not a phylo... attempting to convert.")
    tree <- as.phylo(tree)
  }

  # Check if tree taxa match tip.response table
  treenames <- tree$tip.label
  if (identical(treenames, rownames(tip.resp)[match(treenames, rownames(tip.resp))]) == F){
    stop("Tree names and tip.resp names do not match.")
  }
  # various other checks
  if (is.ultrametric(tree) == F) stop("Tree is not ultrametric.")
  if (class(unc)!="logical") stop("Argument 'unc' must be TRUE or FALSE.")

  ###################
  ### 1: Get data ###
  ###################

  # Sort the response table so that the order of taxa matches that of the tree
  tip.resp <- as.data.frame(tip.resp[match(treenames, rownames(tip.resp)),])

  # Run ape::ace on every column (=trait) of the tip.resp table and output the results in a list
  acelist <- apply(tip.resp, 2, function(x) ace(x, tree, type="continuous", method=ace.method))

  #############################
  ### 2: Create data matrix ###
  #############################

  # Get various information
  ntips  <- Ntip(tree)
  nnodes <- Nnode(tree)
  treenames <- c(treenames, paste0(rep("Node", nnodes), 1:nnodes))
  times <- c(rep(0, ntips), branching.times(tree))
  names(times) <- treenames
  traitnames <- colnames(tip.resp)

  # combine branching times, present-day traits, and nodal traits into a single table
  table <- rbind(tip.resp, do.call(cbind, lapply(acelist, function(x) x$ace)))
  table <- cbind(times, table)
  rownames(table) <- treenames

  # if uncertainty is being considered, add this information
  if (unc==T){

    # get 95% confidence intervals from acelist
    CI95s <- lapply(acelist, function(x) x$CI95)

    # name the CI95 samples
    for (i in 1:length(CI95s)) colnames(CI95s[[i]]) <- paste0(traitnames[i], c("_lCI", "_uCI"))

    # turn into a table
    CI95s <- do.call(cbind, CI95s)

    # generate interleaving index vector
    index <- 1; n <- 0
    for (i in seq(from=2, length.out=length(traitnames))){
      index <- c(index, i, grep(traitnames[i-1], colnames(CI95s))+ncol(table))
    }

    # get filler values for present-day taxa and add to "CI95s"
    fill <- data.frame(table[1:ntips, rep(2:length(table), each=2)])
    colnames(fill) <- colnames(CI95s)
    CI95s <- rbind(fill, CI95s)

    # interleave with master table
    table <- cbind(table, CI95s)[index]
  }

  # the BM ace model is agnostic to what kind of data each variable is so some of them can become nonsensical, i.e. a negative stdev
  # convert all stdevs below zero to their absolute value (usually any negative stdevs will be small so it will still be a low number)
  # clamping at zero unfortunately doesn't make much sense either
  table[,grep("_stdev", colnames(table))] <- abs(table[,grep("_stdev", colnames(table))])
  # note: for certain variables, a mean<0 is also nonsensical (i.e. rainfall), but for some it is not (i.e. temperature).
  # Because the software is agnostic as to what each variable represents, I have to take the risk that a mean can fall below zero when it isn't supposed to

  #####################
  ### 3. Timeslices ###
  #####################

  # if timeslice is provided, take timeslices for each branch-taxon present
  if (is.null(timeslice) == F){

    if (max(timeslice) > max(times)) stop("At least one timeslice is older than the tree.")

    # get edge matrix from ape showing how nodes and tips are connected
    edge <- tree$edge

    # output a list of trait tables, one per timeslice
    table <- lapply(timeslice, function(ts){

      # get each branch subtending that timeslice (each branch (br) is a row in "edge")
      t <- do.call(rbind, apply(edge, 1, function(br){

        # test if that branch subtends the timeslice
        if (table[br[1],1] > ts && ts > table[br[2],1]){

          # get "branch-taxon" name and indices of branch start/end-points for output table
          info <- data.frame(branch_start=br[1],
                             branch_end=br[2])
          rownames(info) <- paste0(rownames(table)[br[1]], "-", rownames(table)[br[2]])

          # for every trait, calculate a, the value of the trait at that timeslice for that branch
          values <- do.call(cbind, apply(table[,-1], 2, function(trait){
            xi <- trait[br[1]]                 # trait value for tip/node i (branch origin)
            xj <- trait[br[2]]                 # trait value for tip/node j (branch end)
            vi <- abs(ts-table[br[1],1])       # distance from tip/node i to timeslice
            vj <- abs(ts-table[br[2],1])       # distance from tip/node j to timeslice
            a <- ((xi/vi)+(xj/vj))/((1/vi)+(1/vj)) # calculate a
            list(a)
          }))

          # combine info and values
          cbind(info, values)
        }
      }))
      if (verbose==T) print(paste0("Processed timeslice ", ts, "."))
      return(t)
    })
    names(table) <- paste0("timeslice_", timeslice)
  } else {
    # if timeslice is not provided, format the original tips & nodes only table as a single-element list and return that
    table <- list(tips_and_nodes=table)
  }
  # write parameters to file if csv.name was provided
  if (is.null(csv.name) == FALSE){
    if (verbose==T) print(paste0("Writing output to ", csv.name, "."))
    output <- data.frame()
    for (i in 1:length(table)){
      add <- data.frame(rep(names(table)[i], nrow(table[[i]])),
                        rownames(table[[i]]),
                        apply(table[[i]], 2, as.numeric))
      output <- rbind(output, add)
    }
    colnames(output)[1:2] <- c("scenario", "taxon")
    write.csv(output, csv.name, row.names=F)
  }
  return(table)
}

#'Load saved ace outputs from .csv file
#'
#'Loads outputs saved with the 'csv.name' argument of machu.2.ace() into a format compatible with machu.3.anc.niche().
#'
#'@param file String specifying the file to load.
#'
#'@details
#'Because of the idiosyncrasies of the machu.2.ace() code, the result of machu.ace.load() is not technically identical to one obtained from the former function. Both should be compatible with machu.3.anc.niche(), but let me know if you encounter bugs.
#'
#'@return A list of one element per scenario, either one or more timeslices or all tips and nodes, each a table of reconstructed climate response parameters.
#'
#'@examples
#'ace <- machu.ace.load("ace.csv")
#'@export
machu.ace.load <- function(file){

  # load file and get different scenarios
  raw <- read.csv(file)
  scenarios <- unique(raw$scenario)

  # convert raw table to proper format
  out <- lapply(scenarios, function(s){

    # create dataframe
    scenario <- subset(raw, scenario==s)
    rownames(scenario) <- scenario$taxon
    scenario <- scenario[,-(1:2)]
  })
  names(out) <- scenarios
  return(out)
}

#'Convert ancestral response values to niche models
#'
#'Create Bioclim niche models for each taxon in each timeslice
#'
#'@param ace Output of machu.2.ace(). Can be subsetted.
#'@param clim Paleoclimatic data. The best format is SpatRaster from the terra package, or a list of SpatRasters (via the function list()). A group of SpatRasters from multiple timeslices joined with c() is acceptable, but requires the use of raster.sets to differentiate. A RasterStack or list of RasterStacks (via list()) from the raster package is also acceptable.
#'@param taxa Optional, specifies which taxa to run. If the input is not timeslices but tips and nodes (names(ace)=="tips_and_nodes"), 'taxa' can be a numeric or character vector (i.e., 1:3, c(1,4), or c("taxon1", "taxon2")). If the input is timeslices (names(ace)!="tips_and_nodes"), only a character vector can be used. If none of the taxa are present in a given timeslice, that timeslice will be skipped.
#'@param resp.curv If TRUE, create ancestral niche models with skew-normal response curves. If FALSE, assume response is a uniform distribution within certain confidence limits specified by 'clip.amt'. This produces a binary niche model. Default = T.
#'@param clip.Q If TRUE, clip the tails of each response curve at certain confidence limits specified by 'clip.amt' (default is 0.95). Produces "cleaner" models. Default = T.
#'@param clip.amt Float value between 0 and 1 specifying confidence limits at which clip.Q=T or resp.curv=F operates. For example, at the default (0.95), the limits are the 0.025\% and 0.975\% quantiles for each response curve.
#'@param clip.samples The number of random samples taken by rsn() while determining confidence limits when clip.Q=T or resp.curv=F. Default = 10000. Decreasing may result in a speed boost at the cost of accuracy.
#'@param ocean A string or numeric specifying null or "ocean" pixels in the input rasters. Default = NA.
#'@param output.folder A string specifying a folder name to write outputs to. Outputs are only written if this argument is specified. Output will be in .tif format (GeoTiff) and given default names in the format "scenario_taxon.tif".
#'@param verbose If TRUE, print to screen certain checks, statuses, and progress updates. Default = F.
#'
#'@details
#'For a single timeslice, the SpatRaster or RasterStack object should have multiple layers, each corresponding to a climate variable. The layers should be consistent in number and type across all timeslices. Within a timeslice, layers should share the same geography (i.e., ocean pixels should be the same).
#'
#'The function associates timeslices with rasters based on the inputs. When the same number of timeslices and rastersets are provided (length(ace)==length(clim)), the first timeslice will be matched to the first rasterset, and so on. When only one rasterset and multiple timeslices are provided, the timeslices will all be matched to the single rasterset, and vice versa. When 'ace' and 'clim' have different lengths, and lengths > 1, only the first 'x' elements of 'ace' and 'clim' will be used, 'x' being the length of the shorter input, with the first timeslice matched to the first rasterset, and so on. When verbose=TRUE, these rasterset-timeslice associations are printed to the screen.
#'
#'Whether uncertainty is included in the ancestral niche models depends on the 'ace' input, i.e. whether unc=T was specified in machu.2.ace(). Uncertainty is encoded by the mean, upper confidence limit, and lower confidence limit for each parameter (mean, standard deviation, and skew). When passed to machu.3.anc.niche(), multiple response curves consisting of every combination of these values (n=3^3=27) are used to construct suitabilities for each raster and then summed. This tends to produce much broader (more conservative) models than when unc=F.
#'
#'When uncertainty is included and resp.curv=F (i.e., binary models are produced), the behavior differs from when uncertainty is not included. A second "clipping" occurs after the multiple suitabilities for that raster (from the 27 combinations of response curve parameters) are summed, because at this point the summed suitability will still not be binary (0 or 1). Suitability values corresponding to raster pixels in the '1-clip.amt' (default=0.05) percentile of suitable values are set to 0, and all values above it are set to 1. This produces a cleaner look more in line with expectations.
#'
#'You can plot the output models on your own, or use machu.plotmap() to do so.
#'
#'The output option only outputs tiffs, and constructs automatic filenames using the raster set and taxon for each layer. This may be too limited for your purposes, in which case you should write your own script to output the file types and/or names that you desire for your project.
#'
#'@return A list of one element per scenario (timeslice or rasterset), each itself a list of one SpatRaster (ancestral niche model) per taxon.
#'
#'@examples
#'## acceptable 'clim' formats
#'# Single SpatRaster (terra)
#'clim1 <- c(rast(list.files(rasterfolder, pattern="T1_", full.names=T)))
#'# List of SpatRasters (preferred)
#'clim2 <- list(c(rast(list.files(rasterfolder, pattern="T1_", full.names=T))), c(rast(list.files(rasterfolder, pattern="T2_", full.names=T))))
#'# Multiple rastersets in a single SpatRaster (requires 'raster.sets' be specified)
#'clim3 <- c(c(rast(list.files(rasterfolder, pattern="T1_", full.names=T))), c(rast(list.files(rasterfolder, pattern="T2_", full.names=T))))
#'# Single RasterStack (raster)
#'clim4 <- stack(rast(list.files(rasterfolder, pattern="T1_", full.names=T)))
#'# List of RasterStacks
#'clim5 <- list(stack(list.files(rasterfolder, pattern="T1_", full.names=T)), stack(list.files(rasterfolder, pattern="T2_", full.names=T)))
#'
#'## various 'ace' formats
#'# Two timeslices, no uncertainty
#'ace_n <- machu.2.ace(resp, tree, timeslice=c(1e6, 2e6))
#'# Two timeslices, uncertainty
#'ace_u <- machu.2.ace(resp, tree, timeslice=c(1e6, 2e6), unc=T)
#'# No timeslices--only tips and nodes, no uncertainty
#'ace_t <- machu.2.ace(resp, tree)
#'
#'## ways of running the actual function
#'# run with default settings, no uncertainty, two timeslices, two rastersets in SpatRaster format
#'output.models <- machu.3.anc.niche(ace_n, clim2)
#'# same as above, but don't clip the models
#'output.models <- machu.3.anc.niche(ace_n, clim2, clip.Q=F)
#'# same as above, but produce binary models
#'output.models <- machu.3.anc.niche(ace_n, clim2, resp.curv=F)
#'# same as above, but make the binary models more stringent
#'output.models <- machu.3.anc.niche(ace_n, clim2, resp.curv=F, clip.amt=0.9)
#'# same as above, but output the models to the present working directory
#'output.models <- machu.3.anc.niche(ace_n, clim2, resp.curv=F, clip.amt=0.9, output.folder=getwd())
#'
#'# run with uncertainty (default settings, two timeslices, two rastersets)
#'output.models <- machu.3.anc.niche(ace_u, clim2)
#'
#'# run with two timeslices, but project them both into the same single rasterset
#'output.models <- machu.3.anc.niche(ace_n, clim1)
#'# run with one timeslice, but project it into two different rastersets
#'output.models <- machu.3.anc.niche(ace_n[1], clim2)
#'
#'# run with only certain taxa
#'output.models <- machu.3.anc.niche(ace_n, clim2, taxa=c("taxon1", "taxon2"))
#'
#'# run with no actual timeslices, only tips and nodes
#'output.models <- machu.3.anc.niche(ace_t, clim1)
#'# same as above, but project all taxa into multiple timeslices
#'output.models <- machu.3.anc.niche(ace_t, clim2)
#'# same as above, but only use the first two taxa (both are equivalent)
#'output.models <- machu.3.anc.niche(ace_t, clim1, taxa=c("taxon1", "taxon2"))
#'output.models <- machu.3.anc.niche(ace_t, clim1, taxa=1:2)
#'
#'# plot the first taxon's ancestral niche model from the first timeslice
#'plot(output.models[[1]][[1]])
#'@import terra
#'@import sn
#'@export
machu.3.anc.niche <- function(ace, clim, taxa=NULL, resp.curv=T, clip.Q=T, clip.amt=0.95, clip.samples=10000, ocean=NA, output.folder=NULL, verbose=F) {

  #############################################
  ### 1: Prune taxon set and perform checks ###
  #############################################

  if (class(ace)!="list") stop("'ace' must be a list. If subsetting, use [] rather than [[]].")

  if (is.null(taxa)==F){
    if (class(taxa)!="character" && names(ace)[1]!="tips_and_nodes") stop("If 'ace' is composed of timeslices, 'taxa' must be a character vector, because the same taxa may have different numerical indices in different timeslices.")
    # if specific taxa are selected, modify ace to only include those taxa
    ace <- lapply(ace, function(ts) if (is.null(taxa)==F) na.omit(ts[taxa,]))
    # get rid of any timeslices where no taxa were retained
    remove <- which(lapply(ace, nrow)==0)
    if (length(remove) > 0) ace <- ace[-remove]
    if (length(ace)==0) stop("No taxa retained from 'ace'. Check 'taxa' for misspellings or other errors.")
    if (verbose==T){
      print("Retained the following scenarios and taxa:")
      for (i in 1:length(ace)) print(paste0(names(ace)[i], ": ", paste(rownames(ace[[i]]), collapse=", ")))
    }
  }

  # detect whether unc was applied to 'ace' in machu.2.ace()
  if (length(grep("_uCI$", colnames(ace[[1]]))) > 0 && length(grep("_lCI$", colnames(ace[[1]]))) > 0){
    unc <- T
    if (verbose==T) print("Detected uncertainty samples in 'ace' ('unc' in machu.2.ace()).")
  } else {
    unc <- F
    if (verbose==T) print("Did not detect uncertainty samples in 'ace' ('unc' in machu.2.ace()).")
  }

  # detect which columns need to be removed later because they include accessory info (depends on whether the form of 'ace' is "tips_and_nodes")
  if (names(ace)[1]!="tips_and_nodes") cr <- c(1,2) else cr <- 1

  # check that clip.amt is between 0 and 1
  if (clip.amt < 0 | clip.amt > 1) stop("'clip.amt' must be between 0 and 1.")

  ##############################
  ### 2: Format climate data ###
  ##############################

  # solo RasterStack: convert to SpatRaster and put into a single-element list
  if (class(clim)=="RasterStack"){
    clim <- list(rast(clim))
    if (verbose==T) print("Solo RasterStack converted to SpatRaster format.")
  }
  # list of RasterStacks: convert to SpatRaster
  if (class(clim)=="list" && class(clim[[1]])=="RasterStack"){
    clim <- lapply(clim, rast)
    if (verbose==T) print("Multiple RasterStacks converted to SpatRaster format.")
  }
  # solo SpatRaster: put into a single-element list
  if (class(clim)=="SpatRaster" && is.null(raster.sets)==T){
    clim <- list(clim)
    print("Warning: solo SpatRaster detected, with no raster.sets specified; treating as a single timeslice.")
  }
  # Multiple SpatRasters, joined with c(): separate into list elements based on raster.sets
  if (class(clim)=="SpatRaster" && is.null(raster.sets)==F){
    clim <- lapply(raster.sets, function(x) clim[[grep(x, names(clim))]])
    print("Warning: solo SpatRaster detected, with raster.sets specified; separated into multiple timeslices.")
  }
  # add names to clim if there aren't any (the names are very basic and agnostic to the corresponding time period)
  if (is.null(names(clim))==T) names(clim) <- paste0("rasterSet", 1:length(clim))
  # final check for proper climate data format
  if (class(clim)!="list" && class(clim[[1]])!="SpatRaster") stop("'clim' is not provided in the proper format (RasterStack/SpatRaster or list thereof, or raster-containing folder).")

  ########################################
  ### 3: reconstruct climate responses ###
  ########################################

  # create association table where each row matches a raster to a timeslice (by index)
  # the first column of 'assoc' is raster indices ('clim'); the second column is timeslice indices ('ace')
  # 'all.out.names' will be the names given to the main elements of the output list
  # 'i' is a column index counter used in reports to the screen
  # 1. clim & ace are same length: raster 1 corresponds to timeslice 1, and so on
  if (length(clim)==length(ace)){
    assoc <- matrix(rep(1:length(ace), 2), ncol=2)
    all.out.names <- names(ace)
    i=2
  }
  # 2. clim of length 1 & multiple timeslices: all timeslices correspond to the same singular raster
  if (length(clim)==1 && length(clim) < length(ace)){
    assoc <- matrix(c(rep(1, length(ace)), 1:length(ace)), ncol=2)
    all.out.names <- names(ace)
    i=2
  }
  # 3. multiple rasters & ace of length 1: all rasters correspond to the same singular timeslice
  if (length(ace)==1 && length(ace) < length(clim)){
    assoc <- matrix(c(1:length(clim), rep(1, length(clim))), ncol=2)
    all.out.names <- names(clim)
    i=1
  }
  # 4. ace & clim are unequal length, & all length>1: raster 1 corresponds to timeslice 1, and so on, until either ace or clim ends (some will be unused)
  if (length(ace)>1 && length(clim)>1 && length(clim)!=length(ace)){
    assoc <- matrix(rep(intersect(1:length(ace), 1:length(clim)), 2), ncol=2)
    all.out.names <- names(ace)
    i=2
  }
  colnames(assoc) <- c("clim", "ace")
  if (verbose==T) print("Associating scenarios ('ace') with rastersets ('clim') (numbers are indices, rows are associations):")
  if (verbose==T) print(assoc)

  # for every timeslice/rasterset combo (i.e., a row in assoc)...
  all.out <- apply(assoc, 1, function(timeset){

    # get all taxa in that timeslice
    timeslice <- ace[[timeset[2]]][,-cr]
    # 'cr' was set earlier and removes columns containing extraneous info (branch start/end for timeslices, time for tips_and_nodes)

    # for every taxon in that timeslice...
    out.bioclims <- lapply(split(timeslice, seq(nrow(timeslice))), function(taxon){

      # create a matrix of indices for each climate response value ('crv')
      # each row corresponds to a climate layer, each number an index of a climate response value for it
      if (unc==F) n.lyr <- length(taxon)/3 else if (unc==T) n.lyr <- length(grep("uCI$", names(taxon)))/3
      crv <- matrix(1:length(taxon), nrow=n.lyr, byrow=T)

      # for every climate layer (by index)...
      lyrs <- lapply(1:n.lyr, function(lyr.ind){

        # get principal climate response values for this layer
        lyr.pars <- names(taxon[crv[lyr.ind,]])
        mean  <- taxon[,lyr.pars[grep("mean",  lyr.pars)]]
        stdev <- taxon[,lyr.pars[grep("stdev", lyr.pars)]]
        skew  <- taxon[,lyr.pars[grep("skew",  lyr.pars)]]

        # this section is included for compatibility with loaded ace files from machu.ace.load()
        # in that case, the above section is enough. For results direct from machu.2.ace(), this section is necessary to unlist climate response values
        if (class(mean)!="numeric")  mean  <- unlist(mean)
        if (class(stdev)!="numeric") stdev <- unlist(stdev)
        if (class(skew)!="numeric")  skew  <- unlist(skew)

        # get raster layer from within the rasterset & convert to vector form
        lyr <- c(as.matrix(clim[[timeset[1]]][[lyr.ind]]))
        lyr[is.na(lyr)] <- 0

        # convert raster to suitability:
        # (what's cool about this part is that it works regardless of if unc==T or not)
        # get every possible combination (n=27) of median, 95% upper CI, and 95$ lower CI for mean, stdev, and skew
        combos <- expand.grid(list(mean=mean, stdev=stdev, skew=skew))

        # project raster into suitability according to each combination of skew-normal parameters
        # 'projVal' is a table of 27 cols where each col is the suitability under a different combo of sn params
        projVal <- apply(combos, 1, function(combo){
          combo <- as.numeric(combo)
          # project raster into suitability using skew-normal dist (xi=mean, omega=sd, alpha=skew)
          y <- dsn(lyr, xi=combo[1], omega=combo[2], alpha=combo[3])

          # if clip.Q is on, clip quantiles from distributions to produce cleaner models
          if (clip.Q==T | resp.curv==F){
            # estimate quantiles for the skew-normal distribution (decrease clip.samples for possible speed boost, but less accuracy)
            Q <- quantile(x=rsn(clip.samples, xi=combo[1], omega=combo[2], alpha=combo[3]),
                           probs=c((1-clip.amt)/2, clip.amt+(1-clip.amt)/2), names=F)
            # remove suitability values beyond those quantiles according to their position in the raster
            y[lyr < Q[1] | lyr > Q[2]] <- 0
            # if resp.curv is off, go further and convert all suitability values between quantiles to 1, producing a binary model
            if (resp.curv==F) y[lyr >= Q[1] & lyr <= Q[2]] <- 1
          }
          return(y)
        })

        # sum the different runs, should provide a more generalized result
        projVal <- rowSums(projVal)

        # if unc=T & resp.curv=F, the output bioclim will not be binary. Convert to binary here
        if (unc==T & resp.curv==F){
          # get a lower quantile specific to the distribution of unique summed suitability values in projVal
          LQ <- quantile(x=unique(projVal), probs=1-clip.amt, names=F)
          # trim the model by converting all values lower than LQ to 0
          projVal[projVal < LQ] <- 0
          # make a binary model by converting all nonzero values to 1
          projVal[projVal >= LQ] <- 1
        }

        # rescale suitability values from 0-1 so that they can be compared
        projVal <- projVal/max(projVal)

        # convert suitability values back into raster form
        r <- matrix(projVal, nrow=dim(clim[[timeset[1]]][[lyr.ind]])[1], ncol=dim(clim[[timeset[1]]][[lyr.ind]])[2], byrow=T)
        r <- rast(r)
        set.ext(r, ext(clim[[timeset[1]]][[lyr.ind]])) # match extent
        crs(r) <- crs(clim[[timeset[1]]][[lyr.ind]]) # match coordinate reference system
        names(r) <- names(clim[[timeset[1]]][[lyr.ind]])

        return(r)
      })

      # convert 'lyrs' to a single SpatRaster object with multiple layers
      # yes I know this is the weirdest possible way to do it but idk how else. I guess the smart thing would have been to use app or lapp from the beginning but I didn't know about it then and now it's too late
      out.bioclim <- eval(parse(text=paste0("c(", paste0("lyrs[[", 1:length(lyrs), "]]", collapse=","),")")))

      # find limiting values from all rasters---this actually creates the Bioclim niche model
      out.bioclim <- app(out.bioclim, min)

      # create a mask that excludes ocean pixels
      mask <- clim[[timeset[[1]]]][[1]] # just get one example of raster from this timeslice (all layers should have same ocean pixels)
      if (is.na(ocean)) mask[!is.na(mask)] <- 1 else mask[mask != ocean] <- 1

      # use mask to exclude ocean pixels from output niche model
      out.bioclim <- mask*out.bioclim

      names(out.bioclim) <- rownames(taxon)
      return(out.bioclim)
    })
    names(out.bioclims) <- rownames(timeslice)
    if (verbose==T) print(paste("Processed scenario", timeset[i]))
    return(out.bioclims)
  })
  names(all.out) <- all.out.names

  # write rasters to files if an output folder was specified
  if (is.null(output.folder)==F){
    if (verbose==T) print(paste0("Writing ", sum(sapply(all.out, length)), " rasters to ", output.folder, ":"))
    for (i in 1:length(all.out)){
      for (j in 1:length(all.out[[i]])){
        outname <- paste0(names(all.out)[i], "_", names(all.out[[i]])[j], ".tif")
        terra::writeRaster(all.out[[i]][[j]], outname, filetype="GTiff", overwrite=T)
        if (verbose==T) print(outname)
      }
    }
  }
  return(all.out)
}

#'Map model output from machu.3.anc.niche()
#'
#'Plot one or more niche models using model output from machu.3.anc.niche(). User can specify color ramps inspired by the viridis package and plot maps in various configurations.
#'
#'@param model Output from machu.3.anc.niche(), a list of lists of SpatRasters, generally one per taxon per timeslice. Subsettable, though this may interfere with the format of titles.
#'@param col Desired color ramp, usually a numeric value from 1 to 6. 1 = base; 2 = plasma; 3 = viridis; 4 = l17 pure white; 5 = r3 pure blue; 6 = black and white. Ramps can also be specified with strings, e.g. "viridis" or "v". Custom color ramps (character vectors with length > 1) can also be supplied, for example by the colorRampPalette function. Default = 1 ("base").
#'@param axes Whether to display axes (lat/long). Default = T.
#'@param title Whether to display titles. These are supposed to take the format "timeslice, taxon", but certain input formats from machu.3.anc.niche, or subsetting the input, may interfere with this. Default = T.
#'@param title.cex Specify size of titles. Default = 1.
#'@param to.scale Whether to color each raster according to the same scale (i.e., the min and max value among all rasters) so that they can be directly compared. Default = F.
#'@param plot How to arrange plots of multiple climate variables. When "separate", plot each in its own window. When "together", plot all in the same window, using n2mfrow to calculate the best arrangement given a specified aspect ratio ('plot.asp'). Default = "separately".
#'@param plot.asp Aspect ratio used to calculate the best arrangement of plots when plot="together". Default = 16/9.
#'
#'@details
#'Titles may not display properly (timeslice, taxon) depending on the input format from machu.3.anc.niche and the amount of subsetting of the input.
#'
#'Implementing scale bars was more trouble than I wanted to deal with so currently they are not supported with this function. Using the plot() function will automatically generate one for your raster though.
#'
#'When to.scale=T, the function calculates the minimum and maximum suitability value among all the input rasters, and scales the color for each one using this range. This enables the rasters to be directly compared; as such, when to.scale=F, comparing the rasters can be misleading.
#'
#'@return Plots of niche models in various configurations.
#'
#'@examples
#'# plot all models, separately, with default color ramp, axes, and titles
#'machu.plotmap(model)
#'# together
#'machu.plotmap(model, plot="t")
#'# no axes or titles
#'machu.plotmap(model, plot="t", axes=F, title=F)
#'# plot all rasters on the same color scale
#'machu.plotmap(model, plot="t", axes=F, title=F, to.scale=T)
#'
#'# plot only the first raster, using the "viridis" color ramp
#'machu.plotmap(model[[1]][[1]], col=3)
#'
#'# plot only rasters from the first timeslice/scenario, using a custom "Mardi Gras" color ramp
#'mardi.gras <- colorRampPalette(c("purple", "forestgreen", "gold"))(256)
#'machu.plotmap(model[[1]], col=mardi.gras, plot="t")
#'@import terra
#'@export
machu.plotmap <- function(model, col=1, axes=T, title=T, title.cex=1, to.scale=F, plot="separately", plot.asp=16/9) {

  # pick color
  colr.o <- list()

  if (length(col)>1 & is.character(col)==T){
    colr.o <- col
  } else {
    # 1. base
    if (col==1 | col=="base" | col=="b") colr.o <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
    # 2. plasma
    if (col==2 | col=="plasma" | col=="p") colr.o<- c("#000000", "#000C7D", "#000D7E", "#000D80", "#000E82", "#000E84", "#000E86", "#000F87", "#000F89", "#00108B", "#00108C", "#00118E", "#001190", "#001191", "#001293", "#001294", "#001296", "#001397", "#001399", "#00139A", "#00149B", "#00149D", "#00149E", "#00149F", "#0015A0", "#0015A1", "#0015A2", "#0015A3", "#0015A4", "#0015A5", "#0016A6", "#0016A7", "#0016A7", "#0016A8", "#0016A9", "#0016A9", "#0016A9", "#0A16AA", "#1516AA", "#1D15AA", "#2315AA", "#2915AA", "#2F15A9", "#3414A9", "#3914A8", "#3E13A7", "#4313A6", "#4712A5", "#4C12A4", "#5011A3", "#5311A2", "#5710A1", "#5A0FA0", "#5E0F9F", "#610E9F", "#640E9E", "#670D9D", "#6A0D9C", "#6C0C9B", "#6F0B9A", "#720B99", "#740A99", "#770A98", "#790997", "#7C0897", "#7E0896", "#800795", "#820795", "#850694", "#870693", "#890693", "#8B0592", "#8D0592", "#8F0491", "#910491", "#930490", "#950390", "#97038F", "#99038F", "#9B028E", "#9D028E", "#9F028E", "#A1018D", "#A3018D", "#A5018C", "#A7018C", "#A9008B", "#AB008B", "#AC008A", "#AE008A", "#B00089", "#B20089", "#B40088", "#B60088", "#B80087", "#B90087", "#BB0087", "#BD0086", "#BF0086", "#C10085", "#C30085", "#C40084", "#C60084", "#C80083", "#CA0083", "#CC0082", "#CE0082", "#CF0081", "#D10081", "#D30080", "#D50080", "#D6007F", "#D8007F", "#DA007E", "#DB007E", "#DD007D", "#DE007C", "#E0017C", "#E2027B", "#E3047B", "#E5067A", "#E6087A", "#E80B79", "#E90D78", "#EA1078", "#EC1277", "#ED1477", "#EE1676", "#F01875", "#F11A75", "#F21C74", "#F41E73", "#F52073", "#F62272", "#F72471", "#F82671", "#F92870", "#FB2A6F", "#FC2C6F", "#FD2E6E", "#FE306D", "#FF326C", "#FF346C", "#FF366B", "#FF386A", "#FF3A6A", "#FF3D69", "#FF3F68", "#FF4167", "#FF4366", "#FF4566", "#FF4765", "#FF4964", "#FF4B63", "#FF4D62", "#FF5062", "#FF5261", "#FF5460", "#FF565F", "#FF585E", "#FF5A5D", "#FF5D5C", "#FF5F5B", "#FF615B", "#FF635A", "#FF6559", "#FF6758", "#FF6A57", "#FF6C56", "#FF6E55", "#FF7054", "#FF7253", "#FF7452", "#FF7651", "#FF7850", "#FF7A4E", "#FF7C4D", "#FF7E4C", "#FF7F4B", "#FF814A", "#FF8349", "#FF8548", "#FF8747", "#FF8845", "#FF8A44", "#FF8C43", "#FF8E42", "#FF8F40", "#FF913F", "#FF933E", "#FF953C", "#FF963B", "#FF9839", "#FF9A38", "#FF9B36", "#FF9D35", "#FF9F33", "#FFA032", "#FFA230", "#FFA32F", "#FFA52E", "#FFA62C", "#FFA82B", "#FFA92A", "#FFAB29", "#FFAC28", "#FFAE27", "#FFAF26", "#FFB126", "#FFB225", "#FFB424", "#FFB523", "#FFB623", "#FFB822", "#FFB922", "#FFBB21", "#FFBC20", "#FFBD20", "#FFBF1F", "#FFC01F", "#FFC21F", "#FFC31E", "#FFC41E", "#FFC61E", "#FFC71D", "#FFC81D", "#FFCA1D", "#FFCB1D", "#FFCC1D", "#FFCE1D", "#FFCF1C", "#FFD01C", "#FFD21C", "#FFD31C", "#FFD41C", "#FFD61C", "#FFD71D", "#FFD81D", "#FFDA1D", "#FFDB1D", "#FFDC1D", "#FFDE1D", "#FFDF1E", "#FFE01E", "#FFE21E", "#FFE31E", "#FFE41F", "#FFE61F", "#FFE71F", "#FFE820", "#FFE920", "#FFEB21", "#FFEC21", "#FFED22", "#FFEF22", "#FFF123")
    # 3. viridis
    if (col==3 | col=="viridis" | col=="v") colr.o <- colorRampPalette(c("#440154", "#440154ff", "#440558ff", "#450a5cff", "#450e60ff", "#451465ff", "#461969ff", "#461d6dff", "#462372ff", "#472775ff", "#472c7aff", "#46307cff", "#45337dff", "#433880ff", "#423c81ff", "#404184ff", "#3f4686ff", "#3d4a88ff", "#3c4f8aff", "#3b518bff", "#39558bff", "#37598cff", "#365c8cff", "#34608cff", "#33638dff", "#31678dff", "#2f6b8dff", "#2d6e8eff", "#2c718eff", "#2b748eff", "#29788eff", "#287c8eff", "#277f8eff", "#25848dff", "#24878dff", "#238b8dff", "#218f8dff", "#21918dff", "#22958bff", "#23988aff", "#239b89ff", "#249f87ff", "#25a186ff", "#25a584ff", "#26a883ff", "#27ab82ff", "#29ae80ff", "#2eb17dff", "#35b479ff", "#3cb875ff", "#42bb72ff", "#49be6eff", "#4ec16bff", "#55c467ff", "#5cc863ff", "#61c960ff", "#6bcc5aff", "#72ce55ff", "#7cd04fff", "#85d349ff", "#8dd544ff", "#97d73eff", "#9ed93aff", "#a8db34ff", "#b0dd31ff", "#b8de30ff", "#c3df2eff", "#cbe02dff", "#d6e22bff", "#e1e329ff", "#eae428ff", "#f5e626ff", "#fde725ff")) (256)
    # 4. l17 pure white
    if (col==4 | col=="l17" | col=="l" | col=="white" | col=="w") colr.o<- colorRampPalette(c("#FFFFFF", "#FFFEFC", "#FEFDF9", "#FEFDF7", "#FDFCF4", "#FDFBF1", "#FCFAEE", "#FCFAEB", "#FBF9E9", "#FAF7E3", "#FAF7E0", "#F9F6DD", "#F8F5DB", "#F8F4D8", "#F7F3D5", "#F7F3D2", "#F6F2CF", "#F6F1CD", "#F4F0C7", "#F4EFC4", "#F3EEC2", "#F2EDBF", "#F2ECBD", "#F2EBBB", "#F2EABA", "#F2E9B8", "#F3E8B6", "#F3E6B3", "#F3E5B1", "#F3E4AF", "#F3E3AE", "#F3E2AC", "#F3E1AA", "#F3E0A9", "#F2DFA7", "#F2DEA5", "#F2DCA2", "#F2DBA0", "#F2DA9F", "#F2D99D", "#F2D79C", "#F2D69A", "#F2D598", "#F2D497", "#F2D396", "#F2D193", "#F3D092", "#F3CF91", "#F3CD90", "#F3CC8E", "#F3CB8D", "#F3CA8C", "#F3C98B", "#F3C688", "#F4C587", "#F4C486", "#F4C385", "#F4C284", "#F4C183", "#F4BF81", "#F4BE80", "#F4BD7F", "#F4BB7D", "#F4BA7C", "#F4B87B", "#F4B77A", "#F4B67A", "#F5B579", "#F5B378", "#F5B277", "#F5B177", "#F5AE75", "#F5AD75", "#F5AC74", "#F5AB73", "#F5AA73", "#F5A872", "#F5A771", "#F5A670", "#F5A470", "#F5A26E", "#F5A16E", "#F59F6D", "#F59E6C", "#F59D6C", "#F59C6C", "#F59A6B", "#F5996B", "#F5986B", "#F5956A", "#F5946A", "#F5936A", "#F5916A", "#F5906A", "#F48F69", "#F48D69", "#F48C69", "#F48A69", "#F48868", "#F48768", "#F48668", "#F48468", "#F38367", "#F38167", "#F38067", "#F37F67", "#F27C67", "#F27B67", "#F27A68", "#F27868", "#F17768", "#F17668", "#F17468", "#F07368", "#F07269", "#EF6F69", "#EF6E69", "#EF6C69", "#EE6B6A", "#EE6A6A", "#EE686A", "#ED676A", "#ED666A", "#ED646A", "#EC616B", "#EC606B", "#EB5F6B", "#EA5D6C", "#EA5C6C", "#E95B6D", "#E85A6D", "#E8586E", "#E7576E", "#E6556F", "#E5536F", "#E55270", "#E45170", "#E34F71", "#E34E71", "#E24D72", "#E14B72", "#E04973", "#E04773", "#DF4674", "#DE4474", "#DE4375", "#DD4175", "#DC4075", "#DB3F76", "#DA3E77", "#D83C78", "#D73B79", "#D63A79", "#D5397A", "#D4377A", "#D3367B", "#D2357C", "#D1347C", "#D0337D", "#CE307E", "#CD2F7F", "#CC2E7F", "#CB2D80", "#CA2B80", "#C92A81", "#C82982", "#C72882", "#C52683", "#C32584", "#C12485", "#C02385", "#BE2386", "#BD2287", "#BC2287", "#BA2188", "#B92089", "#B72089", "#B41F8B", "#B31E8B", "#B11D8C", "#B01D8C", "#AE1C8D", "#AD1C8E", "#AB1B8E", "#AA1B8F", "#A71990", "#A51991", "#A31892", "#A21892", "#A01893", "#9E1993", "#9C1994", "#9A1A94", "#981A95", "#941B96", "#921B97", "#901B97", "#8E1C98", "#8C1C99", "#891C99", "#871D9A", "#851D9A", "#831D9B", "#7F1E9C", "#7C1E9D", "#7A1E9D", "#781E9E", "#751F9E", "#731F9F", "#701F9F", "#6E20A0", "#6B21A0", "#6522A1", "#6223A1", "#5F23A2", "#5C24A2", "#5924A2", "#5525A3", "#5225A3", "#4E26A3", "#4B26A4", "#4327A4", "#3E27A5", "#3A28A5", "#3428A6", "#2F29A6", "#2929A6", "#2129A7", "#182AA7", "#002AA8"))(256)
    # 5. r3 pure blue
    if (col==5 | col=="r3" | col=="r" | col=="blue") colr.o<-colorRampPalette(c("#085CF8", "#0F5FF4", "#1361F1", "#1763ED", "#1965E9", "#1B67E5", "#1C6AE1", "#1D6CDE", "#1D6EDA", "#1D72D2", "#1D74CE", "#1C75CB", "#1B77C7", "#1979C3", "#187BBF", "#167DBB", "#147EB8", "#1380B4", "#1283AC", "#1385A8", "#1486A4", "#1788A0", "#1A899C", "#1D8A98", "#208C93", "#238D8F", "#278E8B", "#2D9082", "#2F927E", "#329379", "#349475", "#359570", "#37966C", "#389767", "#399862", "#3A9A5E", "#3B9C54", "#3B9D50", "#3C9E4B", "#3C9F46", "#3CA042", "#3DA13D", "#3DA239", "#3EA335", "#3FA431", "#42A62B", "#44A728", "#47A826", "#49A824", "#4CA922", "#4EAA21", "#51AA20", "#54AB20", "#5AAC1F", "#5CAD1E", "#5FAE1E", "#62AE1E", "#65AF1E", "#67AF1D", "#6AB01D", "#6CB11D", "#6FB11D", "#74B21D", "#77B31C", "#79B41C", "#7BB41C", "#7EB51C", "#80B51B", "#83B61B", "#85B61B", "#87B71B", "#8CB81A", "#8EB91A", "#91B91A", "#93BA1A", "#95BA19", "#97BB19", "#9ABB19", "#9CBC19", "#9EBD18", "#A3BE18", "#A5BE18", "#A7BF17", "#A9BF17", "#ABC017", "#AEC016", "#B0C116", "#B2C116", "#B4C215", "#B8C315", "#BBC314", "#BDC414", "#BFC514", "#C1C513", "#C3C613", "#C5C613", "#C7C712", "#CCC812", "#CEC811", "#D0C911", "#D2C910", "#D4CA10", "#D6CA0F", "#D8CB0F", "#DBCB0E", "#DDCB0E", "#E1CC0E", "#E3CD0E", "#E5CD0E", "#E7CD0F", "#E9CE10", "#EBCE12", "#EDCE14", "#EFCD16", "#F1CD19", "#F4CC1F", "#F6CB23", "#F7CA26", "#F8C929", "#F9C82D", "#F9C730", "#FAC633", "#FAC536", "#FBC339", "#FCC03F", "#FCBF41", "#FCBE44", "#FCBC46", "#FCBB48", "#FDB94B", "#FDB84D", "#FDB64F", "#FDB551", "#FDB255", "#FEB057", "#FEAF59", "#FEAD5B", "#FEAC5D", "#FEAA5F", "#FEA961", "#FEA763", "#FEA466", "#FEA368", "#FEA16A", "#FEA06C", "#FF9E6D", "#FF9D6F", "#FF9B71", "#FF9A72", "#FF9874", "#FF9577", "#FF9379", "#FF927A", "#FF907C", "#FF8F7D", "#FF8D7F", "#FE8B80", "#FE8A82", "#FE8884", "#FE8587", "#FE8388", "#FE818A", "#FE808B", "#FE7E8C", "#FE7C8E", "#FE7A8F", "#FD7991", "#FD7792", "#FD7395", "#FD7197", "#FD7098", "#FD6E99", "#FC6C9B", "#FC6A9C", "#FC689E", "#FC669F", "#FC64A0", "#FB60A3", "#FB5EA4", "#FB5CA5", "#FA5AA6", "#FA58A7", "#FA56A8", "#FA54A8", "#F952A8", "#F94EA7", "#F84CA6", "#F84AA5", "#F849A3", "#F747A1", "#F7469E", "#F7449C", "#F64399", "#F64296", "#F53F8F", "#F53E8C", "#F43D88", "#F43C85", "#F33A82", "#F3397E", "#F2387B", "#F13778", "#F13674", "#F0336D", "#EF326A", "#EE3167", "#EE2F63", "#ED2E60", "#EC2D5D", "#EB2C59", "#EB2A56", "#EA2953", "#E8264C", "#E82549", "#E72445", "#E62242", "#E5213F", "#E4203B", "#E31E38", "#E21D35", "#E21B31", "#E0182A", "#DF1626", "#DE1423", "#DD131F", "#DC111B", "#DB0F17", "#DA0C12", "#D90A0C", "#D70500"))(256)
    # 6. black and white
    if (col==6 | col=="black and white" | col=="b&w") colr.o <- colorRampPalette(c("#FFFFFF", "#404040", "#000000"))(256)
  }
  # set up plotting space
  mar.bot=0; mar.left=0; mar.top=0; mar.right=0 # default margins
  if (title==T) mar.top=3
  if (axes==T) mar.bot=3.1; mar.left=2.1; mar.right=1.1
  if (plot == "together" | plot == "tog" | plot == "t"){
    arr <- n2mfrow(do.call(sum, lapply(model, length)), asp=plot.asp) # calculate best plot arrangement
    par(mfrow=arr, mar=c(mar.bot, mar.left, mar.top, mar.right))
  } else if (plot == "separate" | plot == "separately" | plot == "sep" | plot == "s"){
    par(mfrow=c(1,1),  mar=c(mar.bot, mar.left, mar.top, mar.right))
  }

  if (to.scale==T) {
    # get mins and maxes for all rasters
    minmaxes <- do.call(c, lapply(model, function(x){
      y=lapply(x, minmax)
      unlist(y)
    }))
    min <- min(minmaxes)
    max <- max(minmaxes)
  }

  # for every scenario...
  invisible(lapply(names(model), function(scenario){
    # and for every taxon in that scenario...
    lapply(model[[scenario]], function(taxon){
      if (axes==F & to.scale==F) image(taxon, col=colr.o, ann=F, axes=F)
      if (axes==F & to.scale==T) image(taxon, col=colr.o, ann=F, axes=F, zlim=c(min, max))
      if (axes==T & to.scale==F) image(taxon, col=colr.o, ann=F, axes=T)
      if (axes==T & to.scale==T) image(taxon, col=colr.o, ann=F, axes=T, zlim=c(min, max))
      # plot title if specified
      if (title == TRUE) title(main=paste0(scenario, ", ", names(taxon)), cex.main=title.cex)
      return(NULL)
    })
  }))
}
