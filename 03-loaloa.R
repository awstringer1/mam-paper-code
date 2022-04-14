### Loa Loa example ###

## Load Packages ----
# install.packages("geostatsp", repos="http://R-Forge.R-project.org")
library(geostatsp)
library(mgcv)
library(ggplot2)

# The mam package
# install.packages('devtools')
# devtools::install_github('awstringer1/mam')
library(mam)

## Setup ----

# Constants #
# Resolution of spatial grid for plotting
ny <- 150 # Cells in the y direction
nx <- ny*2 # Map is about twice as wide as it is tall
br <- c(0,.05,.1,.15,.2,.25,.3,.4,.5,.6,1) # Breaks for plotting probabilities
bU <- c(-4.5,-2,-1,-.5,-.25,0,.25,.5,1,2,4) # Breaks for plotting spatial random effects
MAPINSET <- .05
MAPHEIGHT = 3.5
MAPWIDTH <- 2*MAPHEIGHT
PLOTWIDTH <- PLOTHEIGHT <- 7
GGPLOTTEXTSIZE <- 26
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Paths
globalpath <- tempdir()
datapath <- file.path(globalpath,"data")
plotpath <- file.path(globalpath,"figures")
if (!dir.exists(datapath)) dir.create(datapath)
if (!dir.exists(plotpath)) dir.create(plotpath)

# Data
data(loaloa,package = "geostatsp")
dat <- loaloa@data
rng <- 4.22e04 # Fixed range param from Brown (2011)
# Expand to a dataframe with N rows with y 1's per village
set.seed(748023) # Order is random, but inferences are invariant to this
gen_group <- function(N,y,id) {
  ynew <- rep(0,N)
  idx <- sample.int(N,size=y)
  ynew[idx] <- 1
  data.frame(y=ynew,id=id)
}
datlist <- vector(mode='list',length=length(unique(dat$villageID)))
for (i in 1:length(datlist)) datlist[[i]] <- gen_group(dat[i, ]$N,dat[i, ]$y,dat[i, ]$villageID)
datexpand <- Reduce(rbind,datlist)
stopifnot(nrow(datexpand) == sum(dat$N))
stopifnot(sum(datexpand$y) == sum(dat$y))
stopifnot(all(unique(datexpand$id) == dat$villageID))

# Covariates
# From Brown (2011)
elevationLoa <- elevationLoa - 750
rcl <- rbind(c(9, 8), c(5, 2), c(11, 2), c(12, 14), c(13, 14))
ltLoaRe <- reclassify(ltLoa, rcl)
levels(ltLoaRe) = levels(ltLoa)

# Model data
gamdat <- within(loaloa@data,{
  failures = N-y
  long = loaloa@coords[ ,'long']
  lat = loaloa@coords[ ,'lat']
})
gamdat$elevation <- extract(elevationLoa,loaloa)
gamdat$evi <- extract(eviLoa,loaloa)
gamdat$land <- extract(ltLoaRe,loaloa)
# Centre and scale covariates
zscore <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
# Keep mean/sd for later
covmoments <- with(gamdat,list(
  elevation = c(mean(elevation),sd(elevation)),
  evi = c(mean(evi),sd(evi))
))
gamdat$elevation <- zscore(gamdat$elevation)
gamdat$evi <- zscore(gamdat$evi)
gamdat$land <- factor(gamdat$land)

# "Expanded" Bernoulli data for MAM
gamdatexpand <- merge(datexpand,gamdat,by.x = 'id',by.y = 'villageID')
gamdatexpand$y <- gamdatexpand$y.x
gamdatexpand$y.x <- gamdatexpand$y.y <- gamdatexpand$failures <- gamdatexpand$N <- NULL



# Spatial background data, for plotting results
cameroonBorderLL = raster::getData("GADM", country=c('CMR'), level=2)
nigeriaBorderLL = raster::getData("GADM", country=c('NGA'), level=2)
cameroonBorder = spTransform(cameroonBorderLL, projection(loaloa))
nigeriaBorder = spTransform(nigeriaBorderLL, projection(loaloa))
cameroonBorderouter <- rgeos::gUnaryUnion(cameroonBorder)
nigeriaBorderouter <- rgeos::gUnaryUnion(nigeriaBorder)


fullborder <- raster::bind(cameroonBorder,nigeriaBorder)
fullborderouter <- raster::bind(cameroonBorderouter,nigeriaBorderouter)

fullborder <- crop(fullborder,loaloa)
fullborderouter <- crop(fullborderouter,loaloa)

plot_loaloa <- function(plotraster,breaks) {
  predcols <- mapmisc::colourScale(
    plotraster,
    breaks = breaks,
    style = "fixed",
    col = "Spectral",
    rev = TRUE
  )

  plotraster <- mask(plotraster,fullborderouter)

  mapmisc::map.new(loaloa,legendRight = TRUE)
  plot(plotraster,
       col = predcols$col,
       breaks = predcols$breaks,
       legend=FALSE, add=TRUE)
  points(loaloa,pch = 4)
  plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
  plot(fullborderouter,add = TRUE)
  mapmisc::legendBreaks('right', predcols, cex=1, bty='o',inset = MAPINSET)
}

# Raster and points for plotting
rr <- raster(loaloa,nrows = ny,ncols = nx)
loagridpoints <- rasterToPoints(rr,spatial = TRUE)


## Models ----
gamformula <- cbind(y,failures) ~ s(long,lat,bs='gp',m = c(-3,rng)) + s(evi) + s(elevation)

# Fit the quasibinomial, for dispersion parameter
gammodquasibinomial <- gam(
  gamformula,
  family = quasibinomial(),
  method = 'ML',
  data = gamdat
)
# summary(gammodquasibinomial)
gamquasidisperse <- gammodquasibinomial$scale

# Binomial GAM #
gammodbinomial <- gam(
  gamformula,
  family = binomial(),
  method = 'ML',
  data = gamdat
)

# Bernoulli MAM #
# Need to provide some prediction data (for now)
preddat <- data.frame(
  long = loagridpoints@coords[ ,'x'],
  lat = loagridpoints@coords[ ,'y'],
  elevation = zscore(extract(elevationLoa,loagridpoints)),
  evi = zscore(extract(eviLoa,loagridpoints))
)
narows <- apply(preddat,1,function(x) any(is.na(x)))
preddat <- preddat[!narows, ]
loagridpoints <- loagridpoints[!narows, ]
# Fit the model
mamsmoothlist <- interpret.gam(gamformula)$smooth.spec

themam <- mam(
  smooth = mamsmoothlist,
  re = y ~ (1|id),
  fe = NULL,
  dat = gamdatexpand,
  margdat = gamdatexpand,
  preddat = preddat,
  control = mam_control(
    k = 5,
    retcond = FALSE,
    method = "BFGS",
    varmethod = 1,
    verbose = TRUE,
    centered = TRUE
  )
)


## Write Results ----
## Estimated dispersion parameter ##
readr::write_csv(data.frame(dispersion = gamquasidisperse),file = file.path(plotpath,"dispersion.csv"))

## Variance component estimate ##
logsigmaest <- themam$variance$optresults$solution[1]
logsigmasd <- sqrt(solve(themam$variance$optresults$hessian)[1,1])
mor <- exp(sqrt(2)*exp(logsigmaest)*qnorm(0.75))
morlower <- exp(sqrt(2)*exp(logsigmaest - 2*logsigmasd)*qnorm(0.75))
morupper <- exp(sqrt(2)*exp(logsigmaest + 2*logsigmasd)*qnorm(0.75))

readr::write_csv(
  data.frame(
    varcomp = exp(c(logsigmaest-2*logsigmasd,logsigmaest,logsigmaest+2*logsigmasd)),
    mor = c(morlower,mor,morupper)
  ),
  file = file.path(plotpath,'variance-components.csv')
)

## GAM plots ##
preddat <- data.frame(
  long = loagridpoints@coords[ ,'x'],
  lat = loagridpoints@coords[ ,'y'],
  elevation = zscore(extract(elevationLoa,loagridpoints)),
  evi = zscore(extract(eviLoa,loagridpoints)),
  land = factor(extract(ltLoaRe,loagridpoints))
)
# Spatial effects #
binomspatpred <- predict(gammodbinomial,newdata = preddat,type = 'terms',se.fit=TRUE)
binomspatpred_fit <- binomspatpred$fit[ ,'s(long,lat)']
binomspatpred_lower <- binomspatpred$fit[ ,'s(long,lat)'] - 2*binomspatpred$se.fit[ ,'s(long,lat)']
binomspatpred_upper <- binomspatpred$fit[ ,'s(long,lat)'] + 2*binomspatpred$se.fit[ ,'s(long,lat)']
# Estimate
smoothpreds <- data.frame(pred = binomspatpred_fit)
loagridpointsframe <- SpatialPointsDataFrame(coords = loagridpoints@coords,data = smoothpreds,proj4string = loagridpoints@proj4string,bbox = loagridpoints@bbox)
loagridpixels <- SpatialPixelsDataFrame(loagridpointsframe,data = loagridpointsframe@data)
loagridraster <- raster(loagridpixels)
pdf(file.path(plotpath,"gam-spatialeffect-est.pdf"),width=MAPWIDTH,height=MAPHEIGHT)
plot_loaloa(loagridraster,breaks = bU)
dev.off()

# Lower
smoothpreds <- data.frame(pred = binomspatpred_lower)
loagridpointsframe <- SpatialPointsDataFrame(coords = loagridpoints@coords,data = smoothpreds,proj4string = loagridpoints@proj4string,bbox = loagridpoints@bbox)
loagridpixels <- SpatialPixelsDataFrame(loagridpointsframe,data = loagridpointsframe@data)
loagridraster <- raster(loagridpixels)
pdf(file.path(plotpath,"gam-spatialeffect-lower.pdf"),width=MAPWIDTH,height=MAPHEIGHT)
plot_loaloa(loagridraster,breaks = bU)
dev.off()

# Upper
smoothpreds <- data.frame(pred = binomspatpred_upper)
loagridpointsframe <- SpatialPointsDataFrame(coords = loagridpoints@coords,data = smoothpreds,proj4string = loagridpoints@proj4string,bbox = loagridpoints@bbox)
loagridpixels <- SpatialPixelsDataFrame(loagridpointsframe,data = loagridpointsframe@data)
loagridraster <- raster(loagridpixels)
pdf(file.path(plotpath,"gam-spatialeffect-upper.pdf"),width=MAPWIDTH,height=MAPHEIGHT)
plot_loaloa(loagridraster,breaks = bU)
dev.off()

# Smooth Effects #
evipreddat <- data.frame(evi=with(gamdat,seq(min(evi),max(evi),length.out=1e03)),elevation=0,long=0,lat=0)
evipred <- predict(gammodbinomial,newdata=evipreddat,se.fit=TRUE,type='terms')
evi_gg_data <- data.frame(
  evi = evipreddat$evi,
  pred = evipred$fit[ ,'s(evi)'],
  lower = evipred$fit[ ,'s(evi)'] - 2*evipred$se.fit[ ,'s(evi)'],
  upper = evipred$fit[ ,'s(evi)'] + 2*evipred$se.fit[ ,'s(evi)']
)
myylim = c(-2,2)
evi_gam_plot <- ggplot(data=evi_gg_data,aes(x=evi,y=pred))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="gray",alpha=0.65)+
  geom_line()+
  ylim(myylim)+
  xlab("Evi")+
  ylab(expression(f^M~(evi))) +
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggsave(filename = file.path(plotpath,"evi-gam.pdf"),plot = evi_gam_plot,width=PLOTWIDTH,height=PLOTHEIGHT)


elevpreddat <- data.frame(elevation=with(gamdat,seq(min(elevation),max(elevation),length.out=1e03)),evi=0,long=0,lat=0)
elevpred <- predict(gammodbinomial,newdata=elevpreddat,se.fit=TRUE,type='terms')
elev_gg_data <- data.frame(
  elevation = elevpreddat$elevation,
  pred = elevpred$fit[ ,'s(elevation)'],
  lower = elevpred$fit[ ,'s(elevation)'] - 2*elevpred$se.fit[ ,'s(elevation)'],
  upper = elevpred$fit[ ,'s(elevation)'] + 2*elevpred$se.fit[ ,'s(elevation)']
)
myylim = c(-5,2.5)
elev_gam_plot <- ggplot(data=elev_gg_data,aes(x=elevation,y=pred))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="gray",alpha=0.65)+
  geom_line()+
  ylim(myylim)+
  xlab("Elevation")+
  ylab(expression(f^M~(elevation))) +
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggsave(filename = file.path(plotpath,"elev-gam.pdf"),plot = elev_gam_plot,width=PLOTWIDTH,height=PLOTHEIGHT)

## MAM ##

## Estimates of individual smooths ##
spatidx <- c(4:(4+32-1))
eviidx <- c(2,(max(spatidx)+1):(max(spatidx)+8))
elevationidx <- c(3,(max(eviidx)+1):(max(eviidx)+8))
# Unfortunately, need to recreate the basis stuff, having a different prediction dataset for each smooth...
SS <- lapply(lapply(mamsmoothlist,mgcv::smoothCon,data=gamdatexpand,absorb.cons = TRUE),'[[',1)
numsmooth <- length(smooth) # Number of smooth terms
EE <- lapply(lapply(lapply(SS,'[[','S'),'[[',1),eigen)
p <- sapply(lapply(EE,'[[','vectors'),ncol)
r <- sapply(lapply(EE,'[[','values'),function(x) sum(x>.Machine$double.eps))
m <- p-r
URlist <- mapply(function(x,y) x[ ,1:y],lapply(EE,'[[','vectors'),r,SIMPLIFY = FALSE)
UFlist <- mapply(
  function(x,y,z) {
    if (y<z) return(x[ ,(1+y):z])
    newsparsemat(z,z)
  },lapply(EE,'[[','vectors'),r,p,SIMPLIFY = FALSE)
URlist <- lapply(URlist,cbind) # Ensure they stay matrices
UFlist <- lapply(UFlist,cbind) # Ensure they stay matrices
UR <- Reduce(bdiag,URlist)
UF <- Reduce(bdiag,UFlist)
# if m=1 UF gets coerced to numeric
if (!is.matrix(UF)) UF <- cbind(UF)
Dpilist <- lapply(lapply(EE,'[[','values'),function(x) x[x>.Machine$double.eps])
Dpilist <- lapply(Dpilist,function(x) Diagonal(length(x),1/sqrt(x)))
Dpi <- bdiag(Dpilist)
Xlist <- lapply(SS,'[[','X')
X <- Reduce(cbind,Xlist)
Xr <- as.matrix(X %*% UR %*% Dpi)
Xf <- as.matrix(X %*% UF)
dups <- !duplicated(t(Xf)) & apply(Xf,2,function(x) !all(x==0)) # Remove the duplicated intercepts
if (length(dups) > 1) Xf <- Xf[ ,which(dups)]

eviSSidx <- 2 # Which smooth term it is
evipred <- data.frame(evi=with(gamdat,seq(min(evi),max(evi),length.out=1e03)),elevation=0,long=0,lat=0)
XXpred <- PredictMat(SS[[eviSSidx]],data=evipred)
Xrpred <- as.matrix(XXpred %*% URlist[[eviSSidx]] %*% Dpilist[[eviSSidx]])
Xfpred <- as.matrix(XXpred %*% UFlist[[eviSSidx]])
# Take out intercept (if no sum to zero constraint)
interceptloc <- which(apply(Xfpred,2,function(x) all(x==1)))
if (length(interceptloc)>0) Xfpred <- Xfpred[ ,-interceptloc]
Xpred <- cbind(Xfpred,Xrpred)
eviest_marg <- as.numeric(Xpred %*% themam$mam$coefsmooth[eviidx])
eviest_marg_se <- sqrt(with(themam$variance,colSums((mamvarfactor_cond[ ,eviidx] %*% t(Xpred))^2) + colSums((mamvarfactor_marg[ ,eviidx] %*% t(Xpred))^2)))

evi_gg_data <- data.frame(
  evi = evipred$evi,
  pred = eviest_marg,
  lower = eviest_marg - 2*eviest_marg_se,
  upper = eviest_marg + 2*eviest_marg_se
)
myylim = c(-2,2)
evi_mam_plot <- ggplot(data=evi_gg_data,aes(x=evi,y=pred))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="gray",alpha=0.65)+
  geom_line()+
  ylim(myylim)+
  xlab("Evi")+
  ylab(expression(f^M~(evi))) +
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggsave(filename = file.path(plotpath,"evi-mam.pdf"),plot = evi_mam_plot,width=PLOTWIDTH,height=PLOTHEIGHT)


# Elevation
elevSSidx <- 3 # Which smooth term it is
elevpred <- data.frame(elevation = with(gamdat,seq(min(elevation),max(elevation),length.out=1e03)),evi=0,long=0,lat=0)
XXpred <- PredictMat(SS[[elevSSidx]],data=elevpred)
Xrpred <- as.matrix(XXpred %*% URlist[[elevSSidx]] %*% Dpilist[[elevSSidx]])
Xfpred <- as.matrix(XXpred %*% UFlist[[elevSSidx]])
# Take out intercept
interceptloc <- which(apply(Xfpred,2,function(x) all(x==1)))
if (length(interceptloc)>0) Xfpred <- Xfpred[ ,-interceptloc]
Xpred <- cbind(Xfpred,Xrpred)
elevest_marg <- as.numeric(Xpred %*% themam$mam$coefsmooth[elevationidx])
elevest_marg_se <- sqrt(with(themam$variance,colSums((mamvarfactor_cond[ ,elevationidx] %*% t(Xpred))^2) + colSums((mamvarfactor_marg[ ,elevationidx] %*% t(Xpred))^2)))
myylim <- c(-5,2.5)

elev_gg_data <- data.frame(
  elevation = elevpred$elevation,
  pred = elevest_marg,
  lower = elevest_marg - 2*elevest_marg_se,
  upper = elevest_marg + 2*elevest_marg_se
)
elev_mam_plot <- ggplot(data=elev_gg_data,aes(x=elevation,y=pred))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="gray",alpha=0.65)+
  geom_line()+
  ylim(myylim)+
  xlab("Elevation")+
  ylab(expression(f^M~(elevation))) +
  theme(text = element_text(size = GGPLOTTEXTSIZE))
ggsave(filename = file.path(plotpath,"elev-mam.pdf"),plot = elev_mam_plot,width=PLOTWIDTH,height=PLOTHEIGHT)


## Spatial effect ##

spatSSidx <- 1
spatpred <- preddat
XXpred <- PredictMat(SS[[spatSSidx]],data=spatpred)
Xrpred <- as.matrix(XXpred %*% URlist[[spatSSidx]] %*% Dpilist[[spatSSidx]])
Xpred <- Xrpred # UF = 0 here after constraint "absorbed"
spatest_marg <- as.numeric(Xpred %*% themam$mam$coefsmooth[spatidx])
spatest_marg_se <- sqrt(with(themam$variance,colSums((mamvarfactor_cond[ ,spatidx] %*% t(Xpred))^2) + colSums((mamvarfactor_marg[ ,spatidx] %*% t(Xpred))^2)))
# Estimates
smoothpreds <- data.frame(pred = spatest_marg)
loagridpointsframe <- SpatialPointsDataFrame(coords = loagridpoints@coords,data = smoothpreds,proj4string = loagridpoints@proj4string,bbox = loagridpoints@bbox)
loagridpixels <- SpatialPixelsDataFrame(loagridpointsframe,data = loagridpointsframe@data)
loagridraster <- raster(loagridpixels)
pdf(file.path(plotpath,"marg-spatialfield-est.pdf"),width=MAPWIDTH,height=MAPHEIGHT)
plot_loaloa(loagridraster,breaks = bU)
dev.off()

# Lower CI
smoothpreds <- data.frame(pred = spatest_marg - 2*spatest_marg_se)
loagridpointsframe <- SpatialPointsDataFrame(coords = loagridpoints@coords,data = smoothpreds,proj4string = loagridpoints@proj4string,bbox = loagridpoints@bbox)
loagridpixels <- SpatialPixelsDataFrame(loagridpointsframe,data = loagridpointsframe@data)
loagridraster <- raster(loagridpixels)
pdf(file.path(plotpath,"marg-spatialfield-lower.pdf"),width=MAPWIDTH,height=MAPHEIGHT)
plot_loaloa(loagridraster,breaks = bU)
dev.off()

# Upper CI
smoothpreds <- data.frame(pred = spatest_marg + 2*spatest_marg_se)
loagridpointsframe <- SpatialPointsDataFrame(coords = loagridpoints@coords,data = smoothpreds,proj4string = loagridpoints@proj4string,bbox = loagridpoints@bbox)
loagridpixels <- SpatialPixelsDataFrame(loagridpointsframe,data = loagridpointsframe@data)
loagridraster <- raster(loagridpixels)
pdf(file.path(plotpath,"marg-spatialfield-upper.pdf"),width=MAPWIDTH,height=MAPHEIGHT)
plot_loaloa(loagridraster,breaks = bU)
dev.off()

## The random intercepts ##

Uest <- as.numeric(themam$conditional$predU)
summary(Uest)
pdf(file.path(plotpath,"U-histogram.pdf"),width=7,height=7)
hist(Uest,freq=FALSE,main="",xlab="Predicted Random Intercepts",ylab="Density",breaks=30)
dev.off()

predUpoints <- SpatialPointsDataFrame(coords = loaloa@coords,data = data.frame(predU = Uest),bbox = loaloa@bbox,proj4string = loaloa@proj4string)

# Size by village size, colour by signed value on same scale as maps
sz <- loaloa$N/max(loaloa$N) * 2
cl <- mapmisc::colourScale(Uest,breaks = bU,style='fixed',col='Spectral',rev=TRUE)
plotord <- order(abs(predUpoints$predU),decreasing = FALSE)

pdf(file.path(plotpath,"U-map.pdf"),width=MAPWIDTH,height=MAPHEIGHT)
mapmisc::map.new(loaloa,legendRight = TRUE)
plot(fullborder,add = TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
plot(fullborderouter,add = TRUE)
points(predUpoints[plotord, ],col = cl$plot[plotord],cex = 3*sz[plotord],pch=20)
mapmisc::legendBreaks('right', cl, cex=1, bty='o',inset = MAPINSET)
dev.off()
cat("Done. Output saved at",globalpath,"\n")
