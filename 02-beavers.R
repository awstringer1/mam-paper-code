### Beavers example ###

## Load libraries ----
# The mam package
# install.packages('devtools')
# devtools::install_github('awstringer1/mam')
library(mam)
library(mgcv)
library(ggplot2)
GGPLOTTEXTSIZE <- 26
PLOTWIDTH <- PLOTHEIGHT <- 7

## Set paths ----
globalpath <- tempdir()
plotpath <- file.path(globalpath,"figures")
if (!dir.exists(plotpath)) dir.create(plotpath)

## Data ----
# Data included in mam package
data(beavers)
# Create Sex indicator variables. Change 0's to small nonzero value
# due to internal behaviour of lme4 and Matrix packages which convert data zeroes
# to structural zeroes. See https://github.com/lme4/lme4/issues/671
beavers$sex0 <- as.numeric(beavers$sex=="M");beavers$sex0[beavers$sex0==0] <- 0.0000001
beavers$sex1 <- as.numeric(beavers$sex=="F");beavers$sex1[beavers$sex1==0] <- 0.0000001

## Model ----
# Create the prediction data
beaverspred <- data.frame(time=rep(seq(min(beavers$time),max(beavers$time),length=100),2),
                           timeM=c(seq(min(beavers$time),max(beavers$time),length=100),rep(0,100)),
                           timeF=c(rep(0,100),seq(min(beavers$time),max(beavers$time),length=100)),
                           year=factor(levels(beavers$year)[1],levels=levels(beavers$year)),
                           sex=factor(rep(c("M","F"),each=100),levels=levels(beavers$sex)))

# Fit the model
bv.mam <- mam(
  smooth = list(s(timeF),s(timeM)),
  re = y ~ (sex0-1|id)+(sex1-1|id),
  fe=~year,#+sex
  dat = beavers,
  margdat = beavers,
  preddat = beaverspred,
  control = mam_control(
    centered=TRUE,
    method = 'trust',
    varmethod = 1,
    verbose = TRUE,
    retcond = TRUE
  )
)

## Return Output ----
# Standard deviation estimates
opt <- bv.mam$variance$optresults
est_logsd <- opt$solution[1:2] # log(sigma)
est_logsd_se <- sqrt(diag(solve(opt$hessian)))[1:2]

write.csv(
  data.frame(
    sigma1 = exp(c(est_logsd[1]-2*est_logsd_se[1],est_logsd[1],est_logsd[1]+2*est_logsd_se[1])),
    sigma2 = exp(c(est_logsd[2]-2*est_logsd_se[2],est_logsd[2],est_logsd[2]+2*est_logsd_se[2]))
  ),
  file = file.path(plotpath,"variance-components.csv")
)

# Curves
bv.mam.uci <- bv.mam$mam$fitted+1.96*bv.mam$mam$fitted_se
bv.mam.lci <- bv.mam$mam$fitted-1.96*bv.mam$mam$fitted_se
myrange=c(-2,3.2)
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dfplotF_mam <- data.frame(time=beaverspred$time[101:200],
                          fitted=bv.mam$mam$fitted[101:200],
                          uci=bv.mam.uci[101:200],
                          lci=bv.mam.lci[101:200])
ggF <- ggplot(data=dfplotF_mam,aes(x=time,y=fitted))+
  geom_ribbon(aes(ymin=lci,ymax=uci),fill="gray",alpha=0.65)+
  geom_line()+
  ylim(myrange)+
  xlab("Day of Year")+
  ylab(expression(f^M~(time)))+
  ggtitle("Female") +
  theme(text = element_text(size = GGPLOTTEXTSIZE))
dfplotM_mam <- data.frame(time=beaverspred$time[1:100],
                          fitted=bv.mam$mam$fitted[1:100],
                          uci=bv.mam.uci[1:100],
                          lci=bv.mam.lci[1:100])
ggM <- ggplot(data=dfplotM_mam,aes(x=time,y=fitted))+
  geom_ribbon(aes(ymin=lci,ymax=uci),fill="gray",alpha=0.65)+
  geom_line()+
  ylim(myrange)+
  xlab("Day of Year")+
  ylab(expression(f^M~(time)))+
  ggtitle("Male") +
  theme(text = element_text(size = GGPLOTTEXTSIZE))

ggsave(filename = file.path(plotpath,'beaver-time-female.pdf'),plot = ggF,width=PLOTWIDTH,height=PLOTHEIGHT)
ggsave(filename = file.path(plotpath,'beaver-time-male.pdf'),plot = ggM,width=PLOTWIDTH,height=PLOTHEIGHT)
cat("Done. Output saved at",globalpath,"\n")

