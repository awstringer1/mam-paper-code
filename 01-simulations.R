## Simulations ##
set.seed(456321)

## Load libraries ----
# The mam package
# install.packages('devtools')
# devtools::install_github('awstringer1/mam')
library(mam)
library(mgcv)
library(gamm4)
library(tidyverse)
theme_set(theme_bw())
theme_replace(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
GGPLOTTEXTSIZE <- 26
PLOTWIDTH <- PLOTHEIGHT <- 7

library(parallel)
options(mc.cores = parallel::detectCores())

## Set paths ----
globalpath <- tempdir()
plotpath <- file.path(globalpath,"figures")
if (!dir.exists(plotpath)) dir.create(plotpath)
savepath <- file.path(globalpath,"data")
if (!dir.exists(savepath)) dir.create(savepath)

## Set global parameters ----
f1temp <- function(x){(1.5*sin(pi*x)-2*x)}
f1 <- function(x) f1temp(x)-f1temp(0)
f2temp <- function(x){5*(dnorm(2*(x)))}
f2 <- function(x) f2temp(x)-f2temp(0)
beta3 <- 0

sigma0 <- 2 #2 # random intercept sd
sigma3 <- 1 # random slope sd
rho <- 0.5 # correlation between random intercept and slope

get_data_mam <- function(SSmat,x1,x2,x3,K,Nk){
  if(SSmat[2,2]==0){
    V <- as.data.frame(rnorm(K,0,sqrt(SSmat[1,1])))
    colnames(V) <- "intercepts"
    V$id <- 1:K
  }else{
    V <- mvtnorm::rmvnorm(K,sigma = SSmat) # Intercept and slope
    V <- as.data.frame(V)
    colnames(V) <- c("intercepts","slopes")
    V$id <- 1:K
  }
  dat <- data.frame(id = rep(1:K,each=Nk))
  dat <- merge(dat,V,by="id",all.x=TRUE)
  dat$x1 <- x1
  dat$x2 <- x2
  dat$x3 <- x3
  dat$fx1 <- f1(x1)
  dat$fx2 <- f2(x2)
  
  ## approximation from Hedeker
  if(SSmat[2,2]==0){
    Ui <- dat[,c("intercepts")]
    vars <- SSmat[1,1]+ ##
      +(pi^2)/3 # (pi^2)/3 is approximate due to logistic distn
  }else{
    Z <- cbind(1,x3) # ranefs with design matrix
    Ui <- apply(Z*dat[,c("intercepts","slopes")],1,sum)
    vars <- colSums(((chol(SSmat)%*%t(Z)))^2)+ ## faster than #diag(Z%*%SS%*%t(Z))
      +(pi^2)/3 # (pi^2)/3 is approximate due to logistic distn
  }
  estar <- rlogis(K*Nk) # approx: rnorm(K*Nk,0,sqrt((pi^2)/3))
  dat <- within(dat,{y <- as.numeric(pnorm(estar+Ui,0,sqrt(vars)) < mam::ilogit(fx1 + fx2 +               ## non linear fixed effects
                                                                            beta3*x3 +                ## linear fixed effects
                                                                            0))  })  ## random effects on other side
  return(dat)
}

do_simulation_multiple <- function(K,Nk,sigma0,sigma3,rho,iter=0,verbose=TRUE) {
  # note: sigma3 = 0 ==> only intercept
  slopes <- sigma3 > 0
  # note: rho = 0 ==> independent slope/intercept
  indepU <- rho == 0
  # Form the true covariance matrix
  trueSigma <- matrix(c(sigma0^2,
                        rho*sigma0*sigma3,rho*sigma0*sigma3,
                        sigma3^2),nrow=2,ncol=2)
  ## generate data
  x1 <- round(runif(K*Nk,-1,1),3)
  x3 <- round(runif(K*Nk,-1,1),3)
  x2 <- round(runif(K*Nk,-1,1),3)
  dat <- get_data_mam(trueSigma,x1,x2,x3,K,Nk)
  gridlen=100
  dat2pred <- data.frame(x1 = c(seq(-1,1,length=gridlen),rep(0,2*gridlen)),
                         x2 = c(rep(0,gridlen),seq(-1,1,length=gridlen),rep(0,gridlen)),
                         x3 = c(rep(0,2*gridlen),seq(-1,1,length=gridlen)))
  dat2pred$fx1 <- f1(dat2pred$x1)
  dat2pred$fx2 <- f2(dat2pred$x2)
  
  ## fit models
  # GAM #
  
  tm <- Sys.time()
  gam1 <- gam(y~s(x1)+s(x2)+x3,data=dat,family=binomial(),method="REML")
  gamfit <- predict(gam1,newdata=dat2pred,se.fit=TRUE)
  gamfitted <- gamfit$fit
  gamfitted_se <- gamfit$se.fit
  if (verbose) cat("Fitting gam took",as.numeric(difftime(Sys.time(),tm,units='secs')),"seconds.\n")
  
  # MAM (proposed approach)
  tm <- Sys.time()
  if(!slopes){
    themam <- mam(smooth = list(s(x1),s(x2)),
                  re = y ~ (1|id),
                  fe = ~ x3,
                  dat = dat,
                  margdat = dat,
                  preddat = dat2pred,
                  control = mam_control(
                    # method = 'trust',
                    method = 'BFGS',
                    varmethod = 1,
                    verbose = verbose,
                    retcond = TRUE))
  }else if(indepU){
    themam <- mam(smooth = list(s(x1),s(x2)),
                  re = y ~ (1|id)+(x3-1|id),
                  fe = ~ x3,
                  dat = dat,
                  margdat = dat,
                  preddat = dat2pred ,
                  control = mam_control(
                    # method = 'trust',
                    method = 'BFGS',
                    varmethod = 1,
                    verbose = verbose,
                    retcond = TRUE))
  }else{
    themam <- NULL
    try( ## error handling due to optimizer
      themam <- mam::mam(smooth = list(s(x1),s(x2)),
                         re = y ~ (1+x3|id),
                         fe = ~ x3,
                         dat = dat,
                         margdat = dat,
                         preddat = dat2pred ,
                         control = mam_control(
                           # method = 'trust',
                           method = 'BFGS',
                           varmethod = 1,
                           verbose = verbose,
                           retcond = TRUE))
    )
    if(is.null(themam)){
      themam <- mam::mam(smooth = list(s(x1),s(x2)),
                         re = y ~ (1+x3|id),
                         fe = ~ x3,
                         dat = dat,
                         margdat = dat,
                         preddat = dat2pred ,
                         control = mam_control(
                           method = 'BFGS',
                           varmethod = 1,
                           verbose = verbose,
                           retcond = TRUE))
    }
    
  }
  if (verbose) cat("Fitting mam took",as.numeric(difftime(Sys.time(),tm,units='secs')),"seconds.\n")

  
  # (v) double-penalized mam (i.e. penalty at the conditional and marginal stage)
  tm <- Sys.time()
  dat$gpix <- logit(themam$marginal$prob) #old: #mam$marginalmean$fitted
  mam2 <- gam(gpix~s(x1)+s(x2)+x3,data=dat)
  mam2fit <- predict(mam2,newdata=dat2pred,se.fit=TRUE)
  mam2fitted <- mam2fit$fit
  mam2fitted_se <- mam2fit$se.fit ## these are naive
  if (verbose) cat("Fitting double mam took",as.numeric(difftime(Sys.time(),tm,units='secs')),"seconds.\n")


  # (vi) GAMM (conditional)
  tm <- Sys.time()
  if(slopes==FALSE){
    gamm <- gamm4(y ~ s(x1)+s(x2)+x3,random=~(1|id), data = dat,family = binomial())
  }else if(indepU==TRUE){
    gamm <- gamm4(y ~ s(x1)+s(x2)+x3,random=~(1|id)+(x3-1|id), data = dat,family = binomial())
  }else{
    gamm <- gamm4(y ~ s(x1)+s(x2)+x3,random=~(1+x3|id), data = dat,family = binomial())
  }
  gammfit <- predict(gamm$gam,newdata=dat2pred,se.fit=TRUE)
  gammfitted <- gammfit$fit
  gammfitted_se <- gammfit$se.fit
  gammpred <- ranef(gamm$mer)$id
  if (verbose) cat("Fitting gamm took",as.numeric(difftime(Sys.time(),tm,units='secs')),"seconds.\n")

  
  
  ## get fitted values
  tm <- Sys.time()
  fitteddat <- as_tibble(dat2pred) %>%
    mutate(
      ## Simulation parameters
      K=K,Nk=Nk,sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
      ## gam
      gammargfitted = as.numeric(gamfitted),
      gammargse = as.numeric(gamfitted_se),
      gammarglower = gammargfitted - 1.96*gammargse,
      gammargupper = gammargfitted + 1.96*gammargse,
      gammargcovrx1 = fx1 >= gammarglower & fx1 <= gammargupper,
      gammargcovrx2 = fx2 >= gammarglower & fx2 <= gammargupper,
      gammargbiasx1 = gammargfitted - fx1,
      gammargbiasx2 = gammargfitted - fx2,
      ## mam
      mamcondfitted = themam$conditional$fitted,
      mamcondse = themam$conditional$fitted_se,
      mamcondlower = mamcondfitted - 1.96*mamcondse,
      mamcondupper = mamcondfitted + 1.96*mamcondse,
      mamcondcovrx1 = fx1 >= mamcondlower & fx1 <= mamcondupper,
      mamcondcovrx2 = fx2 >= mamcondlower & fx2 <= mamcondupper,
      mamcondbiasx1 = mamcondfitted - fx1,
      mamcondbiasx2 = mamcondfitted - fx2,
      mammargfitted = as.numeric(themam$mam$fitted),
      mammargse = as.numeric(themam$mam$fitted_se),
      mammarglower = mammargfitted - 1.96*mammargse,
      mammargupper = mammargfitted + 1.96*mammargse,
      mammargcovrx1 = fx1 >= mammarglower & fx1 <= mammargupper,
      mammargcovrx2 = fx2 >= mammarglower & fx2 <= mammargupper,
      mammargbiasx1 = mammargfitted - fx1,
      mammargbiasx2 = mammargfitted - fx2,
      ## doubly penalized mam
      mam2margfitted = as.numeric(mam2fitted),
      mam2margse = as.numeric(mam2fitted_se),
      mam2marglower = mam2margfitted - 1.96*mam2margse,
      mam2margupper = mam2margfitted + 1.96*mam2margse,
      mam2margcovrx1 = fx1 >= mam2marglower & fx1 <= mam2margupper,
      mam2margcovrx2 = fx2 >= mam2marglower & fx2 <= mam2margupper,
      mam2margbiasx1 = mam2margfitted - fx1,
      mam2margbiasx2 = mam2margfitted - fx2,
      ## gamm
      gammcondfitted = gammfitted,
      gammcondse = gammfitted_se,
      gammcondlower = gammcondfitted - 1.96*gammcondse,
      gammcondupper = gammcondfitted + 1.96*gammcondse,
      gammcondcovrx1 = fx1 >= gammcondlower & fx1 <= gammcondupper,
      gammcondcovrx2 = fx2 >= gammcondlower & fx2 <= gammcondupper,
      gammcondbiasx1 = gammcondfitted - fx1,
      gammcondbiasx2 = gammcondfitted - fx2,
      iter=iter
    )
  
  ## get random effects
  b0 <- dat$intercepts[seq(1,nrow(dat),by=Nk)] ## cluster level random intercept
  b1 <- dat$slopes[seq(1,nrow(dat),by=Nk)] ## cluster level random slope
  bki <- dat$intercepts+dat$slopes*dat$x3 ## unit level random effect: b0+x3*b1
  
  ## mam
  mamMSEPb0 <- mean((themam$conditional$predU[,1]-b0)^2)
  if(slopes){
    mamMSEPb1 <- mean((themam$conditional$predU[,2]-b1)^2)
    mamMSEPbki <- mean( ( apply(themam$conditional$predU[dat$id,]*cbind(1,dat$x3),1,sum)-bki )^2)
  }else{
    mamMSEPb1  <- NA
    mamMSEPbki <- NA
  }
  
  ## gamm
  gammMSEPb0 <- mean((gammpred[,1]-b0)^2)
  if(slopes){
    gammMSEPb1 <- mean((gammpred[,2]-b1)^2)
    gammMSEPbki <- mean( c(( apply(gammpred[dat$id,]*cbind(1,dat$x3),1,sum)-bki )^2))
  }else{
    gammMSEPb1  <- NA
    gammMSEPbki <- NA
  }
  
  preddat <- data.frame(
    K=K,Nk=Nk,sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
    mamMSEPb0=mamMSEPb0,  
    mamMSEPb1=mamMSEPb1,  
    mamMSEPbki=mamMSEPbki,
    gammMSEPb0=gammMSEPb0,
    gammMSEPb1=gammMSEPb1,
    gammMSEPbki=gammMSEPbki  
  )
  
  
  ## estimate variance components
  if(slopes){
    Sig <- themam$conditional$theta
    mamtheta <- c(sqrt(Sig[1,1]),Sig[1,2]/sqrt(Sig[1,1]*Sig[2,2]),sqrt(Sig[2,2]))
    if(!indepU){
      gammV <- VarCorr(gamm$mer)$id
      gammtheta <- c(sqrt(gammV[1,1]),gammV[1,2]/sqrt(gammV[1,1]*gammV[2,2]),sqrt(gammV[2,2]))
    }else{
      gammtheta <- c(sqrt(VarCorr(gamm$mer)$id.1[1]),0,sqrt(VarCorr(gamm$mer)$id[1]))
    }
    names(mamtheta) <- names(gammtheta) <- c("sd-Intercept","cor","sd-x3")
  }else{
    mamtheta <-c(sqrt(themam$conditional$theta[1,1]),NA,NA)
    gammtheta <- c(sqrt(VarCorr(gamm$mer)$id[1,1]),NA,NA)
    names(mamtheta) <- names(gammtheta) <- c("sd-Intercept","cor","sd-x3")
  }
  
  varcomp <- data.frame(
    K=K,Nk=Nk,sigma0=sigma0,sigma3=sigma3,rho=rho,iter=iter,
    mam_sd_int=mamtheta['sd-Intercept'],
    mam_sd_slope=mamtheta['sd-x3'],
    mam_cor=mamtheta['cor'],
    gamm_sd_int=gammtheta['sd-Intercept'],
    gamm_sd_slope=gammtheta['sd-x3'],
    gamm_cor=gammtheta['cor']
  )
  if (verbose) cat("Post-processing took",as.numeric(difftime(Sys.time(),tm,units='secs')),"seconds.\n")

  
  return(list(fitteddat, preddat, varcomp))
}

## Simulations ##
B <- 2000

simstodo_intercept <- expand.grid(
  K = c(100,200),
  Nk = c(10,20),
  sigma0 = c(1,2),
  sigma3 = 0,
  rho = 0,
  itr = 1:B
)
simstodo_slope <- expand.grid(
  K = c(100,200),
  Nk = c(10,20),
  sigma0 = 2,
  sigma3 = 1,
  rho = 0.5,
  itr = 1:B
)
simstodo <- rbind(simstodo_intercept,simstodo_slope)
simstodo$simidx <- 1:nrow(simstodo)

simlisttodo <- vector(mode='list',length=nrow(simstodo))
for (i in 1:length(simlisttodo)) simlisttodo[[i]] <- as.numeric(simstodo[i, ])
dosim <- function(lst) {
  cat("Simulation",lst[7],"\n")
  sim <- tryCatch(do_simulation_multiple(K=lst[1],Nk=lst[2],sigma0=lst[3],sigma3=lst[4],rho=lst[5],iter=lst[6]),error=function(e) e)
  # if (inherits(sim,'condition')) return(NULL)
  sim$sim <- lst[6]
  sim
}

cat("Doing",B,"sims.\n")
tm <- Sys.time()
sims2 <- mclapply(simlisttodo,dosim)

which_errs <- which(Reduce(c,lapply(sims2,inherits,what='condition')))
err_messages <- lapply(sims2[which_errs],'[[','message')
err_messages <- t(as.data.frame(err_messages))
colnames(err_messages) <- 'message'
rownames(err_messages) <- NULL




suffix <- gsub('-', '', as.character(Sys.Date()))

save(sims2,file = file.path(savepath,paste0("simsraw-",suffix,".RData")))
# save(simlist2,file = file.path(savepath,paste0("sims-",suffix,".RData")))

tm2 <- round(as.numeric(difftime(Sys.time(),tm,units='secs')),1)
cat("Done. It took ",tm2," seconds.\n",sep="")


# errs <- !sapply(simlist2[1:numsim],inherits,what='tbl_df')
errs <- !sapply(sims2,function(x) is.tbl(x[[1]])) ## changed error handling since simlist2 already ignores the NULL elements and screws up the rest

numsim <- length(sims2[!errs]) ## without errors
simlist2 <- Reduce(rbind,sims2[!errs])
save(simlist2,file = file.path(savepath,paste0("sims-",suffix,".RData")))


simfitted <- bind_rows(simlist2[(0*numsim+(1:numsim))])
simpred <-   bind_rows(simlist2[(1*numsim+(1:numsim))])
simtheta <-  bind_rows(simlist2[(2*numsim+(1:numsim))])


###############################################################################
## Summarize estimation and inference

fitted_x1 <- simfitted %>% filter(x2==0, x3==0)
fitted_x2 <- simfitted %>% filter(x1==0, x3==0)
alph <- 0.2
## Spaghetti Plots
make_spaghetti_plot <- function(K,Nk,sigma0,sigma3,rho,xp) {
  if (xp==1) {
    rangespag_x1 <- c(-4,4)
    
    plotdat <- filter(fitted_x1,K==!!K,Nk==!!Nk,sigma0==!!sigma0,sigma3==!!sigma3,rho==!!rho)
    
    gamplot <- ggplot(plotdat,aes(x1,y=gammargfitted,group=iter)) +
      geom_line(colour="gray",alpha=alph)+
      geom_line(aes(y=fx1,x=x1),size=1)+
      ylim(rangespag_x1)+
      ylab("Fitted Curve")+
      ggtitle(expression(GAM:~f[1](x[1]))) +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
    mamplot <- ggplot(plotdat,aes(x1,y=mammargfitted,group=iter)) +
      geom_line(colour="gray",alpha=alph)+
      geom_line(aes(y=fx1,x=x1),size=1)+
      ylim(rangespag_x1)+
      ylab("Fitted Curve")+
      ggtitle(expression(MAM:~f[1](x[1]))) +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
    
  } else if (xp==2) {
    rangespag_x2 <- c(-4,1)
    plotdat <- filter(fitted_x2,K==!!K,Nk==!!Nk,sigma0==!!sigma0,sigma3==!!sigma3,rho==!!rho)
    
    gamplot <- ggplot(plotdat,aes(x2,y=gammargfitted,group=iter)) +
      geom_line(colour="gray",alpha=alph)+
      geom_line(aes(y=fx2,x=x2),size=1)+
      ylim(rangespag_x2)+
      ylab("Fitted Curve")+
      ggtitle(expression(GAM:~f[2](x[2]))) +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
    mamplot <- ggplot(plotdat,aes(x2,y=mammargfitted,group=iter)) +
      geom_line(colour="gray",alpha=alph)+
      geom_line(aes(y=fx2,x=x2),size=1)+
      ylim(rangespag_x2)+
      ylab("Fitted Curve")+
      ggtitle(expression(MAM:~f[2](x[2]))) +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
  } else {
    stop("Unknown covariate.")
  }
  
  list(gam=gamplot,mam=mamplot)
}
# Save plots to disk
uniquesims <- filter(simstodo,itr==1)
for (i in 1:nrow(uniquesims)) {
  rw <- uniquesims[i, ,drop=FALSE]
  for (xp in c(1,2)) {
    with(rw,{
      tmp <- make_spaghetti_plot(K,Nk,sigma0,sigma3,rho,xp)
      nm <- paste0("-K",K,"-Nk",Nk,"-sl",sigma0,"-ss",sigma3,"-r",10*rho,"-x",xp,".pdf")
      ggsave(file = file.path(plotpath,paste0("gam",nm)),tmp$gam,width=PLOTWIDTH,height=PLOTHEIGHT)
      ggsave(file = file.path(plotpath,paste0("mam",nm)),tmp$mam,width=PLOTWIDTH,height=PLOTHEIGHT)
    })
  }
}

## Summarize Operating Characteristics

covr_x1 <- simfitted %>%
  filter(x2==0, x3==0) %>%
  group_by(x1,K,Nk,sigma0,sigma3,rho) %>%
  summarize(
    gammargcovr=mean(gammargcovrx1),
    mammargcovr=mean(mammargcovrx1),
    gammcondcovr=mean(gammcondcovrx1)
  )

covr_x2 <- simfitted %>%
  filter(x1==0, x3==0) %>%
  group_by(x2,K,Nk,sigma0,sigma3,rho) %>%
  summarize(
    gammargcovr=mean(gammargcovrx2),
    mammargcovr=mean(mammargcovrx2),
    gammcondcovr=mean(gammcondcovrx2)
  )

## Coverage Plots
make_coverage_plot <- function(K,Nk,sigma0,sigma3,rho,xp) {
  if (xp==1) {
    rangecovr_x1 <- c(0.5,1)
    
    plotdat <- filter(covr_x1,K==!!K,Nk==!!Nk,sigma0==!!sigma0,sigma3==!!sigma3,rho==!!rho)
    
    gamplot <- ggplot(data=plotdat,aes(x=x1,y=gammargcovr)) +
      geom_line() +
      ylim(rangecovr_x1) +
      ylab("Coverage") +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
    mamplot <- ggplot(data=plotdat,aes(x=x1,y=mammargcovr)) +
      geom_line() +
      ylim(rangecovr_x1) +
      ylab("Coverage") +
    theme(text = element_text(size = GGPLOTTEXTSIZE))
    
  } else if (xp==2) {
    rangecovr_x2 <- c(0.5,1)
    plotdat <- filter(covr_x2,K==!!K,Nk==!!Nk,sigma0==!!sigma0,sigma3==!!sigma3,rho==!!rho)
    
    gamplot <- ggplot(data=plotdat,aes(x=x2,y=gammargcovr)) +
      geom_line() +
      ylim(rangecovr_x2) +
      ylab("Coverage") +
    theme(text = element_text(size = GGPLOTTEXTSIZE))
    mamplot <- ggplot(data=plotdat,aes(x=x2,y=mammargcovr)) +
      geom_line() +
      ylim(rangecovr_x2) +
      ylab("Coverage") +
    theme(text = element_text(size = GGPLOTTEXTSIZE))
  } else {
    stop("Unknown covariate.")
  }
  list(gam=gamplot,mam=mamplot)
}
# Save plots to disk
uniquesims <- filter(simstodo,itr==1)
for (i in 1:nrow(uniquesims)) {
  rw <- uniquesims[i, ,drop=FALSE]
  for (xp in c(1,2)) {
    with(rw,{
      tmp <- make_coverage_plot(K,Nk,sigma0,sigma3,rho,xp)
      nm <- paste0("-covr-K",K,"-Nk",Nk,"-sl",sigma0,"-ss",sigma3,"-r",10*rho,"-x",xp,".pdf")
      ggsave(file = file.path(plotpath,paste0("gam",nm)),tmp$gam,width=PLOTWIDTH,height=PLOTHEIGHT)
      ggsave(file = file.path(plotpath,paste0("mam",nm)),tmp$mam,width=PLOTWIDTH,height=PLOTHEIGHT)
    })
  }
}

### Relative bias in SEs

sesd_x1 <- simfitted %>%
  filter(x2==0, x3==0) %>%
  group_by(x1,K,Nk,sigma0,sigma3,rho) %>%
  summarize(
    sesd_gam = (mean(gammargse)/sd(gammargfitted))-1,
    sesd_mam = (mean(mammargse)/sd(mammargfitted))-1
  )

sesd_x2 <- simfitted %>%
  filter(x1==0, x3==0) %>%
  group_by(x2,K,Nk,sigma0,sigma3,rho) %>%
  summarize(
    sesd_gam = (mean(gammargse)/sd(gammargfitted))-1,
    sesd_mam = (mean(mammargse)/sd(mammargfitted))-1
  )

make_sesd_plot <- function(K,Nk,sigma0,sigma3,rho,xp) {
  if (xp==1) {
    rangesesd_x1 <- c(-0.5,0.5)
    
    plotdat <- filter(sesd_x1,K==!!K,Nk==!!Nk,sigma0==!!sigma0,sigma3==!!sigma3,rho==!!rho)
    
    gamplot <- ggplot(data=plotdat,aes(x=x1,y=sesd_gam)) +
      geom_line() +
      ylim(rangesesd_x1) +
      ylab("Rel. Bias(SE)") +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
    mamplot <- ggplot(data=plotdat,aes(x=x1,y=sesd_mam)) +
      geom_line() +
      ylim(rangesesd_x1) +
      ylab("Rel. Bias(SE)") +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
    
  } else if (xp==2) {
    rangesesd_x2 <- c(-0.5,0.5)
    plotdat <- filter(sesd_x2,K==!!K,Nk==!!Nk,sigma0==!!sigma0,sigma3==!!sigma3,rho==!!rho)
    
    gamplot <- ggplot(data=plotdat,aes(x=x2,y=sesd_gam)) +
      geom_line() +
      ylim(rangesesd_x2) +
      ylab("Rel. Bias(SE)") +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
    mamplot <- ggplot(data=plotdat,aes(x=x2,y=sesd_mam)) +
      geom_line() +
      ylim(rangesesd_x2) +
      ylab("Rel. Bias(SE)") +
      theme(text = element_text(size = GGPLOTTEXTSIZE))
  } else {
    stop("Unknown covariate.")
  }
  list(gam=gamplot,mam=mamplot)
}



###############################################################################
## Table for paper
estimationtable_x1 <- simfitted %>%
  filter(x2==0,x3==0) %>%
  group_by(K,Nk,sigma0,sigma3,rho,x1) %>%
  summarize(
    bias_x1_gam = mean(gammargbiasx1),
    covr_x1_gam = mean(gammargcovrx1),
    sesd_x1_gam = (mean(gammargse)/sd(gammargfitted)),
    bias_x1_mam = mean(mammargbiasx1),
    covr_x1_mam = mean(mammargcovrx1),
    sesd_x1_mam = (mean(mammargse)/sd(mammargfitted))
  )%>%
  group_by(K,Nk,sigma0,sigma3,rho) %>% ## second grouping because now we collapse over x1 (now that we have ratios of SE to SD)
  summarize_all(mean)
  
estimationtable_x2 <- simfitted %>%
  filter(x1==0,x3==0) %>%
  group_by(K,Nk,sigma0,sigma3,rho,x2) %>%
  summarize(
    bias_x2_gam = mean(gammargbiasx2),
    covr_x2_gam = mean(gammargcovrx2),
    sesd_x2_gam = (mean(gammargse)/sd(gammargfitted)),
    bias_x2_mam = mean(mammargbiasx2),
    covr_x2_mam = mean(mammargcovrx2),
    sesd_x2_mam = (mean(mammargse)/sd(mammargfitted))
  ) %>%
  group_by(K,Nk,sigma0,sigma3,rho) %>% ## second grouping because now we collapse over x2 (now that we have ratios of SE to SD)
  summarize_all(mean)


predinttable <- simpred %>%
  group_by(K,Nk,sigma0,sigma3,rho) %>%
  summarize(
    gammRMSEP_int = sqrt(mean(gammMSEPb0,na.rm=TRUE)),
    gammRMSEP_slope = sqrt(mean(gammMSEPb1,na.rm=TRUE)),
    mamRMSEP_int = sqrt(mean(mamMSEPb0,na.rm=TRUE)),
    mamRMSEP_slope = sqrt(mean(mamMSEPb1,na.rm=TRUE))
  )
varcomptable <- simtheta %>%
  group_by(K,Nk,sigma0,sigma3,rho) %>%
  summarize(
    mam_bias_int = mean(mam_sd_int - sigma0),
    gamm_bias_int = mean(gamm_sd_int - sigma0),
    
    mam_bias_slope = mean(mam_sd_slope - sigma3),
    gamm_bias_slope = mean(gamm_sd_slope - sigma3),
    
    mam_bias_cor = mean(mam_cor - rho),
    gamm_bias_cor = mean(gamm_cor - rho)
  )

joinvars <- c("K","Nk","sigma0","sigma3","rho")
resultstable <- estimationtable_x1 %>%
  left_join(estimationtable_x2,by = joinvars) %>%
  left_join(varcomptable,by = joinvars) %>%
  left_join(predinttable,by = joinvars) 

write_csv(resultstable,file = file.path(plotpath,"numeric-summaries.csv"))

###############################################################################
## Checking the double penalization
pdf(file = file.path(plotpath,"double-pen-intercepts.pdf"),width=PLOTWIDTH,height=PLOTHEIGHT)
plotdat <- filter(simfitted,K==200,Nk==20,sigma0==2,sigma3==0,rho==0)
plot(plotdat$mam2margfitted~plotdat$mammargfitted,cex=0.5,lwd=0.2,
     xlab="MAM",ylab="Doubly Penalized MAM");abline(0,1)
dev.off()
pdf(file = file.path(plotpath,"double-pen-slopes.pdf"),width=PLOTWIDTH,height=PLOTHEIGHT)
plotdat <- filter(simfitted,K==200,Nk==20,sigma0==2,sigma3==1,rho==0.5)
plot(plotdat$mam2margfitted~plotdat$mammargfitted,cex=0.5,lwd=0.2,
     xlab="MAM",ylab="Doubly Penalized MAM");abline(0,1)
dev.off()


cat("Done. Output saved at",globalpath,"\n")

# # Helper code for printing the table, not run
# # Table 1 (manuscript)
# resultstable %>%
#   filter(sigma0==2,sigma3==1,rho==.5) %>%
#   arrange(sigma0,K,Nk) %>%
#   mutate(across(contains('covr'),~round(.x*100))) %>%
#   mutate(across(contains('sesd'),~((.x-1)))) %>%
#   dplyr::select(contains(c("K","Nk",'mam'))) %>%
#   knitr::kable(
#     digits = 2,
#     format = 'latex' # markdown, not latex, easier copying
#   )
# # Table 1 (supplememnt)
# resultstable %>%
#   filter(sigma3==0,rho==0) %>%
#   arrange(sigma0,K,Nk) %>%
#   mutate(across(contains('covr'),~round(.x*100))) %>%
#   mutate(across(contains('sesd'),~((.x-1)))) %>%
#   dplyr::select(-contains(c("slope","cor"))) %>%
#   dplyr::select(contains(c("K","Nk","sigma0",'gam'))) %>% # For printing each variable individually
#   dplyr::select(!contains(c('gamm'))) %>%
#   knitr::kable(
#     digits = 2,
#     format = 'latex' # markdown, not latex, easier copying
#   )





