#///////////////////////////////////////////////////////////////////
#### Bay of Fundy sea scallop SSM: simulate data and fit model  ####
#///////////////////////////////////////////////////////////////////

# Y. Yin, W. H. Aeberhard, S. J. Smith, and J. Mills Flemming (2018).
# Identifiable state-space assessment models, a case study of the Bay of Fundy
# sea scallop fishery. Canadian Journal of Statistics

rm(list=ls())
setwd("Y:/Projects/CANSSI/scripts")
direct <- "Y:/Projects/CANSSI/"
direct <- "D:/Dropbox/TMB/"
direct <- "D:/TMB/"
### load libraries and functions
library(TMB)
library(readxl)
# The function used to simulate the data
source(paste0(direct,'scripts/SS2v4simul.r')) # simulate data according to alternative model
# Load in the data
compile(paste0(direct,"scripts/SS2v4_no_popcorn.cpp")) # alternative model
dyn.load(dynlib(paste0(direct,"scripts/SS2v4_no_popcorn")))

# Are you simulating the data or using real data
simulate<- FALSE
# Do you want to use the catch equation or simply use the model with catch as a known value
use.catch.eq <- TRUE
# Do you want to run the popcorn model or just run the mortalizty random walk
popcorn <- TRUE
# If you are not simulating the data, which area are you getting the actual data from...
# Currently 6 possibilities, "SPA3", "SPA4", and "GBa",BBn, SPA1A, SPA1B, and SPA6.  This won't work with 29 anytime soon.
area <- "SPA4"

if(area %in% c("GBa","BBn")) 
{
  load(paste0(direct,"Data/GBa_model_data.RData")) # This has BBn and GBa in it...
  dat <- as.list(mod.dat[[area]][mod.dat[[area]]$year > 1900,])
  dat$years <- dat$year
  dat$C <- dat$catch
}
  NY <- length(dat$years)
  datalist<- dat[c("I","IR","clappers","C","g","gR","N")]
  names(datalist)[[3]]<- "L"
  # Get this as an instantaneous rate, won't generally make much difference...
  datalist$mobs <- 1-exp(-datalist$L/datalist$N)
  #datalist$I <- datalist$I/1000
  #datalist$IR <- datalist$IR/1000
  #datalist$C <- datalist$C/1000
  datalist<- lapply(datalist,tail,NY)
  yearsvec <- 1:NY
}
# When enable_catcheq is set to 1, the fitted model is the full alternative model; 
# when enable_catcheq is set to 0, equation(9) is ignored and a reduced model is fitted. 
if(use.catch.eq == TRUE) datalist<-c(datalist,list(enable_catcheq = 1))
if(use.catch.eq == FALSE) datalist<-c(datalist,list(enable_catcheq = 0))
# Similarly for the popcorn model, if we want to include it set this to 1, if we want to skip it set it to 0
if(popcorn == TRUE) datalist<-c(datalist,list(popcorn = 1))
if(popcorn == FALSE) datalist<-c(datalist,list(popcorn = 0))

#log_a = 0, log_sigma_C = 0, log_chi = -1
parlist <- list('log_sigma_tau'=0,'log_sigma_phi'=0,'log_sigma_m'=0,
                'log_sigma_epsilon'=0,'log_sigma_upsilon'=0,
                'log_sigma_kappa'=0,'log_sigma_C'=0,
                'log_q_I'=-1,'log_q_R'=-1,
                'log_S'=0,'log_a'=0,'log_chi'=-1,
                # 'log_B'=rep(log(max(dat$It)*10),NY),
                # 'log_R'=rep(log(max(dat$Jt)*10),NY),
                'log_B'=rep(log(max(datalist$I)*10),NY),
                'log_R'=rep(log(max(datalist$IR)*10),NY),
                'log_m'=rep(log(0.05),NY))

#A switch to adjust the model depending on the setting: wheteher enable_catch is equal to 0 or not.
maplist <- list() # This is used if we are using the catch equation and popcorn model.
if(use.catch.eq == FALSE)  maplist <- list(log_a = as.factor(NA),log_sigma_C = as.factor(NA),log_chi = as.factor(NA))
#if(use.catch.eq == FALSE && popcorn == FALSE) maplist <- list(log_a = as.factor(NA),log_sigma_C = as.factor(NA),log_chi = as.factor(NA),
#                                                              log_sigma_kappa = as.factor(NA),log_S = as.factor(NA))
#if(use.catch.eq == TRUE && popcorn == FALSE)  maplist <- list(log_sigma_kappa = as.factor(NA),log_S = as.factor(NA))

obj <- MakeADFun(data=datalist,parameters=parlist,
                 random=c('log_B','log_R','log_m'),
                 map = maplist,
                 DLL="SS2v4_no_popcorn",silent=T)

system.time(opt <- try(nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,
                              control=list(eval.max=1000,iter.max=1000)),T))
# ^ less than 1 sec
# The desirable return codes are 3, 4, 5, and sometimes 6.  See the Bell labs documentation which was saved at
#Y:\Admin\software\NLMINB optimizer - Bell Labs original Usage Summary for Selected Optimization Routines.pdf
opt$message # converged properly?

# Try a different optimizer, it claims it converges but results are basically trash
#system.time(opt <- try(optim(par=obj$par,fn=obj$fn,gr=obj$gr,method = "BFGS",
#                              control=list(maxit=1000))))
#opt$convergence

### look at fixed param estimates
system.time(rep <- sdreport(obj,bias.correct=F))
# ^ less than 1 sec
summ.rep <- summary(rep)
rep
summ.rep[rownames(summ.rep) %in% c("q_I","q_R","S"),]
#summ.rep[(length.theta+3*NY+1):(2*length.theta+3*NY) - 3,]

inst.mort <- mean(summ.rep[dimnames(summ.rep)[[1]]=='m',1])
mean.mort <- 1 -exp(-inst.mort) 
mean.mort
#inst.mort <- -log(-mean.mort+1) # To get back to instaneous mortality I always mess these maths up, never remember to move the 1 to the left side first

# Save the summ.rep object

save(summ.rep,dat,file = paste0(direct,"Results/Model_area_comparisons/TMB_model",area,".RData"))



#######################################
# Now let's implement some model validation for this, first let's see if the Laplace approximation 
# is happy.

# Here is the check of the Laplace approximation
cc <-TMB::checkConsistency(obj, n =100)
summary(cc)

cc
# So based on XXX samples this looks like it is ?? 
summary(cc)$marginal
summary(cc)$joint



# ^ parame estimates and se
## plot predicted randeff with CI as colored envelope
pred.Bt <- summ.rep[dimnames(summ.rep)[[1]]=='B',1]
pred.logBt <- summ.rep[dimnames(summ.rep)[[1]]=='log_B',1]
se.pred.logBt <- summ.rep[dimnames(summ.rep)[[1]]=='log_B',2]
lb.ci.Bt <- exp(pred.logBt-1.96*se.pred.logBt) # 95% CI lower bound
ub.ci.Bt <- exp(pred.logBt+1.96*se.pred.logBt) # 95% CI lower bound

pred.Rt <- summ.rep[dimnames(summ.rep)[[1]]=='R',1]
pred.logRt <- summ.rep[dimnames(summ.rep)[[1]]=='log_R',1]
se.pred.logRt <- summ.rep[dimnames(summ.rep)[[1]]=='log_R',2]
lb.ci.Rt <- exp(pred.logRt-1.96*se.pred.logRt) # 95% CI lower bound
ub.ci.Rt <- exp(pred.logRt+1.96*se.pred.logRt) # 95% CI lower bound

pred.mt <- summ.rep[dimnames(summ.rep)[[1]]=='m',1]
pred.logmt <- summ.rep[dimnames(summ.rep)[[1]]=='log_m',1]
se.pred.logmt <- summ.rep[dimnames(summ.rep)[[1]]=='log_m',2]
lb.ci.mt <- exp(pred.logmt-1.96*se.pred.logmt) # 95% CI lower bound
ub.ci.mt <- exp(pred.logmt+1.96*se.pred.logmt) # 95% CI lower bound


# When running the reduced model (i.e., when enable_catcheq = 0), uncomment the below segement
pred.Bt.Reduced <- pred.Bt
lb.ci.Bt.Reduced <- lb.ci.Bt
ub.ci.Bt.Reduced <- ub.ci.Bt
pred.Rt.Reduced <- pred.Rt
lb.ci.Rt.Reduced <- lb.ci.Rt
ub.ci.Rt.Reduced <- ub.ci.Rt
pred.mt.Reduced <- pred.mt
lb.ci.mt.Reduced <- lb.ci.mt
ub.ci.mt.Reduced <- ub.ci.mt

# When running the full model (i.e., when enable_catcheq = 1), uncomment the below segement
# pred.Bt.Full <- pred.Bt
# lb.ci.Bt.Full <- lb.ci.Bt
# ub.ci.Bt.Full <- ub.ci.Bt
# pred.Rt.Full <- pred.Rt
# lb.ci.Rt.Full <- lb.ci.Rt
# ub.ci.Rt.Full <- ub.ci.Rt
# pred.mt.Full <- pred.mt
# lb.ci.mt.Full <- lb.ci.mt
# ub.ci.mt.Full <- ub.ci.mt

colmed <- c('#2b05ff') # color code for pred
colenv <- paste0(colmed,'30') # color for envelope

par(mfrow=c(3,1))
# Bt
library(scales)

plot(yearsvec,pred.Bt.Reduced,type='o',col=colmed[1],pch=0,
     xlab='Years',ylab=expression(italic(B[t])),
     main=expression('Predicted commercial biomass'~italic(B[t])),
     ylim=c(range(pred.Bt.Reduced)),las=1)
grid(nx=NA,ny=NULL,equilogs=F)
polygon(c(yearsvec,yearsvec[NY:1]),c(lb.ci.Bt.Reduced,ub.ci.Bt.Reduced[NY:1]),
        col=colenv[1],border=NA,xpd=F)

#lines(yearsvec,pred.Bt.Full,type='o', col='red')
#polygon(c(yearsvec,yearsvec[NY:1]),c(lb.ci.Bt.Full,ub.ci.Bt.Full[NY:1]),
#        col=alpha('red',0.4),border=NA,xpd=F)

#legend("topright", col = c(colmed[1], 'red'), lty = 1, bty="o", legend=c("Reduced", "Full"))


# Rt
plot(yearsvec,pred.Rt.Reduced,type='o',col=colmed[1],pch=0,
     xlab='Years',ylab=expression(italic(R[t])),
     main=expression('Predicted recruitment biomass'~italic(R[t])),
     ylim=c(range(pred.Rt.Reduced)),las=1)
grid(nx=NA,ny=NULL,equilogs=F)
polygon(c(yearsvec,yearsvec[NY:1]),c(lb.ci.Rt.Reduced,ub.ci.Rt.Reduced[NY:1]),
        col=colenv[1],border=NA,xpd=F)
#lines(yearsvec,pred.Rt.Full,type='o', col='red')
#polygon(c(yearsvec,yearsvec[NY:1]),c(lb.ci.Rt.Full,ub.ci.Rt.Full[NY:1]),
#        col=alpha('red',0.4),border=NA,xpd=F)

#legend("topright", col = c(colmed[1], 'red'), lty = 1, bty="o", legend=c("Reduced", "Full"))

# mt
plot(yearsvec,pred.mt.Reduced,type='o',col=colmed[1],pch=0,
     xlab='Years',ylab=expression(italic(m[t])),
     main=expression('Predicted natural mortality rate'~italic(m[t])),
     ylim=c(range(pred.mt.Reduced)),las=1)
grid(nx=NA,ny=NULL,equilogs=F)
#polygon(c(yearsvec,yearsvec[NY:1]),c(lb.ci.mt.Reduced,ub.ci.mt.Reduced[NY:1]),
#        col=colenv[1],border=NA,xpd=F)
#lines(yearsvec,pred.mt.Full,type='o', col='red')
#polygon(c(yearsvec,yearsvec[NY:1]),c(lb.ci.mt.Full,ub.ci.mt.Full[NY:1]),
#        col=alpha('red',0.4),border=NA,xpd=F)

# END BoF_SSM


x1_name <- "pred.Bt"
x2_name <- "pred.Rt"
x3_name <- "pred.mt"
y_name <- "yearsvec"


df <- data.frame(yearsvec,pred.Bt,pred.Rt,pred.mt)
colnames(df) <- c(y_name,x1_name,x2_name,x3_name)
print(df)
write.csv(df,'df.csv')
