####################################################################
## Delay difference model.  Necessary inputs are survey estimates for  biomass, abundance, and growth for recruits and commerical sized 
# individuals and the total catch for each year.
###################################################################
# Update history
# January 2016 - Revised by DK 
# April 2016, "Update" script has been overhauled and converted to a function called "Update_function_JAGS"
# May 16, 2016, updated to include options to allow different databases/usernames in the call to the fishery data.
# Sept 20, 2016, updated to allow for different loader files to be called in from the output of Survey_Summary_data.r. Removed spatial figures
# March 14, 2017 - Revision to the Scallop update function, using somewhat simplier DD model...
#####################################  Function Summary ########################################################
####  
##  This function is used within these files:(a.k.a "dependent files") 
##
###############################################################################################################

###############################################################################################################
## This function needs these functions to work (a.k.a. "support files")
##
###############################################################################################################


###############################################################################################################
# Arguments
# 1:  direct:                   The root directory to work from  Default = "d:/r/", 
# 5:  dat                       The data to load, can be  a "csv" or "RData" file.  See previous, default is a "RData" file
#                               located here....paste(direct,"Data/Model_input.RData",sep="")
# 6:  fig:                      Print to 'pdf' or to screen.  default="screen":
# 7:  plot.loc                  If you are saving plots as pdf's where do you want to put them...

###########  Results options.  These options set up how the results from the model will be processed and influence our predictions

#16:  M.priors:                      We can specify our priors for mortality of commercial and recruit "fish".
#                                    default = data.frame(M.a=2,M.b=6,mR.a=2,MR.b=6)
#17:  q.prior                        The catchability prior, this is modeled as a Beta distribution.  Default is NULL which 
#                                    uses a prior of dbeta(20,40), which is a relatively tight prior centered at 0.33 (this is for the scallop model.)
#                                    q.prior = data.frame(a=10,b=50) would give median q around 0.16, see scripts with "A Beta distribution moment" 
#                                    to understand how to best parameterize these beta distribution priors.
########################################################################################################################

###########  Model options.  By and large these options will only be used if run.mod=T
#20:  nchains:                  Number of chains to run.  Default = 8.  When running in parallel each chain gets it's own thread.
#                               The best way to get more saved replicates if the model has converged is to increase the number of chains run.
#21:  niter:                    Number of iterations to run.  Default = 175000
#22:  nburn:                    Number of initial iterations to ignore.  Default = 50000
#23:  nthin:                    Thinning rate of iterations.  Default = 20
#24:  para:                     Run in parallel using jags.parallel? T/F, default =T. Number of processors = nchains
#25:  jags.model:               Get the model to use (same model used in prediction evaluations).  By deFault it looks the the folder
#                               "direct"/Assessment_fns/Model/DDwSE3_jags.bug where direct was specified above.
#26:  parallel:                 Do you want to run JAGS in parallel.  (T/F), F will run JAGS but just using one core.
#27:  seed:                     If running JAGS in parallel you can set a "seed" so that your model results are reproducable.  Default = 123
#28:  parameters:               Model parameters to output.  Default = NULL which will produce all of the priors + the following parameters
#                               'K','P','B','R','mu','Imed','Ipred','Irep', 'IRmed','IRpred',
#                               'IRrep','sIresid','sIRresid','sPresid','Iresid','IRresid','Presid'
#29  export.tables:             Export the Decision tables, should only be done when satisified with results.  T/F, default = F
#30  convergence.check:         Do you want to run the convergence check. T/F, default = T.  The convergence check produces a pdf
#                               which shows the mixing and ACF for each of the parameters/chains in the model.
#######################################################################################################################################

###############################################################################################################

run_model <- function(direct = "d:/Github/Fisheries_module", yr = as.numeric(format(Sys.time(), "%Y"))-1 , fig="screen",
                        dat = dat, 
                        m.priors = data.frame(M.a=2,M.b=6,MR.a=2,MR.b=6,year="all"), q.prior = NULL,
                        # The output options
                        export.tables = F, plot.loc = paste(direct,"Results/",sep=""),
                        # The main model options, these are only used if run.mod = T. (Tho the options starting with "n" could be sent to 
                        # be used with the prediction evaluation model runs if the "pe." options are not specified and run.pre.eval.model=T)
                        nchains = 8,niter = 175000, nburn = 100000, nthin = 20,parallel = T,strt.mod.yr = 1986,
                        jags.model = "delay_difference_model_base",seed = 123,parameters = NULL,convergence.check = T)
  
{
  
  # Load in the functions needed for this function to run.
  
  
  source(paste(direct,"functions/projections.r",sep=""))
  source(paste(direct,"functions/decision.r",sep=""))
  source(paste(direct,"functions/post.plt.R",sep=""))
  source(paste(direct,"functions/exploit.plt.r",sep=""))
  source(paste(direct,"functions/fit.plt.R",sep=""))
  source(paste(direct,"functions/diag.plt.R",sep=""))
  source(paste(direct,"functions/biomass.plt.r",sep=""))
  # The necesary library
  library("R2jags")
  
  
  #############  Section 1 Model#############  Section 1 Model#############  Section 1 Model#############  Section 1 Model ###########
  #############  Section 1 Model#############  Section 1 Model#############  Section 1 Model#############  Section 1 Model ###########
  #############  Section 1 Model#############  Section 1 Model#############  Section 1 Model#############  Section 1 Model ###########
  
  # Set the working directory for figures and tables to be output
  plotsGo <- plot.loc
  
  DD.dat <- mod.dat
  names(DD.dat) <- c("I","IR","g","gR","C","year") # just making catch turn into "C".
  
  
  # Organize the data and set up the model priors/initialization data, then run the model.
  yrs<-min(DD.dat$year):max(DD.dat$year)
  DD.lst<-as.list(DD.dat)
  DD.lst$NY<- length(yrs)
  
  
  # Set up Priors.  First make sure our m priors are lined up properly.
  if(nrow(m.priors) > 1) m.priors <- m.priors[m.priors$year %in% yrs,] # Subset to the correct number of years if you have multile years
  if(nrow(m.priors) == 1) m.priors <- data.frame(M.a=rep(m.priors$M.a,DD.lst$NY),M.b=rep(m.priors$M.b,DD.lst$NY),
                                                 MR.a=rep(m.priors$MR.a,DD.lst$NY),MR.b=rep(m.priors$MR.b,DD.lst$NY)) # If just one year make it a string...
  # Now our catchability priors
  if(is.null(q.prior)) q.pr <- data.frame(a=0,b=1)# If we don't specify a q prior we make it a fairly specific prior centered at 0.33
  if(!is.null(q.prior)) q.pr <- q.prior # if we do specify a q prior just rename it.
  
  
  DDpriors=list(
    # survey catchability fully recruited a= shape1, b=shape2
    q=				    list(a=q.pr$a, 	       b=q.pr$b,		       d="dlnorm",	  l=1),		
    # process error (SD) a = min, b = max
    sigma=			  list(a=0, 		         b=5,		             d="dunif",	  l=1),		
    # measurement error variance survey FR a = shape, b = scale (1/rate)
    I.precision=	list(a=Ip.a,	         b=Ip.b,	           d="dgamma",	l=1),		
    # measurement error variance survey recruits a = shape, b = scale (1/rate)
    IR.precision=	list(a=IRp.a,	         b=IRp.b,            d="dgamma",	l=1),	
    # FR natural mortality, default mostly b/t 0.1-0.4
    M=				    list(a=m.priors$M.a,   b=m.priors$M.b,		 d="dbeta",	  l=1),	
  )
  
  
  
  
  
  #Prepare priors for JAGS
  for(h in 1:length(DDpriors))
  {
    # Get the variances for log-normal and normal converted to precisions, note that in BUGS language the precision is
    # the inverse of the squared standard deviation (which is what you specify in R).  The standard deviation is what
    # was specified in the Prior list (as it is more intuitive)
    if(DDpriors[[h]]$d%in%c("dlnorm","dnorm")) DDpriors[[h]]$b <- 1/DDpriors[[h]]$b^2
    # For a Gamma to convert to precision the precision term is  the inverse of the 'Scale" term in a typical 
    # gamma distribution parameterization, aka this is now knonwn as the rate.
    # Happily this is the same as the parameterization in R dgamma(x,shape,rate) so our b parameter is correct for posterior plots.
    if(DDpriors[[h]]$d=="dgamma")DDpriors[[h]]$b<-1/DDpriors[[h]]$b
  } # end for(h in 1:length(DDpriors))
  # Made a data.frame of the priors, unwrap the list and combine by row.
  prior.dat<- data.frame(par=names(DDpriors),do.call("rbind",lapply(DDpriors,rbind)))
  prior.lst<-list()
  # Now turn this into a list
  for(k in seq(1,nrow(prior.dat)*2,2))
  {
    prior.lst[[k]]<-prior.dat$a[[ceiling(k/2)]]
    prior.lst[[k+1]]<-prior.dat$b[[ceiling(k/2)]]
  } # end for(k in seq(1,nrow(prior.dat)*2,2))
  # And give the list names
  names(prior.lst)<-paste(rep(prior.dat$par,2)[order(rep(1:nrow(prior.dat),2))],rep(c('a','b'),nrow(prior.dat)),sep='.')
  
  # Now if they haven't already been selected grab the parameters you want for the model.
  ifelse(is.null(parameters) == T, parameters <- c(names(DDpriors),'K','P','B','R','mu','Imed','Ipred','Irep', 'IRmed','IRpred','IRrep',
                                                   'sIresid','sIRresid','sPresid','Iresid','m','mR',
                                                   'IRresid','Presid',"G","GR"),parameters)
  # Run the model
  start<-Sys.time()
  ## Call to JAGS, do you want to run in parallel?
  
  if(parallel==F)
  {
    out <- jags(data =  c(prior.lst,DD.lst), inits = NULL,parameters.to.save = parameters,  
                model.file = paste(direct,jags.model,sep=""),n.chains = nchains, n.iter = niter, n.burnin = nburn, 
                n.thin = nthin)
  }
  
  if(parallel==T)
  {
    out <- jags.parallel(data =  c(prior.lst,DD.lst), inits = NULL,parameters.to.save = parameters,  
                         model.file = paste(direct,jags.model,sep=""),n.chains = nchains, n.iter = niter, n.burnin = nburn, 
                         n.thin = 20,jags.seed = seed)
  }
  # How long did that take?
  print(Sys.time()-start)
  
  # Rename the output so I retain the results 
  DD.out <- list(data=c(prior.lst,DD.lst,yrs), sims.list=out$BUGSoutput$sims.list,median=out$BUGSoutput$median,
                 mean=out$BUGSoutput$mean,summary=out$BUGSoutput$summary,priors = prior.lst,parameters=parameters)
  
  # I will also retain the MCMC object produced in case I want it for something.
  mod.out <- out
  
  
  # Save the model results
  save(DD.lst, DDpriors,DD.out,mod.dat,yr,DD.dat,
       file=paste(direct,"Results/Model_results.RData",sep=""))
  
  
  ############# END Section 2 Model#############  END Section 2 Model#############  END Section 2 Model#############  
  ############# END Section 2 Model#############  END Section 2 Model#############  END Section 2 Model#############  
  ############# END Section 2 Model#############  END Section 2 Model#############  END Section 2 Model#############  
  
  
  
  
  ############# Section 3 Some model result summaries and figures############# ############# Section 3 Some model result summaries and figures############# 
  ############# Section 3 Some model result summaries and figures############# ############# Section 3 Some model result summaries and figures############# 
  # If we want the diagnostics  this is the section.  
  
  
  # Some model outputs needed for the Update.  First the mortality
  mort <- 1- exp(-DD.out$mean$m[length(DD.out$mean$m)])
  mort.R <- 1- exp(-DD.out$mean$mR[length(DD.out$mean$mR)])
  # This lines up the column headers with the projected catch...
  TACI <- which(DD.out$data$C.p==(proj.catch))
  # This get us the predicted biomass for next year based on the projected catch
  BM.proj.1yr <- DD.out$median$B.p[TACI]
  # This is only useful for GBa at the moment since BBn doesn't have reference points accepted yet...
  
  # Get the quantiles, this likely would need changed, but which quantile is > our URP (13,284 as of 2015)
  B.quantiles <- quantile(DD.out$sims.list$B[,length(DD.out$sims.list$B[1,])],probs=seq(0,1,0.01))
  # This is the probability (well percentage) that Biomass is below the USR
  prob.below.USR <- names((which(B.quantiles > URP)[1]))
  
  
  # Here we can grab the Fully recruited and recruit biomass for the last 2 years and the median of the time series.
  FR.bm <- DD.out$median$B[(length(DD.out$mean$B)-1):length(DD.out$median$B)]
  # We exclude the current year from the median estimate
  FR.ltm <- median(DD.out$median$B[-length(DD.out$median$B)])
  # Recruit biomass
  rec.bm <- DD.out$median$R[(length(DD.out$median$R)-1):length(DD.out$median$R)]
  # We exclude the current year from the median estimate
  rec.ltm <- median(DD.out$median$R[-length(DD.out$median$R)])
  
  # Get the percent biomass change from the projection. 0 means unchanged, + means % increase, - means % decline
  percent.B.change <- (BM.proj.1yr / DD.out$median$B[length(DD.out$median$B)]) -1
  
  
  save(mort,mort.R,TACI,BM.proj.1yr,B.quantiles,percent.B.change,prob.below.USR,FR.bm,FR.ltm,rec.bm,rec.ltm,
       file=paste(direct,"Results/Model_results_and_diagnostics.RData",sep=""))
  
  #####Plot model diagnostics############## 
  # These plots include the posterior fits, exploitation estimate, Biomass fit to survey and CPUE, residual plot
  # and the model coinvergence plot (which is a 700+ page pdf of the convergence of each parameter + it's ACF.)
  # posterior densities for model parameters
  post.plt(DD.out,DDpriors,years=yrs, graphic=fig,multi=T,path=plotsGo)
  #dev.off()
  ##exploitaiton time series
  exploit.plt(DD.out, years=yrs, plt=c('f','m','mR'),graphic=fig,path=plotsGo)
  #dev.off()
  # model biomass fit to survey
  fit.plt(DD.out, years = yrs, CI=T,graphic=fig,path=plotsGo,CV=T)
  # diagnostic plot
  diag.plt(DD.out, years = yrs,graphic=fig,path=plotsGo)
  # The biomass plot for fully recruited and recruits along with the projection
  biomass.plt(DD.out,years=yrs, graphic=fig,TAC=proj.catch,path=plotsGo,refs = c("LRP","URP","zones"),pred=1,
              URP =URP, LRP=LRP,avg.line=median)
  
  
} # end function
