####################################################################
## Delay difference model.  Necessary inputs are survey estimates for  biomass, abundance, and growth for recruits and commerical sized 
# individuals and the total catch for each year.
###################################################################
# Update history
# Revised function to be used for Fisheries course
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
# 5:  mod.dat                       The data to load, can be  a "csv" or "RData" file.  See previous, default is a "RData" file
#                               located here....paste(direct,"Data/Model_input.RData",sep="")
###########  Results options.  These options set up how the results from the model will be processed and influence our predictions

#16:  M.prior:                      We can specify our priors for mortality of commercial and recruit "fish".
#                                    default = data.frame(a=2,b=6)
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


#######################################################################################################################################

###############################################################################################################

run_model <- function(mod.dat = mod.dat, 
                        m.prior = data.frame(a=2,b=6), 
                        q.prior =  data.frame(a=15,b=35),
                        r.prior = data.frame(a=30000,b=1000),
                        b0.prior= data.frame(a=30000,b=0.5),
                        # The main model options, these are only used if run.mod = T. (Tho the options starting with "n" could be sent to 
                        # be used with the prediction evaluation model runs if the "pe." options are not specified and run.pre.eval.model=T)
                        nchains = 6,niter = 175000, nburn = 100000, nthin = 20,parallel = T,
                        seed = 123,parameters = NULL)
  
{
  
  # Load in the functions needed for this function to run.
  

  # The necessary library
  library("R2jags")
  
  
  #############  Section 1 Model#############  Section 1 Model#############  Section 1 Model#############  Section 1 Model ###########
  #############  Section 1 Model#############  Section 1 Model#############  Section 1 Model#############  Section 1 Model ###########
  #############  Section 1 Model#############  Section 1 Model#############  Section 1 Model#############  Section 1 Model ###########
  
  # Set the working directory for figures and tables to be outpu
  
  DD.dat <- mod.dat
  
  # Organize the data and set up the model priors/initialization data, then run the model.
  yrs<-min(DD.dat$year):max(DD.dat$year)
  DD.lst<-as.list(DD.dat)
  DD.lst$NY<- length(yrs)
  
  DDpriors=list(
    # survey catchability fully recruited a= shape1, b=shape2
    q=				    list(a=q.prior$a, 	       b=q.prior$b,	 d="dbeta",	  l=1),		
    # process error (SD) a = min, b = max
    sigma=			  list(a=0.1, 		         b=5,		         d="dunif",	  l=1),		
    # measurement error variance survey FR a = shape, b = scale (1/rate)
    I.precision=	list(a=0.5,	             b=1,	           d="dgamma",	l=1),		
    # measurement error variance survey recruits a = shape, b = scale (1/rate)
    IR.precision=	list(a=0.5,	             b=1,            d="dgamma",	l=1),	
    # FR natural mortality, default mostly b/t 0.1-0.4
    M=				    list(a=m.prior$a,        b=m.prior$b,		 d="dbeta",	  l=1),
    # scaled recruit biomass, a= meanlog  b = sdlog
    r=				    list(a=log(r.prior$a),   b=r.prior$b,    d="dlnorm",	l=DD.lst$NY),
    # initial biomass estimate
    bs=				    list(a=log(b0.prior$a),  b=b0.prior$b,	 d="dlnorm",	l=1)
    
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
  ifelse(is.null(parameters) == T, parameters <- c(names(DDpriors),'B','R','mu','Imed','Ipred','IRmed','IRpred',
                                                   'sIresid','sIRresid','sPresid','Iresid','IRresid','Presid'),parameters)
  
  # go grab the model...
  tmp <- "https://raw.githubusercontent.com/Dave-Keith/Fisheries_module/master/Scripts/Functions/delay_difference_model_base.bug"
  download.file(tmp,destfile = basename(tmp))
  loc <- paste0(getwd(),"/",basename(tmp))

  
  # Run the model
  start<-Sys.time()
  ## Call to JAGS, do you want to run in parallel?
  
  if(parallel==F)
  {
    out <- jags(data =  c(prior.lst,DD.lst), inits = NULL,parameters.to.save = parameters,  
                model.file = loc,n.chains = nchains, n.iter = niter, n.burnin = nburn, 
                n.thin = nthin)
  }
  
  if(parallel==T)
  {
    out <- jags.parallel(data =  c(prior.lst,DD.lst), inits = NULL,parameters.to.save = parameters,  
                         model.file = loc,n.chains = 8, n.iter = 10000, n.burnin = 1000, 
                         n.thin = 20,jags.seed = seed)
  }
  # How long did that take?
  print(Sys.time()-start)
  
  # Rename the output so I retain the results 
  DD.out <- list(data=c(prior.lst,DD.lst,yrs), sims.list=out$BUGSoutput$sims.list,median=out$BUGSoutput$median,
                 mean=out$BUGSoutput$mean,summary=out$BUGSoutput$summary,priors = prior.lst,parameters=parameters)
  
  # I will also retain the MCMC object produced in case I want it for something.
  mod.out <- out
  
  
  # Return the model results
  assign("DD.lst", DD.lst, pos = 1) 
  assign("DDpriors", DDpriors, pos = 1) 
  assign("DD.out", DD.out, pos = 1) 
  assign("mod.dat", mod.dat, pos = 1) 
  assign("DD.dat", DD.dat, pos = 1) 
  
  ############# END Section 2 Model#############  END Section 2 Model#############  END Section 2 Model#############  
  ############# END Section 2 Model#############  END Section 2 Model#############  END Section 2 Model#############  
  ############# END Section 2 Model#############  END Section 2 Model#############  END Section 2 Model#############  
  
  
} # end function
