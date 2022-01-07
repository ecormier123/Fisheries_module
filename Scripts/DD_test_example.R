####################################################################
## Here I take annual survey results for Redfish and massage them
## So that they can be used within a Delay difference model.
## You will need to figure out the growth terms

# First import the silver hake data, it contains rows with adult and recruit biomass, and their uncertainty, and 
# weights at various sizes
direct = "d:/r/Courses/EAF/WG3/Delay_Difference/"

#dat <- read.csv(paste(direct,"Data/Silver_Hake_raw.csv",sep=""))
dat <- read.csv(paste(direct,"Data/Redfish_raw.csv",sep=""))

######################  Get our model parameters... Biomass ones are easy peasy.
NY <- length(dat$year)
mod.dat <- data.frame(year =  min(dat$year):max(dat$year),I = NA,IR=NA,g=NA,gR=NA)


# Can we say is the recruts in this model are the age 1 individuals so our recruit growth
# should be the growth from age 1 to age 2.  
mod.dat <- data.frame(I = , 
                      IR = , 
                      g = , 
                      GR = ,
                      catch = ,
                      )

###############################################################################################################
###############  Section 2, RUN THE MODEL###############  Section 2, RUN THE MODEL###############  Section 2, RUN THE MODEL
###############  Section 2, RUN THE MODEL###############  Section 2, RUN THE MODEL###############  Section 2, RUN THE MODEL
###############################################################################################################
# Now you are all set to run the model, lucky!!  If you have your own data make sure it is set up like the mod.dat object above
# You will need all the variables that are in the mod.dat file and they need to have the same names as this object does.
# The exception is cpue and cpue.cv if you run DD_base.bug model it does not require any CPUE data.

# Load in the model function
source(paste(direct,"functions/Model_function.R",sep=""))

# Run the model, see the function for an explaination of all the model options.
run_model(dat = mod.dat, direct = direct,fig="pdf"
          jags.model = "models/delay_difference_model_base.bug",
          nchains=6,niter = 1500,nburn = 1000)


## The fully recruited biomass for the previous year and current year
FR.bm
# Long term median fully recruited biomass (not including the most recent year)
FR.ltm
# Recruite Biomass for previous and current year
rec.bm
# Long term median recruitment
rec.ltm
# Range of effective sample size
neff
# Range of rhat values
rhat
# Median projection 1 year out.
BM.proj.1yr
# Natural mortality for the most recent year... 
mort
mort.R

# The probability of currently being below the USR.
prob.below.USR

# The "percent" change in biomass, to be a real percent multiply by 100.
100*percent.B.change

# We can also look at the paramters we used for our projections, hopefully they are as expected!
summary(DD.out$data$projection.parameters)
# How variable are these really...
apply(DD.out$data$projection.parameters,2,function(x) quantile(x,probs=seq(0,1,by=0.1)))


# If you have simulation data you can compare the simulation data with the data you used.
# Probably most importantly what does the simulated vs. real biomass look like?
plot(DD.out$median$B~sim.dat$B,pch=19)
abline(a=0,b=1,col="blue")
plot(DD.out$median$R~sim.dat$r,pch=19)
abline(a=0,b=1,col="blue")

# Another way to visualize this...
plot(DD.out$median$B~DD.out$data$year,pch=19,type="o",lwd=2) # modeled biomass
lines(sim.dat$B~sim.dat$year,col="blue",lty=2) # the actual biomass
# Now the same for the recruits...
plot(DD.out$median$R~DD.out$data$year,pch=19,type="o",lwd=2,ylim=c(min(c(min(sim.dat$r),min(DD.out$median$R))),
                                                                   max(c(max(sim.dat$r),max(DD.out$median$R)))))
lines(sim.dat$r~sim.dat$year,col="blue",lty=2) # the actual recruit biomass


# You can do this for any parameters of interest...
plot(sim.dat$m ~ DD.out$median$m,pch=19)
plot(sim.dat$mR ~ DD.out$median$mR,pch=19)
# The growth terms
plot(sim.dat$g ~ DD.out$median$G,pch=19)
plot(sim.dat$gR ~ DD.out$median$GR,pch=19)

# Catchability
DD.out$median$q
median(sim.dat$q)

# Process errors...
DD.out$median$sigma
median(sim.dat$sigma)

