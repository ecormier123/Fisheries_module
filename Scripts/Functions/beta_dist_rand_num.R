beta.val <- function(mn,sd,fx)
{
library(zipfR)
# This function caluclates a random number from a beta distribution with mean (mn), st. deviation (sd), and
# cumulative distribution function (fx).  This function in MATLAB uses something called the "betainc" function
# will try and use something analogous in R.., see page 277 in book...

if(sd == 0) bb <- mn
    toler <- 1e-9 # how close the CDF value of the answer must be to the input value (Fx)
    var <- sd^2
    (if (var >= (1-mn)* mn) {
     return("You dumb SOB, your variance is too big.. jerk")
     break # End the run
     }
     else
     {
#        return("I like cheese")
    vv <- mn*((mn*(1-mn)/var)-1) # Here we are getting the a and b parameters for the beta distribution
    ww <- (1-mn)*((mn*(1-mn)/var)-1)
     
 up.value <- 1
low.value <- 0

# Start using an initial guess for x, using the random number generator to adds some "wiggle"
# to the start of the search, should help to avoid biases in this search...

     x <- 0.5 + 0.02*runif(1)
   
    i <- Rbeta(x,vv,ww) #This calculates the CDF value of x, THE Rbeta is RIGHT, IBeta is not the same function see R and Matlab help!!
    
# Now the below loop iteratively searches for a value of x, trying to get the value of
# x has a CDF that is within the tolerance set above.  Unless close to 0 or 1, this would cause problems
# and would stop the search.
#  }) # temporary end of the else statement from way above!!
while((toler < abs(i-fx)) && ((x > 1e-6)) && ((1-x) > 1e-6))
    {
        (if (fx > i)
         {
         low.value <- x 
         x <- (up.value+low.value)/2
         }
         else
         {
         up.value <- x
         x <- (up.value+low.value)/2
         })
    i <- Rbeta(x,vv,ww)
    }
  
# Now the below makes values of x kinda random, this will get rid of problems associated with variances that are very small or very large
# This will also truncate x, with small values = tolerance and large values = 1- tolerance.

bbb <- x + toler*0.1*(0.5-runif(1))
if(bbb< toler) bbb <- toler
if(bbb>1) bbb <- 1-toler

bb <- bbb

     }) # end the else statement from way above!!
#return(i)
} # End of function...