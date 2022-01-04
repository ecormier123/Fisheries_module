betaset <- function(vr.types,vr.means,vr.vars,min,max)
{
#this function makes a set of betavalues to be used in to 
# generate correlated vital rate values. 
bbs = zeros(99,1)


if(vr.types != "beta" && vr.types != "st.beta")
{
    if(vr.types != 1 && vr.types != 2)
    {
        message(paste("Dude value '", count, "' of vr.types is '" ,vr.types, "' it gotta be either between 1 and 2 or
        one of \"beta\" or \"st.beta\" "))
    } else
if(vr.types == 1) vr.types <- "beta"
if(vr.types == 2) vr.types <- "st.beta"
} # end vr.types if statement

 	 for(fx99 in 1:99)
         {
  		if(vr.types =="beta") bbs[fx99] = beta.val(vr.means,sqrt(vr.vars),fx99/100)
		if(vr.types =="st.beta") bbs[fx99] = stretch.beta.val(vr.means,sqrt(vr.vars), min, max, fx99/100)
        } #fx100
return(bbs)
}