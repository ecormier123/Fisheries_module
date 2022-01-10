### This function is used to make the posterior historgrams and to compare these with the priors used when applicable.  This also
### plots the prediction posteriors if they are available tho they are far less interesting...

# Update history
# DK revisions April 2016
# updated in Jan 2018 to allow for png's to be saved. for the single posterior plot, for the annuals this needs to stay a pdf due to the number of pages

#####################################  Function Summary ########################################################
####  
##  This function is used within these files:(a.k.a "dependent files") 
##
##  1:  Update_function_JAGS.r
###############################################################################################################

###############################################################################################################
## This function needs these functions to work (a.k.a. "support files")
# 
##
###############################################################################################################


###############################################################################################################
# Arguments
# 1:  model.out:  The model output, works with raw model output and with the projections. Default = missing, must be supplied
# 2:  priors:     The model priors.  Default = missing, must be supplied
# 8:  graphic:    Where to plot the figures.  Options are 'screen' (default) and 'pdf' and 'png'
#11:  direct:     Directory to save if producing a pdf, default = ""

posterior.plt <- function(model.out, priors,  graphic='screen',direct="")
                     
{

  # If graphic is a pdf make the pdf file
	if(graphic=='pdf') pdf(paste0(direct,"/Results/Figures/posterior_plot.pdf"), width = 8.5, height = 11, pointsize = 16)
  if(graphic=='png') png(paste0(direct,"/Results/Figures/posterior_plot.png"), width = 8.5, height = 11,res=920,units="in")
  
  # We don't want the deviance
  model.out$sims.list <- subset(model.out$sims.list,names(model.out$sims.list) != "deviance")
  # This picks out the model parameters for which there is only 1 parameter for the model (rather than the parameters estimated annually)
  posts<- model.out$sims.list[lapply(model.out$sims.list,ncol)==1]
	# If you don't specify both the number of rows (nr) and columns (nc) then use the data to determine them, 
	# if you only specify one of them it will default to calculating from the number of parameters we have...
	# Note that deviance isn't plotted but is part of the "posts"
  # If less than 15 posts we'll make it a 2 column figure, if more than 14 I'll go with 3 columns and hope it works...
  ifelse(length(posts) <= 14,nc <- 2,nc <- 3)
  nr <- ceiling(length(posts)/nc)
	# Set up the plot device
	par(mfrow = c(nr, nc), mar = c(2, 3, 1, 1), omi = c(0.4, 0.6, 0, 0.2))
	
	# Now run the for each parameter calculated once
	for (i in 1:(length(posts)))
	{
		  # If the name in posts[i] is found in the prior list do the below
		  # Set the range from min-max in the data, this is reset for beta, log-normal, and gamma distributions if we have the prior at the next step
		  xl<-range(posts[[i]])
      # If we have a posterior and a prior for a parameter.
		  if(names(posts)[i] %in% names(priors))
			{
			  # If the distribution is beta x goes from 0-1
				if(priors[[names(posts)[i]]]$d =="dbeta") xl=c(0,1)
				# If dlnorm or dgamma it runs from 0 to maximum observed
				if(priors[[names(posts)[i]]]$d%in%c("dlnorm","dgamma")) xl<-c(0,max(posts[[i]]))
				# Make a fake x that goes from the min-max with as many values as there are posterior samples.
				x<-seq(xl[1], xl[2], l = length(posts[[i]]))
				# If we have a log-normal or normal prior our precision term needs to be converted back to a standard deviation
				if(priors[[names(posts)[i]]]$d %in% c("dnorm","dlnorm")) priors[[names(posts)[i]]]$b <- 1/sqrt(priors[[names(posts)[i]]]$b) 
				# Using the prior specified distribution get values across the range of the data.
				p<-get(priors[[names(posts)[i]]]$d)(x,  priors[[names(posts)[i]]]$a,  priors[[names(posts)[i]]]$b)
			} # if(names(posts)[i] %in% names(priors))
			# The histogram of the posteriors
			tmp <- hist(posts[[i]], breaks = 25, main = "", prob = T, ylab = "", las = 1, mgp=c(1,0.4,0), tcl=-0.3, xlab="",cex.axis=1.2,xlim=xl)
			# If we have a uniform prior we need to adjust it so the uniform line shows up at the proper density
			if(names(posts)[i] %in% names(priors) && priors[[names(posts)[i]]]$d == "dunif") 
			{
			  p <- rep((length(posts[[i]])/length(tmp$breaks))/length(posts[[i]])/(tmp$breaks[2]-tmp$breaks[1]),length(posts[[i]]))
			} # end if(names(posts)[i] %in% names(priors) && priors[[names(posts)[i]]]$d == "dunif") 
			
 			# If the posterior had a prior assigned then add the line.
			if(names(posts)[i] %in% names(priors)) lines(x, p, col = 'red')
			# Add the x-axis label.
			mtext(names(posts)[i], 1, 2, cex=1)
			# The y-axis label.
			if(i == 1) mtext("Posterior density", 2, 2, outer = T, adj = 0.5, cex=1.25)
	} # end for (i in 1:(length(posts)))
  # Stick a y axis label on here, hopefully at a nice location...
	mtext("Posterior density", 2, 2, outer = T, adj = 0.5, cex=1.25)
	# If making a pdf shut it down.
	if(graphic %in% c('pdf','png')) dev.off()	
		
	
}# end function
	


