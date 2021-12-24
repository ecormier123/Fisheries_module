# Dalhousie Fisheries Module


## To Do before first class

You will need to have [[R](https://cran.r-project.org/)] installed.


## you will need a wack of packages installed. So copy, paste, and run the below code in to make sure you have the packages you need installed.


```r
req.packages <- c("tidyverse","dplyr","sf","units","cowplot","knitr",'ggthemes',"marmap","RandomFields",
                  "ggplot2","stars","tmaptools","rnaturalearth","rnaturalearthdata","raster","rgdal","RStoolbox",
                  "pals","ggnewscale","ggspatial",'devtools','raster','bookdown','tinytex','gridExtra','scales')
# If you don't have the packages install them + give a heads up that you are
new.packages <- req.packages[!(req.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) 
{
  cat(paste0("Heads up, I have to install these packages for this to work:", new.packages ))
  #wanna.install <- readline(prompt = "If you want to install these package(s) enter 'y': ")
  #if(tolower(wanna.install) == 'y') 
  install.packages(new.packages,repos = "http://cran.us.r-project.org") #else { stop("You didn't want to install the packages so this script does not work.")}
}
```

While not 100% necessary, it will be handy to have tinytext fully installed, so run the below line once you have successfully installed the above packages.

```r

tinytex::install_tinytex()
# to uninstall TinyTeX, run tinytex::uninstall_tinytex() 
```

## R Installation

You can also embed plots, for example:

```r

```


