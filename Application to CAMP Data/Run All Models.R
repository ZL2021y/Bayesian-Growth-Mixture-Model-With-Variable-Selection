
#rm(list = ls())

mydir <- "C:" # directory in which codes are saved
source(paste(mydir,"/Blasso GMM.R",sep=""))
source(paste(mydir,"/Horseshoe GMM.R",sep=""))
source(paste(mydir,"/Noninformative GMM.R",sep=""))
source(paste(mydir,"/SSVS GMM.R",sep=""))


source(paste(mydir,"PG.R",sep=""))
workdir <- "C:" # directory in which data are saved
mydat <- read.csv(paste(workdir,"/CAMP_real_analyze_data.csv",sep=""),header=TRUE)

# pre-process data
mydat$y <- mydat$PREFEV
mydat$time <- mydat$AGE - mean(mydat$AGE)	# centering age
subj <- unique(mydat$NEWCAMP)
N <- length(unique(mydat$NEWCAMP))
n.obs <- NULL; 		# number of observation per subject
id <- NULL
for (i in 1:N) {n.obs <- c(n.obs,length(mydat[mydat$NEWCAMP==subj[i],]$NEWCAMP)); 
		    id  <- c(id , rep(i,length(mydat[mydat$NEWCAMP==subj[i],]$NEWCAMP))) }
mydat$ID <- mydat$id <- id
mydat <- mydat[order(mydat$ID,mydat$time),]  # Sort by ID and time
mydat_uq <- mydat[!duplicated(mydat$id, fromLast=TRUE),] # Keep last observation per ID
covar <- mydat_uq[,c("bPREFEV","race","trt","any_pets","age_current_home","gas_stove","wood_stove","dehumidifier","msmoke","hosp_asthma","age_ast_conf","mother_asthma",
			"white_cell","Hemoglobin")];

#--------------
#- Running ALL Models
#--------------
num.cluster = 4; burn.in=10000;thin=50; per = 1000; max.iter = 60000
res.lasso <- byes.lasso(	mydat = mydat,						   
	 	covar = covar,			
		num.cluster = num.cluster,		
		burn.in=burn.in, 
		thin=thin, 
		per = per , 
		max.iter = max.iter, myseed=8000)
#--------------
set.seed(20200608)
num.cluster = 4; burn.in=10000;thin=50; per = 1000; max.iter = 60000
res.horseshoe <- byes.horseshoe(mydat = mydat,						   
	 	covar = covar,			
		num.cluster = num.cluster,		
		burn.in=burn.in, 
		thin=thin, 
		per = per , 
		max.iter = max.iter, myseed=2020)
#-------------
set.seed(202007)
num.cluster = 4; burn.in=10000;thin=50; per = 1000; max.iter = 60000
res.SSVS1 <- byes.SSVS(mydat = mydat,						   
	 	covar = covar,
		c2 = 100,tau2 = 0.1,				
		num.cluster = num.cluster,		
		burn.in=burn.in, 
		thin=thin, 
		per = per , 
		max.iter = max.iter, myseed=2020)
#--------------
set.seed(1989)
num.cluster = 4; burn.in=10000;thin=50; per = 1000; max.iter = 60000
res.SSVS2 <- byes.SSVS(mydat = mydat,						   
	 	covar = covar,
		c2 = 100,tau2 = 0.04,				
		num.cluster = num.cluster,		
		burn.in=burn.in, 
		thin=thin, 
		per = per , 
		max.iter = max.iter, myseed=2020)


