#!/usr/bin/env Rscript

# simulates lineage-trees with the growth-progression model. 
# Erika Kuchen 2019
 
# how to add a legend in R for the output figures? 
# change parameters in INPUT section for uplaod

# call mat library
library("R.matlab") # this is only necessary if files should be saved as a .mat matlab output file.


################### USER INPUT 

# model parameters
k <- 0.0405 # cell growth rate 
sig_g <- 0.02 # standard deviation of size threshold noise
mu <- 2.76 # mean progression time 
sig_p <- 0.22 # standard deviation of progression time 
alpha <- 0.27 # mother-daughter progression time correlation
gamma <- 0.58  # sibling progression time correlation


# simulation parameters 
growthTerm <- 'exponential' #type of growth function options 'logistic' and 'exponential' 
smax <- 5 # saturation threshold for logistic growth 
corrtype <- "spearman" # type of correlation coefficient calculated, options from R function cor: method = c("pearson", "kendall", "spearman"))

init <- 30  # number of families simulated, up to 1522 cells possible. 
genssimu <- 7 # number of generations simulated per family 
checkstationary <- 0 # options: yes/no (1/0), simulate trees up to xgen generations and then only continue 50% of cells to save memory 
xgen <- 5 # see checkspationary above 

runs <- 100 # number of simulation repeats to generate confidence bounds.
initialRuns <- 100 # number of initialisation runs  


# saving the output 
expname <- 'name_of_experiment' # for saving the output
setwd("my project path") # set path where results should be saved  
plotoutput <- 1 # options: 0,1
saveoutput <- 0 # options: 0,1.  


# parameters that do not need to be adjusted 
snorm <- 1 # normalised threshold cell size 
xd <- seq(0,6,by=0.001) # log(time) over which progression time probability distribution is evaluated 



############## Initialise data frames and arrays 

if(checkstationary==0){
rho <- array(NA,dim=c(runs,9)) # cycle length correlation coefficients 
aX <- array(NA,dim=c(runs,6)) # autocorrelation growth variable 
aC <- array(NA,dim=c(runs,6)) # autocorrelation progression variable 
cX <- array(NA,dim=c(runs,6)) # crosscorrelation growth on progression  
cC <- array(NA,dim=c(runs,6)) # crosscorrelation progression on growth 
mstats <- array(NA,dim=c(runs,4))
sf <- array(NA,dim=c(runs,2))
aX[,1]  <- 1 
aC[,1] <- 1
}
	
if (runs >1){	
	savepop <- data.frame(array(NA,dim=c(1,37))) # data frame to save all repeat runs 
	columnnames <- c("Tdiv","finalSize","cellDivided","progressionTime","limitingProcess","TdivSister","TdivMother","TdivAunt","TdivGrandmother","TdivGreatgrandmother","TdivGreataunt","TdivCousinOnceRemoved1","TdivCousinOnceRemoved2","cellID","TdivCousin1","TdivCousin2","TdivSecondCousin1","TdivSecondCousin2","TdivSecondCousin3","TdivSecondCousin4","progressionTimeMother","progressionTimeGrandmother","progTGreatgrandmother","progTGreatgreatgranny","progTGreatgreatgreatgranny","growthTime","growthTimeMother","growthTimeGrandmother","growthTGreatgranny","growthTGreatgreatgranny","growthTGreatgreatgreatgranny","initialSize","generation","absoluteTimeofDivision","progressionTimeAunt","progTCousin1","progTCousin2") #	
	names(savepop) <- columnnames
}  

# caculate necessary length of output data frame
if (checkstationary){
	xgen <- 5
	pred0 <- init*(2^(1+xgen)-1)
	pred <- pred0*(genssimu-xgen)
} else{
	pred <- init*(2^genssimu-1) 
}

	
# possible ID numbers assigned to a lineage
initnumber <-c(1:10,13:20,23:30,33:40,43:50,53:60,63:70,73:80,83:90,93:100,103:110,113:120,123:130,133:140,143:150,153:160,163:170,173:180,183:190,193:200,203:210,213:220,223:230,233:240,243:250,253:260,263:270,273:280,283:290,293:300,303:310,313:320,323:330,333:340,343:350,353:360,363:370, 373:380,383:390,393:400,403:410, 413:420,423:430,433:440,443:450,453:460,463:470,473:480,483:490,493:500,503:510,513:520,523:530,533:540,543:550,553:560,563:570,573:580,583:590,593:600, 603:610,613:620,623:630,633:640,643:650,653:660,663:670,673:680,683:690,693:700, 703:710,713:720,723:730,733:740,743:750,753:760,763:770,773:780,783:790, 793:800,903:910,913:920,923:930,933:940,943:950,953:960,963:970,973:980,983:990,993:1000,1003:1010,1013:1020,1023:1030,1033:1040,1043:1050,1053:1060, 1063:1070,1073:1080,1083:1090,1093:1100,1103:1110,1113:1120,1123:1130,1133:1140,1143:1150,1153:1160,1163:1170,1173:1180,1183:1190,1193:1200)

episum <- sum(exp(-(xd-mu-alpha*(mu-mu))^2/(2*sig_p^2*(1-alpha^2)))/(sig_p*sqrt((1-alpha^2)*2*pi))) # used for computing the progression time probability distribution.


########## Simulate the lineage trees  
for (reps in 1:runs){ # repeat simulation x times specified by 'runs'
         
   # generate some random numbers for simulating noise                
   condR <- matrix(runif(2*pred),pred,2)
   coI <- matrix(runif(init*initialRuns),init,initialRuns)
   Xt <- matrix(rnorm(2*pred),pred,2)*sig_g + snorm;
   Xtp <- matrix(rnorm(init*initialRuns),init,initialRuns)*sig_g + snorm;

   # data frame to store all lineage tree information    
   pop <- data.frame(array(NA,dim=c(pred,37))) 
   names(pop) <- columnnames 

   pop[1:init,] <- matrix(rep(c(0,1,0,exp(mu),1,array(NA,dim=c(1,27)),1,0,array(NA,dim=c(1,3))),each=init),nrow=init)   
   pop[1:init,14] <- initnumber[1:init] 
   pop[1:init,2] <- rnorm(init)*sig_g + snorm 
 
 
####### Initialise the lineage trees 
   for (su in 1:initialRuns){  # simulate the model without divisions to initialise founder cells
       for (ce in 1:init){ # iterate over all founder cell

            sb <- pop[ce,2]/2 # initial cell size
            sth1 <- Xtp[ce,su] # threshold cell size
                                    
            # determine progression time tp 
            nM <- log(pop[ce,4]) # progression time before 
            ptd <- (exp(-(xd-mu-alpha*(nM-mu))^2/(2*sig_p^2*(1-alpha^2)))/(sig_p*sqrt((1-gamma^2)*2*pi)))/episum # progression time probability distribution dependent on previous progression time 
            ptd <- ptd/sum(ptd)                
            pos <- sum(coI[ce,su] >= cumsum(ptd)) 
            p1 <- exp(xd[pos+1])

			switch(growthTerm,
      			logistic={  
                   tg1 <- log((sth1*(sb-smax))/(sb*(sth1-smax)))/k
                   tg1 <- max(0,tg1) 
                 },
       			exponential={ 
       				tg1 <- log(sth1/sb)/k
       				tg1 <- max(0,tg1)}                             
 			)                                   
                Tdiv1 <- max(p1, tg1) # cycle length 
                fl1 <- (p1<tg1)+1  # fl = 1 if cell was progression limited, else 2 if growth limited 
                  
   			switch(growthTerm,  # calculate the final cell size                                   
                logistic={ se1 <- smax/(1+ ((smax-sb) / sb)*exp(-k*Tdiv1))}, 
                exponential={se1 <- sb*exp(k*Tdiv1)}
            )         
            pop[ce,c(1:2,4,5,26,32)] <- c(Tdiv1,se1,p1,fl1,tg1,sb) # update the lineage tree data frame                 
    	   }
	} 
            
   # distribute birth time of lineage founder cells      
   medpoptot <- median(pop[1:init,1])
   interval = medpoptot/init
   pop[1:init,34] = seq(-medpoptot+interval,0,by=interval) 
 
 
 
######## grow the lineage trees               
   nextc <- init
					
   while (max(pop[,33],na.rm = TRUE)<genssimu){ # run until the maximum generation number specified by genssimu is reached   
                             
     if (checkstationary){ # simulate for longer to check for stability but discontinue randomly selected cells.      
        if (max(pop[,33],na.rm = TRUE)>xgen){ # randomly terminate elongating around 50% of cells. Use to check age distribution and cell size at a particular time.                     
           use <- pop[,33]==max(pop[,33],na.rm = TRUE) & !is.na(pop[,33]) 
           rvals <- sample(2, sum(use,na.rm = TRUE), replace=TRUE)
           pop[use,3] <- rvals-1 # pretend that this cell has divided already and is ignored in the division process below. 
         } 
     }             
                             
	 ind <- which(pop[,3]==0, arr.ind = FALSE, useNames = TRUE) # find cells that have not divided yet and divide these iteratively 
                                
        for (j in ind){
                                    
            pop[j,3] <- 1; # here cell is updated as having divided

            sb <- pop[j,2]/2; # initial cell size 
            sth1 <- Xt[j,1]; # size threshold daughter 1 
            sth2 <- Xt[j,2]; # size threshold daughter 2 

            # determine progression time of daughter 1 dependent on mother progression time 
            nM <- log(pop[j,4]);
            ptd <- (exp(-(xd-mu-alpha*(nM-mu))^2/(2*sig_p^2*(1-alpha^2)))/(sig_p*sqrt((1-alpha^2)*2*pi)))/episum
            ptd <- ptd/sum(ptd) 
            pos <- sum(condR[j,1] >= cumsum(ptd)); 
            p1 <- exp(xd[pos+1]);
                                    
            # determine progression time of daughter 2 dependent on mother and daughter 1                      
            pdgdm <- (exp((((xd[pos+1]-mu)*(alpha^2-sig_p))^2 +((xd-mu)*(alpha^2-1))^2 +((nM-mu)*alpha*(1-gamma))^2-2*(xd-mu)*(1-alpha^2)*((xd[pos+1]-mu)*(gamma-alpha^2) + (1-gamma)*alpha*(nM-mu)) + (xd[pos+1]-mu)*(nM-mu)*2*alpha*(gamma-1)*(alpha^2-gamma))/(2*sig_p^2*(1-alpha^2)*(gamma^2-1+2*alpha^2*(1-gamma))))*(sqrt((1-alpha^2)/(2*pi*(1-gamma^2+2*alpha^2*(gamma-1))))/sig_p))/episum;
            pdgdm <- pdgdm/sum(pdgdm) 
            pos <- sum(condR[j,2] >= cumsum(pdgdm));
            p2 <- exp(xd[pos+1]);    
                              
            # determine growth time of the two daughters                                       
            switch(growthTerm,
             	logistic={  
             	 	tg1 <- log((sth1*(sb-smax))/(sb*(sth1-smax)))/k
                    tg1 <- max(0,tg1) 
                    tg2 <- log((sth2*(sb-smax))/(sb*(sth2-smax)))/k
                    tg2 <- max(0,tg2)}, 
                exponential={
                		tg1<- log(sth1/sb)/k
              		tg2 <- log(sth2/sb)/k
              		tg1 <- max(0,tg1)
              		tg2 <- max(0,tg2)}
            )

			# determine cycle length of the two daughter cells as the maximum of the progression time and the growth time 	
            Tdiv1<- max(p1, tg1) # cycle length daugther 1
            fl1 <- (p1<tg1)+1 # fl = 1 if cell was progression limited, else 2 if growth limited 
            Tdiv2<- max(p2, tg2) 
            fl2 <- (p2<tg2)+1  
              
            # re-evaluate final cell sizes 
            switch(growthTerm,                    
              logistic={
              	se1 <- smax/(1+ ((smax-sb) / sb)*exp(-k*Tdiv1))
              	se2 <- smax/(1+ ((smax-sb) / sb)*exp(-k*Tdiv2))},
            	   exponential={
            	   	se1 <- sb*exp(k*Tdiv1)
            		se2 <- sb*exp(k*Tdiv2)}  
            )
                  
            # give each daughter an ID number derived from the mother ID number                    
            nextc <- nextc+1;
            Id1 <- pop[j,14]*10+1
            Id2 <- pop[j,14]*10+2
            done <- 0; 
            momnum <- toString(pop[j,14])
            n <- nchar(momnum) 
               
            # to save large searches on the lineage tree to find relatives, here the cycle time of all relatives is stored for each cell
            if(checkstationary==0){   # this does not work with checkstationary=1 as it relies on a specific tree structure.                            
               if (pop[j,14]>10 && grepl(paste(c("1","2"),collapse="|"),substr(momnum,n,n))){       
              	  if (grepl("1",substr(momnum,n,n))){   
                     # if mother ends in a 1, the sister has not been done yet                                        
                 	 posshift <- c(2,3)
              	  } else {
                 	 posshift <- c(-2,-1)
                  }
             	  # includes the cousin cycle times and other relatives in the row for each cell.                        
             	  pop[nextc+posshift[1],15:16] <- c(Tdiv1,Tdiv2)
             	  pop[nextc+posshift[2],15:16] <- c(Tdiv1,Tdiv2)
             	  pop[nextc+posshift[1],36:37] <- c(p1,p2)
             	  pop[nextc+posshift[2],36:37] <- c(p1,p2)
               }
      
               if (pop[j,14]>100 && grepl(paste(c("1","2"),collapse="|"),substr(momnum,n,n)) && grepl(paste(c("1","2"),collapse="|"),substr(momnum,n-1,n-1))){
               	   if (grepl("1",substr(momnum,n-1,n-1)) && grepl("1",substr(momnum,n,n))){ 
               		 # if mother ends in a 1, the sister has not been done yet                  	
                 	 pop[nextc+(4:7),17:18] <- matrix(rep(c(Tdiv1,Tdiv2),each=4),nrow=4)
               	   } else if (grepl("1",substr(momnum,n-1,n-1)) && grepl("2",substr(momnum,n,n))){
                 	 pop[nextc+(2:5),19:20] <- matrix(rep(c(Tdiv1,Tdiv2),each=4),nrow=4)        
               	   } else if (grepl("2",substr(momnum,n-1,n-1)) && grepl("1",substr(momnum,n,n))){
                 	 pop[nextc-(1:4),17:18] <- matrix(rep(c(Tdiv1,Tdiv2),each=4),nrow=4)                          
              	   } else if (grepl("2",substr(momnum,n-1,n-1)) && grepl("2",substr(momnum,n,n))){
                 	 pop[nextc-(3:6),19:20] <- matrix(rep(c(Tdiv1,Tdiv2),each=4),nrow=4)
                   }
			   }
 			}
                                    
        		pop[nextc,1:14] <- c(Tdiv1,se1,done,p1,fl1,Tdiv2,pop[j,c(1,6,7,9,8,15,16)],Id1) # upadate the lineage tree data frame, information of first daugther 
       		# now save duration of individual processes for autocorrelation and crosscorrelation calculations
        		pop[nextc,21:25] <- c(pop[j,c(4,21:24)]) #progression time mother, granny, etc 
        		pop[nextc,26:31] <- c(tg1,pop[j,26:30]) #tg mother, granny, etc
			pop[nextc,32:34] <- c(sb,pop[j,33]+1,pop[j,34]+pop[j,1]) # initial size, generation, absolut time of division
			if (j>init){ #progression time of other relatives saved here
		       if (grepl("1",substr(momnum,n,n))){
				  pop[nextc,35] <- pop[j+1,4]# aunt p1 from mother sibling (aunt)
			   }else{
				  pop[nextc,35] <- pop[j-1,4]# aunt p1 from mother sibling (aunt)
			   }
			}else{
			   pop[nextc,35] <- NA
			}
		    # pop column 36,37 information of first cousins
            
            nextc <- nextc+1
     		pop[nextc,1:14] <- c(Tdiv2,se2,done,p2,fl2,Tdiv1,pop[j,c(1,6,7,9,8,15,16)],Id2) # upadate the lineage tree data frame, information of second daugther                             
     		pop[nextc,21:25] <- c(pop[j,c(4,21:24)]) #progression time mother, granny, gg and more, check
     		pop[nextc,26:31] <- c(tg2,pop[j,26:30]) #tg mother, granny, more #aT
     		pop[nextc,32:34] <- c(sb,pop[j,33]+1,pop[j,34]+pop[j,1]) # initial size, generation, absolut time of division
     		if (j>init){
     			if (grepl("1",substr(momnum,n,n))){
					pop[nextc,35] <- pop[j+1,4]# aunt p1 from mother sibling (aunt)
				}else{
					pop[nextc,35] <- pop[j-1,4]# aunt p1 from mother sibling (aunt)
				}
			}else{
				pop[nextc,35] <- NA
			}                                
		}
	}



###### do stats on the simulated lineage trees 

	if (checkstationary==0){ # information on relatives etc will not be saved under checkstationary 

		sf[reps,] <- c(sum(pop[,5]==1),sum(pop[,5]==2)) # number of progression-limited and growth-limited cells 
		mstats[reps,1] <- mean(pop[,1],na.rm=TRUE) # mean cycle length
    		mstats[reps,2] <- median(pop[,1],na.rm=TRUE) # median cycle length
    		mstats[reps,3] <- quantile(pop[,1],0.25,na.rm=TRUE) # quartiles cycle length
    		mstats[reps,4] <- quantile(pop[,1],0.75,na.rm=TRUE)

 		# calculate the correlations in cycle times between relatives 
 		n <- dim(pop)[1]
    		val1 <- pop[(init+1):n,1] # mother-daughter
    		val2 <- pop[(init+1):n,7]
    		rho[reps,2] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype) 
  		val2 <- pop[(init+1):n,6] # siblings
   		rho[reps,1] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
   		val1 <- pop[(3*init+1):n,1]
   		val2 <- pop[(3*init+1):n,9]
   		rho[reps,3] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
   		val2 <- pop[(3*init+1):n,8]; # aunt
   		rho[reps,5] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)				
   		val1 <- c(pop[(3*init+1):n,1],pop[(3*init+1):n,1])
   		val2 <- c(pop[(3*init+1):n,15],pop[(3*init+1):n,16])
   		rho[reps,6] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype) # first cousins
   		val1 <- pop[(7*init+1):n,1]
   		val2 <- pop[(7*init+1):n,10]
   		rho[reps,4] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val2 <- pop[(7*init+1):n,11]
 		rho[reps,7] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
    		val1 <- c(pop[(7*init+1):n,1],pop[(7*init+1):n,1])
    		val2 <- c(pop[(7*init+1):n,12],pop[(7*init+1):n,13])
    		rho[reps,8] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
    		val1 <- c(pop[(7*init+1):n,1],pop[(7*init+1):n,1],pop[(7*init+1):n,1],pop[(7*init+1):n,1])
    		val2 <- c(pop[(7*init+1):n,17],pop[(7*init+1):n,18],pop[(7*init+1):n,19],pop[(7*init+1):n,20])
		rho[reps,9] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)

	
		#auto and cross-correlations
    		val1 <- pop[(init+1):n,4] # mother-daughter
    		val2 <- pop[(init+1):n,21]
		aC[reps,2] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(3*init+1):n,4] # grandmother
    		val2 <- pop[(3*init+1):n,22]
		aC[reps,3] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
	    val1 <- pop[(7*init+1):n,4]
   		val2 <- pop[(7*init+1):n,23] # greatgranny 
		aC[reps,4] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
	    val1 <- pop[(15*init+1):n,4] # 31, 15, 63
   		val2 <- pop[(15*init+1):n,24] # greatgreatgranny 
		aC[reps,5] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(63*init+1):n,4] # 31, 15, 63
   		val2 <- pop[(63*init+1):n,25] # greatgreatgreatgranny  
		aC[reps,6] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)

		val1 <- pop[(init+1):n,26] # mother-daughter
    		val2 <- pop[(init+1):n,27]
		aX[reps,2] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype) 
		val1 <- pop[(3*init+1):n,26]
    		val2 <- pop[(3*init+1):n,28] # granny
		aX[reps,3] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(7*init+1):n,26]
   		val2 <- pop[(7*init+1):n,29] # greatgranny 
		aX[reps,4] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
	    val1 <- pop[(15*init+1):n,26] # 31, 15, 63
   		val2 <- pop[(15*init+1):n,30] # greatgreatgranny	
		aX[reps,5] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(63*init+1):n,26] # 31, 15, 63
   		val2 <- pop[(63*init+1):n,31] # greatgreatgreatgranny  	
		aX[reps,6] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)

		val1 <- pop[,4] # within cell
    		val2 <- pop[,26]
		cX[reps,1] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype) 
		val1 <- pop[(init+1):n,4] # mother-daughter
    		val2 <- pop[(init+1):n,27]
		cX[reps,2] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype) 
		val1 <- pop[(3*init+1):n,4]
    		val2 <- pop[(3*init+1):n,28] # granny
		cX[reps,3] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(7*init+1):n,4]
   		val2 <- pop[(7*init+1):n,29] # greatgranny 
		cX[reps,4] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
	    val1 <- pop[(15*init+1):n,4] # 31, 15, 63
   		val2 <- pop[(15*init+1):n,30] # greatgreatgranny	
		cX[reps,5] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(63*init+1):n,4] # 31, 15, 63
   		val2 <- pop[(63*init+1):n,31] # greatgreatgreatgranny  	
		cX[reps,6] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)

    		val1 <- pop[,26] # within cell 
    		val2 <- pop[,4]
		cC[reps,1] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
    		val1 <- pop[(init+1):n,26] # mother-daughter
    		val2 <- pop[(init+1):n,21]
		cC[reps,2] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(3*init+1):n,26]
    		val2 <- pop[(3*init+1):n,22]
		cC[reps,3] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
	    val1 <- pop[(7*init+1):n,26]
   		val2 <- pop[(7*init+1):n,23] # greatgranny 
		cC[reps,4] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
	    val1 <- pop[(15*init+1):n,26] # 31, 15, 63
   		val2 <- pop[(15*init+1):n,24] # greatgreatgranny 
		cC[reps,5] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)
		val1 <- pop[(63*init+1):n,26] # 31, 15, 63
   		val2 <- pop[(63*init+1):n,25] # greatgreatgreatgranny  
		cC[reps,6] <- cor(val1,val2, use="pairwise.complete.obs", method=corrtype)

	}	

#### save the population structure for every repeat run 
	if (runs >1 && checkstationary==0){
		savepop <- rbind(savepop,pop) # concatenate all repeat runs into one large array
	}
	if(checkstationary){
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(reps),'reps_stationary.mat',sep=""), pop=pop)
	} 
}

if(runs >1 && checkstationary==0){ # cut out the first row in savepop, which is NAs
	n <- dim(savepop)
	savepop <- savepop[2:n[1],]	
}


###### generate some output files
# currently results are only saved as matlab files 
if (saveoutput){ 
	if(checkstationary==0){
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_mstats.mat',sep=""), mstats=mstats)
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_rho.mat',sep=""), rho=rho)
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_savepop.mat',sep=""), savepop=savepop) 
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_aX.mat',sep=""), aX=aX)   
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_aC.mat',sep=""), aC=aC)   
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_cC.mat',sep=""), cC=cC)   
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_cX.mat',sep=""), cX=cX)   
		writeMat(paste(expname,growthTerm,toString(genssimu),'gen',toString(runs),'runs_sf.mat',sep=""), sf=sf)
	}
}


######## output figures 
if (plotoutput){

	par(mfrow=c(2,2))

	# subplot 1 - plot the correlations between relatives
	plot(c(1,9),c(0,0), type='c',xlim=c(0.5, 9.5), ylim=c(-0.2, 0.8), , ann=FALSE
, xaxt = "n", ylab="Spearman rank correlation coefficient")
	title(paste(expname,growthTerm)) 
	for (j in 1:9){
		rl <- quantile(rho[,j],0.025)
   	 	ru <- quantile(rho[,j],0.975)
   	 	lines(c(j,j),c(rl, ru),col='red')	
	}
	points((1:9),colMeans(rho), pch=17, col='red') # mean not working yet
	axis(1,1:9,labels=FALSE)
	text(1:9,labels=as.vector(c("sister","mother","grandmother","greatgrandmother","aunt","first cousin","greataunt","cousin once removed","second cousin")), par("usr")[3]- 0.2, srt = 45, pos = 1, xpd = TRUE)


	# subplot 2 - proportion of growth-limited and progression-limited cells 
	pp <- sum(savepop[,5]==2)/dim(savepop)[1]
	barplot(c(pp,1-pp),horiz=FALSE, names.arg=c("growth", "progression"))

	#subplot 3 - distribution of cycle lengths
	hist(savepop[,1],100, xlab="cycle length (h)", ylab="relative frequency", 		main=NULL)

	# subplot 4 - auto, crosscorrelation
	plot(1:6,colMeans(aX),pch=16, ylim=c(-1, 1), xlim=c(0.8,7), col=c("black"), xlab="ancestral generation", ylab="corerelation (auto, cross)" )
	lines(1:6,colMeans(aX),col=c("black"),lty=3)
	lines(1:6,colMeans(aC),col=c("green"),lty=3)
	#lines(1:6,colMeans(cX),lty=3)
	lines(1:6,colMeans(cC),lty=3,col=c("red"))
	lines(1:4,colMeans(rho[,1:4]),col=c("yellow"),lty=3)
	points(1:6,colMeans(aC),pch=15,col=c("green"))
	points(1:6,colMeans(aX),pch=19,col=c("black"))
	points(1:4,colMeans(rho[,1:4]),pch=17,col=c("yellow"))
	points(1:6,colMeans(cC),pch=18,col=c("red"))
	legend(5.3,1.1, legend = c(expression(paste(tau)),
                              expression(paste(tau["p"])),
                              expression(paste(tau["g"])),
                              expression(paste(tau["p"], tau["g"]))), 			 			 			col=c("yellow","green","black","red"),pch=c(17,15,19,18))
}

