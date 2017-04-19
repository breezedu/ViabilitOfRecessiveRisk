

###############################################
##
## 03-30-2017
## Viability Simulate Under the Power beta!=0
###############################################


###########################################################################################
## Step one
## Read in a vector of allele frequencies.
## The af list was pulled down from ExAC, filtered by CCDS.r14, rare threhold set at 10%
##
############################################################################################

TTN_af <- read.table("/work/AndrewGroup/ViabilitySimulation/QualifyTTN_variants_OnExons.txt", header = T, sep = ",")

## loptop file
## TTN_af <- read.table("D:/PhD/QualifyTTN_variants_OnExons.txt", header = T, sep = ",")

TTN_af <- TTN_af$Allele.Frequency

variants.count <- length( TTN_af )

# head(TTN_af)

###########################################################################################
## 
## Step two
## Pass the allele frequencies we got in step-one to a function simulateGenotypes()
## Return n1 (the counts of heterozygous variants), n2 (the counts of homozygous variants)
##        Pai2g (the sum of Pai2g * Rho) from Equation #8 in the proposal
########################################################################################### 

simulateGenotypes<-function(af.list, alpha, beta){
  
  len <- length(af.list)   # the length of allele frequencies vector
  af.sim1 <- runif(n=len, min=0, max=1)               #generate a vector of random double numbers [0,1]
  af.sim2 <- runif(n=len, min=0, max=1)
  
  ## initialize a vector of zeros to store the simulated allele frequencies
  aflist.sim <- rep(0, len)
  
  n1 <- 0                     ## counts of heterozygous variants
  n2 <- 0                     ## counts of homozygous variants
  
  for(i in 1:len){
    
    ### get counts for the homozygous variants
    if(af.list[i] > af.sim1[i] & af.list[i] > af.sim2[i]){
      
      n2 <- n2 + 1
    } 
    
    ### get counts for the heterozygous
    ## check copy1
    if(af.list[i] > af.sim2[i]){
      n1 <- n1 + 1
      
      aflist.sim[i] <- aflist.sim[i] + 1
    } 
    
    ## check copy2
    if( af.list[i] > af.sim1[i]){
      n1 <- n1 + 1
      
      aflist.sim[i] <- aflist.sim[i] + 1
    }
    
    
  } # end for i in 1:len loop; 
  
  ### calculate Pai2g for current person/genotype
  Pai2g <- 0
  
  if(n2 > 0){
    Pai2g <- 1
    
  } else {
    
    if( n1 < 2) {
      Pai2g <- 0
      
    } else {
      
      Pai2g <- 1 - (0.5)^(n1 - 1)  
      
    }
    
  }
  
  ## check one case when n1>1
  if(n1 > 1){
    print( c(n1, n2, Pai2g))
    #print( aflist.sim )
  }  
  #print( length( c(aflist.sim, Pai2g)))
  
  
  #################################################
  ## Simulate the viability under the pwer, beta!=0
  ## alpha <- -1.92
  ## beta <- 1.80
  
  ## check the Indicator
  # I <- 0
  # if(Pai2g > 0) {
  #  I <- 1
  # }
  
  ## or, we could use I <- Pai2g directly
  I <- Pai2g
  
  viability <- 1 / ( 1 + exp(alpha + beta * I))
  
  random.via <- runif(1)
  
  if( random.via < viability){
    
    return( c(Pai2g, aflist.sim))
  
  } else {
      
      if(Pai2g > 0)
        print(paste(random.via, '>', viability))
      return( simulateGenotypes(af.list, alpha, beta) )
    
  }
  
  ## return a list of variants count on each know allele site, plus the calculated Pai2g for this 'person'
  ## return( c(Pai2g, aflist.sim) )
  
}
  
## test the simulateGenotypes() function
# for(i in 1:100000){
#   simulateGenotypes(TTN_af)
# }

##################################################################
##
## step 2.5
## Calculate the Sum( Pai2g * Rho) for the simulated genotypes
################################################################## 
CalPai2gRho <- function( af.list ){
  
  ## print( summary(af.list) )
  
  Pai2gRho <- 1
  for(af in af.list){
    
    Pai2gRho <- Pai2gRho * ( 1 - af)
  }
  
  Pai2gRho <- ( 1 - Pai2gRho )^2
  
  return(Pai2gRho) 
  
} ## end of CalPai2gRho() function

## Pass the population allele frequencies list to the function
## CalPai2gRho(TTN_af), got 3.632859e-05, the same result Java returned



###########################################################################################
##
## step three
## 
## Simulate genotypes for 100k times, get n1, n2, Pai2g Matrix
## Perform Rao's Score Test on this set of data simulated
## Print the Score and P-value
###########################################################################################

simu100kGenotypes <- function(TTN_af, sample.size, variants.count, alpha, beta){
  
  ## initialize a vector of variants af starts at 0
  variants.sim <- rep(0, variants.count)
  
  ## initialize pai2g.sim
  pai2g.sim <- NULL
  
  ## simulate sample.size (100000) genotypes, return the Pai2gs and a vectore of allele counts for each site
  for(i in 1:sample.size){
    
    genotype <- simulateGenotypes(TTN_af, alpha, beta)
    
    pai2g.sim <- c(pai2g.sim, genotype[1])
    genotype <- genotype[-1]
    
    variants.sim <- variants.sim + genotype
    
  }
  
  
  variants.sim <- variants.sim / (sample.size * 2)
  
  ## ttn_pai2g_exp <- 3.6328594930727866E-5  
  ## Pass the simulated variants alleles list to CalPai2gRho() to calculate the Sum{ Pai2g*Rho }
  ttn_pai2g_exp <- CalPai2gRho( variants.sim )
  
  TTN_pai2g.sim <- pai2g.sim 
  Si.sim <- TTN_pai2g.sim - ttn_pai2g_exp
  
  n.sim <- length(Si.sim)
  
  ## Rao's Score Test
  #Score.sim = ( sum(Si.sim))^2 / ( n.sim * var(Si.sim) )
   I_beta <- TTN_pai2g.sim^2 - ttn_pai2g_exp^2
  
   Score.sim <- ( sum(Si.sim))^2 / (sum(I_beta) )
  
  ## Calculate p-values 
  p.value <- (1 - pchisq(Score.sim, df=1)) 
  
  print(paste('Pai2gRho', ttn_pai2g_exp,'score: ', Score.sim, 'Pvalue:', p.value))
  return( p.value )
}





###########################################################################################
##
## step four
## Perform the simulation for 1000 times, see the PValues returned
###########################################################################################


###########################################################################################
##
## non parallel
#
# PValues <- NULL
# sample.size <- 200000
# print(PValues)
#
# pdf(file = "histPvalues0330_breaksample2ksize200k.pdf")
#
#
# for(i in 1:20){
#  print(c('simulateiong #', i) )
#  PValues <- c(PValues, simu100kGenotypes(TTN_af, sample.size, variants.count))
#  
#  if(i%%5 == 0){
#    
#    hist(PValues, breaks = 20, main = paste('circle', i, 'Hist of P-values'), xlab = paste('samples:', length(PValues)) )
#  }
#  
#}
#
# hist(PValues, breaks = 40, main = 'Hist of P-values', xlab = paste('samples:', length(PValues)) )
#
# dev.off()
# mean(PValues < 0.05)
###########################################################################################




###########################################################################################
## when doing parallel, have to use list<-foreach() to catch all the returned values from all loops!


## Try parellel 
library(foreach)
library(doMC)

## set seed
set.seed(2017)

## Create Alpha and Beta
alpha <- -1.922

beta.vector <- (25:35) / 10
beta.vector


sample.size <- 200000

######################
for(beta in beta.vector){

  PValues <- NULL	  
  registerDoMC(64)

  ##############
  ## the first 1000 samples
  list <- foreach( i = 1:1000) %dopar% {
    
    print(c('simulating: ', i))
    PValues <- c(PValues, simu100kGenotypes(TTN_af, sample.size, variants.count, alpha, beta) )
    
  } ## end parallel foreach() loop
  
  PValues <- c(PValues, unlist(list) )
  
  ## plot Hist into a PDF document
  pdf(file = paste('histPvalues0330Power_sample1ksize200k_beta', beta, '.pdf', sep = '') )
  
  hist(	PValues, 
	breaks = 40, 	
	main = paste('Hist of P-values, breaks 40','alpha=', alpha,'beta=', beta), 
	xlab = paste('samples:', length(PValues), 'Pvalues<0.05', mean(PValues<0.05) )	
	)
	
  hist(PValues, breaks = 20, main = 'Hist of P-values, breaks 20', xlab = paste('samples:', length(PValues)) ) 
  
  
  ## plot the qq uniform fit
  y <- qunif(ppoints( length(PValues) ) )
  qqplot(PValues, y, 
         main = "qqPlot of uniform distribution")
  
  ## finish the pdf ploting
  dev.off() 
  
  ## Print a 'table' of Pvalues we just simulated
  print(PValues)
  
  mean(PValues < 0.05) 
  
  
  
}

########################################################################################################
########################################################################################################

