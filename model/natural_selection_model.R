## Agent-based simulation of natural selection on QTL
## Model assumptions: free recombination, random pairing
rm(list=ls())
## === Settings ===
nAgent <- 300
nLoci <- 10
maxgen <- 500
optimal <- 8
alleles <- c("A", "a")
aPairs <- c("AA","Aa","aa")

## === Custom Functions ===
phenotype <- function(x){
  # x is a matrix of genotypes. Each row is an agent, each column is a locus
  nAA <- rowSums(x=="AA") # yes=1, no=0. rowSums() = number of AA in each row
  nAa <- rowSums(x=="Aa")
  naa <- rowSums(x=="aa")
  nAA*1 + nAa*0.5 + naa*0
}
test <- matrix(c("AA","AA","AA", # 1st row, 3 AA -> 3
                 "AA","Aa","AA", # 2nd row, 2 AA & 1 Aa -> 2.5
                 "aa","aa","aa"),nrow=3,byrow=T) # 3rd row, only aa -> 0
phenotype(test) # should return 3, 2.5, 0

fitness <- function(x){
  # x is a vector of phenotypes
  1 - abs(x - optimal)/nLoci
}
plot(x=0:nLoci, y=fitness(0:nLoci), type='b')

meiosis <- function(x){
  # x is a matrix of genotypes
  y <- x
  y[y=="AA"] <- "A" # replace all "AA" with "A"
  y[y=="aa"] <- "a"
  nAa <- sum(y=="Aa")
  y[y=="Aa"] <- sample(alleles,nAa,replace=T)
  y
}
meiosis(rep("AA",10)) # should return 10 "A"
meiosis(rep("aa",10)) # should return 10 "a"
meiosis(rep("Aa",20)) # should return both types at similar rate

## === Initialize agents ===
# each row is one agent; each column is one diploid locus
genotypes <- matrix(sample(aPairs, nAgent*nLoci, replace=T), nrow=nAgent, ncol=nLoci, byrow=T)

## place to record data
data <- data.frame(mTrait=NA,sdTrait=NA,meanFit=NA,nSurvive=NA,het=NA)
## Iterate through life cycles
for(i in 1:maxgen){
  ## === Phenotypes ===
  trait <- phenotype(genotypes)
  
  ## === Natural selection ===
  fit <- fitness(trait)
  survive <- runif(nAgent) < fit
  # check histogram to see if actual survival matches trait value
  #hist(trait)
  #hist(trait[survive],col=2,add=T)
  
  survivor <- genotypes[survive,]
  nSurvive <- nrow(survivor)
  
  ## === Reproduction ===
  mom <- survivor[1:floor(nSurvive/2),]
  dad <- survivor[ceiling(nSurvive/2):nSurvive,]
  
  # if nSurvive = 5 -> mom = 1:2.5 = 1,2 ; dad = 5:3.5 = 5,4
  # if nSurvive = 6 -> mom = 1:3 ; dad = 6:4
  # one individual in the middle is ignored if odd number.
  
  ## -- shuffle parents and pair up --
  # to keep population size constant, we sample "nAgent" pairs of parents
  x <- sample(nrow(mom),nAgent,replace=T)
  y <- sample(nrow(dad),nAgent,replace=T)
  mom <- mom[x,]
  dad <- dad[y,]
  
  ## -- make gametes --
  egg <- meiosis(mom)
  sperm <- meiosis(dad)
  
  offspring <- paste(egg, sperm, sep='')
  offspring <- matrix(offspring,nrow=nAgent)
  offspring[offspring=="aA"] <- "Aa"
  
  ## record data
  nAa <- sum(x == "Aa")
  y[x == "Aa"] <- sample(alleles, nAa * 2, replace=T)
  
  ## go to next iteration
  genotypes <- offspring
}

## === Plot results ===
layout(1:3)
## mean trait value with SD
par(mar=c(2,2,3,1))
plot(data$mTrait,type='l',ylim=c(0,nLoci),main='Trait (population mean +- sd)')
lines(data$mTrait + data$sdTrait, lty=3)
lines(data$mTrait - data$sdTrait, lty=3)
abline(h=optimal,lty=3,col="red")
## mean fitness and actual survive count
plot(data$nSurvive/nAgent,type='l',col=3,ylim=c(0,1),main="Mean fitness")
lines(data$meanFit)
legend("bottomright",c("mean fitness","actual survival rate"),col=c(1,3),lty=1)

## heterozygosity, amount of "Aa" in the population
plot(data$het,type='l',ylim=c(0,1),main="Heterozygosity (amount of Aa)")

