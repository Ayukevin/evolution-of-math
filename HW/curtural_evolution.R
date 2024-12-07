# parameter
p = c(0.25,0.25,0.25,0.25) # pHA, pHB, pSA ,pSB

a = 0.5 # 0 ≤ a ≤ 1
rA = 0.8
rB = 0.3# 0 ≤ rB < rA ≤ 1
d = 0.2 # 0 ≤ d ≤ 1
c = 0.4 # 0 ≤ c ≤ 1

tMAX = 100

#pA & pH transform
pA_pH = function(p){
  
  pH = p[1]+p[2]
  pA = p[1]+p[3]
  D = p[1]*p[4] - p[2]*p[3]
  
  c(pH,pA,D)
}

#生病和康復
sick_recover = function(p){
  
  pH = p[1]
  pA = p[2]
  D = p[3]
  
  pSA = (1-pH)*pA-D
  pSB = (1-pH)*(1-pA)+D
  
  next_pH = (1-a)*pH+ rA*pSA+ rB*pSB
  next_pA = pA
  next_D  = (1-a)*D+ rA*(1-pA)*pSA +rB*pH*(1-pH)*(1-pA)+D
  
  c(next_pH, next_pA, next_D)  
}

#存活或死亡
live_death =function(p){
  pH = p[1]
  pA = p[2]
  D = p[3]
  
  pSA = (1-pH)*pA-D
  
  next_pH = pH/ (pH + (1-d)*pH)
  next_pA = (pA - d*pSA)/ (pH + (1-d)*pH)
  next_D  = (1-d)*D/ (pH + (1-d)*pH)
  
  c(next_pH, next_pA, next_D) 
} 

#學習病人用的療法
treatment_change = function(p){
  pH = p[1]
  pA = p[2]
  D = p[3]
  
  next_pH = pH- 2*D*c
  next_pA = pA- D*c
  next_D  = D - D*c*(1-pH)
  
  c(next_pH, next_pA, next_D) 
}



#main
data = matrix(0,nrow = tMAX+1, ncol = 3, #num here
               dimnames = list(rep(1:(tMAX+1)), c('pH','pA','D')))
p = pA_pH(p)
data[1,] = pA_pH(p)

for (t in 1:tMAX){
  patient1 = sick_recover(p)
  patient2 = live_death(patient1)
  p = treatment_change(patient2)
  
  data[t+1,] = p
}


data[tMAX+1,]

matplot(data,type = 'l',xlab = 'tMAX' ) # ,main= paste('var'=var())
legend('topright',c('pH','pA','D'),lty=1:3,col=1:3)