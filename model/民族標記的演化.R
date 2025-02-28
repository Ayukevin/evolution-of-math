dev.off() 
#payoff-biased 
b =0.7 
a =0.9
r =0.1
m =0.01
tMAX = 1000
Place1 = c(0.55,0.9,0) #pa,px,D
Place2 = c(0.46,0.05,0) #pa,px,D

decode = function(Place){
  p = Place[1]
  q = Place[2]
  D = Place[3]
  
  gene = c(p*q,p*(1-q),(1-p)*q,(1-p)*(1-q)) + c(1,-1,-1,1)*D
  return(gene) # pax, pay, pbx, pby
}


migration = function(p1,p2){ #place1,place2
  place1 = (1-m)*p1 + m*p2
  place2 = (1-m)*p2 + m*p1
  
  c(place1,place2)
}

payoff_biased = function(p){
  pa = p[1]+p[3]
  px = p[1]+p[2]
  D = p[1]*p[4] - p[3]*p[2]
  
  w  = 1+ b*(c(px,1-px,px,1-px)+a*D/c(pa,-pa,-(1-pa),(1-pa)))
  w[is.na(w)] =1
  p*w /sum(p*w)
}

reorder = function(p,r){
  D = p[1]*p[4] - p[3]*p[2]
  p + c(-1,1,1,-1)*r*D
}

encode =function(p){
  pa = p[1]+p[3]
  px = p[1]+p[2]
  D = p[1]*p[4] - p[3]*p[2]
  
  c(pa,px,D)
}


#main
data1 = matrix(0,nrow = tMAX+1, ncol = length(Place1),
               dimnames = list(rep(1:(tMAX+1)),c('pa','px','D')))
data2 = matrix(0,nrow = tMAX+1, ncol = length(Place2),
               dimnames = list(rep(1:(tMAX+1)),c('pa','px','D')))


data1[1,] = Place1
data2[1,] = Place2

for (t in 1:tMAX){
  Place1 = decode(Place1)
  Place2 = decode(Place2)
  place = migration(Place1,Place2)
  
  temp1 = payoff_biased(place[1:4])
  temp2 = reorder(temp1,r)
  Place1 = encode(temp2)
  data1[t+1,] = Place1
  
  temp1 = payoff_biased(place[5:8])
  temp2 = reorder(temp1,r)
  Place2 = encode(temp2)
  data2[t+1,] = Place2
}

#plot
data1[tMAX+1,]
data2[tMAX+1,]

par(mfrow = c(1, 2))
matplot(data1,type = 'l',xlab = 'tMAX')
legend('topright',c('pa','px','D'),lty=1:3,col=1:3)

matplot(data2,type = 'l',xlab = 'tMAX')
legend('topright',c('pa','px','D'),lty=1:3,col=1:3)
