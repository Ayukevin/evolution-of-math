## ==============================
## =  種化 Speciation 課堂練習  =
## ==============================
## ==============================
## = 一、合子後隔離、單基因座、 =
## ==============================
## = 二、合子後隔離、雙基因座、 =
## ==============================
## =  三、既存的完整合子後隔離  =
## ==============================
## =       四、少數方優勢       =
## ==============================
## =        五、地理分隔        =
## ==============================
## =     六、適應各自棲息地     =
## ===================================
## =七、部分合子後隔離+部分交配前隔離=
## ===================================


## ==============================
## = 一、合子後隔離、單基因座、 =
## ==============================
## =    突變、天擇、隨機交配    =
## ==============================
rm(list = ls())
## == parameters == 
s = 1
mu = 0.01
gene = c(1,0,0) #AA, Aa, aa
fit = c(1,0,1)
tMAX =100

## == function == 
mutation =  function(gene,mu){
  AA = (1-mu)*gene[1]+(mu/2)*gene[2]
  Aa = (1-mu)*gene[2]+mu*(gene[1]+gene[3])
  aa = (1-mu)*gene[3]+(mu/2)*gene[2]
  return(c(AA, Aa, aa))
}

natural_selection =function(gene){
  gene*fit / sum(gene*fit)
}

random_selection =function(gene){
  A = gene[1] + gene[2]/2
  a = gene[3] + gene[2]/2
  
  AA = A^2
  Aa = 2*A*a
  aa = a^2
  
  c(AA,Aa,aa)
}

## == main == 
data = matrix(0,nrow = tMAX+1, ncol = length(gene))
data[1,] = gene
for (t in 1:tMAX){
  a = mutation(gene,mu)
  b = natural_selection(a,fit)
  gene = random_selection(b)
  data[t+1,] =gene
}


## == plot == 
data[tMAX+1,]
matplot(data,type = 'l')
legend('topright',c(AB, Ab, aB, ab),lty=1:4,col=1:4)

## ==============================
## = 二、合子後隔離、雙基因座、 =
## ==============================
## =    突變、天擇、隨機交配    =
## ==============================
rm(list = ls())
## == parameters == 
s = 1
mu = 0.01
gene = c(1,0,0,0) #AB, Ab, aB, ab
fit = c(1,1-s,1,1)
tMAX =100

## == function ==
mutation =  function(gene,mu){
  AB = (1-mu)*gene[1]+(mu/2)*(gene[2]+gene[3])
  Ab = (1-mu)*gene[2]+(mu/2)*(gene[1]+gene[4])
  aB = (1-mu)*gene[3]+(mu/2)*(gene[1]+gene[4])
  ab = (1-mu)*gene[4]+(mu/2)*(gene[2]+gene[3])
  return(c(AB, Ab, aB, ab))
}

natural_selection = function(gene,fit){
  gene*fit / sum(gene*fit)
}

random_selection = function(gene){
  
  D = gene[1]*gene[4] +gene[2]*gene[3]
  AB = gene[1] - D/2
  Ab = gene[2] + D/2
  aB = gene[3] + D/2
  ab = gene[4] - D/2
  c(AB, Ab, aB, ab)
}


## == main == 
data = matrix(0,nrow = tMAX+1, ncol = length(gene))
data[1,] = gene
for (t in 1:tMAX){
  a = mutation(gene,mu)
  b = natural_selection(a,fit)
  gene = random_selection(b)
  data[t+1,] =gene
}


## == plot ==  
data[tMAX+1,]
matplot(data,type = 'l')
legend('topright',c(AB, Ab, aB, ab),lty=1:4,col=1:4)


## ==============================
## =  三、既存的完整合子後隔離  =
## ==============================
## =       天擇、隨機交配       =
## ==============================
## =      比較基因起始狀態      =
## ==============================
rm(list = ls())
## == parameters == 
s = 1
mu = 0.01
gene = c(0.51,0,0,0.49) #AB, Ab, aB, ab
fit = c(1,1-s,1-s,1)
tMAX =100

## == function ==
natural_selection = function(gene,fit){
  gene*fit / sum(gene*fit)
}

random_selection = function(gene){
  D2 = (gene[1]*gene[4] +gene[2]*gene[3])/2
  gene + c(-1,1,1,-1)*D2
}


## == main == 
data = matrix(0,nrow = tMAX+1, ncol = length(gene),
              dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))
data[1,] = gene
for (t in 1:tMAX){
  a = natural_selection(gene,fit)
  gene = random_selection(a)
  data[t+1,] = gene
}

## == plot ==  
data[tMAX+1,]
matplot(data,type = 'l')
legend('topright',c('AB', 'Ab', 'aB', 'ab'),lty=1:4,col=1:4)

## ==============================
## =       四、少數方優勢       =
## ==============================
## =       天擇、隨機交配       =
## ==============================
## =       比較適應度參數       =
## ==============================
rm(list = ls())
## == parameters == 
Sh = 1
Sf = 0.6 #0.3 #天擇的頻率強度
mu = 0.01
gene = c(0.51,0,0,0.49) #AB, Ab, aB, ab
fit = c(1-gene[1]*Sf,1-Sh,1-Sh,1-gene[4]*Sf)
tMAX =100

## == function ==
natural_selection = function(gene,fit){
  gene*fit / sum(gene*fit)
}

random_selection = function(gene){
  D2 = (gene[1]*gene[4] +gene[2]*gene[3])/2
  gene + c(-1,1,1,-1)*D2
}

## == main == 
data = matrix(0,nrow = tMAX+1, ncol = length(gene),
              dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))
data[1,] = gene
for (t in 1:tMAX){
  a = natural_selection(gene,fit)
  gene = random_selection(a)
  data[t+1,] = gene
  fit = c(1-gene[1]*Sf,1-Sh,1-Sh,1-gene[4]*Sf)
}


## == plot ==  
data[tMAX+1,]

matplot(data,type = 'l')
legend('topright',c('AB', 'Ab', 'aB', 'ab'),lty=1:4,col=1:4)

## ==============================
## =        五、地理分隔        =
## ==============================
## =  少量遷徙、天擇、隨機交配  =
## ==============================
## =         比較遷徙比例       =
## ==============================
rm(list = ls())
## == parameters == 
Sh =1
m =0.1 #0.05
Place1 = c(0.99,0,0,0.01) #AB, Ab, aB, ab
Place2 = c(0.02,0,0,0.98) #AB, Ab, aB, ab
fit = c(1,1-Sh,1-Sh,1)
tMAX =100

## == function ==
migration = function(p1,p2){ #place1,placee2
  place1 = (1-m)*p1 + m*p2
  place2 = (1-m)*p2 + m*p1
  
  c(place1,place2)
}

natural_selection = function(gene,fit){
  gene*fit / sum(gene*fit)
}

random_selection = function(gene){
  D2 = (gene[1]*gene[4] +gene[2]*gene[3])/2
  gene + c(-1,1,1,-1)*D2
}

## == main == 
data1 = matrix(0,nrow = tMAX+1, ncol = length(Place1),
               dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))
data2 = matrix(0,nrow = tMAX+1, ncol = length(Place2),
               dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))

data1[1,] = Place1
data2[1,] = Place2


for (t in 1:tMAX){
  place = migration(Place1,Place2)
  
  a = natural_selection(place[1:4],fit)
  Place1 = random_selection(a)
  data1[t+1,] = Place1
  
  a = natural_selection(place[5:8],fit)
  Place2 = random_selection(a)
  data2[t+1,] = Place2
}


## == plot ==  
data1[tMAX+1,]
data2[tMAX+1,]

par(mfrow = c(1, 2))
matplot(data1,type = 'l')
legend('topright',c('AB', 'Ab', 'aB', 'ab'),lty=1:4,col=1:4)

matplot(data2,type = 'l')
legend('topright',c('AB', 'Ab', 'aB', 'ab'),lty=1:4,col=1:4)

## ==============================
## =     六、適應各自棲息地     =
## ==============================
## =  少量遷徙、天擇、隨機交配  =
## ==============================
## =       比較適應度參數       =
## ==============================
rm(list = ls())
## == parameters == 
Sh =1
Sd = 0.01 #0.01
m =0.1 
Place1 = c(0.99,0,0,0.01) #AB, Ab, aB, ab
Place2 = c(0.02,0,0,0.98) #AB, Ab, aB, ab
fit1 = 1+c(1,-1,-1,-1)*c(Sd,Sh,Sh,Sd)
fit2 = 1+c(-1,-1,-1,1)*c(Sd,Sh,Sh,Sd)
tMAX =100

## == function ==
migration = function(p1,p2){ #place1,placee2
  place1 = (1-m)*p1 + m*p2
  place2 = (1-m)*p2 + m*p1
  
  c(place1,place2)
}

natural_selection = function(gene,fit){
  gene*fit / sum(gene*fit)
}

random_selection = function(gene){
  D2 = (gene[1]*gene[4] +gene[2]*gene[3])/2
  gene + c(-1,1,1,-1)*D2
}

## == main == 
data1 = matrix(0,nrow = tMAX+1, ncol = length(Place1),
               dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))
data2 = matrix(0,nrow = tMAX+1, ncol = length(Place2),
               dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))

data1[1,] = Place1
data2[1,] = Place2

for (t in 1:tMAX){
  place = migration(Place1,Place2)
  
  a = natural_selection(place[1:4],fit1)
  Place1 = random_selection(a)
  data1[t+1,] = Place1
  
  a = natural_selection(place[5:8],fit2)
  Place2 = random_selection(a)
  data2[t+1,] = Place2
}

## == plot ==  
data1[tMAX+1,]
data2[tMAX+1,]

par(mfrow = c(1, 2))
matplot(data1,type = 'l')
legend('right',c('AB', 'Ab', 'aB', 'ab'),lty=1:4,col=1:4)

matplot(data2,type = 'l')
legend('right',c('AB', 'Ab', 'aB', 'ab'),lty=1:4,col=1:4)

## ===================================
## =七、部分合子後隔離+部分交配前隔離=
## ===================================
## =  遷徙、天擇、各自棲息地內擇偶   =
## ===================================
## =        比較擇偶偏好度參數       =
## ===================================
rm(list = ls())
## == parameters == 
Sh =0.3
Sd = 0.01
m =0.05 
alpha = 0 
r = 0.05
Place1 = c(0.99,0,0,0.01) #AB, Ab, aB, ab
Place2 = c(0.02,0,0,0.98) #AB, Ab, aB, ab
fit1 = 1+c(1,-1,-1,-1)*c(Sd,Sh,Sh,Sd)
fit2 = 1+c(-1,-1,-1,1)*c(Sd,Sh,Sh,Sd)
preference_matrix = matrix(c(1+alpha,1+alpha,1,1,
                             1+alpha,1+alpha,1,1,
                             1,1,1+alpha,1+alpha,
                             1,1,1+alpha,1+alpha),nrow = 4, ncol=4)
tMAX =1000

## == function ==
migration = function(p1,p2){ #place1,placee2
  place1 = (1-m)*p1 + m*p2
  place2 = (1-m)*p2 + m*p1
  
  c(place1,place2)
}

natural_selection = function(gene,fit){
  gene*fit / sum(gene*fit)
}

sexual_selection = function(gene,preference_matrix){
  
  pair_matrix = outer(gene, gene) * preference_matrix
  pair_matrix = pair_matrix / sum(pair_matrix) 
  
  children <- rep(0, 4)
  
  for (mom in 1:4) {
    for (dad in 1:4) {
      if (mom == dad) { # 父母基因型相同
        children[mom] <- children[mom] + pair_matrix[mom, dad]
      } else if (mom + dad != 5) { # 父母基因型部分相同
        children[c(mom, dad)] <- children[c(mom, dad)] + pair_matrix[mom, dad] / 2
      } else { # 父母基因型完全相反
        children[c(mom, dad)] <- children[c(mom, dad)] + pair_matrix[mom, dad] * (1 - r) / 2
        other <- if (min(mom, dad) == 1) c(2, 3) else c(1, 4)
        children[other] <- children[other] + pair_matrix[mom, dad] * r / 2
      }
    }
  }
  
  children / sum(children)
}

## == main == 
data1 = matrix(0,nrow = tMAX+1, ncol = length(Place1),
               dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))
data2 = matrix(0,nrow = tMAX+1, ncol = length(Place2),
               dimnames = list(rep(1:(tMAX+1)),c('AB', 'Ab', 'aB', 'ab')))

data1[1,] = Place1
data2[1,] = Place2

for (t in 1:tMAX){
  place = migration(Place1,Place2)
  
  a = natural_selection(place[1:4],fit1)
  Place1 = sexual_selection(a,preference_matrix)
  data1[t+1,] = Place1
  
  a = natural_selection(place[5:8],fit2)
  Place2 = sexual_selection(a,preference_matrix)
  data2[t+1,] = Place2
}

## == plot ==  
data1[tMAX+1,]
data2[tMAX+1,]

par(mfrow = c(1,2))
matplot(data1,type = 'l',xlab = 'tMAX', main = paste('Alpha =', alpha))
legend('topright',c('AB', 'Ab', 'aB', 'ab'),lty=1:4,col=1:4)

matplot(data2,type = 'l',xlab = 'tMAX')