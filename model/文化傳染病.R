library(deSolve)

a=1/3
b=1/2
c=3/4
f=60
tmax = 50
N = rep(a,tmax+1)


data = c(0,0,f,
         a,0,0,
         0,b,c)
A = matrix(data = data,nrow = 3)

x = c(10,10,10)
result = matrix(0,nrow = tmax+1,ncol = 3)
result[1,] = x


for (i in 1:tmax){
  result[i+1,] = A %*% result[i,]
}

matplot(c(1:51),result[,1:3],type='l')
matplot(c(1:51),result[,1:3],log='y')


ev = eigen(A)
lamda = ev$values[1]
vector = as.numeric(ev$vectors[,1])

layout(matrix(1:2,1))
barplot(n)
barplot(comp)



