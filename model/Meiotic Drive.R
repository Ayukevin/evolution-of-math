dd = 1000
Dd = 0 
DD = 2
tmax= 50
fertility = c(dd=2,Dd=2,DD=0)
n_dd = c(dd)
n_Dd = c(Dd)
n_DD = c(DD)
total = c(dd+Dd+DD)

for (gen in 1:tmax){
  DD_male = DD/2
  Dd_male = Dd/2
  dd_male = dd/2
  
  DD_female = DD/2
  Dd_female = Dd/2
  dd_female = dd/2
  
  egg_D = (Dd_female/2)*fertility['Dd']
  egg_d = dd_female * fertility["dd"] + (Dd_female/2)*fertility['Dd']
    
  sperm_D = (DD_male + Dd_male)/(DD_male + Dd_male+dd_male)
  sperm_d = dd_male /(DD_male + Dd_male+dd_male)
  
  dd = egg_d * sperm_d
  Dd = egg_d * sperm_D + egg_D * sperm_d
  DD = egg_D * sperm_D
  
  n_dd = c(n_dd,dd)
  n_Dd = c(n_Dd,Dd)
  n_DD = c(n_DD,DD)
  total = c(total,dd+Dd+DD)
}

Numbers = cbind(n_dd,n_Dd,n_DD,total)
colnames(Numbers) <- c("dd", "Dd", "DD", "Total")
time = c(0:tmax)

matplot(time,Numbers,type='l')
legend('topright',legend = colnames(Numbers),col=1:4,lty=1)

matplot(data[,-1],data)