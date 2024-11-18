library(deSolve)

## ==setting==
a= 0.8
b= 1.1
Kx = 80
Ky = 80
Rx = 0.9
Ry = 1.2


## ==grid==
x0 = seq(0,100,length.out = 20)
y0 = seq(0,100,length.out = 20)
x0s = rep(x0,length(y0))
y0s = rep(y0,each=length(x0))

## vector
dxdt = (1-(x0s+a*y0s)/Kx)*Rx*x0s 
dydt = (1-(y0s+b*x0s)/Ky)*Ry*y0s

## convert 
x1s = x0s + dxdt*0.01
y1s = y0s + dydt*0.01

##draw the field
plot(NA, xlim = range(x0), ylim = range(y0), main = "Vector Field")
arrows(x0s, y0s, x1s, y1s, length = 0.1) 

## formula
competition = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dxdt = (1-(x0s+a*y0s)/Kx)*Rx*x0s 
    dydt = (1-(y0s+b*x0s)/Ky)*Ry*y0s
    list(c(dxdt,dydt))
  })
}

state = c(x0s,y0s)
parameter = c(a,b,Rx,Ry,Kx,Ky)

t = seq(0,500,0.1)

solution = ode(y=state,times=t,
               func=competition,
               parms=parameter)

# time series
matplot(t,solution[,-1])
plot(solution[,-1])

