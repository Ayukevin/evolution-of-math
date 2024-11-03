library(deSolve)

SIR = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSdt = -beta*I*S
    dIdt = (beta*S-gamma)*I
    dRdt = gamma*I
    list(c(dSdt, dIdt, dRdt))
  })
}

state=c(S=0.99,I=0.01,R=0)
parameters = c(beta=0.3,gamma=0.1)
time = seq(0,500,0.1)

solution = ode(y=state,times = time,
               func=SIR, parms = parameters)

matplot(time,solution[,-1],type='l')
legend('topright',colnames(solution)[-1],
       lty =1:3,col = 1:3)
    
