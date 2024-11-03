library(deSolve)

SIR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSdt <- -beta * I * S
    dIdt <- (beta * S - gamma) * I
    D <- 1 - S - I
    list(c(dSdt, dIdt, D))
  })
}


state <- c(S = 0.99, I = 0.01, D = 0)
beta <- 0.3
time <- seq(0, 100, 0.1)

#drawing plot
run_simulation <- function(gamma_values) {
  plot_colors <- rainbow(length(gamma_values))
  
  #setting 
  plot(NULL, xlim = range(time), ylim = c(0, 5),
       xlab = "Time", ylab = "D",
       main = "Relation between gamma and D")
  
  #adding lines
  for (i in seq_along(gamma_values)) {
    parameters <- c(beta = beta, gamma = gamma_values[i])
    
    solution <- ode(y = state, times = time,
                    func = SIR, parms = parameters)
    lines(time, solution[,"D"], col = plot_colors[i], lty = 1)
  }
  
  legend("bottomright", legend = paste("gamma =", gamma_values),
         col = plot_colors, lty = 1)
}


#main
gamma_values <- seq(0.1, 0.9, by = 0.1)
run_simulation(gamma_values)