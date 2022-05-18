## Simple homeless system model
## Using stochastic differential equations

# Model contains four compartments
# U = unsheltered
# S = temporary shelter (e.g., emergency shelter, transitional housing)
# P = permanent housing (e.g., rapid rehousing, permanent supportive housing)
# O = other housing in the community
# Because shelter and housing have fixed capacities, their inflows/outflows cancel out
# As a result, only the system inflow and system outflow need to be modeled
# Note that outflow is proportional to the conditions of U, S, and P

# However, in the real world exit rates to permanent housing may be note, 
# but not specifically to other housing in the community.
# To address this, the inflows/outflow from S and P can be incorporated into estimates

rm(list = ls(all = TRUE))  # resets R to fresh

library(ggplot2)

set.seed(78324)
tend <- 10         # years to model
delta_t <- 1/12    # time steps
niter <- 100       # iterations of stochastic model
stochastic <- TRUE # TRUE for stochastic, FALSE for deterministic


#### Initial conditions ####
U_0 <- 1000
S_0 <- 150
P_0 <- 350
O_0 <- 0

#### Model parameters ####

## Inflow
alpha <- 300 #fixed inflow rate per year

## Outflow 
# Rates of flow in years^{-1} (inverse of avg LOS)
u_out <- 365/180
s_out <- 365/90
p_out <- 365/720

# Percent of exits to other community housing
u_pct_o <- 0.1
s_pct_o <- 0.2
p_pct_o <- 0.7

# Final outflow parameters
beta_u <- u_out * u_pct_o
beta_s <- s_out * s_pct_o
beta_p <- p_out * p_pct_o


#### Deterministic model ####


# library(deSolve)
# ## Function for deterministic model
# model_func=function(t, x, vparameters){
#   U = x[1]  
#   S = x[2]  
#   P = x[3]  
#   O = x[4]  
#   
#   with(as.list(vparameters),{
#     dU = alpha - beta_u*U - beta_s*S - beta_p*P            
#     dS = 0 #fixed capacity
#     dP = 0 #fixed capacity
#     dO = beta_u*U + beta_s*S + beta_p*P
#     out = c(dU, dS, dP, dO)
#     list(out)
#   })
# }
# 
# ## Run deterministic model
# vt <- seq(0, tend, delta_t)
# vparameters = c(alpha = alpha, beta_u = beta_u, beta_s = beta_s, beta_p = beta_p)
# inits <- c(U = U_0, S = S_0, P = P_0, O = O_0)
# dmodel <- as.data.frame(lsoda(inits, vt, model_func, vparameters))



#### Stochastic Model ####
# Credit for approach and significant chunks of code to: Sherry Towers
# http://sherrytowers.com/2016/02/06/stochastic-compartmental-modelling-with-stochastic-differential-equations-2/

niter2 <- if(stochastic) niter else 1 # FALSE = only one iteration

zstate <- list()
i <- 1L
for (iter in 1:niter2){
  time <- 0
  vstate <- c(U_0, S_0, P_0, O_0)
  
  # State changes: http://sherrytowers.com/2016/01/02/stochastic-compartmental-modelling/#step2
  K <- length(vstate)  # number of compartments
  J <- 4               # number of possible state changes
  lambda <- matrix(0,nrow=J,ncol=length(vstate))
  lambda[1,] <- c(1,0,0,0)
  lambda[2,] <- c(-1,0,0,1)
  lambda[3,] <- c(-1,0,0,1)
  lambda[4,] <- c(-1,0,0,1)
  
  while(vstate[1] > 0 & time < (tend + delta_t)) {
    zstate[[i]] <- c(vstate, time, iter)
    U <- vstate[1]
    S <- vstate[2]
    P <- vstate[3]
    O <- vstate[4]
    vec_p <- c(alpha, 
              beta_u*U,
              beta_s*S,
              beta_p*P)
    
    G_t <- lambda
    for (irow in 1:nrow(G_t)){
      G_t[irow,] <- G_t[irow,]*sqrt(vec_p[irow])
    }
    G <- t(G_t)
    
    W <- rnorm(ncol(G)) * stochastic # TRUE coerced to 1, FALSE coerced to 0
    delta_U <- delta_t*(alpha - beta_u*U - beta_s*S - beta_p*P) + sqrt(delta_t)*sum(G[1,]*W)
    delta_S <- 0 #delta_t*0                                     + sqrt(delta_t)*sum(G[2,]*W)
    delta_P <- 0 #delta_t*0                                     + sqrt(delta_t)*sum(G[3,]*W)
    delta_O <- delta_t*(beta_u*U + beta_s*S + beta_p*P)         + sqrt(delta_t)*sum(G[3,]*W)
    
    vstate[1] <- vstate[1] + delta_U
    vstate[2] <- vstate[2] + delta_S 
    vstate[3] <- vstate[3] + delta_P 
    vstate[4] <- vstate[4] + delta_O 
    
    vstate[vstate<0] <- 0
    i <- i + 1
    time <- time + delta_t
  }
  cat("Doing SDE realisation:",iter,niter," ",time,vstate,"\n")
}

smodel <- data.frame(matrix(unlist(zstate), ncol = 6, byrow = TRUE))
colnames(smodel) <- c("U", "S", "P", "O", "time", "iter")

#### Plot results ####
subtitle <- if(stochastic) {
  paste("Results of", niter, "runs of stochastic model with time steps of", round(delta_t,3), "year")
} else { paste("Results of deterministic model with time steps of", round(delta_t,3), "year")}
alpha <- if(stochastic) 0.1 else 1

ggplot() +
  geom_line(data = smodel, aes(x = time, y = U, group = iter), alpha = alpha) +
  # geom_line(data = dmodel, aes(x = time, y = U), color = "indianred") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  expand_limits(y = 0) +
  labs(title = "Projection of Unsheltered Homelessness",
       subtitle = subtitle,
       x = "Years",
       y = "Unsheltered Homelessness") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank())
