## Simple homeless system model
## Using stochastic differential equations

## Model contains four compartments
## U = unsheltered
## S = temporary shelter (e.g., emergency shelter, transitional housing)
## P = permanent housing (e.g., rapid rehousing, permanent supportive housing)
## O = other housing in the community

rm(list = ls(all = TRUE))  # resets R to fresh

library(ggplot2)

set.seed(126189)

niter <- 50

#### Deterministic model function ####
model_func=function(t, x, vparameters){
  U = x[1]  
  S = x[2]  
  P = x[3]  
  O = x[4]  
  
  with(as.list(vparameters),{
    dU = alpha - beta_u*U - beta_s*S - beta_p*P            
    dS = 0 #fixed capacity
    dP = 0 #fixed capacity
    dO = beta_u*U + beta_s*S + beta_p*P
    out = c(dU, dS, dP, dO)
    list(out)
  })
}


#### Initial conditions ####
U_0 <- 1000
S_0 <- 150
P_0 <- 350
O_0 <- 0

#### Parameters ####
# Because shelter and housing have fixed capacities, their inflows/outflows cancel out
# As a result, only the system inflow and system outflow need to be modeled
# Note, however, that outflow is proportional to the conditions of U, S, and P

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

#### Run deterministic model ####
tend <- 10
vt <- seq(0,tend,1/12)
vparameters = c(alpha = alpha, beta_u = beta_u, beta_s = beta_s, beta_p = beta_p)
inits <- c(U = U_0, S = S_0, P = P_0, O = O_0)
dmodel <- as.data.frame(lsoda(inits, vt, model_func, vparameters))

ggplot(dmodel, aes(x = time, y = U)) + 
  geom_line() +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  expand_limits(y = 0)
