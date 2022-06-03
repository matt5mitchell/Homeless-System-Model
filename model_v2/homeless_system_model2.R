## Simple homeless system model
## Using stochastic differential equations

# Model contains four compartments
# U = unsheltered
# S = temporary shelter (e.g., emergency shelter, transitional housing)
# P = permanent housing (e.g., rapid rehousing, permanent supportive housing)
# O = other housing in the community
# Because S and P have fixed capacities, their inflows/outflows cancel out 
# (except for white noise, which is scaled according to the outflows from S and P).
# As a result, only the dynamics of system inflow and system outflow are modeled,
# with system outflow being proportional to the conditions of U, S, and P.
# White noise is added to these dynamics, scaled according to their flows.

# Model makes simplified assumptions about adding permanent housing units:
# - Lease up is happens within a single time step
# - Added permanent housing reduces unsheltered homelessness immediately
# - Utilization is 100%

rm(list = ls(all = TRUE))  # resets R to fresh

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(78324)
tend <- 5          # years to model
delta_t <- 1/12    # time steps
niter <- 100       # iterations of stochastic model
stochastic <- TRUE # TRUE for stochastic, FALSE for deterministic


#### Initial conditions ####
U_0 <- 1000 #unsheltered households
S_0 <- 150  #households in shelter/temporary housing
P_0 <- 350  #households in permanent housing within homelessness response system
O_0 <- 0    #households placed in other permanent housing (outside of system)

#### Additional units ####
# Added units are net changes at a point in time, not cumulative
s_new <- data.frame(units = 50, #new shelter units added
                    year = 2)
p_new <- data.frame(units = c(50,50), #new permanent housing units added
                    year = c(1,7)) #year new units are added

#### Model parameters ####

## Inflow
u_in <- 300 #fixed inflow rate per year

## Outflow 
# Rates of flow in years^{-1} (inverse of avg LOS)
u_out <- 365/180
s_out <- 365/90
p_out <- 365/720

# Percent of exits to other community housing
u_pct_o <- 0.1
s_pct_o <- 0.2
p_pct_o <- 0.6

# Final outflow rates
rate_u_o <- u_out * u_pct_o
rate_s_o <- s_out * s_pct_o
rate_p_o <- p_out * p_pct_o

## Fixed capacity compartments
# Inflow = Outflow, except for white noise
# Note that outflow is split between terms above and below
rate_s_in <- s_out
rate_p_in <- p_out
rate_s_x <- s_out - rate_s_o # portion of flow not already accounted for
rate_p_x <- p_out - rate_p_o # portion of flow not already accounted for


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
  # First four state changes are the net inflow and outflow of the system
  # Second four state changes are for the white noise for fixed capacity resources
  K <- length(vstate)  # number of compartments
  J <- 8               # number of possible state changes
  lambda <- matrix(0,nrow=J,ncol=K)
  lambda[1,] <- c( 1, 0, 0, 0) # inflow into U
  lambda[2,] <- c(-1, 0, 0, 1) # (partial) outflow from U, drives net flow to O
  lambda[3,] <- c(-1,-1, 0, 1) # (partial) outflow from S, drives net flow to O
  lambda[4,] <- c(-1, 0,-1, 1) # (partial) outflow from P, drives net flow to O
  lambda[5,] <- c( 0, 1, 0, 0) # inflow into S
  lambda[6,] <- c( 0,-1, 0, 0) # (partial) outflow from S
  lambda[7,] <- c( 0, 0, 1, 0) # inflow into P
  lambda[8,] <- c( 0, 0,-1, 0) # (partial) outflow from P
  
  while(vstate[1] > 0 & time <= tend) {
    
    time <- round(time / delta_t, 0) * delta_t # fix rounding errors
    
    zstate[[i]] <- c(vstate, time, iter) # current state of system
    
    U <- vstate[1]
    S <- vstate[2]
    P <- vstate[3]
    O <- vstate[4]
    
    vec_p <- c(u_in, 
              rate_u_o*U,
              rate_s_o*S,
              rate_p_o*P,
              rate_s_in*S,
              rate_s_x*S,
              rate_p_in*P,
              rate_p_x*P)
    
    # G matrix for scaling white noise terms
    G_t <- lambda * sqrt(vec_p) # state changes * sqrt of flow rates
    G <- t(G_t)
    
    W <- rnorm(ncol(G)) * stochastic # TRUE coerced to 1, FALSE coerced to 0
    
    # Stochastic changes
    ## G contains square roots of flow rates (with signs representing state changes)
    ## W contains draws from Wiener process with mean = 0 and sd = 1
    ## Multiplying G * W scales the variance to approximate poisson process with lambda = flow rates
    ## sqrt(lambda) * norm(0, 1) ~= pois(lambda)
    
    delta_U <- delta_t*(u_in - (rate_u_o*U + rate_s_o*S + rate_p_o*P)) + sqrt(delta_t)*sum(G[1,]*W)
    delta_S <- 0                                                       + sqrt(delta_t)*sum(G[2,]*W)
    delta_P <- 0                                                       + sqrt(delta_t)*sum(G[3,]*W)
    delta_O <- delta_t*(rate_u_o*U + rate_s_o*S + rate_p_o*P)          + sqrt(delta_t)*sum(G[4,]*W)
    
    # Capacity change
    delta_S_cap <- max(0, s_new$units[s_new$year == time])
    delta_P_cap <- max(0, p_new$units[p_new$year == time])
    
    vstate[1] <- vstate[1] + delta_U - delta_S_cap - delta_P_cap
    vstate[2] <- vstate[2] + delta_S + delta_S_cap
    vstate[3] <- vstate[3] + delta_P + delta_P_cap 
    vstate[4] <- vstate[4] + delta_O 
    
    vstate[vstate<0] <- 0
    i <- i + 1
    time <- time + delta_t
  }
  cat("Iteration:", iter, niter, " ", time, vstate, "\n")
}

smodel <- data.frame(matrix(unlist(zstate), ncol = 6, byrow = TRUE))
colnames(smodel) <- c("U", "S", "P", "O", "time", "iter")

#### Plot results ####
subtitle <- if(stochastic) {
  paste("Results of", niter, "runs of stochastic model with time steps of", round(delta_t,3), "year")
} else { paste("Results of deterministic model with time steps of", round(delta_t,3), "year")}
alpha <- if(stochastic) 0.1 else 1

ci <- smodel %>%
  group_by(time) %>%
  summarize(lwr = quantile(U, .025)[[1]],
            med = median(U),
            upr = quantile(U, .975)[[1]])

p <- ggplot(smodel, aes(x = time, y = U, group = iter)) +
  geom_line(alpha = alpha) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  expand_limits(y = 0) +
  labs(title = "Projections of Unsheltered Homelessness",
       subtitle = subtitle,
       x = "Years",
       y = "Households") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank())

if(stochastic) {
  p +
    geom_ribbon(data = ci, aes(ymin = lwr, ymax = upr), alpha = alpha) +
    geom_line(data = ci, aes(y = med), color = "indianred")
} else {p}
