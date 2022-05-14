library(simecol)
library(tidyverse)

#Initial values
yrs <- 10 #Simulation time period
risk_init <- 60000 #At risk of homelessness
hmlss_init <- 3000 #Non-chronic homeless
chron_init <- 1000 #Chronically homeless
trans_init <- 2000 #Initial number of STRA/transitional housing units
perm_init <- 5000 #Initial number of PSH units
trans_add <- 2500 #Number of STRA/transitional housing units to add over time period
perm_add <- 500 #Number of PSH units to add over time period

#Model parameters
a <- 1.015 #Annual population growth
v <- 0.02 #Vacancy rate (between 0 and 0.1)
b1 <- 0.8 #Multiply with (0.1 - vacancy rate) to get effect on entering homelessness
b2 <- 0.5 #Multiply with vacancy rate to get effect on exiting homelessness
c <- 0.1 #Percent of STRA/Trans set aside for chronically homeless
d <- 0.8 #Percent of PSH set aside for chronically homeless

h <- new("odeModel",
  main = function (time, init, parms, ... ) {
    t <- time
    x <- init
    p <- parms
    
    # Addition of housing over time
    # Modeled as outflows from Trans and Perm, but not added to Risk
    t_input <- approxTime1(inputs, time, rule=2)["t_input"]
    p_input <- approxTime1(inputs, time, rule=2)["p_input"]
    
    e1 <- risk_init * p["a"]^t - risk_init
    e2 <- x[1] * p["b1"] * (.1 - p["v"]) - x[2] * p["b2"] * p["v"]
    e8 <- x[4] * p["b2"] + t_input #Addition of housing exerts "pull" on homelessness
    e9 <- x[5] * p["b2"] + p_input #Addition of housing exerts "pull" on homelessness
    e7 <- min(c(e9 * p["d"], x[3]))
    e6 <- min(c(e8 * p["c"], x[3] - e7))
    e4 <- e8 - e6
    e5 <- e9 - e7
    avglos <- 1/(e4 + e5) #Inverse of outflow is average length of stay
    scale <- 4 #Scale parameter for gamma distribution
    e3 <- x[2] * (1 - pgamma(1, shape = avglos/scale, scale = scale)) #Approximate flow over 1 year of homelessness
    
    Risk <- e1 + e8 + e9 - e2 - t_input - p_input #Housing inputs only "pull" on homelessness; do not add to Risk
    Homeless <- e2 - e3 - e4 - e5
    Chronic <- e3 - e6 - e7
    Trans <- e4 + e6 - e8
    Perm <- e5 + e7 - e9
    
    list(c(Risk, Homeless, Chronic, Trans, Perm))
  },
  parms = c(a=a, #Annual population growth
            v=v, #Vacancy rate (between 0 and 0.1)
            b1=b1, #Multiply with 0.1 - vacancy rate to get effect on entering homelessness
            b2=b2, #Multiply with vacancy rate to get effect on exiting homelessness
            c=c, #Percent of STRA/Trans set aside for chronically homeless
            d=d), #Percent of PSH set aside for chronically homeless
  times = c(from=0, to=yrs, by=0.01),
  init = c(Risk=risk_init, Homeless=hmlss_init, Chronic=chron_init, Trans=trans_init, Perm=perm_init),
  solver = "rk4"
)

#Additional housing to add over time
inputs(h) <- as.matrix(data.frame(
  time = 0:yrs,
  t_input = seq(0, trans_add, length.out=length(0:yrs)), #Transitional housing
  p_input = seq(0, perm_add, length.out=length(0:yrs)) #Permanent housing
))

#Run simulation
h <- sim(h)

#Plot homelessness values
out(h) %>%
  gather("type", "value", -time) %>%
  filter(type %in% c("Chronic", "Homeless")) %>%
  ggplot(aes(x=time, y=value, color=type)) +
  geom_line() +
  scale_y_continuous(limits=c(0, 10000))
