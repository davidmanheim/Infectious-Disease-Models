# From https://rpubs.com/srijana/110753

library (deSolve) 

seir_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  L = state_values [2]        # latent
  I = state_values [3]        # infectious
  K = state_values [4]        # killed (immediately)
  R = state_values [5]        # recovered
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      dS = (-beta * S * I)
      dL = (beta * S * I) - (delta * L)
      dI = (delta * L) - (mu * I) -(gamma * I)
      dK = (mu * I)
      dR = (gamma * I)
      
      # combine results
      results = c (dS, dL, dI, dK, dR)
      list (results)
    }
  )
}

contact_rate = 10                     # number of contacts per day
transmission_probability = 0.05       # transmission probability
infectious_period = 3.3                 # infectious period
latent_period = 2                     # latent period
case_fatality = .95

beta_value = contact_rate * transmission_probability
gamma_value = (1 / infectious_period) - case_fatality/infectious_period
delta_value = 1 / latent_period
mu_value = case_fatality/infectious_period

Ro = beta_value / (gamma_value+mu_value)

#Disease dynamics parameters.
parameter_list = c (beta = beta_value, gamma = gamma_value, delta = delta_value, mu=mu_value)
#Initial values for sub-populations.
W = 9990        # susceptible hosts
X = 10           # infectious hosts
Y = 0           # recovered hosts
Z = 0           # exposed hosts
ZED = 0         # pre-killed hosts

#Compute total population.

N = W + X + Y + Z + ZED
#Initial state values for the differential equations.

initial_values = c (S = W/N, E = X/N, I = Y/N, K = ZED/N, R = Z/N)
#Output timepoints.

timepoints = seq (0, 200, by=0.1)
#Simulate the SEIR epidemic.

output = lsoda (initial_values, timepoints, seir_model, parameter_list)
#Plot dynamics of all sub-populations.

# susceptible hosts over time
plot (S ~ time, data = output, type='l', ylim = c(0,1), col = 'blue', ylab = 'S, E, I, R', main = 'SEIR epidemic') 
# remain on same frame
par (new = TRUE)    
# exposed hosts over time
plot (E ~ time, data = output, type='l', ylim = c(0,1), col = 'grey', ylab = '', axes = FALSE)
# remain on same frame
par (new = TRUE) 
# killed hosts over time
plot (K ~ time, data = output, type='l', ylim = c(0,1), col = 'red', ylab = '', axes = FALSE) 
# remain on same frame
par (new = TRUE)
# infectious hosts over time
plot (I ~ time, data = output, type='l', ylim = c(0,1), col = 'black', ylab = '', axes = FALSE) 
# remain on same frame
par (new = TRUE)  
# recovered hosts over time
plot (R ~ time, data = output, type='l', ylim = c(0,1), col = 'green', ylab = '', axes = FALSE)

legend(x="topright", legend=c("S","L","I","K","R"), fill=c('blue','grey', 'black', 'red', 'green'))

#Find totals: