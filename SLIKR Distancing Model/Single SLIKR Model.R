Single_SLIKR_Model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  L = state_values [2]        # latent
  I = state_values [3]        # infectious
  K = state_values [4]        # killed (immediately)
  R = state_values [5]        # recovered
  Tot = sum(state_values) #Total
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      dead_pct = (K/Tot)
      if(dead_pct<max_scary){
        distancing_pct = distancing*min(max(((dead_pct-min_scary)/(max_scary-min_scary)),0),1) #Linear between min_scared and scared % 
      } else {
        if(continue_distancing){
          distancing_pct = distancing + (1-distancing)*(dead_pct-max_scary)/(1-max_scary) #Linear up to 100% 
        } else {distancing_pct = distancing} #Or stay at max.
      }
      adjusted_contact_rate = contact_rate * (1-distancing_pct)
      beta = adjusted_contact_rate * transmission_probability
      gamma = (1 / infectious_period) - case_fatality / infectious_period
      delta = 1 / latent_period
      mu = case_fatality / infectious_period
      
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