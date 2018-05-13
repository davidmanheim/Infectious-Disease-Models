Single_Seir_Model = function (current_timepoint, state_values, parameters)
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

#Speed: 10k runs in 40 minutes, 1/4 second per run.
#Should I bother speeding it up? Probably not.


#The multipandemic model needs some helper functions:
source('Multipandemic_Helper_Functions.R')

MultiPandemic_ode_model = function (current_timepoint, state_values, parameters){
  
  #First, create a local array of state values.
  #The order of State values can be lexicographic, but that seems harder...
  Population = state_values # Flat. (Works well with helper functions.)
  
  # I can recreate the list from before, or pass the structured set in as a variable. 
  # I don't end up using it except to flatten it. I should just use numbers and get names with CNumToName(i) 
  
  Compartments = append(unlist(C_List), "K")
  CompartmentKNum = length(Compartments) #Useful Constant
  dead_pct = Population[length(Population)]
  
  dCompartments = rep(0,length(Compartments))
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    { 
      if(dead_pct<max_scary){
        distancing_pct = distancing*min(max(((dead_pct-min_scary)/(max_scary-min_scary)),0),1) #Linear between min_scared and scared % 
      } else {
        if(continue_distancing){
          distancing_pct = distancing + (1-distancing)*(dead_pct-max_scary)/(1-max_scary) #Linear up to 100% 
        } else {distancing_pct = distancing} #Or stay at max.
      }
      
      adjusted_contact_rate = contact_rate * (1-distancing_pct)
      
      #beta = adjusted_contact_rate * transmission_probability
      #gamma = (1 / infectious_period) - case_fatality / infectious_period
      #delta = 1 / latent_period
      #mu = case_fatality / infectious_period
      greeks = c('beta', 'gamma', 'delta', 'mu')
      
      Disease_Greeks = as.list(rep(0.0,D))
      #names(Disease_Greeks) = greeks
      
      
      for (d in 1:D){
        Disease_Greeks[[d]] = c(
          beta = adjusted_contact_rate * p_s[[d]]$transmission_probability,
          gamma = (1 / p_s[[d]]$infectious_period) - p_s[[d]]$case_fatality / p_s[[d]]$infectious_period,
          delta = 1 / p_s[[d]]$latent_period,
          mu = p_s[[d]]$case_fatality / p_s[[d]]$infectious_period
        )
      }
      
      #The parameters are: D, C, Population, contact_rate, min_scary, max_scary, max_distancing, d[i]
      
      
      # To calculate infectiousness, we need the sum of infectious people across compartments for each disease.
      Infectious=rep(0,D)
      for (compartmentname in Compartments){
        for (place in 1:nchar(compartmentname)){ #i.e. 1:D
          if(substr(compartmentname,place,place)=="I"){
            Infectious[place] = Infectious[place] + Population[CNameToNum(compartmentname)]
          }
        } 
      }
      
      #print(paste("Pop",paste(Population)))
      # We need to loop through compartments and add / subtract components where needed.
      for (compartmentname in Compartments){
        compartmentnumber=CNameToNum(compartmentname)
        compartmentvector=CNameToVec(compartmentname)
        #OK, we have the vector for each. 
        #print(compartmentname)
        #print(compartmentnumber)
        #print(compartmentvector)
        #The derivative is the sum of inputs from each compartment 
        #with hamming distance 1. Now, we populate that:
        #print(compartmentvector)
        for (d in 1:D){
          #Components for each compartment after this one:
          if (compartmentvector[d] == 1){ #S
            dCompartments[compartmentnumber] = dCompartments[compartmentnumber] - (Disease_Greeks[[d]]['beta'] * Population[compartmentnumber] * Infectious[d])
            destinationvector = compartmentvector
            destinationvector[d] = 2
            destinationnumber = CVecToNum(destinationvector)
            dCompartments[destinationnumber] = dCompartments[destinationnumber] + (Disease_Greeks[[d]]['beta'] * Population[compartmentnumber] * Infectious[d])
          }
          if (compartmentvector[d] == 2){ #L
            dCompartments[compartmentnumber] = dCompartments[compartmentnumber] - Disease_Greeks[[d]]['delta'] * Population[compartmentnumber]
            destinationvector = compartmentvector
            destinationvector[d] = 3
            destinationnumber = CVecToNum(destinationvector)
            dCompartments[destinationnumber] = dCompartments[destinationnumber] + Disease_Greeks[[d]]['delta'] * Population[compartmentnumber]
          }
          if (compartmentvector[d] == 3){ #I: - (mu * I) - (gamma * I)
            dCompartments[compartmentnumber] = dCompartments[compartmentnumber] - (Disease_Greeks[[d]]['mu'] * Population[compartmentnumber]) - (Disease_Greeks[[d]]['gamma'] * Population[compartmentnumber])
            #First, add to R
            destinationvector = compartmentvector
            destinationvector[d] = 4
            destinationnumber = CVecToNum(destinationvector)
            dCompartments[destinationnumber] = dCompartments[destinationnumber] + (Disease_Greeks[[d]]['gamma'] * Population[compartmentnumber])
            #Now, add to "K"
            dCompartments[CompartmentKNum] = dCompartments[CompartmentKNum] + (Disease_Greeks[[d]]['mu'] * Population[compartmentnumber])
          }
          
        }
        
      }
      #print(dCompartments)
      return(list(dCompartments))
    }
    
    #dS = (-beta * S * I)
    #dL = (beta * S * I) - (delta * L)
    #dI = (delta * L) - (mu * I) -(gamma * I)
    #dK = (mu * I)
    #dR = (gamma * I)
    
    
    #results = c(x) #List of D_compartment values
  )
}

source('SLIKR_Multimodel_Unit_Tests.R')

SLIKR_Unit_Tests = function(){
  Pass_All=TRUE
  if(SLIKR_Unit_Test_Multi_vs_Single()){
    print("PASS Single Model and Multimodel with 1 Disease outcomes match")
  } else{
    print("FAIL - Single Model and Multimodel with 1 Disease outcomes differ")
    Pass_All=FALSE
  }
  
  if(Conservation_of_Population("s")){
    print("PASS Single Model preserves population")
  } else{
    print("FAIL - Single Model does not preserve population")
    Pass_All=FALSE
  }
  
  if(Conservation_of_Population("m")){
    print("PASS Multi Model preserves population")
  } else{
    print("FAIL - Multi Model does not preserve population")
    Pass_All=FALSE
  }
  
  return(Pass_All)
  }
  
# SLIKR_Unit_Tests()
