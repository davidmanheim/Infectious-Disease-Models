# These functions are needed to set up the runs of the model. 

# They are also run inside this script, for program flow reasons.

# Each compartment is a length D string from alphabet C. All are valid.
# These form a tree, where each node is increment-1, with an order. 
# Origin: SSS.
# Children: SSI, SIS, ISS
# Connections between nodes have transition prob's given in the matrix
# [dS_0,dE_0...]
# [dS_1, ...]
# [dS_D, dE_D, ...]
# (subscripts are the string position number)
# Each overall transition probability is the sum of derivatives in, minus sum of those out.

# The entries are based on several features, which are universal or differ per disease.
# Diseases are characterized by array Disease_Characteristics[i] = c(transmission_probability, latent_period, infectious_period, case_fatality)
# state_values are the values of the compartments, in order. (Ordering matters!)
# Plus an extra compartment for those killed.
# parameters is: list(D, C, Population, contact_rate, min_scary, max_scary, max_distancing, Disease_Characteristics)


setup_model = function(C,D){
  
  Origin = paste(rep(C[1],D),sep="",collapse="") #Start with Susceptibles
  Total_Layers=1+(length(C)-1)*(D)
  C_List=as.list(c(rep('',Total_Layers)))
  C_List[1] = Origin
  
  #C_List[[2]]=c("SL","LS")
  
  for(layer in 1:Total_Layers){
    #In layer N, we look at layer N-1 and move each compartment up 1.
    for (string in C_List[[layer]]){
      #Get the list of places
      stringvector = rep(0,nchar(string))
      for (place in 1:nchar(string)){
        stringvector[place] = match(substr(string,place,place),C)
      }
      #Modify each character up by 1, if possible.
      for (place in 1:nchar(string)){
        new_vector = stringvector
        if (new_vector[place]<length(C)){ #Don't increment recovered people.
          new_vector[place] = new_vector[place] + 1
          New = paste0(C[new_vector], collapse ="")
          #We check to see if that is already in the list, and if not we add it.  
          if(!(New %in% C_List[[layer+1]])){
            if(C_List[[layer+1]][[1]]=='') {
              C_List[[layer+1]][[1]] = New
            }else {C_List[[layer+1]] = append(C_List[[layer+1]],New)}
          }
        }
      }
    }
  }
  
  
  # We want latent infections in each of the second-layer compartments, 
  # and everyone else starts in I or LS, etc,
  Population = rep(0, length(unlist(C_List))+1) #Add 1 for "killed" compartment.
  #Population[(1+length(C_List[[2]])):(1+length(C_List[[2]])+length(C_List[[3]]))+1] = 0.001 #All of the first I compartments
  #Population[1]= 1 - (length(C_List[[3]])*0.001)
  Population[2:(1+length(C_List[[2]]))] = 0.001 #All of the first I compartments
  Population[1]= 1 - (length(C_List[[2]])*0.001)
  
  
  #Compartments = append(unlist(C_List), "K")
  return(list(C_List=C_List,Compartments = append(unlist(C_List), "K"),Population=Population))
}


# ATTACH IS BAD PRACTICE.
# Too bad I can't figure out how to do this another way easily, since I need it to stay there until the end of the model runs...
# Hmmm... Maybe I should split helper functions into multiple files?

attach(setup_model(C,D))

#Some Helper Functions:
CNameToNum = function(compartmentname){
  return(match(compartmentname,Compartments))
}
CNumToName = function(compartmentnum){
  return(Compartments[compartmentnum])
}
CVecToNum = function(compartmentvector){
  if (compartmentvector[1]==666){return(length(Compartments))}
  else return(CNameToNum(paste0(C[compartmentvector],collapse="")))
}
CVecToName = function(compartmentvector){
  if (compartmentvector[1]==666){return("K")}
  else return(paste0(C[compartmentvector],collapse=""))
}
CNameToVec = function(compartmentname){
  if (compartmentname=='K'){compartmentvector=rep(666,D)} else{
    compartmentvector = rep(0,nchar(compartmentname))
    for (place in 1:nchar(compartmentname)){
      compartmentvector[place] = match(substr(compartmentname,place,place),C)
    }}
  return(compartmentvector)
}
CNumToVec= function(compartmentnum){
  return(CNameToVec(CNumToName(compartmentnum)))
}

generate_pandemics = function(D){
  Disease_Characteristics = as.list(rep(0,D))
  
  for (i in 1:D){
    #Each of these are drawn from a semi-reasonable prior:
    Disease_Characteristics[[i]] = list(
      transmission_probability = rbeta(1,2,20),
      latent_period = rlnorm(1,2,1),
      infectious_period = rlnorm(1,1.25,.5),
      case_fatality = rbeta(1,4,4)
    )
  }
  return(Disease_Characteristics)
}

p_s = generate_pandemics(D)
