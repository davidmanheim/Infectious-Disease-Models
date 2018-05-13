SLIKR_Unit_Test_Multi_vs_Single = function(){
  
  D=1
  C=c('S','L','I','R')
  #Plus compartment 'K'
  source('Multipandemic_Helper_Functions.R')
  # That gives me the models and all. Now to run and compare them;
  
  #Don't do distancing, so it will work!
  Overall_variables = as.list(c(min_scary=.9901,max_scary=0.991,contact_rate=30,distancing=.9,continue_distancing=TRUE))
  
  Single_Model_Out = function(Pandemic) { 
    output= lsoda(c(S=.999,L=.001,I=0,K=0,R=0), seq(0,200,0.1), Single_Seir_Model, c(Pandemic,Overall_variables))
    return(output)
  }
  
  SModel_Outcome = Single_Model_Out(p_s[[1]])
  SModel_Outcome = SModel_Outcome[,c(1,2,3,4,6,5)] #rearrange to match multimodel.
  
  Param_Combo=list(D=D, min_scary=Overall_variables$min_scary, max_scary=Overall_variables$max_scary, 
                   contact_rate=Overall_variables$contact_rate, distancing=Overall_variables$distancing,
                   continue_distancing=Overall_variables$continue_distancing,
                   C=as.array(C), p_s=p_s, C_List=as.list(C_List))
  
  MModel_Outcome = lsoda(Population, seq(0,200,0.1),MultiPandemic_ode_model, Param_Combo)
  
  return(sum(abs(SModel_Outcome-MModel_Outcome))==0)
}


Conservation_of_Population  = function(model){
  Overall_variables = as.list(c(min_scary=.025,max_scary=0.075,contact_rate=30,distancing=.9,continue_distancing=TRUE))
  if (model=="s"){
 #Set up single model run.
    Model_fx = Single_Seir_Model
    parameters = c(Pandemic,Overall_variables)
    population = c(S=.999,L=.001,I=0,K=0,R=0)
   } else {if (model=="m"){
     Model_fx = MultiPandemic_ode_model
     parameters = list(D=D, min_scary=Overall_variables$min_scary, max_scary=Overall_variables$max_scary, 
                       contact_rate=Overall_variables$contact_rate, distancing=Overall_variables$distancing,
                       continue_distancing=Overall_variables$continue_distancing,
                       C=as.array(C), p_s=p_s, C_List=as.list(C_List))
     population=c(0.99,rep(0.01/(length(Population)-1), length(Population)-1))
     
     
  }} else return(NULL)

  output= lsoda(population, seq(0,200,0.1), Model_fx, parameters)
  
  #Row picked randomly, exclude time which is the first column.
  row_to_check = sample(1:2000, 1)
  
  if ( sum(output[row_to_check,2:(length(population)+1)])==1 ){
    return(TRUE)
  } else return(FALSE)
  
}