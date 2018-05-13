# Base model from https://rpubs.com/srijana/110753
#Go to my directory, so I can load sources and save pictures.
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # If run from Rstudio,
  #require(rstudioapi)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
} else {setwd(getSrcDirectory()[1])} #if not in Rstudio.

# Outline:

#There are multiple pandemics. 
#Someone can be in any number of different SEI buckets.
#This means we need a bucket for each subset of disease/Status combination
# Compartments are generated programmatically
#For example, D1S-D2E-D3S-D4I-D5E would be exposed to D1 and D5, Infectious with D4, and susceptible to the others.
#Transition probabilities need to treat each differently - combinatorial explosion...

library (deSolve) 

timepoints = seq (0, 200, by=0.1)

# Number of Diseases
D=3
# Set of Compartments per disease:
C=c('S','L','I','R')
#Plus compartment 'K'



#SLIKR DiffEq Model
source('SLIKR_Models.R')

#Helper Functions moved to separate file.
#It is called from the separate SLIKR_Models.R file, which is included above.
# Helper functions include:
# generate_pandemics() function
# setup_model(C,D) for compartment and population creation (returned as list of those variables)
# All the convert from compartment id to name to list functions for the pandemic model. 
# The generation and model setup is done there.



# First, generate D pandemics: (function in Helper File.)
Overall_variables = as.list(c(min_scary=.025,max_scary=0.075,contact_rate=30,distancing=.9,continue_distancing=TRUE))

#p_s = generate_pandemics(D) # Now in helper function file.
#Keep_p_s=p_s
#Yay! It works!

#p_s[[2]]=p_s[[1]]
#p_s[[3]]=p_s[[1]]

#I should probably simulate a couple thousand of these and find the distribution of outcomes.


Single_Model_Out = function(Pandemic) { 
    output= lsoda(c(S=.999,L=.001,I=0,K=0,R=0), seq(0,200,0.1), Single_Seir_Model, c(Pandemic,Overall_variables))
    return(output)
    }

#SModel_Outcome = Single_Model_Out(p_s[[1]])
#SModel_Outcome = SModel_Outcome[,c(1,2,3,4,6,5)] #rearrange to match multimodel.
# The above lets us check that a multimodel with only 1 diseasae matches.

# Now we can look at the empirical distribution of the death toll given these input distributions:
simulate_histogram = function(n) {
t_simstart = Sys.time()
simulation_test_pandemics = generate_pandemics(n)
Death_Toll = rep(0,n)
for (i in 1:n){
Death_Toll[i] = max(Single_Model_Out(simulation_test_pandemics[[i]])[,'K'])
}
t_simend = Sys.time()
print(t_simend - t_simstart)

return(hist(Death_Toll))
}

# histogram = simulate_histogram(10) #YAY! It woks. (3 seconds for 10.)

#histogram = simulate_histogram(10000) # 40 minutes for 10k. Too slow.
#plot(histogram)

#p_s_backup <- p_s


# This actually builds the Population, Compartments, and C_List variables: 
# (From helper functions file)


Param_Combo=list(D=D, min_scary=Overall_variables$min_scary, max_scary=Overall_variables$max_scary, 
                 contact_rate=Overall_variables$contact_rate, distancing=Overall_variables$distancing,
                 continue_distancing=Overall_variables$continue_distancing,
                 C=as.array(C), p_s=p_s, C_List=as.list(C_List))

#Diagnostic info for the runs:



#########################
#    Run the Model!     #
#########################
t1 = Sys.time()
result_solution = lsoda(Population, seq(0,200,0.1),MultiPandemic_ode_model, Param_Combo)
Result <- as.data.frame(result_solution)
names(Result)<-c("Time",Compartments)
t2 = Sys.time()
#result_solution_hi = lsoda(Population, seq(0,200,0.1),MultiPandemic_ode_model, Param_Combo,rtol = 1e-7, atol = 1e-7,hmax=1e-2,maxsteps=5e7)
#Result_hi <- as.data.frame(result_solution)
#names(Result_hi)<-c("Time",Compartments)
#t3 = Sys.time()
t2-t1
# 3 Comparments - 30 seconds.
# 4 Compartments - 2.1 minutes. 
# 5 Compartments - 9.05 minutes. 

# If comparing a single model, we check the following:
# sum(abs(SModel_Outcome - result_solution)) == 0



#Sum across compartments for each disease alone:
entries = length(Result[,'K'])
Aggregate_Result = data.frame(matrix(vector(mode = 'numeric',length = (length(C)*D+1)*entries), nrow = entries, ncol = (length(C)*D+1)))
#Do clever stuff to name everything in form D#C, plus K
names(Aggregate_Result) = append(paste0(rep(paste0(rep("D",D),1:D), each=length(C)),rep(C,D)),"K")

#This is a horribly inefficient solution. Bad Programmer, no cookie!
#OK, it's better, but not great. Good enough.
for (compartment in Compartments){
  for (d in 1:D){
    vec_item = CNameToVec(compartment)[d]
    if (vec_item<(length(C)+1)){ #Not 'K'
      eval(parse(text=paste0("Aggregate_Result[,","\'D",d,C[vec_item],"\'] <- Aggregate_Result[,","\'D",d,C[vec_item],"\'] + ", "Result[,\'",compartment,"\']")))
    } else Aggregate_Result[,"K"] <- Result[,'K']
  }
}
#Multi-counting susceptibles...?
#Yes, kind of, but not a problem.

# How long do the graphs go out? (Problems exist with showing longer-term, slower burn diseases on the same graph as quick ones.)
end_time = 600

pdf(paste0("plots ",format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))


#Insert a plot that just has the diagnostic info.
  Diag_String = ""
  for (i in 1:D) {
    Diag_String = paste(Diag_String, "/n", "Disease", i, "characteristics:",paste(labels(p_s[[i]]), p_s[[i]], collapse = ' '))
  }
  Diag_String = paste(Diag_String, "/n",  "Overall characteristics:",paste(labels(Overall_variables), Overall_variables, collapse = ' '))
  s=strsplit(paste(Diag_String, collapse = " "), "/n")
  par(mar = c(0,0,0,0)) #Text needs no borders!
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, gsub('(.{1,90})(\\s|$)', '\\1\n', s), #Split to have decent line breaks
     cex = 0.6, col = "black")
  par(mar = c(5, 4, 4, 2) + 0.1)

#t0=Sys.time()

# Now plot the multipandemic model:
  plot(0,0,xlim = c(0,end_time),ylim = c(0,1),type = "n")
  compartmentcolors = c('grey','darkgoldenrod2',"red","darkblue")
  
  for (i in 1:(length(names(Aggregate_Result))-1)){
    compartment=names(Aggregate_Result)[i]
    lines(Aggregate_Result[,i],col = compartmentcolors[((i-1)%%4+1)], type = 'l', lty=floor((i-0.01)/length(C))+2, lwd=2)
  }
  lines(Aggregate_Result[,"K"],col = 'black', type = 'l', lwd=2)
  #cl <- rainbow(15)
  legend("topright",legend=paste("D",1:D, sep=""), lty=2:(D+1))
  legend("topleft",legend = c("Susceptible ","Latent","Infectious","Killed","Recovered"), fill=c('grey','darkgoldenrod2',"red","black","darkblue"))
  legend(x="top", legend=c(paste("Final death toll: ", paste(round(100*max(Aggregate_Result[,"K"]), 2), "%", sep=""))))
  #dim(Result)

# some_test_stuff = function(){
# i=2
# t1 = Sys.time()
# Single_Model_Outcome_hiprec = lsoda(c(S=.997,L=.003,I=0,K=0,R=0), seq(0,200,0.1), Single_Seir_Model, c(p_s[[i]],Overall_variables),rtol = 1e-7, atol = 1e-7,hmax=1e-2,maxsteps=5e6)
# t2 = Sys.time()
# #Single_Model_Outcome_lowprec = lsoda(c(S=.997,L=.003,I=0,K=0,R=0), seq(0,200,0.1), Single_Seir_Model, c(p_s[[i]],Overall_variables),rtol = 1e-6, atol = 1e-6)
# #t3 = Sys.time()
# 
# t2-t1
# #t3-t2
# 
# sum(abs(Single_Model_Outcome_hiprec - Single_Model_Outcome_lowprec))
# max(Single_Model_Outcome_hiprec[,'K'])
# max(Single_Model_Outcome_lowprec[,'K'])
# Single_Model_Outcome_hiprec - Single_Model_Outcome_lowprec
# }


#Now Plot each on its own:
for (i in 1:D){
  Single_Model_Outcome = lsoda(c(S=.997,L=.003,I=0,K=0,R=0), seq(0,200,0.1), Single_Seir_Model, c(p_s[[i]],Overall_variables))  
  plot(0,0,xlim = c(0,end_time),ylim = c(0,1),type = "n")
  
  for (c in 1:length(C)){
    compartment=names(Single_Model_Outcome)[c]
    lines(Single_Model_Outcome[,C[c]],
          col = compartmentcolors[c], 
          type = 'l', lty=i+1, lwd=2)
  }
  
  lines(Single_Model_Outcome[,"K"],col = 'black', type = 'l')
  #legend("topright",legend=c("D1","D2","D2"), lty=c(2,3,4))
  legend("topright",legend=c(paste("Disease",i)), lty=(i+1))
  legend("topleft",legend = c("Susceptible ","Latent","Infectious","Killed","Recovered"), fill=c('grey','darkgoldenrod2',"red","black","darkblue"))
  legend(x="top", legend=c(paste("Final death toll: ", paste(round(100*max(Single_Model_Outcome[,"K"]), 2), "%", sep=""))),bty="n")
  legend(x="left", legend=c(paste("CFR: ", paste(round(100*p_s[[i]][['case_fatality']], 2), "%
p(t):", round(100*p_s[[i]][['transmission_probability']], 2), "%
Latent:", round(p_s[[i]][['latent_period']], 1), " days
Infectious:", round(p_s[[i]][['infectious_period']], 1), " days", sep=""))),bty="n")
}
#Done Plotting. 
dev.off() #Save pdf
