#Shiny Epidemic Model
library(shiny)
library(deSolve)

#SLIKR DiffEq Model
seir_model = function (current_timepoint, state_values, parameters)
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
      if(dead_pct<scared){
        distancing_pct = distancing*max(min(((dead_pct-min_scary)/(scared-min_scary)),1),0) #Linear between min_scared and scared % 
      } else {
        if(continue_distancing){
          distancing_pct = distancing + (1-distancing)*(dead_pct-scared)/(1-scared) #Linear up to 100% 
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
      return(list(results))
    }
  )
}


#Disease dynamics parameters.
parameter_list = c (
  contact_rate = 30, # number of pre-distancing contacts per day
  transmission_probability = 0.05,  # transmission probability
  infectious_period = 3.3, # infectious period
  latent_period = 0.001, # latent period
  case_fatality = .95,
  distancing = .75, #How mcuh will people stop interacting, maximum
  min_scary = 0.001, #What percentage of people need to die before people start distancing
  scared = 0.05, #What percentage of people need to die before people are about maximally distancing
  continue_distancing=TRUE #Do people keep distancing?
)

#Population and Setup:

#Initial values for sub-populations.
W = 9990        # susceptible hosts
X = 10           # infectious hosts
Y = 0           # recovered hosts
Z = 0           # exposed hosts
ZED = 0         # pre-killed hosts
#Compute total population.
N = W + X + Y + Z + ZED


#Initial state values for the differential equations.
initial_values = c (S = W/N, L = X/N, I = Y/N, K = ZED/N, R = Z/N)

tick_size = 1/10
timepoints = seq (0, 200, by=tick_size)

ui <- fixedPage(
  title="SLIKR Epidemic with Distancing",
  # App title ----
  titlePanel("Epidemic Model with Pandemic Distancing"),
  # Sidebar layout with input and output definitions ----
      column(4,
             # Input: Slider for the number of bins ----
             sliderInput(inputId = "contact_rate",
                         label = "Number of Daily Contacts Per Person:",
                         min = 0.1,
                         max = 50,
                         value = 20),
             sliderInput(inputId = "transmission_probability",
                         label = "P(transmission)/contact:",
                         min = 0,
                         max = 1,
                         value = .05),
             sliderInput(inputId = "infectious_period",
                         label = "Days Infectious:",
                         min = 0.00001,
                         max = 10,
                         value = 3)
      ),
      column(4,
             sliderInput(inputId = "latent_period",
                         label = "Days Latent:",
                         min = 0.0001,
                         max = 10,
                         value = 2),
             sliderInput(inputId = "case_fatality",
                         label = "P(death):",
                         min = 0,
                         max = 1,
                         value = 0.75)
      ),
  column(4,
         sliderInput(inputId = "distancing",
                     label = "Initial decrease in contacts due to fear:",
                     min = 0,
                     max = 1,
                     value = .75),
  sliderInput("scary_range", label = "Decrease occurs when deaths are between:", 
              min = 0, 
              max = 1, 
              value = c(0.005, 0.025)
  ),
  checkboxInput("continue_distancing", 
                label = "Does distancing continue as population thins?", 
                value = TRUE)
  ),
    # Main panel for displaying outputs ----
    column(12,
      # Output: Histogram ----
      plotOutput(outputId = "distPlot"),
      textOutput("Deaths")
    )
  )


server <- function(input, output) {

  output$distPlot <- renderPlot({
    parameter_list = c(
      contact_rate = input$contact_rate, # number of pre-distancing contacts per day
      transmission_probability = input$transmission_probability,  # transmission probability
      infectious_period = input$infectious_period, # infectious period
      latent_period = input$latent_period, # latent period
      case_fatality = input$case_fatality,
      distancing = input$distancing, #How mcuh will people stop interacting, maximum
      min_scary = input$scary_range[1], #What percentage of people need to die before people start distancing
      scared = input$scary_range[2], #What percentage of people need to die before people are about maximally distancing
      continue_distancing = input$continue_distancing #Do people thin out?
    )
    
    #Simulate
    output = lsoda (initial_values, timepoints, seir_model, parameter_list)
    
    #Find where it is stable?
    #i=1;abs_diff=1
    #while (abs_diff>1e-7){
    #  abs_diff=sum(abs(output[i]-output[i+1])[2:6])
    #}
    #
    #i = i*1.1
    #xmax = round(i / tick_size, 0)
        
    #Plot dynamics of all sub-populations.
    # susceptible hosts over time
    plot (S ~ time, data = output, type='l', ylim = c(0,1), col = 'blue', lwd=2, ylab = 'Pct in Compartment', main = 'SLIKR Epidemic with Distancing') 
    # remain on same frame
    par (new = TRUE)    
    # exposed hosts over time
    plot (L ~ time, data = output, type='l', ylim = c(0,1), col = 'grey', lwd=2, xlab='', ylab = '', axes = FALSE)
    # remain on same frame
    par (new = TRUE) 
    # killed hosts over time
    plot (K ~ time, data = output, type='l', ylim = c(0,1), col = 'red', lwd=2, xlab='', ylab = '', axes = FALSE) 
    # remain on same frame
    par (new = TRUE)
    # infectious hosts over time
    plot (I ~ time, data = output, type='l', ylim = c(0,1), col = 'black', lwd=2, xlab='', ylab = '', axes = FALSE) 
    # remain on same frame
    par (new = TRUE)  
    # recovered hosts over time
    plot (R ~ time, data = output, type='l', ylim = c(0,1), col = 'green', lwd=2, xlab='', ylab = '', axes = FALSE)
    legend(x="top", legend=c(paste("Final death toll: ", paste(round(100*max(output[,'K']), 2), "%", sep=""))
                             ,paste("Percentage Unexposed: ", paste(round(100*output[2000,'S'], 2), "%", sep=""))))

    abline(v=min(which(output[,'K']>input$scary_range[1])*tick_size))
    abline(v=min(which(output[,'K']>input$scary_range[2])*tick_size))
    legend(x="topright", legend=c("Susceptible ","Latent","Infectious","Dead","Recovered"), fill=c('blue','grey', 'black', 'red', 'green'))
    
 
      })

}
#Run!
shinyApp(ui, server)