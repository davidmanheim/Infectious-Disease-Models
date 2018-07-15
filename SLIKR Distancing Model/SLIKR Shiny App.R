#Shiny Epidemic Model
library(shiny)
library(deSolve)

source('Single SLIKR Model.R')

#Disease dynamics parameters.
parameter_list = c (
  contact_rate = 30, # number of pre-distancing contacts per day
  transmission_probability = 0.05,  # transmission probability
  infectious_period = 3.3, # infectious period
  latent_period = 0.01, # latent period
  case_fatality = .95,
  distancing = .75, #How mcuh will people stop interacting, maximum
  min_scary = 0.025, #What percentage of people need to die before people start distancing
  max_scary = 0.075, #What percentage of people need to die before people are about maximally distancing
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
ui <- fluidPage(
  title="SLIKR Epidemic with Distancing",
  # App title ----
  titlePanel("Epidemic Model with Pandemic Distancing"),
  # Sidebar layout with input and output definitions ----
      column(width=4,
             # Input: Slider for the number of bins ----
             sliderInput(inputId = "contact_rate",
                         label = "Initial Daily Contacts Per Person:",
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
      column(width=4,
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
  column(width=4,
         sliderInput(inputId = "distancing",
                     label = "Percentage decrease in contacts due to fear",
                     min = 0,
                     max = 1,
                     value = .75),
  sliderInput("scary_range", label = "Decrease occurs when deaths are between:", 
              min = 0.01, 
              max = 1, 
              value = c(0.01, 0.03)
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
      max_scary = input$scary_range[2], #What percentage of people need to die before people are about maximally distancing
      continue_distancing = input$continue_distancing #Do people thin out?
    )
    
    #Simulate
    output = lsoda (initial_values, timepoints, Single_SLIKR_Model, parameter_list)
    
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
    if (max(output[,'K'])>input$scary_range[1]){
      abline(v=min(which(output[,'K']>input$scary_range[1])*tick_size))
      if (max(output[,'K'])>input$scary_range[2]){  
        abline(v=min(which(output[,'K']>input$scary_range[2])*tick_size))
    }
    }
    legend(x="topright", legend=c("Susceptible ","Latent","Infectious","Killed","Recovered"), fill=c('blue','grey', 'black', 'red', 'green'))
    
 
      })

}
#Run!
shinyApp(ui, server)