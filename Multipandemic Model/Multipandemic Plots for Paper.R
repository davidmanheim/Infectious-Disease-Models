D=3
#Before running the script, load the dataset from the relevant run!
if(D==3){N=1000; load("~/Pandemic Models/SIR Models/3_disease_1k_Runs.RData.gz")}
if(D==4){N=500; load("~/Pandemic Models/SIR Models/4_disease_500_Runs.RData.gz")}
if(D==5){N=100; load("~/Pandemic Models/SIR Models/5_disease_100_Runs.RData.gz")}

DeathRates=data.frame(as.list(c(1:(D+1))))
for (i in 1:N){
  DeathRates[i,1:D]=unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 5))
  DeathRates[i,D+1]=Multipandemic_List[[i]][[D+1]]
}

# I'm thinking for each multipandemic I want a ellipse with width = max(DeathRate)-min(DeathRate),
# height = max(Latency)-min(Latency), and an arrow from (avg(DeathRate),avg(Latency)) to (Multipandemic_Deathrate, avg(Latency))

# I need a data.frame with the correct Points.

Ellipses_DF=data.frame(as.list(c(1:8)))
names(Ellipses_DF)<-c('Max_fatal','Min_fatal','Max_latent','Min_latent','Avg_fatal','Avg_latent','Multi_fatal', 'Fastest-fatal')
for (i in 1:N){
  Ellipses_DF[i,1]=max(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 5)))
  Ellipses_DF[i,2]=min(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 5)))
  Ellipses_DF[i,3]=max(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 2)))
  Ellipses_DF[i,4]=min(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 2)))
  Ellipses_DF[i,5]=mean(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 5)))
  Ellipses_DF[i,6]=mean(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 2)))
  Ellipses_DF[i,7]=Multipandemic_List[[i]][[D+1]]
  Ellipses_DF[i,8]=Multipandemic_List[[i]][[which(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 2))==min(unlist(lapply(Multipandemic_List[[i]][1:D],`[`, 2))))]]$population_fatality_rate
}

#Now make the plot:

p=plot(1, type="n", xlab="Epidemic Fatality Rate", ylab="Disease Latency Periods", xlim=c(min(Ellipses_DF[,2]*0.9), max(Ellipses_DF[,1]*1.1)), 
       ylim=c(max(mean(Ellipses_DF[,4]-sd(Ellipses_DF[,3])),0), mean(Ellipses_DF[,4]+2*sd(Ellipses_DF[,3]))))
# Reduce Y-axis plot area to be reasonable.

#And fill it in:
num=10 #N
for (i in 1:num){
  #plotellipse(rx=(Ellipses_DF[i,1]-Ellipses_DF[i,2])/2, ry=(Ellipses_DF[i,3]-Ellipses_DF[i,4])/2, 
  #            mid=c((Ellipses_DF[i,1]+Ellipses_DF[i,2])/2, (Ellipses_DF[i,3]+Ellipses_DF[i,4])/2), lwd=0.25)
  arrows(x0=(Ellipses_DF[i,1]),y0=(Ellipses_DF[i,3]+Ellipses_DF[i,4])/2,x1=(Ellipses_DF[i,7]),y1=(Ellipses_DF[i,3]+Ellipses_DF[i,4])/2, length=0.05, col="red") #From Max to Multi
  
  segments(x1=(Ellipses_DF[i,1]),y1=(Ellipses_DF[i,3]+Ellipses_DF[i,4])/2,x0=(Ellipses_DF[i,5]),y0=(Ellipses_DF[i,3]+Ellipses_DF[i,4])/2, lty=3, col="blue", lwd=1) #Average vs. Max
}
legend(x="topright", legend=c("Mean Rate to Max Rate","Max Rate to Multipandemic Rate"), fill=c("blue","red"))

if (CreatePlots==TRUE) {
  pdf(paste0(D,"-Disease Syndemic - Histograms ",format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
  
  par(mfrow = c(3,2))
  
  hist(unlist(DeathRates[,1:3]), main="Single Pandemic Death Rate",xlab="% Death Rate",breaks=20)
  
  hist(Ellipses_DF[,5], main="Average Single Rate / Syndemic",xlab="% Death Rate",breaks=20)
  
  hist(Ellipses_DF[,1], main="Syndemic Death Rate",xlab="% Death Rate",breaks=20)
  
  hist(Ellipses_DF[,7], main="Worst Single Rate / Syndemic",xlab="% Death Rate",breaks=20)
  
  hist((Ellipses_DF[,1]-Ellipses_DF[,7]), main="Syndemic Rate minus Worst Rate",xlab="% Death Rate",breaks=20)
  
  hist((Ellipses_DF[,5]-Ellipses_DF[,7]), main="Syndemic Rate minus Average Rate",xlab="% Death Rate",breaks=20)
  
  dev.off()
  
  pdf(paste0(D,"-Disease Syndemic - Simulation ",format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
  par(mfrow = c(1,1))
  
  num=100 #Try just having them in a row
  p=plot(1, type="n", xlab="Epidemic Fatality Rate", ylab="Trial Number", xlim=c(min(Ellipses_DF[,2]*0.9), max(Ellipses_DF[,1]*1.1)), 
         ylim=c(0,num+0.6),yaxs="i")
  par(mar=c(1.1, 4.1, 1.1, 2.1), xpd=TRUE)
  plot_range = par("usr")
  
  for (i in 1:num){
    if(i%%2==0){rect(xleft=plot_range[1]+0.25, xright=plot_range[2]-0.25, ybottom=i-(25/num), ytop=i+(25/num), density=100, col='gray90')}
    points(x=Multipandemic_List[[i]][D+1], y=i, col='red',bg='red', cex=.5, pch=4)
    for (d in 1:D){ 
      points(x=Multipandemic_List[[i]][[d]][5], y=i, col='blue',bg='blue', cex=.25)
    }
  }
  
  legend(x="top", legend=c("Per-Cluster Pandemics Death Rates","Syndemic Death Rate"), fill=c("blue","red"), inset=c(0,-0.12))
  
  dev.off()
  
}

# We pick the cases with least and most extreme disparities between syndemic and worst-case disease, and show a population graph of each.
#Do this only for D=3

# This is done using the Multipandemic model script, since I need to rerun it to generate the population vs. time graphs.