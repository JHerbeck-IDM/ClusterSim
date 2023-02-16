## Branching process model of epidemic spread

# This code is adapted from the study of [Riou and Althaus](https://www.biorxiv.org/content/10.1101/2020.01.23.917351v1), using the [code from GitHub](https://github.com/jriou/wcov), stripped down and rewritten for clarity.

# https://gist.github.com/sdwfrost/b3d0c4cbff7ff460562162affceffb17

library(ggplot2)

# Set random number seed.

set.seed(1234)

# Define Bellman-Harris branching process model.

bhbp <- function(R0,k,shape,scale, index_cases,max_cases,max_time){
  t <- rep(0, index_cases)
  times <- t
  tmax <- 0
  cases <- index_cases
  while(cases > 0 & length(times) < max_cases) {
    secondary = rnbinom(cases, size=k, mu=R0)
    t.new = numeric()
    for(j in 1:length(secondary)) {
      t.new = c(t.new, t[j] + rgamma(secondary[j], shape = shape, scale = scale))
    }
    t.new = t.new[t.new < max_time]
    cases = length(t.new)
    t = t.new
    times = c(times, t.new)
  }
  #times <- sort(times)
  return(times)
}


# Set parameter values, using values from [Imai et al.](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-2019-nCoV-transmissibility.pdf) and assuming a gamma distribution of generation times with mean `mu` 8.4 days and standard deviation `sd` of 3.8 days, taken from [Lipsitch et al. (2003](https://science.sciencemag.org/content/300/5627/1966.full).

R0 <- 2.6
k <- 0.16
mu <- 8.4
stdev <- 3.8
shape <- (mu/stdev)^2
scale <- (stdev^2)/mu


# Plot the generation time distribution.

t <- seq(0,30,by=0.01)
g <- dgamma(t,shape=shape,scale=scale)
ggplot(data.frame(GenerationTime=t,Probability=g))+geom_line(aes(x=GenerationTime,y=Probability))

# Plot the offspring distribution.

i <- seq(0,10)
d <- dnbinom(i, size=k, mu=R0)
ggplot(data.frame(Number=i,Probability=d))+geom_bar(aes(x=Number,y=Probability),stat="identity")

# Initial and stopping conditions.

index_cases <- 40
max_cases <- 5e4
max_time <- 90


# Set the number of simulations (note - Imai et al. used 5000).

nsims <- 50


# Run simulations.

l <- list()
for(i in 1:nsims){
  times <- bhbp(R0,k,shape,scale,index_cases,max_cases,max_time)
  # Generates cumulative counts per day
  # Note that this includes the index cases
  l[[i]] <- cumsum(hist(times, breaks = 0:max_time,plot=FALSE)$counts)
}

# Combine individual runs into a dataframe.

results <- as.data.frame(do.call(cbind,l))

median_cases <- apply(results,1,median)
lq_cases <- apply(results,1,quantile,0.025)
uq_cases <- apply(results,1,quantile,0.975)
summary_results <- data.frame(Day=seq(1,max_time),
                              Date=as.Date("2019-12-01")+seq(1,max_time),
                              Median=median_cases,
                              Lower=lq_cases,Upper=uq_cases)

# Add day/dates with day 0 corresponding to 2019-12-01.

results$Day <- seq(1,max_time)
results$Date <- as.Date("2019-12-01")+results$Day

results$Day <- seq(1,max_time)
results$Date <- as.Date("2019-12-01")+results$Day

# Reshape results into 'long' format.

df <- reshape(results,varying=paste("V",1:nsims,sep=""),direction="long",sep="",idvar="Day",timevar="Run",v.names=c("Cases"))

# Plot trajectories over time, highlighting 4000 cases on 2020-01-18.

ntraj <- 20
ggplot(summary_results)+
  geom_line(aes(x=Date,y=Median),color="red")+
  geom_line(aes(x=Date,y=Upper),color="blue")+
  geom_line(aes(x=Date,y=Lower),color="blue")+
  geom_ribbon(aes(x=Date,ymin=Lower,ymax=Upper, linetype = NA), fill="blue", alpha = 0.1)+
  geom_line(aes(x=Date,y=Cases,group=factor(Run)),data=df[df$Run<=ntraj,],color="gray")+
  coord_cartesian(xlim=c(as.Date("2019-12-01"),as.Date("2020-01-30")),ylim=c(1,20000))+
  geom_vline(xintercept=as.Date("2020-01-18"))+
  geom_hline(yintercept=4000)+
  theme_classic()
