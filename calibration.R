### Model calibration with ABC


#### This is not done



# The goals of this modeling experiment require us to vary biological and technical 
# processes across several different simulations, all of which have similar average 
# epidemic patterns, i.e. all simulations have R-eff of ~1, or mean viral generation 
# times of X years, etc., and to evaluate if phylogenetic clustering patterns are
# affected by different processes.

# So we have to set target stats for the calibration, and we need to choose a parameter, or
# parameters, to vary. 

### Target stat

So first specify the target stat.

```{r target stats}

time <- c(180,360,540,720,900,1080)  # every 6 months for 3 years
mean.epidemic.pattern.stat <- rep(0.3, 6)    # 
target.stats <- data.frame(time, placebo.incidence, VE.target)
```

### Function to get output target stat for each simulation replicate

```{r function.for.target.stat}

target.stat.from.sim <- function(ClusterSim.output){

# pull in the output from each simulation (from the population_summary file)
# estimate the target stat

return(output df)

}
```

### Calibration input model

Function to wrap the ClusterSim model run and make the output file.  

```{r easyABC.input.model}

f <- function(x) {
  param <- param.dcm(beta = 0.004, 
                     c = 90/365,
                     prev = 0.10,    
                     #lambda = beta*c*prev,
                     #lambda = 0.000096,
                     lambda = 8e-06,
                     epsilon = x[1]) 
                     #risk = 10)
                     #risk = x[2])
  init <- init.dcm(Sp = 5000, Ip = 0,
                 Sv = 5000, Iv = 0,
                 Sph = 500, Iph = 0,    #placebo, high risk
                 Spm = 4000, Ipm = 0,   #placebo, medium risk
                 Spl = 500, Ipl = 0,    #placebo, low risk
                 Svh = 500, Ivh = 0,    #vaccine
                 Svm = 4000, Ivm = 0,   #vaccine
                 Svl = 500, Ivl = 0,    #vaccine
                 SIp.flow = 0, SIv.flow = 0, 
                 SIph.flow = 0, SIpm.flow = 0, SIpl.flow = 0,
                 SIvh.flow = 0, SIvm.flow = 0, SIvl.flow = 0)
  control <- control.dcm(nsteps = 365*3, new.mod = si_ode)
  mod <- dcm(param, init, control)
  
  #mod <- mutate_epi(mod, rate.Placebo = (SIp.flow/Sp)*365*100) #Incidence, placebo
  #mod <- mutate_epi(mod, rate.Placebo.het = (SIp.flow/Sp)*365*100) #Incidence, placebo
  
  #sim.ve <- VE.from.sim(mod)
  sim.ve <- mod.manipulate(mod)
  sim.ve <- as.data.frame(sim.ve)
  matchedTimes <- sim.ve$time %in% target.stats$time
  #out <- sim.ve$VE.inst[matchedTimes]
  out <- sim.ve$rate.Placebo.het[matchedTimes]
  return(out)
}
```

### Specify priors for parameters (or initial conditions)

The `beta` and `c` parameters are the transmission rate per contact and the contact rate, respectively. For now we keep these together in the `lambda` parameter, due to non-identifiability. Below we just set priors for `lambda`. We also will want to set priors for some of the initial conditions eventually. 

```{r priors}
# In order of "x" in the above function.
priors  <- list(#c("unif", 0.003, 0.008),      # beta
                #c("unif", 60/365, 120/365))    # c 
                #c("unif", 0.000005, 0.0001))  # lambda
                c("unif", 0.01, 0.5))         # epsilon
                #c("unif", 1, 20))            # risk multiplier
```


```{r eval=FALSE, include=FALSE}

fit.seq <- ABC_sequential(method = "Lenormand",
                       model = f,
                       prior = priors,
                       summary_stat_target = target.stats$placebo.incidence,
#                       summary_stat_target = target.stats$VE.target,
                       nb_simul = 10,
                       p_acc = 0.10,
                       alpha = 0.5,
                       progress_bar = TRUE)

```


```{r eval=FALSE, include=TRUE}

fit.rej <- ABC_rejection(model = f,
                     prior = priors,
                     nb_simul = 1000,
                     summary_stat_target = target.stats$placebo.incidence,
                     tol = 0.25,
                     progress_bar = TRUE)

```


```{r eval=FALSE, include=TRUE}

par(mfrow = c(1, 2))

fit <- fit.rej

plot(density(fit$param[, 1], from = 0.0000001,  to = 0.0001),
     #0.0000001, 0.00009
     main = "epsilon", 
     xlim = c(0.0000001, 0.0001),
     #ylim = c(0, 10),
     col=2)
#lines(density(fit$param[, 1], from = 0.01,  to = 0.5), col = 2)
abline(v = VE.target[1], lty = 2, col = 1)
legend("topright", legend = c("Truth", "Posterior"),
      lty = c(1, 2), col = 1:2, lwd = 2)

plot(density(fit$param[, 2], from = 1, to = 20),
     main = "risk multiplier", ylim = c(0, 0.5))
#lines(density(fit3$param[, 2], from = 1, to = 20), col = 2)
#abline(v = 1.5, lty = 2)
legend("topright", legend = c("Posterior"),
      lty = 1, col = 1, lwd = 2)
```

Other examples from Sam

### Use the mean of the parameters for model selection

```{r}

ms <- colMeans(fit$param)
param <- param.dcm(tau = ms[1], c = ms[2], D = 7)
init <- init.dcm(s.num = 9999, i.num = 1, si.flow = 0, is.flow = 0)
control <- control.dcm(nsteps = 104, new.mod = SISmod)
sim <- dcm(param, init, control)

par(mfrow = c(1,1))
plot(sim, y = "P")
arrows(myDat$time, myDat$uci, myDat$time, myDat$lci,
       col = "grey", len = 0.025, angle = 90, code = 3)
points(myDat$time, myDat$sampPrev, col = "black", pch = 16, cex = 1)

### Or use the full posterior distribution, as before

param <- param.dcm(tau = fit3$param[, 1], c = fit3$param[, 2], D = 7)
sim <- dcm(param, init, control)

plot(sim, y = "P")
arrows(myDat$time, myDat$uci, myDat$time, myDat$lci,
       col = "grey", len = 0.025, angle = 90, code = 3)
points(myDat$time, myDat$sampPrev, col = "black", pch = 16, cex = 1)
```
