
pars <- c((1/20),1/3,1/4,1/5,1/10,74/100000,164/100000,0.15,0.08,0.02,0.01,0.60,0.4,0.14,0.35,0.62,0.20,6.4/100000,0.35,0.9*0.9,0.5*1,-log(0.99)/10,-log(0.96)/10)

pars_month <- as.vector(mode="numeric",length(pars))

for(i in 1:length(pars)){
  pars_month[i] <- -log(1-pars[i])/12
}

pars_month

plot(x=pars, y=pars_month)

