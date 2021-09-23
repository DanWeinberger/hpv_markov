#Specify some paramters
 time.step='month'
 WidthAgeClassMonth <- c(24,24,60) #how many time steps are in each age category
 DurationMatImmunityDays=60
 dur.days1=365
 dur.days2=180
 dur.days3=750
 N.people.age <- c( 1e4, 1e4, 1e5)
 
 #Initialize matrix y.init is the number of people in each state at t=1
 N.states <- 4 #How many compartments?
 N.ages <- 3 #How many age groups?
 N.times <- 120
 y.init <- matrix(0, nrow=N.ages, ncol=N.states)
 colnames(y.init) <-c('S0','I1','S1','I2')
 rownames(y.init) <-c('age1','age2','age3')
 y.init[,'I1'] <- 10 #10 infected people in each age group initially
y.init[,'S0'] <- N.people.age - y.init[,'I1'] #how many uninfected initially?
 
  States<-array(NA, dim=c(dim(y.init),N.times) )
  dimnames(States) <- c(dimnames(y.init),list((1:N.times)))
  
  States[,,1] <- y.init
  
  if(time.step=='month'){
    length.step=30.44 #days
  }else if(time.step=='week'){
    length.step=7 #days
  }
  
  omega = 1/(DurationMatImmunityDays/length.step)
  
  mu= 1/WidthAgeClassMonth
  if(time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)
  }
  
  p_1_0= 1/(dur.days1/length.step)  #converts 1/days to 1/lenth.step
  p_1_0= 1/(dur.days2/length.step)  
  gamma3= 1/(dur.days3/length.step)  
  gamma4= gamma3  
  
  #Pull out the states  for the model as matrices
  S0 <-  States[,'S0',]
  
  I1 <-  States[,'I1',]
  S1 <-  States[,'S1',]
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #B=annual birth rate
  
  #https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1004591.s016&type=supplementary
  #mu represents aging to the next class
  #um is death rate
  
  Aging.Prop <- mu
  
 for( i in 2:ncol(S0)){
  #Need to change to add in new people to S0 at some rate
  States[,'S0',i] <- States[,'S0',i-1] +
    (mu + um)*States[,'S0',i-1] + #
    Aging.Prop*States[,'S0',i-1] 
  
  States[,'I1',i] <-   lambda*States[,'I1',i-1] - 
    (gamma1 + mu + um)*States[,'I1',i-1] 

  States[,'S1',i] <- gamma1*States[,'S1',i-1] - 
    sigma1*lambda*States[,'I1',i-1] - 
    (mu+um)*States[,'I1',i-1] + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  
  
 }
  
  res <- list(States)
