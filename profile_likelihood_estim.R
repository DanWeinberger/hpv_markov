### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Dec 30, 2021

### This model is for the undetected compartments ONLY

### Profile likelihood estimation

# Starting parameters
t = 1000
N.states = 5
Pop_size = 154566548 

# aging into cohort
p_age_in <- 1/20
p_age_out <- p_age_in
# Dying
p_die <- 0.01

# Create years label
prefix = "Year"
suffix = seq(1:t)
years = paste(prefix, suffix, sep=" ")

# Create empty array
mat1 = matrix(NA, t, N.states, dimnames=list(years, c("Normal","LSIL","HSIL","Cancer","Cancer Death")))

# Assign starting prevalence of each state
prev1 = 0.8
prev2 = 0.1
prev3 = 0.1
prev4 = 0
prev5 = 0

set.seed(123)


# Profile likelihood estimation function

# Evaluatte likelihood within function; only output log likelihood

profile <- function(p_1_2,p_2_1,p_2_3,p_3_2,p_3_4,p_4_5){

  # Starting parameters
  t = 1000
  N.states = 5
  Pop_size = 154566548 
  
  # aging into cohort
  p_age_in <- 1/20
  p_age_out <- p_age_in
  # Dying
  p_die <- 0.01
  
  # Create years label
  prefix = "Year"
  suffix = seq(1:t)
  years = paste(prefix, suffix, sep=" ")
  
  # Create empty array
  mat1 = matrix(NA, t, N.states, dimnames=list(years, c("Normal","LSIL","HSIL","Cancer","Cancer Death")))
  
  # Assign starting prevalence of each state
  prev1 = 0.8
  prev2 = 0.1
  prev3 = 0.1
  prev4 = 0
  prev5 = 0
  
  set.seed(123)
  
  # Assign starting states
  mat1[1,] <- rmultinom(1, Pop_size, prob=c(prev1,prev2,prev3,prev4,prev5)) 

  # Run model
  for(i in 2:t){
    # Normal
    mat1[i,1] <- (mat1[(i-1),1]) +  
      (Pop_size*p_age_in) +   
      (mat1[(i-1),2])*p_2_1 -   
      (mat1[(i-1),1])*p_1_2 - 
      (mat1[(i-1),1])*p_age_out - 
      (mat1[(i-1),1])*p_die
    # LSIL  
    mat1[i,2] <-   (mat1[(i-1),2]) + 
      (mat1[(i-1),1])*p_1_2 +
      (mat1[(i-1),3])*p_3_2 - 
      (mat1[(i-1),2])*p_2_1 - 
      (mat1[(i-1),2])*p_2_3 - 
      (mat1[(i-1),2])*p_age_out -
      (mat1[(i-1),2])*p_die
    # HSIL  
    mat1[i,3] <-   (mat1[(i-1),3]) + 
      (mat1[(i-1),2])*p_2_3 -   
      (mat1[(i-1),3])*p_3_2 - 
      (mat1[(i-1),3])*p_3_4 - 
      (mat1[(i-1),3])*p_age_out -
      (mat1[(i-1),3])*p_die
    # Cancer
    mat1[i,4] <- (mat1[(i-1),4]) +
      (mat1[(i-1),3])*p_3_4 -   
      (mat1[(i-1),4])*p_4_5 - 
      (mat1[(i-1),4])*p_age_out -
      (mat1[(i-1),4])*p_die
    # Cancer deaths
    mat1[i,5] <-
      (mat1[(i-1),4])*p_4_5        #progress from cancer
  
  }

  mat1 <- round(mat1,0)
  
  # Vector of HSIL, cancer, and cancer death at steady state
  result <- c(mat1[1000,3],mat1[1000,4],mat1[1000,5])
  
  # Log likelihood of result vs predicted
  dpois(result, c(216000,11778,3663), log=TRUE)
  
}


################### Run function with many inputs #######################

## Would like to use seq by 0.01, but uses too much memory

p_1_2_seq=seq(0.05,0.5,0.05) # Normal > LSIL
p_2_1_seq=seq(0.3,0.7,0.05) # LSIL > normal

p_2_3_seq=0.2 # LSIL > HSIL
p_3_2_seq=0.2 # HSIL > LSIL
p_3_4_seq=0.04 # HSIL > Cancer
p_4_5_seq=0.01 # Cancer > cancer death

#p_2_3_seq=seq(0.1,0.4,0.05) # LSIL > HSIL
#p_3_2_seq=seq(0.1,0.4,0.05) # HSIL > LSIL
#p_3_4_seq=seq(0.04,0.4,0.05) # HSIL > Cancer
#p_4_5_seq=seq(0.01,0.5,0.05) # Cancer > cancer death


# Create empty vector here (or array - 5 dimensional, i, j, k, l, m)
arr1 <- array(NA, dim=c(length(p_1_2_seq),length(p_2_1_seq),length(p_2_3_seq),length(p_3_2_seq),length(p_3_4_seq),length(p_4_5_seq)))
arr1 <- array()

# Fill in vector/array

# Time start
ptm <- proc.time()
# Run for loop
for(i in 1:length(p_1_2_seq)){
  for(j in 1:length(p_2_1_seq)){
    for(k in 1:length(p_2_3_seq)){
      for(l in 1:length(p_3_2_seq)){
        for(m in 1:length(p_3_4_seq)){
          for(n in 1:length(p_4_5_seq)){
            LL <- profile(p_1_2=i,p_2_1=j,p_2_3=k,p_3_2=l,p_3_4=m,p_4_5=n)
            LL_sum <- sum(LL)
          }
        }
      }
    }
  }
}
#Time end
proc.time() - ptm




# apply function?
params <- expand.grid(p_1_2_seq,p_2_1_seq,p_2_3_seq,p_3_2_seq,p_3_4_seq,p_4_5_seq)
# Need to change the names for this to work
library(purrr)
pmap(params, profile(p_1_2=..1,p_2_1=..2,p_2_3=..3,p_3_2=..4,p_3_4=..5,p_4_5=..6))
pmap(params, profile(p_1_2=p_1_2_seq,
                     p_2_1=p_2_1_seq,
                     p_2_3=p_2_3_seq,
                     p_3_2=p_3_2_seq,
                     p_3_4=p_3_4_seq,
                     p_4_5=p_4_5_seq))



