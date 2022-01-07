### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Jan 7, 2021

###### Undetected and detected combined ######

#### Starting parameters
t = 1000
N.states = 8
Pop_size = 1000 

#### Transition probabilities

# aging into/out of cohort
age_in <- 1/20
age_out <- p_age_in

# dying of unrelated cause
die <- 0.01

# undetected
norm_ulsil <- 0.1
ulsil_norm <- 0.5
ulsil_uhsil <- 0.18
uhsil_ulsil <- 0.2
uhsil_ucan <- 0.04
  
# detected
norm_dlsil <- 0.1
dlsil_norm <- 0.5
dlsil_dhsil <- 0.18
dhsil_dlsil <- 0.2
dhsil_dcan <- 0.04
dcan_dcandeath <- 0.1
# treatment
dhsil_norm <- 0.9*0.9
dcan_norm <- 0.5*1

# BRFSS rates!!
# screening
ulsil_dlsil <-
uhsil_dhsil <-
ucan_dcan <-

# loss to follow up
dlsil_uhsil <-
dhsil_ucan <-



# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:t)
years = paste(prefix, suffix, sep=" ")

# Create empty array
mat1 = matrix(NA, t, N.states, dimnames=list(years, c("Normal","Undet_LSIL","Det_LSIL","Undet_HSIL","Det_HSIL","Undet_Cancer","Det_Cancer","Cancer Death")))

# Assign starting prevalence of each state
prev1 = 0.8
prev2 = 0.1
prev3 = 0
prev4 = 0.1
prev5 = 0
prev6 = 0
prev7 = 0
prev8 = 0

# Assign starting states
set.seed(123)
mat1[1,] <- rmultinom(1, Pop_size, prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 

## Need to update this
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

# For steady state need (Pop_size*p_age_in) + (mat1[(i-1),1])*p_1_1 + (mat1[(i-1),2])*p_2_1 =
#                        (mat1[(i-1),2])*p_2_1 + (mat1[(i-1),1])*p_1_2 + (mat1[(i-1),1])*p_age_out 

mat1 <- round(mat1,0)

#Check Progression over time

States <- apply(mat1, c(1,2), sum)

plot(States[,'Normal'], type='l')

plot(States[,'LSIL'], type='l')

plot(States[,'HSIL'], type='l')

plot(States[,'Cancer'], type='l')

plot(States[,'Cancer Death'], type='l')

pop_size_t <- apply(mat1,1, sum)

plot(pop_size_t)
