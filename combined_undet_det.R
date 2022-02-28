### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Jan 23, 2021

###### Undetected and detected combined ######

#### Starting parameters
t = 1000
N.states = 8
Pop_size = 154566548

#### Transition probabilities

# aging into/out of cohort
age_in <- 1/20
age_out <- age_in

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
ulsil_dlsil <- 0.83/3
uhsil_dhsil <- 0.83/3
ucan_dcan <- 0.83/3

# loss to follow up
dlsil_uhsil <- 0.17/3
dhsil_ucan <- 0.17/3


profile.likelihood
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

# Run model
for(i in 2:t){
  # Normal (1)
  mat1[i,1] <- (mat1[(i-1),1]) +  
    (Pop_size*age_in) +        # age in (need to adjust this)
    (mat1[(i-1),2])*ulsil_norm + # undetected LSIL regress to normal
    (mat1[(i-1),5])*dhsil_norm + # detected HSIL treated to normal
    (mat1[(i-1),7])*dcan_norm +  # detected cancer treated to normal
    (mat1[(i-1),3])*dlsil_norm - # detected LSIL regress to normal
    (mat1[(i-1),1])*norm_ulsil - # normal progress to undetected LSIL
    (mat1[(i-1),1])*norm_dlsil - # normal progress to detected LSIL
    (mat1[(i-1),1])*age_out -  # age out
    (mat1[(i-1),1])*die        # die (unrelated)
  # LSIL undetected (2)
  mat1[i,2] <- (mat1[(i-1),2]) + 
    # Age in (need to add this)
    (mat1[(i-1),1])*norm_ulsil +  # normal progress to undetected LSIL
    (mat1[(i-1),4])*uhsil_ulsil - # undetected HSIL regress to undetected LSIL
    (mat1[(i-1),2])*ulsil_norm -  # undetected LSIL regress to normal
    (mat1[(i-1),2])*ulsil_uhsil - # undetected LSIL progress to undetected HSIL
    (mat1[(i-1),2])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (mat1[(i-1),2])*age_out -   # age out
    (mat1[(i-1),2])*die         # die (unrelated)
  # LSIL detected (3)
  mat1[i,3] <- (mat1[(i-1),3]) + 
    # Age in (need to add this)
    (mat1[(i-1),1])*norm_dlsil +  # normal progress to detected LSIL
    (mat1[(i-1),5])*dhsil_dlsil + # detected HSIL regress to detected LSIL
    (mat1[(i-1),2])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (mat1[(i-1),3])*dlsil_norm -  # detected LSIL regress to normal
    (mat1[(i-1),3])*dlsil_dhsil - # detected LSIL progress to detected HSIL
    (mat1[(i-1),3])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
    (mat1[(i-1),3])*age_out -   # age out
    (mat1[(i-1),3])*die         # die (unrelated)
  # HSIL undetected (4) 
  mat1[i,4] <- (mat1[(i-1),4]) + 
    (mat1[(i-1),2])*ulsil_uhsil + # undetected LSIL progress to undetected HSIL 
    (mat1[(i-1),3])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (mat1[(i-1),4])*uhsil_ulsil - # undetected HSIL regress to undetected LSIL
    (mat1[(i-1),4])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (mat1[(i-1),4])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (mat1[(i-1),4])*age_out -   # age out
    (mat1[(i-1),4])*die         # die (unrelated)
  # HSIL detected (5)
  mat1[i,5] <- (mat1[(i-1),5]) + 
    (mat1[(i-1),3])*dlsil_dhsil + # detected LSIL progress to detected HSIL
    (mat1[(i-1),4])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (mat1[(i-1),5])*dhsil_dlsil - # detected HSIL regress to detected LSIL
    (mat1[(i-1),5])*dhsil_ucan -  # detected HSIL LTFU to undetected cancer
    (mat1[(i-1),5])*dhsil_dcan -  # detected HSIL progress to detected cancer
    (mat1[(i-1),5])*dhsil_norm -  # detected HSIL treated to normal
    (mat1[(i-1),5])*age_out -   # age out
    (mat1[(i-1),5])*die         # die (unrelated)
  # Cancer undetected (6)
  mat1[i,6] <- (mat1[(i-1),6]) +
    (mat1[(i-1),4])*uhsil_ucan +  # undetected HSIL progress to undetected cancer
    (mat1[(i-1),5])*dhsil_ucan -  # detected HSIL progress to undetected cancer
    (mat1[(i-1),6])*ucan_dcan -   # undetected cancer screened to detected cancer
    (mat1[(i-1),6])*age_out -   # age out
    (mat1[(i-1),6])*die         # die (unrelated)
  # Cancer detected (7)
  mat1[i,7] <- (mat1[(i-1),7]) +
    (mat1[(i-1),5])*dhsil_dcan +     # detected HSIL progress to detected cancer
    (mat1[(i-1),6])*ucan_dcan -      # undetected cancer screened to detected cancer
    (mat1[(i-1),7])*dcan_dcandeath - # detected cancer progress to cancer death
    (mat1[(i-1),7])*dcan_norm -      # detected cancer treated to normal
    (mat1[(i-1),7])*age_out -      # age out
    (mat1[(i-1),7])*die            # die (unrelated)
  # Cancer deaths (8)
  mat1[i,8] <-
    (mat1[(i-1),7])*dcan_dcandeath # detected cancer progress to cancer death
  
}

mat1 <- round(mat1,0)

#Check Progression over time

States <- apply(mat1, c(1,2), sum)

plot(States[,'Normal'], type='l')

plot(States[,'Undet_LSIL'], type='l')

plot(States[,'Det_LSIL'], type='l')

plot(States[,'Undet_HSIL'], type='l')

plot(States[,'Det_HSIL'], type='l')

plot(States[,'Undet_Cancer'], type='l')

plot(States[,'Det_Cancer'], type='l')

plot(States[,'Cancer Death'], type='l')

pop_size_t <- apply(mat1,1, sum)

plot(pop_size_t)
