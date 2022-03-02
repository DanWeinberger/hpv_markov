## FINAL MODEL 2018 PARAMETERS ##

library(ggplot2)

### UNVACCINATED ###

#### Starting parameters
t = 1030
N.states = 8
#Pop_size = 154566548 #all women
#Pop_size = 44978865*(1-0.032) # women 18-39*adjustment for hysterectomies
Pop_size_1 = 14683822*0.4
Pop_size_2 = 10201392*0.4
Pop_size_3 = 19939084*0.4

# aging into cohort
age_in <- (1/20) # multiply by pop under 18
# aging up within cohort
age_up_0 <- 1/3 # proportion turning 21
age_up_1 <- 1/4 # proportion turning 25
age_up_2 <- 1/5 # proportion turning 30
# aging out of cohort
age_out <- 1/10 # proportion turning 40
# Dying
p_die_1 <- 74/100000 # proportion dying age 18-24
p_die_2_3 <- 164/100000 # proportion dying age 25-39

# undetected
norm_ulsil_0 <- 0.15
norm_ulsil_1 <- 0.08
norm_ulsil_2 <- 0.02
norm_ulsil_3 <- 0.01
ulsil_norm_1 <- 0.60
ulsil_norm_2_3 <- 0.4
ulsil_uhsil_1 <- 0.14
ulsil_uhsil_2_3 <- 0.35
uhsil_ulsil_1 <- 0.62
uhsil_ulsil_2_3 <- 0.20
uhsil_ucan <- 6.4/100000

# detected
norm_dlsil_1 <- norm_ulsil_1
norm_dlsil_2 <- norm_ulsil_2
norm_dlsil_3 <- norm_ulsil_3
dlsil_norm_1 <- ulsil_norm_1
dlsil_norm_2_3 <- ulsil_norm_2_3
dlsil_dhsil_1 <- ulsil_uhsil_1
dlsil_dhsil_2_3 <- ulsil_uhsil_2_3
dhsil_dlsil_1 <-  uhsil_ulsil_1
dhsil_dlsil_2_3 <-  uhsil_ulsil_2_3
dhsil_dcan <- uhsil_ucan
dcan_dcandeath <- 0.35

# treatment
dhsil_norm <- 0.9*0.9 # % success * % treated
dcan_norm <- 0.5*1 # % success (based on ~5 year survival) * % treated

# hysterectomies
hyst_1 <- -log(0.99)/10
hyst_2 <- -log(0.96)/10


# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:t)
years = paste(prefix, suffix, sep=" ")

# Create empty array
arr1 = array(NA, dim=c(t, N.states, 4), dimnames=list(years, c("Normal","Undet_LSIL","Det_LSIL","Undet_HSIL","Det_HSIL","Undet_Cancer","Det_Cancer","Cancer_Death"), c("18-20","21-24","25-29","30-39")))

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
arr1[1,,] <- rmultinom(1, Pop_size, prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 

# Run model
for(i in 2:t){
  
  ########################### Variable rates ########################
  
  t.index <- i
  
  # screening
  ulsil_dlsil <- 0.83/3
  uhsil_dhsil <- 0.83/3
  ucan_dcan <- 0.83/3
  
  # loss to follow up
  dlsil_uhsil <- 0.17/3
  dhsil_ucan <- 0.17/3
  
  
  ########################### 18-20 ########################
  
  # Normal 18-20 (1)
  arr1[i,1,1] <- (arr1[(i-1),1,1]) +  
    (Pop_size_1*age_in) +        # age in (need to adjust this)
    (arr1[(i-1),2,1])*ulsil_norm_1 - # undetected LSIL regress to normal
    #(arr1[(i-1),5,1])*dhsil_norm + # detected HSIL treated to normal
    #(arr1[(i-1),3,1])*dlsil_norm_1 - # detected LSIL regress to normal
    (arr1[(i-1),1,1])*norm_ulsil_1 - # normal progress to undetected LSIL
    #(arr1[(i-1),1,1])*norm_dlsil_1 - # normal progress to detected LSIL
    (arr1[(i-1),1,1])*age_up_0 -   # age up
    (arr1[(i-1),1,1])*p_die_1        # die (unrelated)
  # LSIL undetected 18-20 (2)
  arr1[i,2,1] <- (arr1[(i-1),2,1]) + 
    # Age in (need to add this)
    (arr1[(i-1),1,1])*norm_ulsil_0 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    #(arr1[(i-1),2,1])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,1])*age_up_0 -    # age up
    (arr1[(i-1),2,1])*p_die_1         # die (unrelated)
  # LSIL detected 18-20 (3) = 0
  arr1[i,3,1] <- (arr1[(i-1),3,1]) 
  # Age in (need to add this)
  #(arr1[(i-1),1,1])*norm_dlsil_1 +  # normal progress to detected LSIL
  #(arr1[(i-1),5,1])*dhsil_dlsil_1 + # detected HSIL regress to detected LSIL
  #(arr1[(i-1),2,1])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
  #(arr1[(i-1),3,1])*dlsil_norm_1 -  # detected LSIL regress to normal
  #(arr1[(i-1),3,1])*dlsil_dhsil_1 - # detected LSIL progress to detected HSIL
  #(arr1[(i-1),3,1])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
  #(arr1[(i-1),3,1])*age_up_1 -    # age up
  #(arr1[(i-1),3,1])*p_die_1         # die (unrelated)
  # HSIL undetected 18-20 (4) 
  arr1[i,4,1] <- (arr1[(i-1),4,1]) + 
    (arr1[(i-1),2,1])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,1])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    #(arr1[(i-1),4,1])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,1])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,1])*age_up_0 -    # age up
    (arr1[(i-1),4,1])*p_die_1         # die (unrelated)
  # HSIL detected 18-20 (5) = 0
  arr1[i,5,1] <- (arr1[(i-1),5,1]) 
  #(arr1[(i-1),3,1])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
  #(arr1[(i-1),4,1])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
  #(arr1[(i-1),5,1])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
  #(arr1[(i-1),5,1])*dhsil_norm -  # detected HSIL treated to normal
  #(arr1[(i-1),5,1])*age_up_1 -    # age up
  #(arr1[(i-1),5,1])*p_die_1         # die (unrelated)
  # Cancer undetected 18-20 (6) = 0
  arr1[i,6,1] <- (arr1[(i-1),6,1])
  # Cancer detected 18-20 (7) = 0
  arr1[i,7,1] <- (arr1[(i-1),7,1])
  # Cancer deaths 18-20 (8) = 0
  arr1[i,8,1] <- (arr1[(i-1),8,1])
  
  ########################### 21-24 ########################
  
  # Normal 21-24 (1)
  arr1[i,1,2] <- (arr1[(i-1),1,2]) +  
    (arr1[(i-1),1,1])*age_up_0 +   # age up
    (arr1[(i-1),2,2])*ulsil_norm_1 + # undetected LSIL regress to normal
    (arr1[(i-1),5,2])*dhsil_norm + # detected HSIL treated to normal
    (arr1[(i-1),3,2])*dlsil_norm_1 - # detected LSIL regress to normal
    (arr1[(i-1),1,2])*norm_ulsil_1 - # normal progress to undetected LSIL
    (arr1[(i-1),1,2])*norm_dlsil_1 - # normal progress to detected LSIL
    (arr1[(i-1),1,2])*age_up_1 -   # age up
    (arr1[(i-1),1,2])*hyst_1-        # hysterectomy
    (arr1[(i-1),1,2])*p_die_1        # die (unrelated)
  # LSIL undetected 21-24 (2)
  arr1[i,2,2] <- (arr1[(i-1),2,2]) + 
    (arr1[(i-1),2,1])*age_up_0 +   # age up
    (arr1[(i-1),1,2])*norm_ulsil_1 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,2])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,2])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,2])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,2])*age_up_1 -    # age up
    (arr1[(i-1),2,2])*hyst_1-        # hysterectomy
    (arr1[(i-1),2,2])*p_die_1         # die (unrelated)
  # LSIL detected 21-24 (3)
  arr1[i,3,2] <- (arr1[(i-1),3,2]) + 
    (arr1[(i-1),3,1])*age_up_0 +   # age up
    (arr1[(i-1),1,2])*norm_dlsil_1 +  # normal progress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,2])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,2])*dlsil_norm_1 -  # detected LSIL regress to normal
    (arr1[(i-1),3,2])*dlsil_dhsil_1 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,2])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,2])*age_up_1 -    # age up
    (arr1[(i-1),3,2])*hyst_1-        # hysterectomy
    (arr1[(i-1),3,2])*p_die_1         # die (unrelated)
  # HSIL undetected 21-24 (4) 
  arr1[i,4,2] <- (arr1[(i-1),4,2]) + 
    (arr1[(i-1),4,1])*age_up_0 +   # age up
    (arr1[(i-1),2,2])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,2])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    #(arr1[(i-1),4,1])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,2])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,2])*age_up_1 -    # age up
    (arr1[(i-1),4,2])*hyst_1-        # hysterectomy
    (arr1[(i-1),4,2])*p_die_1         # die (unrelated)
  # HSIL detected 21-24 (5)
  arr1[i,5,2] <- (arr1[(i-1),5,2]) + 
    (arr1[(i-1),5,1])*age_up_0 +   # age up
    (arr1[(i-1),3,2])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,2])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_norm -  # detected HSIL treated to normal
    (arr1[(i-1),5,2])*age_up_1 -    # age up
    (arr1[(i-1),5,2])*hyst_1-        # hysterectomy
    (arr1[(i-1),5,2])*p_die_1         # die (unrelated)
  # Cancer undetected 18-24 (6) = 0
  arr1[i,6,2] <- (arr1[(i-1),6,2])
  # Cancer detected 18-24 (7) = 0
  arr1[i,7,2] <- (arr1[(i-1),7,2])
  # Cancer deaths 18-24 (8) = 0
  arr1[i,8,2] <- (arr1[(i-1),8,2])
  
  ########################### 25-29 ########################
  
  # Normal 25-29 (1)
  arr1[i,1,3] <- (arr1[(i-1),1,3]) +  
    (arr1[(i-1),1,2])*age_up_1 +   # age up
    (arr1[(i-1),2,3])*ulsil_norm_2_3 + # undetected LSIL regress to normal
    (arr1[(i-1),5,3])*dhsil_norm + # detected HSIL treated to normal
    (arr1[(i-1),7,3])*dcan_norm +  # detected cancer treated to normal
    (arr1[(i-1),3,3])*dlsil_norm_2_3 - # detected LSIL regress to normal
    (arr1[(i-1),1,3])*norm_ulsil_2 - # normal progress to undetected LSIL
    (arr1[(i-1),1,3])*norm_dlsil_2 - # normal progress to detected LSIL
    (arr1[(i-1),1,3])*age_up_2 -   # age up 
    (arr1[(i-1),1,3])*hyst_1-        # hysterectomy
    (arr1[(i-1),1,3])*p_die_2_3      # die (unrelated)
  # LSIL undetected 25-29 (2)
  arr1[i,2,3] <- (arr1[(i-1),2,3]) + 
    (arr1[(i-1),2,2])*age_up_1 +   # age up
    (arr1[(i-1),1,3])*norm_ulsil_2 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,3])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,3])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,3])*age_up_2 -     # age up
    (arr1[(i-1),2,3])*hyst_1-        # hysterectomy
    (arr1[(i-1),2,3])*p_die_2_3       # die (unrelated)
  # LSIL detected 25-29 (3)
  arr1[i,3,3] <- (arr1[(i-1),3,3]) + 
    (arr1[(i-1),3,2])*age_up_1 +   # age up
    (arr1[(i-1),1,3])*norm_dlsil_2 +  # normal progress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,3])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,3])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,3])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,3])*age_up_2 -     # age up
    (arr1[(i-1),3,3])*hyst_1-        # hysterectomy
    (arr1[(i-1),3,3])*p_die_2_3       # die (unrelated)
  # HSIL undetected 25-29 (4) 
  arr1[i,4,3] <- (arr1[(i-1),4,3]) + 
    (arr1[(i-1),4,2])*age_up_1 +   # age up
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,3])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,3])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,3])*age_up_2 -     # age up
    (arr1[(i-1),4,3])*hyst_1-        # hysterectomy
    (arr1[(i-1),4,3])*p_die_2_3       # die (unrelated)
  # HSIL detected 25-29 (5)
  arr1[i,5,3] <- (arr1[(i-1),5,3]) + 
    (arr1[(i-1),5,2])*age_up_1 +   # age up
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,3])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_ucan -  # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,3])*dhsil_dcan -  # detected HSIL progress to detected cancer
    (arr1[(i-1),5,3])*dhsil_norm -  # detected HSIL treated to normal
    (arr1[(i-1),5,3])*age_up_2 -     # age up
    (arr1[(i-1),5,3])*hyst_1-        # hysterectomy
    (arr1[(i-1),5,3])*p_die_2_3       # die (unrelated)
  # Cancer undetected 25-29 (6)
  arr1[i,6,3] <- (arr1[(i-1),6,3]) +
    (arr1[(i-1),4,3])*uhsil_ucan +  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,3])*dhsil_ucan -  # detected HSIL progress to undetected cancer
    (arr1[(i-1),6,3])*ucan_dcan -   # undetected cancer screened to detected cancer
    (arr1[(i-1),6,3])*age_up_2 -     # age up
    (arr1[(i-1),6,3])*hyst_1-        # hysterectomy
    (arr1[(i-1),6,3])*p_die_2_3       # die (unrelated)
  # Cancer detected 25-29 (7)
  arr1[i,7,3] <- (arr1[(i-1),7,3]) +
    (arr1[(i-1),5,3])*dhsil_dcan +     # detected HSIL progress to detected cancer
    (arr1[(i-1),6,3])*ucan_dcan -      # undetected cancer screened to detected cancer
    (arr1[(i-1),7,3])*dcan_dcandeath - # detected cancer progress to cancer death
    (arr1[(i-1),7,3])*dcan_norm -      # detected cancer treated to normal
    (arr1[(i-1),7,3])*age_up_2 -        # age out
    (arr1[(i-1),7,3])*p_die_2_3          # die (unrelated)
  # Cancer deaths 25-29 (8)
  arr1[i,8,3] <-
    (arr1[(i-1),7,3])*dcan_dcandeath # detected cancer progress to cancer death
  
  ########################### 30-39 ########################
  
  # Normal 30-39 (1)
  arr1[i,1,4] <- (arr1[(i-1),1,4]) +  
    (arr1[(i-1),1,3])*age_up_2 +   # age up
    (arr1[(i-1),2,4])*ulsil_norm_2_3 + # undetected LSIL regress to normal
    (arr1[(i-1),5,4])*dhsil_norm + # detected HSIL treated to normal
    (arr1[(i-1),7,4])*dcan_norm +  # detected cancer treated to normal
    (arr1[(i-1),3,4])*dlsil_norm_2_3 - # detected LSIL regress to normal
    (arr1[(i-1),1,4])*norm_ulsil_3 - # normal progress to undetected LSIL
    (arr1[(i-1),1,4])*norm_dlsil_3 - # normal progress to detected LSIL
    (arr1[(i-1),1,4])*age_out -    # age out
    (arr1[(i-1),1,4])*hyst_2-        # hysterectomy
    (arr1[(i-1),1,4])*p_die_2_3      # die (unrelated)
  # LSIL undetected 25-39 (2)
  arr1[i,2,4] <- (arr1[(i-1),2,4]) + 
    (arr1[(i-1),2,3])*age_up_2 +   # age up
    (arr1[(i-1),1,4])*norm_ulsil_3 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,4])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,4])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,4])*age_out -     # age out
    (arr1[(i-1),2,4])*hyst_2-        # hysterectomy
    (arr1[(i-1),2,4])*p_die_2_3       # die (unrelated)
  # LSIL detected 25-39 (3)
  arr1[i,3,4] <- (arr1[(i-1),3,4]) + 
    (arr1[(i-1),3,3])*age_up_2 +   # age up
    (arr1[(i-1),1,4])*norm_dlsil_3 +  # normal progress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,4])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,4])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,4])*age_out -     # age out
    (arr1[(i-1),3,4])*hyst_2-        # hysterectomy
    (arr1[(i-1),3,4])*p_die_2_3       # die (unrelated)
  # HSIL undetected 25-39 (4) 
  arr1[i,4,4] <- (arr1[(i-1),4,4]) + 
    (arr1[(i-1),4,3])*age_up_2 +   # age up
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,4])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,4])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,4])*age_out -     # age out
    (arr1[(i-1),4,4])*hyst_2-        # hysterectomy
    (arr1[(i-1),4,4])*p_die_2_3       # die (unrelated)
  # HSIL detected 25-39 (5)
  arr1[i,5,4] <- (arr1[(i-1),5,4]) + 
    (arr1[(i-1),5,3])*age_up_2 +   # age up
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,4])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_ucan -  # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,4])*dhsil_dcan -  # detected HSIL progress to detected cancer
    (arr1[(i-1),5,4])*dhsil_norm -  # detected HSIL treated to normal
    (arr1[(i-1),5,4])*age_out -     # age out
    (arr1[(i-1),5,4])*hyst_2-        # hysterectomy
    (arr1[(i-1),5,4])*p_die_2_3       # die (unrelated)
  # Cancer undetected 25-39 (6)
  arr1[i,6,4] <- (arr1[(i-1),6,4]) +
    (arr1[(i-1),6,4])*age_up_2 +   # age up
    (arr1[(i-1),4,4])*uhsil_ucan +  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,4])*dhsil_ucan -  # detected HSIL progress to undetected cancer
    (arr1[(i-1),6,4])*ucan_dcan -   # undetected cancer screened to detected cancer
    (arr1[(i-1),6,4])*age_out -     # age out
    (arr1[(i-1),6,4])*hyst_2-        # hysterectomy
    (arr1[(i-1),6,4])*p_die_2_3       # die (unrelated)
  # Cancer detected 25-39 (7)
  arr1[i,7,4] <- (arr1[(i-1),7,4]) +
    (arr1[(i-1),7,3])*age_up_2 +   # age up
    (arr1[(i-1),5,4])*dhsil_dcan +     # detected HSIL progress to detected cancer
    (arr1[(i-1),6,4])*ucan_dcan -      # undetected cancer screened to detected cancer
    (arr1[(i-1),7,4])*dcan_dcandeath - # detected cancer progress to cancer death
    (arr1[(i-1),7,4])*dcan_norm -      # detected cancer treated to normal
    (arr1[(i-1),7,4])*age_out -        # age out
    (arr1[(i-1),7,4])*p_die_2_3          # die (unrelated)
  # Cancer deaths 25-39 (8)
  arr1[i,8,4] <-
    (arr1[(i-1),7,4])*dcan_dcandeath # detected cancer progress to cancer death
}

arr1 <- round(arr1,0)

# Vector of HSIL, cancer, and cancer death at steady state for iteration
result_noCov <- apply(arr1, 2L, rowSums)
final_result_noCov <- as.data.frame(result_noCov[nrow(result_noCov),])

# Log likelihood of result vs predicted
LL <- sum(dpois(final_result_noCov[c(5,7,8),], c(196000*0.4,10510*0.4,3400*0.4), log=TRUE))

Results_noCov <- final_result_noCov[c(5,7,8),]

Results_3_noCov <- final_result_noCov[c(5,6,8),]

# Plot results

result_2_noCov <- as.data.frame(result_noCov)
result_3_noCov <- as.data.frame(tibble::rownames_to_column(result_2_noCov, "Year"))
result_3_noCov$Year.No <- seq(1:t)



### VACCINATED ###

#### Starting parameters
t = 1030
N.states = 8
#Pop_size = 154566548 #all women
#Pop_size = 44978865*(1-0.032) # women 18-39*adjustment for hysterectomies
Pop_size_1 = 14683822*0.6
Pop_size_2 = 10201392*0.6
Pop_size_3 = 19939084*0.6

# aging into cohort
age_in <- (1/20) # multiply by pop under 18
# aging up within cohort
age_up_0 <- 1/3 # proportion turning 21
age_up_1 <- 1/4 # proportion turning 25
age_up_2 <- 1/5 # proportion turning 30
# aging out of cohort
age_out <- 1/10 # proportion turning 40
# Dying
p_die_1 <- 74/100000 # proportion dying age 18-24
p_die_2_3 <- 164/100000 # proportion dying age 25-39

# undetected
#ifelse(t<990, norm_ulsil_0 <- 0.15, norm_ulsil_0 <- 0.15*0.2)
#ifelse(t<990, norm_ulsil_1 <- 0.08, norm_ulsil_1 <- 0.08*0.2)
#ifelse(t<990, norm_ulsil_2 <- 0.02, norm_ulsil_2 <- 0.02*0.2)
norm_ulsil_3 <- 0.01*0.1
ulsil_norm_1 <- 0.60
ulsil_norm_2_3 <- 0.4
ulsil_uhsil_1 <- 0.14
ulsil_uhsil_2_3 <- 0.35
uhsil_ulsil_1 <- 0.62
uhsil_ulsil_2_3 <- 0.20
uhsil_ucan <- 6.4/100000

# detected
norm_dlsil_1 <- norm_ulsil_1
norm_dlsil_2 <- norm_ulsil_2
norm_dlsil_3 <- norm_ulsil_3
dlsil_norm_1 <- ulsil_norm_1
dlsil_norm_2_3 <- ulsil_norm_2_3
dlsil_dhsil_1 <- ulsil_uhsil_1
dlsil_dhsil_2_3 <- ulsil_uhsil_2_3
dhsil_dlsil_1 <-  uhsil_ulsil_1
dhsil_dlsil_2_3 <-  uhsil_ulsil_2_3
dhsil_dcan <- uhsil_ucan
dcan_dcandeath <- 0.35

# treatment
dhsil_norm <- 0.9*0.9
dcan_norm <- 0.5*1

# hysterectomies
hyst_1 <- -log(0.99)/10
hyst_2 <- -log(0.96)/10

# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:t)
years = paste(prefix, suffix, sep=" ")

# Create empty array
arr1 = array(NA, dim=c(t, N.states, 4), dimnames=list(years, c("Normal","Undet_LSIL","Det_LSIL","Undet_HSIL","Det_HSIL","Undet_Cancer","Det_Cancer","Cancer Death"), c("18-20","21-24","25-29","30-39")))

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
arr1[1,,] <- rmultinom(1, Pop_size, prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 

# Run model
for(i in 2:t){
  
  ########################### Variable rates ########################
  
  t.index <- i
  
  norm_ulsil_0 <- ifelse(t.index<988, 0.15, 0.15*0.2)
  norm_ulsil_1 <- ifelse(t.index<988, 0.08, 0.08*0.2)
  norm_ulsil_2 <- ifelse(t.index<988, 0.02, 0.02*0.2)
  
  # screening
  ulsil_dlsil <- 0.83/3
  uhsil_dhsil <- 0.83/3
  ucan_dcan <- 0.83/3
  
  # loss to follow up
  dlsil_uhsil <- 0.17/3
  dhsil_ucan <- 0.17/3
  
  
  ########################### 18-20 ########################
  
  # Normal 18-20 (1)
  arr1[i,1,1] <- (arr1[(i-1),1,1]) +  
    (Pop_size_1*age_in) +        # age in (need to adjust this)
    (arr1[(i-1),2,1])*ulsil_norm_1 - # undetected LSIL regress to normal
    #(arr1[(i-1),5,1])*dhsil_norm + # detected HSIL treated to normal
    #(arr1[(i-1),3,1])*dlsil_norm_1 - # detected LSIL regress to normal
    (arr1[(i-1),1,1])*norm_ulsil_1 - # normal progress to undetected LSIL
    #(arr1[(i-1),1,1])*norm_dlsil_1 - # normal progress to detected LSIL
    (arr1[(i-1),1,1])*age_up_0 -   # age up
    (arr1[(i-1),1,1])*p_die_1        # die (unrelated)
  # LSIL undetected 18-20 (2)
  arr1[i,2,1] <- (arr1[(i-1),2,1]) + 
    # Age in (need to add this)
    (arr1[(i-1),1,1])*norm_ulsil_0 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    #(arr1[(i-1),2,1])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,1])*age_up_0 -    # age up
    (arr1[(i-1),2,1])*p_die_1         # die (unrelated)
  # LSIL detected 18-20 (3) = 0
  arr1[i,3,1] <- (arr1[(i-1),3,1]) 
  # Age in (need to add this)
  #(arr1[(i-1),1,1])*norm_dlsil_1 +  # normal progress to detected LSIL
  #(arr1[(i-1),5,1])*dhsil_dlsil_1 + # detected HSIL regress to detected LSIL
  #(arr1[(i-1),2,1])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
  #(arr1[(i-1),3,1])*dlsil_norm_1 -  # detected LSIL regress to normal
  #(arr1[(i-1),3,1])*dlsil_dhsil_1 - # detected LSIL progress to detected HSIL
  #(arr1[(i-1),3,1])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
  #(arr1[(i-1),3,1])*age_up_1 -    # age up
  #(arr1[(i-1),3,1])*p_die_1         # die (unrelated)
  # HSIL undetected 18-20 (4) 
  arr1[i,4,1] <- (arr1[(i-1),4,1]) + 
    (arr1[(i-1),2,1])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,1])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    #(arr1[(i-1),4,1])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,1])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,1])*age_up_0 -    # age up
    (arr1[(i-1),4,1])*p_die_1         # die (unrelated)
  # HSIL detected 18-20 (5) = 0
  arr1[i,5,1] <- (arr1[(i-1),5,1]) 
  #(arr1[(i-1),3,1])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
  #(arr1[(i-1),4,1])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
  #(arr1[(i-1),5,1])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
  #(arr1[(i-1),5,1])*dhsil_norm -  # detected HSIL treated to normal
  #(arr1[(i-1),5,1])*age_up_1 -    # age up
  #(arr1[(i-1),5,1])*p_die_1         # die (unrelated)
  # Cancer undetected 18-20 (6) = 0
  arr1[i,6,1] <- (arr1[(i-1),6,1])
  # Cancer detected 18-20 (7) = 0
  arr1[i,7,1] <- (arr1[(i-1),7,1])
  # Cancer deaths 18-20 (8) = 0
  arr1[i,8,1] <- (arr1[(i-1),8,1])
  
  ########################### 21-24 ########################
  
  # Normal 21-24 (1)
  arr1[i,1,2] <- (arr1[(i-1),1,2]) +  
    (arr1[(i-1),1,1])*age_up_0 +   # age up
    (arr1[(i-1),2,2])*ulsil_norm_1 + # undetected LSIL regress to normal
    (arr1[(i-1),5,2])*dhsil_norm + # detected HSIL treated to normal
    (arr1[(i-1),3,2])*dlsil_norm_1 - # detected LSIL regress to normal
    (arr1[(i-1),1,2])*norm_ulsil_1 - # normal progress to undetected LSIL
    (arr1[(i-1),1,2])*norm_dlsil_1 - # normal progress to detected LSIL
    (arr1[(i-1),1,2])*age_up_1 -   # age up
    (arr1[(i-1),1,2])*p_die_1        # die (unrelated)
  # LSIL undetected 21-24 (2)
  arr1[i,2,2] <- (arr1[(i-1),2,2]) + 
    (arr1[(i-1),2,1])*age_up_0 +   # age up
    (arr1[(i-1),1,2])*norm_ulsil_1 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,2])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,2])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,2])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,2])*age_up_1 -    # age up
    (arr1[(i-1),2,2])*p_die_1         # die (unrelated)
  # LSIL detected 21-24 (3)
  arr1[i,3,2] <- (arr1[(i-1),3,2]) + 
    (arr1[(i-1),3,1])*age_up_0 +   # age up
    (arr1[(i-1),1,2])*norm_dlsil_1 +  # normal progress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,2])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,2])*dlsil_norm_1 -  # detected LSIL regress to normal
    (arr1[(i-1),3,2])*dlsil_dhsil_1 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,2])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,2])*age_up_1 -    # age up
    (arr1[(i-1),3,2])*p_die_1         # die (unrelated)
  # HSIL undetected 21-24 (4) 
  arr1[i,4,2] <- (arr1[(i-1),4,2]) + 
    (arr1[(i-1),4,1])*age_up_0 +   # age up
    (arr1[(i-1),2,2])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,2])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    #(arr1[(i-1),4,1])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,2])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,2])*age_up_1 -    # age up
    (arr1[(i-1),4,2])*p_die_1         # die (unrelated)
  # HSIL detected 21-24 (5)
  arr1[i,5,2] <- (arr1[(i-1),5,2]) + 
    (arr1[(i-1),5,1])*age_up_0 +   # age up
    (arr1[(i-1),3,2])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,2])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_norm -  # detected HSIL treated to normal
    (arr1[(i-1),5,2])*age_up_1 -    # age up
    (arr1[(i-1),5,2])*p_die_1         # die (unrelated)
  # Cancer undetected 18-24 (6) = 0
  arr1[i,6,2] <- (arr1[(i-1),6,2])
  # Cancer detected 18-24 (7) = 0
  arr1[i,7,2] <- (arr1[(i-1),7,2])
  # Cancer deaths 18-24 (8) = 0
  arr1[i,8,2] <- (arr1[(i-1),8,2])
  
  ########################### 25-29 ########################
  
  # Normal 25-29 (1)
  arr1[i,1,3] <- (arr1[(i-1),1,3]) +  
    (arr1[(i-1),1,2])*age_up_1 +   # age up
    (arr1[(i-1),2,3])*ulsil_norm_2_3 + # undetected LSIL regress to normal
    (arr1[(i-1),5,3])*dhsil_norm + # detected HSIL treated to normal
    (arr1[(i-1),7,3])*dcan_norm +  # detected cancer treated to normal
    (arr1[(i-1),3,3])*dlsil_norm_2_3 - # detected LSIL regress to normal
    (arr1[(i-1),1,3])*norm_ulsil_2 - # normal progress to undetected LSIL
    (arr1[(i-1),1,3])*norm_dlsil_2 - # normal progress to detected LSIL
    (arr1[(i-1),1,3])*age_up_2 -   # age up 
    (arr1[(i-1),1,3])*p_die_2_3      # die (unrelated)
  # LSIL undetected 25-29 (2)
  arr1[i,2,3] <- (arr1[(i-1),2,3]) + 
    (arr1[(i-1),2,2])*age_up_1 +   # age up
    (arr1[(i-1),1,3])*norm_ulsil_2 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,3])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,3])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,3])*age_up_2 -     # age out
    (arr1[(i-1),2,3])*p_die_2_3       # die (unrelated)
  # LSIL detected 25-29 (3)
  arr1[i,3,3] <- (arr1[(i-1),3,3]) + 
    (arr1[(i-1),3,2])*age_up_1 +   # age up
    (arr1[(i-1),1,3])*norm_dlsil_2 +  # normal progress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,3])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,3])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,3])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,3])*age_up_2 -     # age out
    (arr1[(i-1),3,3])*p_die_2_3       # die (unrelated)
  # HSIL undetected 25-29 (4) 
  arr1[i,4,3] <- (arr1[(i-1),4,3]) + 
    (arr1[(i-1),4,2])*age_up_1 +   # age up
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,3])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,3])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,3])*age_up_2 -     # age out
    (arr1[(i-1),4,3])*p_die_2_3       # die (unrelated)
  # HSIL detected 25-29 (5)
  arr1[i,5,3] <- (arr1[(i-1),5,3]) + 
    (arr1[(i-1),5,2])*age_up_1 +   # age up
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,3])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_ucan -  # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,3])*dhsil_dcan -  # detected HSIL progress to detected cancer
    (arr1[(i-1),5,3])*dhsil_norm -  # detected HSIL treated to normal
    (arr1[(i-1),5,3])*age_up_2 -     # age out
    (arr1[(i-1),5,3])*p_die_2_3       # die (unrelated)
  # Cancer undetected 25-29 (6)
  arr1[i,6,3] <- (arr1[(i-1),6,3]) +
    (arr1[(i-1),4,3])*uhsil_ucan +  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,3])*dhsil_ucan -  # detected HSIL progress to undetected cancer
    (arr1[(i-1),6,3])*ucan_dcan -   # undetected cancer screened to detected cancer
    (arr1[(i-1),6,3])*age_up_2 -     # age out
    (arr1[(i-1),6,3])*p_die_2_3       # die (unrelated)
  # Cancer detected 25-29 (7)
  arr1[i,7,3] <- (arr1[(i-1),7,3]) +
    (arr1[(i-1),5,3])*dhsil_dcan +     # detected HSIL progress to detected cancer
    (arr1[(i-1),6,3])*ucan_dcan -      # undetected cancer screened to detected cancer
    (arr1[(i-1),7,3])*dcan_dcandeath - # detected cancer progress to cancer death
    (arr1[(i-1),7,3])*dcan_norm -      # detected cancer treated to normal
    (arr1[(i-1),7,3])*age_up_2 -        # age out
    (arr1[(i-1),7,3])*p_die_2_3          # die (unrelated)
  # Cancer deaths 25-29 (8)
  arr1[i,8,3] <-
    (arr1[(i-1),7,3])*dcan_dcandeath # detected cancer progress to cancer death
  
  ########################### 30-39 ########################
  
  # Normal 30-39 (1)
  arr1[i,1,4] <- (arr1[(i-1),1,4]) +  
    (arr1[(i-1),1,3])*age_up_2 +   # age up
    (arr1[(i-1),2,4])*ulsil_norm_2_3 + # undetected LSIL regress to normal
    (arr1[(i-1),5,4])*dhsil_norm + # detected HSIL treated to normal
    (arr1[(i-1),7,4])*dcan_norm +  # detected cancer treated to normal
    (arr1[(i-1),3,4])*dlsil_norm_2_3 - # detected LSIL regress to normal
    (arr1[(i-1),1,4])*norm_ulsil_3 - # normal progress to undetected LSIL
    (arr1[(i-1),1,4])*norm_dlsil_3 - # normal progress to detected LSIL
    (arr1[(i-1),1,4])*age_out -    # age out
    (arr1[(i-1),1,4])*p_die_2_3      # die (unrelated)
  # LSIL undetected 25-39 (2)
  arr1[i,2,4] <- (arr1[(i-1),2,4]) + 
    (arr1[(i-1),2,3])*age_up_2 +   # age up
    (arr1[(i-1),1,4])*norm_ulsil_3 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,4])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,4])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,4])*age_out -     # age out
    (arr1[(i-1),2,4])*p_die_2_3       # die (unrelated)
  # LSIL detected 25-39 (3)
  arr1[i,3,4] <- (arr1[(i-1),3,4]) + 
    (arr1[(i-1),3,3])*age_up_2 +   # age up
    (arr1[(i-1),1,4])*norm_dlsil_3 +  # normal progress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,4])*ulsil_dlsil - # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,4])*dlsil_uhsil - # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,4])*age_out -     # age out
    (arr1[(i-1),3,4])*p_die_2_3       # die (unrelated)
  # HSIL undetected 25-39 (4) 
  arr1[i,4,4] <- (arr1[(i-1),4,4]) + 
    (arr1[(i-1),4,3])*age_up_2 +   # age up
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,4])*dlsil_uhsil - # detected LSIL progress to undetected HSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ucan -  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,4])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,4])*age_out -     # age out
    (arr1[(i-1),4,4])*p_die_2_3       # die (unrelated)
  # HSIL detected 25-39 (5)
  arr1[i,5,4] <- (arr1[(i-1),5,4]) + 
    (arr1[(i-1),5,3])*age_up_2 +   # age up
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,4])*uhsil_dhsil - # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_ucan -  # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,4])*dhsil_dcan -  # detected HSIL progress to detected cancer
    (arr1[(i-1),5,4])*dhsil_norm -  # detected HSIL treated to normal
    (arr1[(i-1),5,4])*age_out -     # age out
    (arr1[(i-1),5,4])*p_die_2_3       # die (unrelated)
  # Cancer undetected 25-39 (6)
  arr1[i,6,4] <- (arr1[(i-1),6,4]) +
    (arr1[(i-1),6,4])*age_up_2 +   # age up
    (arr1[(i-1),4,4])*uhsil_ucan +  # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,4])*dhsil_ucan -  # detected HSIL progress to undetected cancer
    (arr1[(i-1),6,4])*ucan_dcan -   # undetected cancer screened to detected cancer
    (arr1[(i-1),6,4])*age_out -     # age out
    (arr1[(i-1),6,4])*p_die_2_3       # die (unrelated)
  # Cancer detected 25-39 (7)
  arr1[i,7,4] <- (arr1[(i-1),7,4]) +
    (arr1[(i-1),7,3])*age_up_2 +   # age up
    (arr1[(i-1),5,4])*dhsil_dcan +     # detected HSIL progress to detected cancer
    (arr1[(i-1),6,4])*ucan_dcan -      # undetected cancer screened to detected cancer
    (arr1[(i-1),7,4])*dcan_dcandeath - # detected cancer progress to cancer death
    (arr1[(i-1),7,4])*dcan_norm -      # detected cancer treated to normal
    (arr1[(i-1),7,4])*age_out -        # age out
    (arr1[(i-1),7,4])*p_die_2_3          # die (unrelated)
  # Cancer deaths 25-39 (8)
  arr1[i,8,4] <-
    (arr1[(i-1),7,4])*dcan_dcandeath # detected cancer progress to cancer death
}

arr1 <- round(arr1,0)

# Vector of HSIL, cancer, and cancer death at steady state for iteration
result_5_noCov <- apply(arr1, 2L, rowSums)
final_result_2_noCov <- as.data.frame(result_5_noCov[nrow(result_5_noCov),])

# Log likelihood of result vs predicted
LL_2 <- sum(dpois(final_result_2_noCov[c(5,7,8),], c(196000*0.6,10510*0.6,3400*0.6), log=TRUE))

Results_2_noCov <- final_result_2_noCov[c(5,7,8),]

Results_4_noCov <- final_result_2_noCov[c(4,6,8),]


# Plot results

result_6_noCov <- as.data.frame(result_5_noCov)
result_7_noCov <- as.data.frame(tibble::rownames_to_column(result_6_noCov, "Year"))
result_7_noCov$Year.No <- seq(1:t)




################### COMBINED RESULTS ###################

Results_all_noCov = Results_noCov+Results_2_noCov # final results

result_tot_noCov <- result_3_noCov[,-1] + result_7_noCov[,-1]
#View(result_tot)

result_tot_noCov$Year.No <- result_tot_noCov$Year.No/2
result_tot_noCov$Year <- result_tot_noCov$Year.No+(2008-987)


# Normal
p_Norm<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Normal), colour="red", alpha=0.5) + 
  geom_line(data=result_tot, aes(y=Normal), colour="red") + 
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,20000000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3)
p_Norm


# LSIL
p_LSIL<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Undet_LSIL), colour="red", alpha=0.5) + 
  geom_line(data=result_tot, aes(y=Undet_LSIL), colour="red") + 
  geom_line(data=result_tot_noCov, aes(y=Det_LSIL), colour="blue", alpha=0.5) +
  geom_line(data=result_tot, aes(y=Det_LSIL), colour="blue") +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,1000000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3)
p_LSIL


# HSIL
p_HSIL<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Undet_HSIL), colour="red", alpha=0.5) + 
  geom_line(data=result_tot, aes(y=Undet_HSIL), colour="red") + 
  geom_line(data=result_tot_noCov, aes(y=Det_HSIL), colour="blue", alpha=0.5) +
  geom_line(data=result_tot, aes(y=Det_HSIL), colour="blue") +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,500000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3)
p_HSIL


p_Cancer<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Undet_Cancer), colour="red", alpha=0.5) + 
  geom_line(data=result_tot, aes(y=Undet_Cancer), colour="red") + 
  geom_line(data=result_tot_noCov, aes(y=Det_Cancer), colour="blue", alpha=0.5) +
  geom_line(data=result_tot, aes(y=Det_Cancer), colour="blue") +
  geom_line(data=result_tot_noCov, aes(y=Cancer_Death), colour="black", alpha=0.5) +
  geom_line(data=result_tot, aes(y=Cancer_Death), colour="black") +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,40000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3)
p_Cancer


## Integration attempt

test <- integrate(approxfun(result_tot$Year, result_tot$Undet_HSIL), 2020,2050)

test <- integrate(splinefun(result_tot$Year, result_tot$Undet_HSIL), 2020,2050)
# This gives lower absolute error



