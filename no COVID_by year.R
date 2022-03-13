## FINAL MODEL 2018 PARAMETERS ##

library(ggplot2)

#################### UNVACCINATED ####################

#### Starting parameters ####
t = 1030 # 1030 years
N.states = 8

# Population size
# 154566548 = all women
Pop_size_0 = 35704873*0.4 # <18*proportion unvaxxed
Pop_size_1 = 14683822*0.4 # 18-24*proportion unvaxxed
Pop_size_2 = 10201392*0.4 # 25-29*proportion unvaxxed
Pop_size_3 = 19939084*0.4 # 30-39*proportion unvaxxed

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
norm_ulsil_0 <- 0.15 # normal > undetected LSIL (18-20)
norm_ulsil_1 <- 0.08 # normal > undetected LSIL (21-24)
norm_ulsil_2 <- 0.02 # normal > undetected LSIL (25-29)
norm_ulsil_3 <- 0.01 # normal > undetected LSIL (30-39)
ulsil_norm_1 <- 0.60 # undetected LSIL > normal (18-24)
ulsil_norm_2_3 <- 0.4 # undetected LSIL > normal (25-39)
ulsil_uhsil_1 <- 0.14 # undetected LSIL > undetected HSIL (18-24)
ulsil_uhsil_2_3 <- 0.35 # undetected LSIL > undetected HSIL (25-39)
uhsil_ulsil_1 <- 0.62 # undetected HSIL > undetected LSIL (18-24)
uhsil_ulsil_2_3 <- 0.20 # undetected HSIL > undetected LSIL (25-39)
uhsil_ucan <- 6.4/100000 # undetected HSIL > undetected cancer

# detected
norm_dlsil_1 <- norm_ulsil_1 # normal > detected LSIL (21-24)
norm_dlsil_2 <- norm_ulsil_2 # normal > detected LSIL (24-29)
norm_dlsil_3 <- norm_ulsil_3 # normal > detected LSIL (30-39)
dlsil_norm_1 <- ulsil_norm_1 # detected LSIL > normal (21-24)
dlsil_norm_2_3 <- ulsil_norm_2_3 # detected LSIL > normal (25-39)
dlsil_dhsil_1 <- ulsil_uhsil_1 # detected LSIL > detected HSIL (21-24)
dlsil_dhsil_2_3 <- ulsil_uhsil_2_3 # detected LSIL > detected HSIL (25-39)
dhsil_dlsil_1 <-  uhsil_ulsil_1 # detected HSIL > detected LSIL (21-24)
dhsil_dlsil_2_3 <-  uhsil_ulsil_2_3 # detected HSIL > detected LSIL (25-39)
dhsil_dcan <- uhsil_ucan # detected HSIL > detected cancer
dcan_dcandeath <- 0.35 # detected cancer > detected cancer death

# treatment
dhsil_norm <- 0.9*0.9 # detected HSIL > normal (% success * % treated)
dcan_norm <- 0.5*1 # detected cancer > normal (% success (based on ~5 year survival) * % treated)

# hysterectomies
hyst_1 <- -log(0.99)/10 # hysterectomy 21-29 (normalized from 10 year rate > 1 year rate)
hyst_2 <- -log(0.96)/10 # hysterectomy 30-39 (normalized from 10 year rate > 1 year rate)


#### Create array ####

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
arr1[1,,] <- rmultinom(1, (Pop_size_1+Pop_size_2+Pop_size_3), prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 

### Run model ###

for(i in 2:t){
  
  ########################### Variable rates ########################
  
  t.index <- i
  
  # screening
  ulsil_dlsil <- 0.83/3 # no change in screening
  uhsil_dhsil <- 0.83/3 # no change in screening
  ucan_dcan <- 0.83/3 # no change in screening
  
  # loss to follow up
  dlsil_uhsil <- 0.17/3 # no change in LTFU
  dhsil_ucan <- 0.17/3 # no change in LTFU
  
  ########################### 18-20 ########################
  
  # Normal 18-20 (1)
  arr1[i,1,1] <- (arr1[(i-1),1,1]) +  
    (Pop_size_1*age_in) +             # age in
    (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),1,1])*norm_ulsil_1 -  # normal progress to undetected LSIL
    (arr1[(i-1),1,1])*age_up_0 -      # age up
    (arr1[(i-1),1,1])*p_die_1         # die (unrelated)
  # LSIL undetected 18-20 (2)
  arr1[i,2,1] <- (arr1[(i-1),2,1]) +  
    (arr1[(i-1),1,1])*norm_ulsil_0 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,1])*age_up_0 -      # age up
    (arr1[(i-1),2,1])*p_die_1         # die (unrelated)
  # LSIL detected 18-20 (3) = 0
  arr1[i,3,1] <- (arr1[(i-1),3,1])    # stay the same (0)
  # HSIL undetected 18-20 (4) 
  arr1[i,4,1] <- (arr1[(i-1),4,1]) + 
    (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,1])*age_up_0 -      # age up
    (arr1[(i-1),4,1])*p_die_1         # die (unrelated)
  # HSIL detected 18-20 (5) = 0
  arr1[i,5,1] <- (arr1[(i-1),5,1])    # stay the same (0)
  # Cancer undetected 18-20 (6) = 0
  arr1[i,6,1] <- (arr1[(i-1),6,1])    # stay the same (0)
  # Cancer detected 18-20 (7) = 0
  arr1[i,7,1] <- (arr1[(i-1),7,1])    # stay the same (0)
  # Cancer deaths 18-20 (8) = 0
  arr1[i,8,1] <- (arr1[(i-1),8,1])    # stay the same (0)
  
  ########################### 21-24 ########################
  
  # Normal 21-24 (1)
  arr1[i,1,2] <- (arr1[(i-1),1,2]) +  
    (arr1[(i-1),1,1])*age_up_0 +      # age up
    (arr1[(i-1),2,2])*ulsil_norm_1 +  # undetected LSIL regress to normal
    (arr1[(i-1),5,2])*dhsil_norm +    # detected HSIL treated to normal
    (arr1[(i-1),3,2])*dlsil_norm_1 -  # detected LSIL regress to normal
    (arr1[(i-1),1,2])*norm_ulsil_1 -  # normal progress to undetected LSIL
    (arr1[(i-1),1,2])*norm_dlsil_1 -  # normal progress to detected LSIL
    (arr1[(i-1),1,2])*age_up_1 -      # age up
    (arr1[(i-1),1,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),1,2])*p_die_1         # die (unrelated)
  # LSIL undetected 21-24 (2)
  arr1[i,2,2] <- (arr1[(i-1),2,2]) + 
    (arr1[(i-1),2,1])*age_up_0 +      # age up
    (arr1[(i-1),1,2])*norm_ulsil_1 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,2])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,2])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,2])*ulsil_dlsil -   # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,2])*age_up_1 -      # age up
    (arr1[(i-1),2,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),2,2])*p_die_1         # die (unrelated)
  # LSIL detected 21-24 (3)
  arr1[i,3,2] <- (arr1[(i-1),3,2]) + 
    (arr1[(i-1),3,1])*age_up_0 +      # age up
    (arr1[(i-1),1,2])*norm_dlsil_1 +  # normal progress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,2])*ulsil_dlsil -   # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,2])*dlsil_norm_1 -  # detected LSIL regress to normal
    (arr1[(i-1),3,2])*dlsil_dhsil_1 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,2])*dlsil_uhsil -   # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,2])*age_up_1 -      # age up
    (arr1[(i-1),3,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),3,2])*p_die_1         # die (unrelated)
  # HSIL undetected 21-24 (4) 
  arr1[i,4,2] <- (arr1[(i-1),4,2]) + 
    (arr1[(i-1),4,1])*age_up_0 +      # age up
    (arr1[(i-1),2,2])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,2])*dlsil_uhsil -   # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,2])*age_up_1 -      # age up
    (arr1[(i-1),4,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),4,2])*p_die_1         # die (unrelated)
  # HSIL detected 21-24 (5)
  arr1[i,5,2] <- (arr1[(i-1),5,2]) + 
    (arr1[(i-1),5,1])*age_up_0 +      # age up
    (arr1[(i-1),3,2])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_norm -    # detected HSIL treated to normal
    (arr1[(i-1),5,2])*age_up_1 -      # age up
    (arr1[(i-1),5,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),5,2])*p_die_1         # die (unrelated)
  # Cancer undetected 18-24 (6) = 0
  arr1[i,6,2] <- (arr1[(i-1),6,2])    # stay the same (0)
  # Cancer detected 18-24 (7) = 0
  arr1[i,7,2] <- (arr1[(i-1),7,2])    # stay the same (0)
  # Cancer deaths 18-24 (8) = 0
  arr1[i,8,2] <- (arr1[(i-1),8,2])    # stay the same (0)
  
  ########################### 25-29 ########################
  
  # Normal 25-29 (1)
  arr1[i,1,3] <- (arr1[(i-1),1,3]) +  
    (arr1[(i-1),1,2])*age_up_1 +        # age up
    (arr1[(i-1),2,3])*ulsil_norm_2_3 +  # undetected LSIL regress to normal
    (arr1[(i-1),5,3])*dhsil_norm +      # detected HSIL treated to normal
    (arr1[(i-1),7,3])*dcan_norm +       # detected cancer treated to normal
    (arr1[(i-1),3,3])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),1,3])*norm_ulsil_2 -    # normal progress to undetected LSIL
    (arr1[(i-1),1,3])*norm_dlsil_2 -    # normal progress to detected LSIL
    (arr1[(i-1),1,3])*age_up_2 -        # age up 
    (arr1[(i-1),1,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),1,3])*p_die_2_3         # die (unrelated)
  # LSIL undetected 25-29 (2)
  arr1[i,2,3] <- (arr1[(i-1),2,3]) + 
    (arr1[(i-1),2,2])*age_up_1 +        # age up
    (arr1[(i-1),1,3])*norm_ulsil_2 +    # normal progress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,3])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,3])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,3])*age_up_2 -        # age up
    (arr1[(i-1),2,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),2,3])*p_die_2_3         # die (unrelated)
  # LSIL detected 25-29 (3)
  arr1[i,3,3] <- (arr1[(i-1),3,3]) + 
    (arr1[(i-1),3,2])*age_up_1 +        # age up
    (arr1[(i-1),1,3])*norm_dlsil_2 +    # normal progress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,3])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,3])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,3])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,3])*age_up_2 -        # age up
    (arr1[(i-1),3,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),3,3])*p_die_2_3         # die (unrelated)
  # HSIL undetected 25-29 (4) 
  arr1[i,4,3] <- (arr1[(i-1),4,3]) + 
    (arr1[(i-1),4,2])*age_up_1 +        # age up
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,3])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ucan -      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,3])*age_up_2 -        # age up
    (arr1[(i-1),4,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),4,3])*p_die_2_3         # die (unrelated)
  # HSIL detected 25-29 (5)
  arr1[i,5,3] <- (arr1[(i-1),5,3]) + 
    (arr1[(i-1),5,2])*age_up_1 +        # age up
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,3])*dhsil_dcan -      # detected HSIL progress to detected cancer
    (arr1[(i-1),5,3])*dhsil_norm -      # detected HSIL treated to normal
    (arr1[(i-1),5,3])*age_up_2 -        # age up
    (arr1[(i-1),5,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),5,3])*p_die_2_3         # die (unrelated)
  # Cancer undetected 25-29 (6)
  arr1[i,6,3] <- (arr1[(i-1),6,3]) +
    (arr1[(i-1),4,3])*uhsil_ucan +      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),6,3])*age_up_2 -        # age up
    (arr1[(i-1),6,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),6,3])*p_die_2_3         # die (unrelated)
  # Cancer detected 25-29 (7)
  arr1[i,7,3] <- (arr1[(i-1),7,3]) +
    (arr1[(i-1),5,3])*dhsil_dcan +      # detected HSIL progress to detected cancer
    (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),7,3])*dcan_dcandeath -  # detected cancer progress to cancer death
    (arr1[(i-1),7,3])*dcan_norm -       # detected cancer treated to normal
    (arr1[(i-1),7,3])*age_up_2 -        # age out
    (arr1[(i-1),7,3])*p_die_2_3         # die (unrelated)
  # Cancer deaths 25-29 (8)
  arr1[i,8,3] <-
    (arr1[(i-1),7,3])*dcan_dcandeath    # detected cancer progress to cancer death
  
  ########################### 30-39 ########################
  
  # Normal 30-39 (1)
  arr1[i,1,4] <- (arr1[(i-1),1,4]) +  
    (arr1[(i-1),1,3])*age_up_2 +        # age up
    (arr1[(i-1),2,4])*ulsil_norm_2_3 +  # undetected LSIL regress to normal
    (arr1[(i-1),5,4])*dhsil_norm +      # detected HSIL treated to normal
    (arr1[(i-1),7,4])*dcan_norm +       # detected cancer treated to normal
    (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),1,4])*norm_ulsil_3 -    # normal progress to undetected LSIL
    (arr1[(i-1),1,4])*norm_dlsil_3 -    # normal progress to detected LSIL
    (arr1[(i-1),1,4])*age_out -         # age out
    (arr1[(i-1),1,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),1,4])*p_die_2_3         # die (unrelated)
  # LSIL undetected 25-39 (2)
  arr1[i,2,4] <- (arr1[(i-1),2,4]) + 
    (arr1[(i-1),2,3])*age_up_2 +        # age up
    (arr1[(i-1),1,4])*norm_ulsil_3 +    # normal progress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,4])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,4])*age_out -         # age out
    (arr1[(i-1),2,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),2,4])*p_die_2_3         # die (unrelated)
  # LSIL detected 25-39 (3)
  arr1[i,3,4] <- (arr1[(i-1),3,4]) + 
    (arr1[(i-1),3,3])*age_up_2 +        # age up
    (arr1[(i-1),1,4])*norm_dlsil_3 +    # normal progress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,4])*age_out -         # age out
    (arr1[(i-1),3,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),3,4])*p_die_2_3         # die (unrelated)
  # HSIL undetected 25-39 (4) 
  arr1[i,4,4] <- (arr1[(i-1),4,4]) + 
    (arr1[(i-1),4,3])*age_up_2 +        # age up
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ucan -      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,4])*age_out -         # age out
    (arr1[(i-1),4,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),4,4])*p_die_2_3         # die (unrelated)
  # HSIL detected 25-39 (5)
  arr1[i,5,4] <- (arr1[(i-1),5,4]) + 
    (arr1[(i-1),5,3])*age_up_2 +        # age up
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,4])*dhsil_dcan -      # detected HSIL progress to detected cancer
    (arr1[(i-1),5,4])*dhsil_norm -      # detected HSIL treated to normal
    (arr1[(i-1),5,4])*age_out -         # age out
    (arr1[(i-1),5,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),5,4])*p_die_2_3         # die (unrelated)
  # Cancer undetected 25-39 (6)
  arr1[i,6,4] <- (arr1[(i-1),6,4]) +
    (arr1[(i-1),6,4])*age_up_2 +        # age up
    (arr1[(i-1),4,4])*uhsil_ucan +      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),6,4])*age_out -         # age out
    (arr1[(i-1),6,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),6,4])*p_die_2_3         # die (unrelated)
  # Cancer detected 25-39 (7)
  arr1[i,7,4] <- (arr1[(i-1),7,4]) +
    (arr1[(i-1),7,3])*age_up_2 +        # age up
    (arr1[(i-1),5,4])*dhsil_dcan +      # detected HSIL progress to detected cancer
    (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),7,4])*dcan_dcandeath -  # detected cancer progress to cancer death
    (arr1[(i-1),7,4])*dcan_norm -       # detected cancer treated to normal
    (arr1[(i-1),7,4])*age_out -         # age out
    (arr1[(i-1),7,4])*p_die_2_3         # die (unrelated)
  # Cancer deaths 25-39 (8)
  arr1[i,8,4] <-
    (arr1[(i-1),7,4])*dcan_dcandeath    # detected cancer progress to cancer death
}

arr1 <- round(arr1,0)

# Collapse age stratification
result_noCov <- apply(arr1, 2L, rowSums)
# pull out last row of resultant df
final_result_noCov <- as.data.frame(result_noCov[nrow(result_noCov),])

# Log likelihood of result vs predicted
LL <- sum(dpois(final_result_noCov[c(5,7,8),], c(196000*0.4,10510*0.4,3400*0.4), log=TRUE))

# Pull out detected HSIL, detected cancer, and cancer death
Results_noCov <- final_result_noCov[c(5,7,8),]
# Pull out detected HSIL, undetected cancer, and cancer death
Results_3_noCov <- final_result_noCov[c(5,6,8),]

# Make resultant array into df
result_2_noCov <- as.data.frame(result_noCov)
# Create column name for Years column
result_3_noCov <- as.data.frame(tibble::rownames_to_column(result_2_noCov, "Year"))
# Add year # column
result_3_noCov$Year.No <- seq(1:t)



#################### UNVACCINATED ####################

# Population size
# 154566548 = all women
Pop_size_0 = 35704873*0.6 # <18*proportion vaxxed
Pop_size_1 = 14683822*0.6 # 18-24*proportion vaxxed
Pop_size_2 = 10201392*0.6 # 25-29*proportion vaxxed
Pop_size_3 = 19939084*0.6 # 30-39*proportion vaxxed

#### Create array ####

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
arr1[1,,] <- rmultinom(1, (Pop_size_1+Pop_size_2+Pop_size_3), prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 


### Run model ###

for(i in 2:t){
  
  ########################### Variable rates ########################
  
  t.index <- i
  
  # vaccination
  norm_ulsil_0 <- ifelse(t.index<988, 0.15, 0.15*0.2) # starting in year 988 (2008), vaccination drops incidence by 80%
  norm_ulsil_1 <- ifelse(t.index<988, 0.08, 0.08*0.2) # starting in year 988 (2008), vaccination drops incidence by 80%
  norm_ulsil_2 <- ifelse(t.index<988, 0.02, 0.02*0.2) # starting in year 988 (2008), vaccination drops incidence by 80%
  
  # screening
  ulsil_dlsil <- 0.83/3 # no change in screening
  uhsil_dhsil <- 0.83/3 # no change in screening
  ucan_dcan <- 0.83/3 # no change in screening
  
  # loss to follow up
  dlsil_uhsil <- 0.17/3 # no change in LTFU
  dhsil_ucan <- 0.17/3 # no change in LTFU

  ########################### 18-20 ########################
  
  # Normal 18-20 (1)
  arr1[i,1,1] <- (arr1[(i-1),1,1]) +  
    (Pop_size_1*age_in) +             # age in
    (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),1,1])*norm_ulsil_1 -  # normal progress to undetected LSIL
    (arr1[(i-1),1,1])*age_up_0 -      # age up
    (arr1[(i-1),1,1])*p_die_1         # die (unrelated)
  # LSIL undetected 18-20 (2)
  arr1[i,2,1] <- (arr1[(i-1),2,1]) +  
    (arr1[(i-1),1,1])*norm_ulsil_0 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,1])*age_up_0 -      # age up
    (arr1[(i-1),2,1])*p_die_1         # die (unrelated)
  # LSIL detected 18-20 (3) = 0
  arr1[i,3,1] <- (arr1[(i-1),3,1])    # stay the same (0)
  # HSIL undetected 18-20 (4) 
  arr1[i,4,1] <- (arr1[(i-1),4,1]) + 
    (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,1])*age_up_0 -      # age up
    (arr1[(i-1),4,1])*p_die_1         # die (unrelated)
  # HSIL detected 18-20 (5) = 0
  arr1[i,5,1] <- (arr1[(i-1),5,1])    # stay the same (0)
  # Cancer undetected 18-20 (6) = 0
  arr1[i,6,1] <- (arr1[(i-1),6,1])    # stay the same (0)
  # Cancer detected 18-20 (7) = 0
  arr1[i,7,1] <- (arr1[(i-1),7,1])    # stay the same (0)
  # Cancer deaths 18-20 (8) = 0
  arr1[i,8,1] <- (arr1[(i-1),8,1])    # stay the same (0)
  
  ########################### 21-24 ########################
  
  # Normal 21-24 (1)
  arr1[i,1,2] <- (arr1[(i-1),1,2]) +  
    (arr1[(i-1),1,1])*age_up_0 +      # age up
    (arr1[(i-1),2,2])*ulsil_norm_1 +  # undetected LSIL regress to normal
    (arr1[(i-1),5,2])*dhsil_norm +    # detected HSIL treated to normal
    (arr1[(i-1),3,2])*dlsil_norm_1 -  # detected LSIL regress to normal
    (arr1[(i-1),1,2])*norm_ulsil_1 -  # normal progress to undetected LSIL
    (arr1[(i-1),1,2])*norm_dlsil_1 -  # normal progress to detected LSIL
    (arr1[(i-1),1,2])*age_up_1 -      # age up
    (arr1[(i-1),1,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),1,2])*p_die_1         # die (unrelated)
  # LSIL undetected 21-24 (2)
  arr1[i,2,2] <- (arr1[(i-1),2,2]) + 
    (arr1[(i-1),2,1])*age_up_0 +      # age up
    (arr1[(i-1),1,2])*norm_ulsil_1 +  # normal progress to undetected LSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,2])*ulsil_norm_1 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,2])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,2])*ulsil_dlsil -   # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,2])*age_up_1 -      # age up
    (arr1[(i-1),2,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),2,2])*p_die_1         # die (unrelated)
  # LSIL detected 21-24 (3)
  arr1[i,3,2] <- (arr1[(i-1),3,2]) + 
    (arr1[(i-1),3,1])*age_up_0 +      # age up
    (arr1[(i-1),1,2])*norm_dlsil_1 +  # normal progress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,2])*ulsil_dlsil -   # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,2])*dlsil_norm_1 -  # detected LSIL regress to normal
    (arr1[(i-1),3,2])*dlsil_dhsil_1 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,2])*dlsil_uhsil -   # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,2])*age_up_1 -      # age up
    (arr1[(i-1),3,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),3,2])*p_die_1         # die (unrelated)
  # HSIL undetected 21-24 (4) 
  arr1[i,4,2] <- (arr1[(i-1),4,2]) + 
    (arr1[(i-1),4,1])*age_up_0 +      # age up
    (arr1[(i-1),2,2])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,2])*dlsil_uhsil -   # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,2])*age_up_1 -      # age up
    (arr1[(i-1),4,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),4,2])*p_die_1         # die (unrelated)
  # HSIL detected 21-24 (5)
  arr1[i,5,2] <- (arr1[(i-1),5,2]) + 
    (arr1[(i-1),5,1])*age_up_0 +      # age up
    (arr1[(i-1),3,2])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,2])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,2])*dhsil_norm -    # detected HSIL treated to normal
    (arr1[(i-1),5,2])*age_up_1 -      # age up
    (arr1[(i-1),5,2])*hyst_1-         # hysterectomy
    (arr1[(i-1),5,2])*p_die_1         # die (unrelated)
  # Cancer undetected 18-24 (6) = 0
  arr1[i,6,2] <- (arr1[(i-1),6,2])    # stay the same (0)
  # Cancer detected 18-24 (7) = 0
  arr1[i,7,2] <- (arr1[(i-1),7,2])    # stay the same (0)
  # Cancer deaths 18-24 (8) = 0
  arr1[i,8,2] <- (arr1[(i-1),8,2])    # stay the same (0)
  
  ########################### 25-29 ########################
  
  # Normal 25-29 (1)
  arr1[i,1,3] <- (arr1[(i-1),1,3]) +  
    (arr1[(i-1),1,2])*age_up_1 +        # age up
    (arr1[(i-1),2,3])*ulsil_norm_2_3 +  # undetected LSIL regress to normal
    (arr1[(i-1),5,3])*dhsil_norm +      # detected HSIL treated to normal
    (arr1[(i-1),7,3])*dcan_norm +       # detected cancer treated to normal
    (arr1[(i-1),3,3])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),1,3])*norm_ulsil_2 -    # normal progress to undetected LSIL
    (arr1[(i-1),1,3])*norm_dlsil_2 -    # normal progress to detected LSIL
    (arr1[(i-1),1,3])*age_up_2 -        # age up 
    (arr1[(i-1),1,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),1,3])*p_die_2_3         # die (unrelated)
  # LSIL undetected 25-29 (2)
  arr1[i,2,3] <- (arr1[(i-1),2,3]) + 
    (arr1[(i-1),2,2])*age_up_1 +        # age up
    (arr1[(i-1),1,3])*norm_ulsil_2 +    # normal progress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,3])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,3])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,3])*age_up_2 -        # age up
    (arr1[(i-1),2,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),2,3])*p_die_2_3         # die (unrelated)
  # LSIL detected 25-29 (3)
  arr1[i,3,3] <- (arr1[(i-1),3,3]) + 
    (arr1[(i-1),3,2])*age_up_1 +        # age up
    (arr1[(i-1),1,3])*norm_dlsil_2 +    # normal progress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,3])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,3])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,3])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,3])*age_up_2 -        # age up
    (arr1[(i-1),3,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),3,3])*p_die_2_3         # die (unrelated)
  # HSIL undetected 25-29 (4) 
  arr1[i,4,3] <- (arr1[(i-1),4,3]) + 
    (arr1[(i-1),4,2])*age_up_1 +        # age up
    (arr1[(i-1),2,3])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,3])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,3])*uhsil_ucan -      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,3])*age_up_2 -        # age up
    (arr1[(i-1),4,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),4,3])*p_die_2_3         # die (unrelated)
  # HSIL detected 25-29 (5)
  arr1[i,5,3] <- (arr1[(i-1),5,3]) + 
    (arr1[(i-1),5,2])*age_up_1 +        # age up
    (arr1[(i-1),3,3])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,3])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,3])*dhsil_dcan -      # detected HSIL progress to detected cancer
    (arr1[(i-1),5,3])*dhsil_norm -      # detected HSIL treated to normal
    (arr1[(i-1),5,3])*age_up_2 -        # age up
    (arr1[(i-1),5,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),5,3])*p_die_2_3         # die (unrelated)
  # Cancer undetected 25-29 (6)
  arr1[i,6,3] <- (arr1[(i-1),6,3]) +
    (arr1[(i-1),4,3])*uhsil_ucan +      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),6,3])*age_up_2 -        # age up
    (arr1[(i-1),6,3])*hyst_1-           # hysterectomy
    (arr1[(i-1),6,3])*p_die_2_3         # die (unrelated)
  # Cancer detected 25-29 (7)
  arr1[i,7,3] <- (arr1[(i-1),7,3]) +
    (arr1[(i-1),5,3])*dhsil_dcan +      # detected HSIL progress to detected cancer
    (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),7,3])*dcan_dcandeath -  # detected cancer progress to cancer death
    (arr1[(i-1),7,3])*dcan_norm -       # detected cancer treated to normal
    (arr1[(i-1),7,3])*age_up_2 -        # age out
    (arr1[(i-1),7,3])*p_die_2_3         # die (unrelated)
  # Cancer deaths 25-29 (8)
  arr1[i,8,3] <-
    (arr1[(i-1),7,3])*dcan_dcandeath    # detected cancer progress to cancer death
  
  ########################### 30-39 ########################
  
  # Normal 30-39 (1)
  arr1[i,1,4] <- (arr1[(i-1),1,4]) +  
    (arr1[(i-1),1,3])*age_up_2 +        # age up
    (arr1[(i-1),2,4])*ulsil_norm_2_3 +  # undetected LSIL regress to normal
    (arr1[(i-1),5,4])*dhsil_norm +      # detected HSIL treated to normal
    (arr1[(i-1),7,4])*dcan_norm +       # detected cancer treated to normal
    (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),1,4])*norm_ulsil_3 -    # normal progress to undetected LSIL
    (arr1[(i-1),1,4])*norm_dlsil_3 -    # normal progress to detected LSIL
    (arr1[(i-1),1,4])*age_out -         # age out
    (arr1[(i-1),1,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),1,4])*p_die_2_3         # die (unrelated)
  # LSIL undetected 25-39 (2)
  arr1[i,2,4] <- (arr1[(i-1),2,4]) + 
    (arr1[(i-1),2,3])*age_up_2 +        # age up
    (arr1[(i-1),1,4])*norm_ulsil_3 +    # normal progress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),2,4])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
    (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),2,4])*age_out -         # age out
    (arr1[(i-1),2,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),2,4])*p_die_2_3         # die (unrelated)
  # LSIL detected 25-39 (3)
  arr1[i,3,4] <- (arr1[(i-1),3,4]) + 
    (arr1[(i-1),3,3])*age_up_2 +        # age up
    (arr1[(i-1),1,4])*norm_dlsil_3 +    # normal progress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
    (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
    (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
    (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),3,4])*age_out -         # age out
    (arr1[(i-1),3,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),3,4])*p_die_2_3         # die (unrelated)
  # HSIL undetected 25-39 (4) 
  arr1[i,4,4] <- (arr1[(i-1),4,4]) + 
    (arr1[(i-1),4,3])*age_up_2 +        # age up
    (arr1[(i-1),2,4])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
    (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
    (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
    (arr1[(i-1),4,4])*uhsil_ucan -      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),4,4])*age_out -         # age out
    (arr1[(i-1),4,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),4,4])*p_die_2_3         # die (unrelated)
  # HSIL detected 25-39 (5)
  arr1[i,5,4] <- (arr1[(i-1),5,4]) + 
    (arr1[(i-1),5,3])*age_up_2 +        # age up
    (arr1[(i-1),3,4])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
    (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
    (arr1[(i-1),5,4])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
    (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),5,4])*dhsil_dcan -      # detected HSIL progress to detected cancer
    (arr1[(i-1),5,4])*dhsil_norm -      # detected HSIL treated to normal
    (arr1[(i-1),5,4])*age_out -         # age out
    (arr1[(i-1),5,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),5,4])*p_die_2_3         # die (unrelated)
  # Cancer undetected 25-39 (6)
  arr1[i,6,4] <- (arr1[(i-1),6,4]) +
    (arr1[(i-1),6,4])*age_up_2 +        # age up
    (arr1[(i-1),4,4])*uhsil_ucan +      # undetected HSIL progress to undetected cancer
    (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
    (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),6,4])*age_out -         # age out
    (arr1[(i-1),6,4])*hyst_2-           # hysterectomy
    (arr1[(i-1),6,4])*p_die_2_3         # die (unrelated)
  # Cancer detected 25-39 (7)
  arr1[i,7,4] <- (arr1[(i-1),7,4]) +
    (arr1[(i-1),7,3])*age_up_2 +        # age up
    (arr1[(i-1),5,4])*dhsil_dcan +      # detected HSIL progress to detected cancer
    (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
    (arr1[(i-1),7,4])*dcan_dcandeath -  # detected cancer progress to cancer death
    (arr1[(i-1),7,4])*dcan_norm -       # detected cancer treated to normal
    (arr1[(i-1),7,4])*age_out -         # age out
    (arr1[(i-1),7,4])*p_die_2_3         # die (unrelated)
  # Cancer deaths 25-39 (8)
  arr1[i,8,4] <-
    (arr1[(i-1),7,4])*dcan_dcandeath    # detected cancer progress to cancer death
}

arr1 <- round(arr1,0)

# Collapse array - remove age stratification
result_5_noCov <- apply(arr1, 2L, rowSums)
# Vector of HSIL, cancer, and cancer death at steady state for iteration
final_result_2_noCov <- as.data.frame(result_5_noCov[nrow(result_5_noCov),])

# Log likelihood of result vs predicted
LL_2 <- sum(dpois(final_result_2_noCov[c(5,7,8),], c(196000*0.6,10510*0.6,3400*0.6), log=TRUE))

# Vector of detected HSIL, detected cancer, and cancer death
Results_2_noCov <- final_result_2_noCov[c(5,7,8),]
# Vector of detected HSIL, detected cancer, and cancer death
Results_4_noCov <- final_result_2_noCov[c(4,6,8),]

# Create df from resultant array
result_6_noCov <- as.data.frame(result_5_noCov)
# Add column name to Years column
result_7_noCov <- as.data.frame(tibble::rownames_to_column(result_6_noCov, "Year"))
# Add year # column
result_7_noCov$Year.No <- seq(1:t)


################### COMBINED RESULTS ###################

# Combine unvaccinated and vaccinated into one vector (final results only)
Results_all_noCov = Results_noCov+Results_2_noCov # final results

# Combine vaccinated and unvaccinated into one df (entire df)
result_tot_noCov <- result_3_noCov[,-1] + result_7_noCov[,-1]
#View(result_tot)

# After adding, divide by two to fix Year #
result_tot_noCov$Year.No <- result_tot_noCov$Year.No/2
# Subtract to match year # to calendar year
result_tot_noCov$Year <- result_tot_noCov$Year.No+(2008-987)

# Combined LSIL undetected and detected
result_tot_noCov$All_LSIL <- result_tot_noCov$Undet_LSIL + result_tot_noCov$Det_LSIL
# Combine HSIL undetected and detected
result_tot_noCov$All_HSIL <- result_tot_noCov$Undet_HSIL + result_tot_noCov$Det_HSIL
# Combined Cancer undetected and detected
result_tot_noCov$All_Cancer <- result_tot_noCov$Undet_Cancer + result_tot_noCov$Det_Cancer



# Expected: 216000,11778,3663
result_tot_noCov[which(result_tot$Year==2008),]
# Expected: 196000,10510,3400
result_tot_noCov[which(result_tot$Year==2016),]
# 2050
result_tot[which(result_tot$Year==2008),]
result_tot[which(result_tot$Year==2016),]

## Integration

HSIL_Cov <- integrate(splinefun(result_tot$Year, result_tot$All_HSIL), 2020,2050)
HSIL_noCov <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_HSIL), 2020,2050)
HSIL_Cov
HSIL_noCov

Undet_HSIL_Cov <- integrate(splinefun(result_tot$Year, result_tot$Undet_HSIL), 2020,2050)
Undet_HSIL_noCov <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Undet_HSIL), 2020,2050)
Undet_HSIL_Cov
Undet_HSIL_noCov

cancer_Cov <- integrate(splinefun(result_tot$Year, result_tot$Det_Cancer), 2020,2050)
cancer_noCov <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_Cancer), 2020,2050)

cancerdeath_Cov <- integrate(splinefun(result_tot$Year, result_tot$Cancer_Death), 2020,2050)
cancerdeath_noCov <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Cancer_Death), 2020,2050)
cancerdeath_Cov
cancerdeath_noCov



