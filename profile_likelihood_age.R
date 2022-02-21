### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Feb 4, 2021

### COMBINED MODEL ###

### Profile likelihood estimation with age stratification

# Profile likelihood estimation function

profile <- function(
  norm_ulsil_1,norm_ulsil_2,norm_ulsil_3,
  ulsil_norm_1,ulsil_norm_2_3,
  ulsil_uhsil_1,ulsil_uhsil_2_3,
  uhsil_ulsil_1,uhsil_ulsil_2_3,
  dcan_dcandeath){
  
  #### Starting parameters
  t = 1000
  N.states = 8
  #Pop_size = 154566548 #all women
  #Pop_size = 44978865*(1-0.032) # women 18-39*adjustment for hysterectomies
  Pop_size_1 = 14522300
  Pop_size_2 = 10089177
  Pop_size_3 = 18882313
  
  # aging into cohort
  age_in <- 1/20 # multiply by pop under 18
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
  #norm_ulsil_1 <- 0.05
  #norm_ulsil_2 <- 0.03
  #norm_ulsil_3 <- 0.01
  #ulsil_norm_1 <- 0.60
  #ulsil_norm_2_3 <- 0.4
  #ulsil_uhsil_1 <- 0.14
  #ulsil_uhsil_2_3 <- 0.20
  #uhsil_ulsil_1 <- 0.62
  #uhsil_ulsil_2_3 <- 0.30
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
  #dcan_dcandeath <- 0.1
  
  # treatment
  dhsil_norm <- 0.9*0.9
  dcan_norm <- 0.5*1
  
  # BRFSS rates!!
  # screening
  ulsil_dlsil <- 0.41/3
  uhsil_dhsil <- 0.83/3
  ucan_dcan <- 0.83/3
  
  # loss to follow up
  dlsil_uhsil <- 0.17/3
  dhsil_ucan <- 0.17/3
  
  
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
  result <- apply(arr1, 2L, rowSums)
  final_result <- as.data.frame(result[1000,])
  
  # Log likelihood of result vs predicted
  LL <- sum(dpois(final_result[c(5,7,8),], c(216000,11778,3663), log=TRUE))
  
  output <- list("Results"=final_result[c(5,7,8),], "LL"=LL)
  
  return(output)
  
}



################### Run function with many inputs #######################

# Inputs
norm_ulsil_1_seq <- 0.08
#norm_ulsil_1_seq <- seq(0.03,0.10,0.01)
norm_ulsil_2_seq <- 0.02
#norm_ulsil_2_seq <- seq(0.02,0.05,0.01)
norm_ulsil_3_seq <- 0.01
#norm_ulsil_3_seq <- seq(0.01,0.03,0.01)
ulsil_norm_1_seq <- 0.60
ulsil_norm_2_3_seq <- 0.4
ulsil_uhsil_1_seq <- 0.14
ulsil_uhsil_2_3_seq <- 0.35
#ulsil_uhsil_2_3_seq <- seq(0.10,0.35,0.05)
uhsil_ulsil_1_seq <- 0.62
uhsil_ulsil_2_3_seq <- 0.20
#uhsil_ulsil_2_3_seq <- seq(0.20,0.50,0.05)
#dcan_dcandeath <- 0.1
#dcan_dcandeath_seq <- seq(0.05,0.40,0.05)
dcan_dcandeath_seq <- 0.35

# Empty vector for function output
LL <- array(NA, dim=c(length(norm_ulsil_1_seq),length(norm_ulsil_2_seq),length(norm_ulsil_3_seq),length(ulsil_norm_1_seq),length(ulsil_norm_2_3_seq),length(ulsil_uhsil_1_seq),length(ulsil_uhsil_2_3_seq),length(uhsil_ulsil_1_seq),length(uhsil_ulsil_2_3_seq),length(dcan_dcandeath_seq)))

Results <- array(NA, dim=c(3,length(norm_ulsil_1_seq),length(norm_ulsil_2_seq),length(norm_ulsil_3_seq),length(ulsil_norm_1_seq),length(ulsil_norm_2_3_seq),length(ulsil_uhsil_1_seq),length(ulsil_uhsil_2_3_seq),length(uhsil_ulsil_1_seq),length(uhsil_ulsil_2_3_seq),length(dcan_dcandeath_seq)))

# Time start
ptm <- proc.time()
# Run for loop

##### START HERE #####
for(i in 1:length(norm_ulsil_1_seq)){
  for(j in 1:length(norm_ulsil_2_seq)){
    for(k in 1:length(norm_ulsil_3_seq)){
      for(l in 1:length(ulsil_norm_1_seq)){
        for(m in 1:length(ulsil_norm_2_3_seq)){
          for(n in 1:length(ulsil_uhsil_1_seq)){
            for(o in 1:length(ulsil_uhsil_2_3_seq)){
              for(p in 1:length(uhsil_ulsil_1_seq)){
                for(q in 1:length(uhsil_ulsil_2_3_seq)){
                  for(r in 1:length(dcan_dcandeath_seq)){
                  
                  LL[i,j,k,l,m,n,o,p,q,r] <- profile(norm_ulsil_1=norm_ulsil_1_seq[i],norm_ulsil_2=norm_ulsil_2_seq[j],norm_ulsil_3=norm_ulsil_3_seq[k],ulsil_norm_1=ulsil_norm_1_seq[l],ulsil_norm_2_3=ulsil_norm_2_3_seq[m],ulsil_uhsil_1=ulsil_uhsil_1_seq[n],ulsil_uhsil_2_3=ulsil_uhsil_2_3_seq[o],uhsil_ulsil_1=uhsil_ulsil_1_seq[p],uhsil_ulsil_2_3=uhsil_ulsil_2_3_seq[q],dcan_dcandeath=dcan_dcandeath_seq[r])$LL
                  Results[,i,j,k,l,m,n,o,p,q,r] <- profile(norm_ulsil_1=norm_ulsil_1_seq[i],norm_ulsil_2=norm_ulsil_2_seq[j],norm_ulsil_3=norm_ulsil_3_seq[k],ulsil_norm_1=ulsil_norm_1_seq[l],ulsil_norm_2_3=ulsil_norm_2_3_seq[m],ulsil_uhsil_1=ulsil_uhsil_1_seq[n],ulsil_uhsil_2_3=ulsil_uhsil_2_3_seq[o],uhsil_ulsil_1=uhsil_ulsil_1_seq[p],uhsil_ulsil_2_3=uhsil_ulsil_2_3_seq[q],dcan_dcandeath=dcan_dcandeath_seq[r])$Results
                  
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
#Time end
proc.time() - ptm

# Determine which parameters maximize likelihood
which(LL==max(LL), arr.ind=TRUE)

 # Plot results
plot(norm_ulsil_1_seq, LL[,1,1,1,1,1,1,1,1])
plot(ulsil_uhsil_2_3_seq, LL[1,1,1,1,1,1,,1,1])
plot(dcan_dcandeath_seq, LL[1,1,1,1,1,1,1,1,1,])


