
############################################################################################################

######### Analysis of impact of COVID-19 on LSIL, HSIL, cervical cancer, and cervical cancer death 
######### Guinevere Oliver
######### April 11, 2022

############################################################################################################


# Function to assess impact based on decrease in screening and increase in LTFU

screen_decr <- function(decrease){

  #################### UNVACCINATED ####################
  
  t = 1030 # 1030 years
  N.states = 8 # 8 states
  
  # Population size
  # 154566548 = all women
  Pop_size_0 = 35704873*0.5 # <18*proportion unvaxxed
  Pop_size_1 = 14683822*0.5 # 18-24*proportion unvaxxed
  Pop_size_2 = 10201392*0.5 # 25-29*proportion unvaxxed
  Pop_size_3 = 19939084*0.5 # 30-39*proportion unvaxxed
  Pop_size_4 = (154566548-Pop_size_0-Pop_size_1-Pop_size_2-Pop_size_3)*0.5 # (40+)*proportion unvaxxed
  
  # aging into cohort
  age_in <- (1/20)
  # aging up within cohort
  age_up_0 <- 1/3 # proportion turning 21
  age_up_1 <- 1/4 # proportion turning 25
  age_up_2 <- 1/5 # proportion turning 30
  age_up_3 <- 1/10 # proportion turning 40
  # dying
  die_1 <- 74/100000 # proportion dying age 18-24
  die_2_3 <- 164/100000 # proportion dying age 25-39
  die_4 <- 1500/100000 # proportion dying age 40+
  
  # undetected
  norm_ulsil_0 <- 0.15 # normal > undetected LSIL (18-20)
  norm_ulsil_1 <- 0.08 # normal > undetected LSIL (21-24)
  norm_ulsil_2 <- 0.02 # normal > undetected LSIL (25-29)
  norm_ulsil_3 <- 0.01 # normal > undetected LSIL (30-39)
  norm_ulsil_4 <- 0 # normal > undetected LSIL (40+)
  ulsil_norm_1 <- 0.60 # undetected LSIL > normal (18-24)
  ulsil_norm_2_3 <- 0.4 # undetected LSIL > normal (25-39)
  ulsil_uhsil_1 <- 0.14 # undetected LSIL > undetected HSIL (18-24)
  ulsil_uhsil_2_3 <- 0.30 # undetected LSIL > undetected HSIL (25-39)
  ulsil_uhsil_4 <- 0.30 # undetected LSIL > undetected HSIL (40+)
  uhsil_ulsil_1 <- 0.62 # undetected HSIL > undetected LSIL (18-24)
  uhsil_ulsil_2_3 <- 0.35 # undetected HSIL > undetected LSIL (25-39)
  uhsil_ulsil_4 <- 0.30 # undetected HSIL > undetected LSIL (40+)
  uhsil_ucan_2_3 <- 6/100000 # undetected HSIL > undetected cancer (21-39)
  uhsil_ucan_4 <- 12/100000 # undetected HSIL > undetected cancer (40+)
  
  # detected
  norm_dlsil_1 <- norm_ulsil_1 # normal > detected LSIL (21-24)
  norm_dlsil_2 <- norm_ulsil_2 # normal > detected LSIL (24-29)
  norm_dlsil_3 <- norm_ulsil_3 # normal > detected LSIL (30-39)
  norm_dlsil_4 <- norm_ulsil_4 # normal > detected LSIL (40+)
  dlsil_norm_1 <- ulsil_norm_1 # detected LSIL > normal (21-24)
  dlsil_norm_2_3 <- ulsil_norm_2_3 # detected LSIL > normal (25-39)
  dlsil_dhsil_1 <- ulsil_uhsil_1 # detected LSIL > detected HSIL (21-24)
  dlsil_dhsil_2_3 <- ulsil_uhsil_2_3 # detected LSIL > detected HSIL (25-39)
  dlsil_dhsil_4 <- ulsil_uhsil_4 # detected LSIL > detected HSIL (40+)
  dhsil_dlsil_1 <-  uhsil_ulsil_1 # detected HSIL > detected LSIL (21-24)
  dhsil_dlsil_2_3 <-  uhsil_ulsil_2_3 # detected HSIL > detected LSIL (25-39)
  dhsil_dlsil_4 <-  uhsil_ulsil_4 # detected HSIL > detected LSIL (40+)
  dhsil_dcan_2_3 <- uhsil_ucan_2_3 # detected HSIL > detected cancer
  dhsil_dcan_4 <- uhsil_ucan_4 # detected HSIL > detected cancer
  dcan_dcandeath <- 0.35 # detected cancer > detected cancer death
  
  # treatment
  dhsil_norm <- 0.9*0.9 # detected HSIL > normal (% success * % treated)
  dcan_norm <- 0.5*1 # detected cancer > normal (% success (based on ~5 year survival) * % treated)
  
  # hysterectomies
  hyst_1 <- -log(0.99)/10 # hysterectomy 21-29 (normalized from 10 year rate > 1 year rate)
  hyst_2 <- -log(0.96)/10 # hysterectomy 30-39 (normalized from 10 year rate > 1 year rate)
  hyst_3 <- -log(0.70)/10
  
  
  #### Create array ####
  
  # Create years label
  prefix = "Year"
  suffix = seq(1:t)
  years = paste(prefix, suffix, sep=" ")
  
  # Create empty array
  arr1 = array(NA, dim=c(t, N.states, 5), dimnames=list(years, c("Normal","Undet_LSIL","Det_LSIL","Undet_HSIL","Det_HSIL","Undet_Cancer","Det_Cancer","Cancer_Death"), c("18-20","21-24","25-29","30-39","40+")))
  
  # Assign starting prevalence of each state
  prev1 = 0.8 # normal
  prev2 = 0.1 # LSIL undetected
  prev3 = 0 # LSIL detected
  prev4 = 0.1 # HSIL undetected
  prev5 = 0 # HSIL detected
  prev6 = 0 # cancer undetected
  prev7 = 0 # cancer detected
  prev8 = 0 # cancer death
  
  # Assign starting states
  set.seed(123)
  arr1[1,,] <- rmultinom(1, (Pop_size_1+Pop_size_2+Pop_size_3+Pop_size_4), prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 
  
  ### Run model ###
  
  for(i in 2:t){
    
    ########################### Variable rates ########################
    
    t.index <- i
    
    # screening
    ulsil_dlsil <- ifelse(t.index<1000|t.index>1000,0.83/3,(0.83/3)*(1-decrease)) # in the year 1000 (2020), screening drops by "decrease"%
    uhsil_dhsil <- ifelse(t.index<1000|t.index>1000,0.83/3,(0.83/3)*(1-decrease)) # in the year 1000 (2020), screening drops by "decrease"%
    ucan_dcan <- 0.83/3 # cancer screening does not change (symptomatic)
  
    # loss to follow up
    dlsil_uhsil <- ifelse(t.index<1000|t.index>1000,0.17/3,(0.17/3)*(1+decrease)) # in the year 1000 (2020), LTFU increases by "1+decrease"%
    dhsil_ucan <- ifelse(t.index<1000|t.index>1000,0.17/3,(0.17/3)*(1+decrease)) # in the year 1000 (2020), LTFU increases by "1+decrease"%

    ########################### 18-20 ########################
    
    # Normal 18-20 (1)
    arr1[i,1,1] <- (arr1[(i-1),1,1]) +  
      (Pop_size_1*age_in) +             # age in
      (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
      (arr1[(i-1),1,1])*norm_ulsil_1 -  # normal progress to undetected LSIL
      (arr1[(i-1),1,1])*age_up_0 -      # age up
      (arr1[(i-1),1,1])*die_1           # die (unrelated)
    # LSIL undetected 18-20 (2)
    arr1[i,2,1] <- (arr1[(i-1),2,1]) +  
      (arr1[(i-1),1,1])*norm_ulsil_0 +  # normal progress to undetected LSIL
      (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
      (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
      (arr1[(i-1),2,1])*age_up_0 -      # age up
      (arr1[(i-1),2,1])*die_1           # die (unrelated)
    # LSIL detected 18-20 (3) = 0
    arr1[i,3,1] <- (arr1[(i-1),3,1])    # stay the same (0)
    # HSIL undetected 18-20 (4) 
    arr1[i,4,1] <- (arr1[(i-1),4,1]) + 
      (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,1])*age_up_0 -      # age up
      (arr1[(i-1),4,1])*die_1           # die (unrelated)
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
      (arr1[(i-1),1,2])*die_1           # die (unrelated)
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
      (arr1[(i-1),2,2])*die_1           # die (unrelated)
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
      (arr1[(i-1),3,2])*die_1           # die (unrelated)
    # HSIL undetected 21-24 (4) 
    arr1[i,4,2] <- (arr1[(i-1),4,2]) + 
      (arr1[(i-1),4,1])*age_up_0 +      # age up
      (arr1[(i-1),2,2])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,2])*dlsil_uhsil -   # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,2])*age_up_1 -      # age up
      (arr1[(i-1),4,2])*hyst_1-         # hysterectomy
      (arr1[(i-1),4,2])*die_1           # die (unrelated)
    # HSIL detected 21-24 (5)
    arr1[i,5,2] <- (arr1[(i-1),5,2]) + 
      (arr1[(i-1),5,1])*age_up_0 +      # age up
      (arr1[(i-1),3,2])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,2])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,2])*dhsil_norm -    # detected HSIL treated to normal
      (arr1[(i-1),5,2])*age_up_1 -      # age up
      (arr1[(i-1),5,2])*hyst_1-         # hysterectomy
      (arr1[(i-1),5,2])*die_1           # die (unrelated)
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
      (arr1[(i-1),1,3])*die_2_3           # die (unrelated)
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
      (arr1[(i-1),2,3])*die_2_3           # die (unrelated)
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
      (arr1[(i-1),3,3])*die_2_3           # die (unrelated)
    # HSIL undetected 25-29 (4) 
    arr1[i,4,3] <- (arr1[(i-1),4,3]) + 
      (arr1[(i-1),4,2])*age_up_1 +        # age up
      (arr1[(i-1),2,3])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,3])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,3])*uhsil_ucan_2_3 -  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,3])*age_up_2 -        # age up
      (arr1[(i-1),4,3])*hyst_1-           # hysterectomy
      (arr1[(i-1),4,3])*die_2_3           # die (unrelated)
    # HSIL detected 25-29 (5)
    arr1[i,5,3] <- (arr1[(i-1),5,3]) + 
      (arr1[(i-1),5,2])*age_up_1 +        # age up
      (arr1[(i-1),3,3])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,3])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),5,3])*dhsil_dcan_2_3 -  # detected HSIL progress to detected cancer
      (arr1[(i-1),5,3])*dhsil_norm -      # detected HSIL treated to normal
      (arr1[(i-1),5,3])*age_up_2 -        # age up
      (arr1[(i-1),5,3])*hyst_1-           # hysterectomy
      (arr1[(i-1),5,3])*die_2_3           # die (unrelated)
    # Cancer undetected 25-29 (6)
    arr1[i,6,3] <- (arr1[(i-1),6,3]) +
      (arr1[(i-1),4,3])*uhsil_ucan_2_3 +  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),6,3])*age_up_2 -        # age up
      (arr1[(i-1),6,3])*hyst_1-           # hysterectomy
      (arr1[(i-1),6,3])*die_2_3           # die (unrelated)
    # Cancer detected 25-29 (7)
    arr1[i,7,3] <- (arr1[(i-1),7,3]) +
      (arr1[(i-1),5,3])*dhsil_dcan_2_3 +  # detected HSIL progress to detected cancer
      (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),7,3])*dcan_dcandeath -  # detected cancer progress to cancer death
      (arr1[(i-1),7,3])*dcan_norm -       # detected cancer treated to normal
      (arr1[(i-1),7,3])*age_up_2 -        # age up
      (arr1[(i-1),7,3])*die_2_3           # die (unrelated)
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
      (arr1[(i-1),1,4])*age_up_3 -        # age up
      (arr1[(i-1),1,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),1,4])*die_2_3           # die (unrelated)
    # LSIL undetected 30-39 (2)
    arr1[i,2,4] <- (arr1[(i-1),2,4]) + 
      (arr1[(i-1),2,3])*age_up_2 +        # age up
      (arr1[(i-1),1,4])*norm_ulsil_3 +    # normal progress to undetected LSIL
      (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),2,4])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
      (arr1[(i-1),2,4])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
      (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),2,4])*age_up_3 -        # age up
      (arr1[(i-1),2,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),2,4])*die_2_3           # die (unrelated)
    # LSIL detected 30-39 (3)
    arr1[i,3,4] <- (arr1[(i-1),3,4]) + 
      (arr1[(i-1),3,3])*age_up_2 +        # age up
      (arr1[(i-1),1,4])*norm_dlsil_3 +    # normal progress to detected LSIL
      (arr1[(i-1),5,4])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
      (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
      (arr1[(i-1),3,4])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
      (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),3,4])*age_up_3 -        # age up
      (arr1[(i-1),3,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),3,4])*die_2_3           # die (unrelated)
    # HSIL undetected 30-39 (4) 
    arr1[i,4,4] <- (arr1[(i-1),4,4]) + 
      (arr1[(i-1),4,3])*age_up_2 +        # age up
      (arr1[(i-1),2,4])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,4])*uhsil_ucan_2_3 -  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,4])*age_up_3 -        # age up
      (arr1[(i-1),4,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),4,4])*die_2_3           # die (unrelated)
    # HSIL detected 30-39 (5)
    arr1[i,5,4] <- (arr1[(i-1),5,4]) + 
      (arr1[(i-1),5,3])*age_up_2 +        # age up
      (arr1[(i-1),3,4])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,4])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),5,4])*dhsil_dcan_2_3 -  # detected HSIL progress to detected cancer
      (arr1[(i-1),5,4])*dhsil_norm -      # detected HSIL treated to normal
      (arr1[(i-1),5,4])*age_up_3 -        # age up
      (arr1[(i-1),5,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),5,4])*die_2_3           # die (unrelated)
    # Cancer undetected 25-39 (6)
    arr1[i,6,4] <- (arr1[(i-1),6,4]) +
      (arr1[(i-1),6,3])*age_up_2 +        # age up
      (arr1[(i-1),4,4])*uhsil_ucan_2_3 +  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),6,4])*age_up_3 -        # age up
      (arr1[(i-1),6,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),6,4])*die_2_3           # die (unrelated)
    # Cancer detected 30-39 (7)
    arr1[i,7,4] <- (arr1[(i-1),7,4]) +
      (arr1[(i-1),7,3])*age_up_2 +        # age up
      (arr1[(i-1),5,4])*dhsil_dcan_2_3 +  # detected HSIL progress to detected cancer
      (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),7,4])*dcan_dcandeath -  # detected cancer progress to cancer death
      (arr1[(i-1),7,4])*dcan_norm -       # detected cancer treated to normal
      (arr1[(i-1),7,4])*age_up_3 -        # age up
      (arr1[(i-1),7,4])*die_2_3           # die (unrelated)
    # Cancer deaths 25-39 (8)
    arr1[i,8,4] <-
      (arr1[(i-1),7,4])*dcan_dcandeath    # detected cancer progress to cancer death
    
    ########################### 40+ ########################
    
    # Normal 40+ (1)
    arr1[i,1,5] <- (arr1[(i-1),1,5]) +  
      (arr1[(i-1),1,4])*age_up_3 +        # age up
      (arr1[(i-1),2,5])*ulsil_norm_2_3 +  # undetected LSIL regress to normal
      (arr1[(i-1),5,5])*dhsil_norm +      # detected HSIL treated to normal
      (arr1[(i-1),7,5])*dcan_norm +       # detected cancer treated to normal
      (arr1[(i-1),3,5])*dlsil_norm_2_3 -  # detected LSIL regress to normal
      (arr1[(i-1),1,5])*norm_ulsil_4 -    # normal progress to undetected LSIL
      (arr1[(i-1),1,5])*norm_dlsil_4 -    # normal progress to detected LSIL
      (arr1[(i-1),1,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),1,5])*die_4             # die (unrelated)
    # LSIL undetected 40+ (2)
    arr1[i,2,5] <- (arr1[(i-1),2,5]) + 
      (arr1[(i-1),2,4])*age_up_3 +        # age up
      (arr1[(i-1),1,5])*norm_ulsil_4 +    # normal progress to undetected LSIL
      (arr1[(i-1),4,5])*uhsil_ulsil_4 -   # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),2,5])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
      (arr1[(i-1),2,5])*ulsil_uhsil_4 -   # undetected LSIL progress to undetected HSIL
      (arr1[(i-1),2,5])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),2,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),2,5])*die_4             # die (unrelated)
    # LSIL detected 40+ (3)
    arr1[i,3,5] <- (arr1[(i-1),3,5]) + 
      (arr1[(i-1),3,4])*age_up_3 +        # age up
      (arr1[(i-1),1,5])*norm_dlsil_4 +    # normal progress to detected LSIL
      (arr1[(i-1),5,5])*dhsil_dlsil_4 +   # detected HSIL regress to detected LSIL
      (arr1[(i-1),2,5])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),3,5])*dlsil_norm_2_3 -  # detected LSIL regress to normal
      (arr1[(i-1),3,5])*dlsil_dhsil_4 -   # detected LSIL progress to detected HSIL
      (arr1[(i-1),3,5])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),3,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),3,5])*die_4             # die (unrelated)
    # HSIL undetected 40+ (4) 
    arr1[i,4,5] <- (arr1[(i-1),4,5]) + 
      (arr1[(i-1),4,4])*age_up_3 +        # age up
      (arr1[(i-1),2,5])*ulsil_uhsil_4 +   # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,5])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,5])*uhsil_ulsil_4 -   # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,5])*uhsil_ucan_2_3 -  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),4,5])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),4,5])*die_4             # die (unrelated)
    # HSIL detected 40+ (5)
    arr1[i,5,5] <- (arr1[(i-1),5,5]) + 
      (arr1[(i-1),5,4])*age_up_3 +        # age up
      (arr1[(i-1),3,5])*dlsil_dhsil_4 +   # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,5])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,5])*dhsil_dlsil_4 -   # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,5])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),5,5])*dhsil_dcan_4 -    # detected HSIL progress to detected cancer
      (arr1[(i-1),5,5])*dhsil_norm -      # detected HSIL treated to normal
      (arr1[(i-1),5,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),5,5])*die_4             # die (unrelated)
    # Cancer undetected 40+ (6)
    arr1[i,6,5] <- (arr1[(i-1),6,5]) +
      (arr1[(i-1),6,4])*age_up_3 +        # age up
      (arr1[(i-1),4,5])*uhsil_ucan_4 +    # undetected HSIL progress to undetected cancer
      (arr1[(i-1),5,5])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),6,5])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),6,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),6,5])*die_4             # die (unrelated)
    # Cancer detected 40+ (7)
    arr1[i,7,5] <- (arr1[(i-1),7,5]) +
      (arr1[(i-1),7,4])*age_up_3 +        # age up
      (arr1[(i-1),5,5])*dhsil_dcan_4 +    # detected HSIL progress to detected cancer
      (arr1[(i-1),6,5])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),7,5])*dcan_dcandeath -  # detected cancer progress to cancer death
      (arr1[(i-1),7,5])*dcan_norm -       # detected cancer treated to normal
      (arr1[(i-1),7,5])*die_4             # die (unrelated)
    # Cancer deaths 25-39 (8)
    arr1[i,8,5] <-
      (arr1[(i-1),7,5])*dcan_dcandeath    # detected cancer progress to cancer death
    
  }
  
  # Round numbers to nearest whole number
  arr1 <- round(arr1,0)
  
  # Collapse age stratification
  result <- apply(arr1, 2L, rowSums)

  # Make resultant array into df
  result_2 <- as.data.frame(result)
  # Create column name for Years column
  result_3 <- as.data.frame(tibble::rownames_to_column(result_2, "Year"))
  # Add year # column
  result_3$Year.No <- seq(1:t)


  #################### VACCINATED ####################
  
  # Population size
  # 154566548 = all women
  Pop_size_0 = 35704873*0.5 # <18*proportion vaxxed
  Pop_size_1 = 14683822*0.5 # 18-24*proportion vaxxed
  Pop_size_2 = 10201392*0.5 # 25-29*proportion vaxxed
  Pop_size_3 = 19939084*0.5 # 30-39*proportion vaxxed
  Pop_size_4 = (154566548-Pop_size_0-Pop_size_1-Pop_size_2-Pop_size_3)*0.5 # (40+)*proportion vaxxed
  
  #### Create array ####
  
  arr1 = array(NA, dim=c(t, N.states, 5), dimnames=list(years, c("Normal","Undet_LSIL","Det_LSIL","Undet_HSIL","Det_HSIL","Undet_Cancer","Det_Cancer","Cancer_Death"), c("18-20","21-24","25-29","30-39","40+")))
  
  # Assign starting prevalence of each state
  prev1 = 0.8 # normal
  prev2 = 0.1 # LSIL undetected
  prev3 = 0 # LSIL detected
  prev4 = 0.1 # LSIL detected
  prev5 = 0 # HSIL detected
  prev6 = 0 # cancer undetected
  prev7 = 0 # cancer detected
  prev8 = 0 # cancer death
  
  # Assign starting states
  set.seed(123)
  arr1[1,,] <- rmultinom(1, (Pop_size_1+Pop_size_2+Pop_size_3+Pop_size_4), prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 
  
  
  ### Run model ###
  
  for(i in 2:t){
    
    ########################### Variable rates ########################
    
    t.index <- i
    
    # vaccination
    norm_ulsil_0 <- ifelse(t.index<988, 0.15, 0.15*0.2) # starting in year 988 (2008), vaccination drops incidence by 80%
    norm_ulsil_1 <- ifelse(t.index<988, 0.08, 0.08*0.2) # starting in year 988 (2008), vaccination drops incidence by 80%
    norm_ulsil_2 <- ifelse(t.index<988, 0.02, 0.02*0.2) # starting in year 988 (2008), vaccination drops incidence by 80%
    
    # screening
    ulsil_dlsil <- ifelse(t.index<1000|t.index>1000,0.83/3,(0.83/3)*(1-decrease)) # in the year 1000 (2020), screening drops "decrease"%
    uhsil_dhsil <- ifelse(t.index<1000|t.index>1000,0.83/3,(0.83/3)*(1-decrease)) # in the year 1000 (2020), screening drops "decrease"%
    ucan_dcan <- 0.83/3 # cancer screening does not change (symptomatic)
  
    # loss to follow up
    dlsil_uhsil <- ifelse(t.index<1000|t.index>1000,0.17/3,(0.17/3)*(1+decrease)) # in the year 1000 (2020), LTFU increases "1+decrease"%
    dhsil_ucan <- ifelse(t.index<1000|t.index>1000,0.17/3,(0.17/3)*(1+decrease)) # in the year 1000 (2020), LTFU increases "1+decrease"%

    ########################### 18-20 ########################
    
    # Normal 18-20 (1)
    arr1[i,1,1] <- (arr1[(i-1),1,1]) +  
      (Pop_size_1*age_in) +             # age in
      (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
      (arr1[(i-1),1,1])*norm_ulsil_1 -  # normal progress to undetected LSIL
      (arr1[(i-1),1,1])*age_up_0 -      # age up
      (arr1[(i-1),1,1])*die_1           # die (unrelated)
    # LSIL undetected 18-20 (2)
    arr1[i,2,1] <- (arr1[(i-1),2,1]) +  
      (arr1[(i-1),1,1])*norm_ulsil_0 +  # normal progress to undetected LSIL
      (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),2,1])*ulsil_norm_1 -  # undetected LSIL regress to normal
      (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL
      (arr1[(i-1),2,1])*age_up_0 -      # age up
      (arr1[(i-1),2,1])*die_1           # die (unrelated)
    # LSIL detected 18-20 (3) = 0
    arr1[i,3,1] <- (arr1[(i-1),3,1])    # stay the same (0)
    # HSIL undetected 18-20 (4) 
    arr1[i,4,1] <- (arr1[(i-1),4,1]) + 
      (arr1[(i-1),2,1])*ulsil_uhsil_1 - # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),4,1])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,1])*age_up_0 -      # age up
      (arr1[(i-1),4,1])*die_1           # die (unrelated)
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
      (arr1[(i-1),1,2])*die_1           # die (unrelated)
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
      (arr1[(i-1),2,2])*die_1           # die (unrelated)
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
      (arr1[(i-1),3,2])*die_1           # die (unrelated)
    # HSIL undetected 21-24 (4) 
    arr1[i,4,2] <- (arr1[(i-1),4,2]) + 
      (arr1[(i-1),4,1])*age_up_0 +      # age up
      (arr1[(i-1),2,2])*ulsil_uhsil_1 + # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,2])*dlsil_uhsil -   # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,2])*uhsil_ulsil_1 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,2])*age_up_1 -      # age up
      (arr1[(i-1),4,2])*hyst_1-         # hysterectomy
      (arr1[(i-1),4,2])*die_1           # die (unrelated)
    # HSIL detected 21-24 (5)
    arr1[i,5,2] <- (arr1[(i-1),5,2]) + 
      (arr1[(i-1),5,1])*age_up_0 +      # age up
      (arr1[(i-1),3,2])*dlsil_dhsil_1 + # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,2])*uhsil_dhsil -   # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,2])*dhsil_dlsil_1 - # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,2])*dhsil_norm -    # detected HSIL treated to normal
      (arr1[(i-1),5,2])*age_up_1 -      # age up
      (arr1[(i-1),5,2])*hyst_1-         # hysterectomy
      (arr1[(i-1),5,2])*die_1           # die (unrelated)
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
      (arr1[(i-1),1,3])*die_2_3           # die (unrelated)
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
      (arr1[(i-1),2,3])*die_2_3           # die (unrelated)
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
      (arr1[(i-1),3,3])*die_2_3           # die (unrelated)
    # HSIL undetected 25-29 (4) 
    arr1[i,4,3] <- (arr1[(i-1),4,3]) + 
      (arr1[(i-1),4,2])*age_up_1 +        # age up
      (arr1[(i-1),2,3])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,3])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,3])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,3])*uhsil_ucan_2_3 -  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,3])*age_up_2 -        # age up
      (arr1[(i-1),4,3])*hyst_1-           # hysterectomy
      (arr1[(i-1),4,3])*die_2_3           # die (unrelated)
    # HSIL detected 25-29 (5)
    arr1[i,5,3] <- (arr1[(i-1),5,3]) + 
      (arr1[(i-1),5,2])*age_up_1 +        # age up
      (arr1[(i-1),3,3])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,3])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,3])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),5,3])*dhsil_dcan_2_3 -  # detected HSIL progress to detected cancer
      (arr1[(i-1),5,3])*dhsil_norm -      # detected HSIL treated to normal
      (arr1[(i-1),5,3])*age_up_2 -        # age up
      (arr1[(i-1),5,3])*hyst_1-           # hysterectomy
      (arr1[(i-1),5,3])*die_2_3           # die (unrelated)
    # Cancer undetected 25-29 (6)
    arr1[i,6,3] <- (arr1[(i-1),6,3]) +
      (arr1[(i-1),4,3])*uhsil_ucan_2_3 +  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),5,3])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),6,3])*age_up_2 -        # age up
      (arr1[(i-1),6,3])*hyst_1-           # hysterectomy
      (arr1[(i-1),6,3])*die_2_3           # die (unrelated)
    # Cancer detected 25-29 (7)
    arr1[i,7,3] <- (arr1[(i-1),7,3]) +
      (arr1[(i-1),5,3])*dhsil_dcan_2_3 +  # detected HSIL progress to detected cancer
      (arr1[(i-1),6,3])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),7,3])*dcan_dcandeath -  # detected cancer progress to cancer death
      (arr1[(i-1),7,3])*dcan_norm -       # detected cancer treated to normal
      (arr1[(i-1),7,3])*age_up_2 -        # age up
      (arr1[(i-1),7,3])*die_2_3           # die (unrelated)
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
      (arr1[(i-1),1,4])*age_up_3 -        # age up
      (arr1[(i-1),1,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),1,4])*die_2_3           # die (unrelated)
    # LSIL undetected 30-39 (2)
    arr1[i,2,4] <- (arr1[(i-1),2,4]) + 
      (arr1[(i-1),2,3])*age_up_2 +        # age up
      (arr1[(i-1),1,4])*norm_ulsil_3 +    # normal progress to undetected LSIL
      (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),2,4])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
      (arr1[(i-1),2,4])*ulsil_uhsil_2_3 - # undetected LSIL progress to undetected HSIL
      (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),2,4])*age_up_3 -        # age up
      (arr1[(i-1),2,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),2,4])*die_2_3           # die (unrelated)
    # LSIL detected 30-39 (3)
    arr1[i,3,4] <- (arr1[(i-1),3,4]) + 
      (arr1[(i-1),3,3])*age_up_2 +        # age up
      (arr1[(i-1),1,4])*norm_dlsil_3 +    # normal progress to detected LSIL
      (arr1[(i-1),5,4])*dhsil_dlsil_2_3 + # detected HSIL regress to detected LSIL
      (arr1[(i-1),2,4])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),3,4])*dlsil_norm_2_3 -  # detected LSIL regress to normal
      (arr1[(i-1),3,4])*dlsil_dhsil_2_3 - # detected LSIL progress to detected HSIL
      (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),3,4])*age_up_3 -        # age up
      (arr1[(i-1),3,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),3,4])*die_2_3           # die (unrelated)
    # HSIL undetected 30-39 (4) 
    arr1[i,4,4] <- (arr1[(i-1),4,4]) + 
      (arr1[(i-1),4,3])*age_up_2 +        # age up
      (arr1[(i-1),2,4])*ulsil_uhsil_2_3 + # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,4])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,4])*uhsil_ulsil_2_3 - # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,4])*uhsil_ucan_2_3 -  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,4])*age_up_3 -        # age up
      (arr1[(i-1),4,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),4,4])*die_2_3           # die (unrelated)
    # HSIL detected 30-39 (5)
    arr1[i,5,4] <- (arr1[(i-1),5,4]) + 
      (arr1[(i-1),5,3])*age_up_2 +        # age up
      (arr1[(i-1),3,4])*dlsil_dhsil_2_3 + # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,4])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,4])*dhsil_dlsil_2_3 - # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),5,4])*dhsil_dcan_2_3 -  # detected HSIL progress to detected cancer
      (arr1[(i-1),5,4])*dhsil_norm -      # detected HSIL treated to normal
      (arr1[(i-1),5,4])*age_up_3 -        # age up
      (arr1[(i-1),5,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),5,4])*die_2_3           # die (unrelated)
    # Cancer undetected 25-39 (6)
    arr1[i,6,4] <- (arr1[(i-1),6,4]) +
      (arr1[(i-1),6,3])*age_up_2 +        # age up
      (arr1[(i-1),4,4])*uhsil_ucan_2_3 +  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),5,4])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),6,4])*age_up_3 -        # age up
      (arr1[(i-1),6,4])*hyst_2-           # hysterectomy
      (arr1[(i-1),6,4])*die_2_3           # die (unrelated)
    # Cancer detected 30-39 (7)
    arr1[i,7,4] <- (arr1[(i-1),7,4]) +
      (arr1[(i-1),7,3])*age_up_2 +        # age up
      (arr1[(i-1),5,4])*dhsil_dcan_2_3 +  # detected HSIL progress to detected cancer
      (arr1[(i-1),6,4])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),7,4])*dcan_dcandeath -  # detected cancer progress to cancer death
      (arr1[(i-1),7,4])*dcan_norm -       # detected cancer treated to normal
      (arr1[(i-1),7,4])*age_up_3 -        # age up
      (arr1[(i-1),7,4])*die_2_3           # die (unrelated)
    # Cancer deaths 25-39 (8)
    arr1[i,8,4] <-
      (arr1[(i-1),7,4])*dcan_dcandeath    # detected cancer progress to cancer death
    
    ########################### 40+ ########################
    
    # Normal 40+ (1)
    arr1[i,1,5] <- (arr1[(i-1),1,5]) +  
      (arr1[(i-1),1,4])*age_up_3 +        # age up
      (arr1[(i-1),2,5])*ulsil_norm_2_3 +  # undetected LSIL regress to normal
      (arr1[(i-1),5,5])*dhsil_norm +      # detected HSIL treated to normal
      (arr1[(i-1),7,5])*dcan_norm +       # detected cancer treated to normal
      (arr1[(i-1),3,5])*dlsil_norm_2_3 -  # detected LSIL regress to normal
      (arr1[(i-1),1,5])*norm_ulsil_4 -    # normal progress to undetected LSIL
      (arr1[(i-1),1,5])*norm_dlsil_4 -    # normal progress to detected LSIL
      (arr1[(i-1),1,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),1,5])*die_4             # die (unrelated)
    # LSIL undetected 40+ (2)
    arr1[i,2,5] <- (arr1[(i-1),2,5]) + 
      (arr1[(i-1),2,4])*age_up_3 +        # age up
      (arr1[(i-1),1,5])*norm_ulsil_4 +    # normal progress to undetected LSIL
      (arr1[(i-1),4,5])*uhsil_ulsil_4 -   # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),2,5])*ulsil_norm_2_3 -  # undetected LSIL regress to normal
      (arr1[(i-1),2,5])*ulsil_uhsil_4 -   # undetected LSIL progress to undetected HSIL
      (arr1[(i-1),2,5])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),2,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),2,5])*die_4             # die (unrelated)
    # LSIL detected 40+ (3)
    arr1[i,3,5] <- (arr1[(i-1),3,5]) + 
      (arr1[(i-1),3,4])*age_up_3 +        # age up
      (arr1[(i-1),1,5])*norm_dlsil_4 +    # normal progress to detected LSIL
      (arr1[(i-1),5,5])*dhsil_dlsil_4 +   # detected HSIL regress to detected LSIL
      (arr1[(i-1),2,5])*ulsil_dlsil -     # undetected LSIL screened to detected LSIL
      (arr1[(i-1),3,5])*dlsil_norm_2_3 -  # detected LSIL regress to normal
      (arr1[(i-1),3,5])*dlsil_dhsil_4 -   # detected LSIL progress to detected HSIL
      (arr1[(i-1),3,5])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),3,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),3,5])*die_4             # die (unrelated)
    # HSIL undetected 40+ (4) 
    arr1[i,4,5] <- (arr1[(i-1),4,5]) + 
      (arr1[(i-1),4,4])*age_up_3 +        # age up
      (arr1[(i-1),2,5])*ulsil_uhsil_4 +   # undetected LSIL progress to undetected HSIL 
      (arr1[(i-1),3,5])*dlsil_uhsil -     # detected LSIL LTFU to undetected HSIL
      (arr1[(i-1),4,5])*uhsil_ulsil_4 -   # undetected HSIL regress to undetected LSIL
      (arr1[(i-1),4,5])*uhsil_ucan_2_3 -  # undetected HSIL progress to undetected cancer
      (arr1[(i-1),4,5])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),4,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),4,5])*die_4             # die (unrelated)
    # HSIL detected 40+ (5)
    arr1[i,5,5] <- (arr1[(i-1),5,5]) + 
      (arr1[(i-1),5,4])*age_up_3 +        # age up
      (arr1[(i-1),3,5])*dlsil_dhsil_4 +   # detected LSIL progress to detected HSIL
      (arr1[(i-1),4,5])*uhsil_dhsil -     # undetected HSIL screened to detected HSIL
      (arr1[(i-1),5,5])*dhsil_dlsil_4 -   # detected HSIL regress to detected LSIL
      (arr1[(i-1),5,5])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),5,5])*dhsil_dcan_4 -    # detected HSIL progress to detected cancer
      (arr1[(i-1),5,5])*dhsil_norm -      # detected HSIL treated to normal
      (arr1[(i-1),5,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),5,5])*die_4             # die (unrelated)
    # Cancer undetected 40+ (6)
    arr1[i,6,5] <- (arr1[(i-1),6,5]) +
      (arr1[(i-1),6,4])*age_up_3 +        # age up
      (arr1[(i-1),4,5])*uhsil_ucan_4 +    # undetected HSIL progress to undetected cancer
      (arr1[(i-1),5,5])*dhsil_ucan -      # detected HSIL LTFU to undetected cancer
      (arr1[(i-1),6,5])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),6,5])*hyst_3-           # hysterectomy
      (arr1[(i-1),6,5])*die_4             # die (unrelated)
    # Cancer detected 40+ (7)
    arr1[i,7,5] <- (arr1[(i-1),7,5]) +
      (arr1[(i-1),7,4])*age_up_3 +        # age up
      (arr1[(i-1),5,5])*dhsil_dcan_4 +    # detected HSIL progress to detected cancer
      (arr1[(i-1),6,5])*ucan_dcan -       # undetected cancer screened to detected cancer
      (arr1[(i-1),7,5])*dcan_dcandeath -  # detected cancer progress to cancer death
      (arr1[(i-1),7,5])*dcan_norm -       # detected cancer treated to normal
      (arr1[(i-1),7,5])*die_4             # die (unrelated)
    # Cancer deaths 25-39 (8)
    arr1[i,8,5] <-
      (arr1[(i-1),7,5])*dcan_dcandeath    # detected cancer progress to cancer death

  }
  
  # Round to whole numbers
  arr1 <- round(arr1,0)
    
  # Collapse array - remove age stratification
  result_5 <- apply(arr1, 2L, rowSums)

  # Create df from resultant array
  result_6 <- as.data.frame(result_5)
  # Add column name to Years column
  result_7 <- as.data.frame(tibble::rownames_to_column(result_6, "Year"))
  # Add year # column
  result_7$Year.No <- seq(1:t)



  ################### COMBINED RESULTS ###################

  # Combine vaccinated and unvaccinated into one df (entire df)
  result_tot <- result_3[,-1] + result_7[,-1]
  #View(result_tot)

  # After adding, divide by two to fix Year #
  result_tot$Year.No <- result_tot$Year.No/2
  # Subtract to match year # to calendar year
  result_tot$Year <- result_tot$Year.No+(2008-987)

  # Combined LSIL undetected and detected
  result_tot$All_LSIL <- result_tot$Undet_LSIL + result_tot$Det_LSIL
  # Combine HSIL undetected and detected
  result_tot$All_HSIL <- result_tot$Undet_HSIL + result_tot$Det_HSIL
  # Combined Cancer undetected and detected 
  result_tot$All_Cancer <- result_tot$Undet_Cancer + result_tot$Det_Cancer

  return(result_tot)
}


# Run model
# 60% decrease in screening, 60% increase in LFTU
screen_0.6 <- screen_decr(0.6)
# 40% decrease in screening, 40% increase in LFTU
screen_0.4 <- screen_decr(0.4)
# 80% decrease in screening, 80% increase in LFTU
screen_0.8 <- screen_decr(0.8)
# rename
result_tot <- screen_0.6


