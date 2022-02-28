
## FINAL MODEL 2018 PARAMETERS ##

## Time step = MONTH ##


### INCLUDE LTFU


### UNVACCINATED ###

#### Starting parameters
t = 12*500
  # 12*500 = burn-in (1-6000)
  # 12*30 = COVID + 30ish years (6001-6360)

N.states = 8

Pop_size_1 = 14683822*0.4
Pop_size_2 = 10201392*0.4
Pop_size_3 = 19939084*0.4

# aging into cohort
age_in <- -log(1-1/20)/12 # multiply by pop under 18
# aging up within cohort
age_up_0 <- -log(1-1/3)/12 # proportion turning 21
age_up_1 <- -log(1-1/4)/12 # proportion turning 25
age_up_2 <- -log(1-1/5)/12 # proportion turning 30
# aging out of cohort
age_out <- -log(1-1/10)/12 # proportion turning 40
# Dying
p_die_1 <- -log(1-74/100000)/12 # proportion dying age 18-24
p_die_2_3 <- -log(1-164/100000)/12 # proportion dying age 25-39

# undetected
norm_ulsil_0 <- -log(1-0.15)/12
norm_ulsil_1 <- -log(1-0.08)/12
norm_ulsil_2 <- -log(1-0.02)/12
norm_ulsil_3 <- -log(1-0.01)/12
ulsil_norm_1 <- -log(1-0.60)/12
ulsil_norm_2_3 <- -log(1-0.4)/12
ulsil_uhsil_1 <- -log(1-0.14)/12
ulsil_uhsil_2_3 <- -log(1-0.35)/12
uhsil_ulsil_1 <- -log(1-0.62)/12
uhsil_ulsil_2_3 <- -log(1-0.20)/12
uhsil_ucan <- -log(1-6.4/100000)/12

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
dcan_dcandeath <- -log(1-0.35)/12

# treatment
dhsil_norm <- -log(1-0.9*0.9)/12
dcan_norm <- -log(1-0.5*1)/12

# screening
ulsil_dlsil <- -log(1-0.83/3)/12
uhsil_dhsil <- -log(1-0.83/3)/12
ucan_dcan <- -log(1-0.83/3)/12

#ulsil_dlsil <- ifelse(t<6000|t>6036,-log(1-0.41/3),0)
#uhsil_dhsil <- ifelse(t<6000|t>6036,-log(1-0.83/3),0)
#ucan_dcan <- ifelse(t<6000|t>6036,-log(1-0.83/3),0)

#ulsil_dlsil <- ifelse(t<6000,-log(1-0.83/3),1)
#uhsil_dhsil <- ifelse(t<6000,-log(1-0.83/3),1)
#ucan_dcan <- ifelse(t<6000,-log(1-0.83/3),1)


# loss to follow up
dlsil_uhsil <- -log(1-0.17/3)/12
dhsil_ucan <- -log(1-0.17/3)/12

# hysterectomies
hyst_1 <- -log(1+log(0.99)/10)/12
hyst_2 <- -log(1+log(0.96)/10)/12


# Starting number in each state

# Create years label
prefix = "Month"
suffix = seq(1:t)
months = paste(prefix, suffix, sep=" ")

# Create empty array
arr1 = array(NA, dim=c(t, N.states, 4), dimnames=list(months, c("Normal","Undet_LSIL","Det_LSIL","Undet_HSIL","Det_HSIL","Undet_Cancer","Det_Cancer","Cancer Death"), c("18-20","21-24","25-29","30-39")))

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
arr1[1,,] <- rmultinom(1, Pop_size_1+Pop_size_2+Pop_size_3, prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 

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
result <- apply(arr1, 2L, rowSums)
final_result <- as.data.frame(result[max(t),])

# Log likelihood of result vs predicted
LL <- sum(dpois(final_result[c(5,7,8),], c(196000*0.4,10510*0.4,3400*0.4), log=TRUE))

Results <- final_result[c(5,7,8),]

result_2 <- as.data.frame(result)
result_3 <- as.data.frame(tibble::rownames_to_column(result_2, "Year"))
result_3$Year.No <- seq(1:t)


p_Norm<-ggplot(result_3, aes(x=Year.No)) +
  geom_line(aes(y=Normal), colour="red") + 
  coord_cartesian(
    xlim = NULL,
    ylim = c(0,10000000))

p_Norm


p_LSIL<-ggplot(result_3, aes(x=Year.No)) +
  geom_line(aes(y=Undet_LSIL), colour="red") + 
  geom_line(aes(y=Det_LSIL), colour="blue") +
  coord_cartesian(
    xlim = NULL,
    ylim = c(0,500000))

p_LSIL

p_HSIL<-ggplot(result_3, aes(x=Year.No)) +
  geom_line(aes(y=Undet_HSIL), colour="red") + 
  geom_line(aes(y=Det_HSIL), colour="blue") +
  coord_cartesian(
    xlim = NULL,
    ylim = c(0,200000))

p_HSIL

# Make vector of paramters from monthly vs annual model - should form line
# Make dataframe - do operation vector 



### VACCINATED ###

#### Starting parameters
t = 12*500
  # 12*500 = burn-in (1-6000)
  # 12*30 = COVID + 30ish years (6001-6360)

N.states = 8
#Pop_size = 154566548 #all women
#Pop_size = 44978865*(1-0.032) # women 18-39*adjustment for hysterectomies
Pop_size_1 = 14683822*0.6
Pop_size_2 = 10201392*0.6
Pop_size_3 = 19939084*0.6

# aging into cohort
age_in <- -log(1-(1/20))/12 # multiply by pop under 18
# aging up within cohort
age_up_0 <- -log(1-1/3)/12 # proportion turning 21
age_up_1 <- -log(1-1/4)/12 # proportion turning 25
age_up_2 <- -log(1-1/5)/12 # proportion turning 30
# aging out of cohort
age_out <- -log(1-1/10)/12 # proportion turning 40
# Dying
p_die_1 <- -log(1-74/100000)/12 # proportion dying age 18-24
p_die_2_3 <- -log(1-164/100000)/12 # proportion dying age 25-39

# undetected
ifelse(t<(490*12), norm_ulsil_0 <- -log(1-0.15)/12, norm_ulsil_0 <- -log(1-0.15*0.2)/12)
ifelse(t<(490*12), norm_ulsil_1 <- -log(1-0.08)/12, norm_ulsil_1 <- -log(1-0.08*0.2)/12)
ifelse(t<(490*12), norm_ulsil_2 <- -log(1-0.02)/12, norm_ulsil_2 <- -log(1-0.02*0.2)/12)
norm_ulsil_3 <- -log(1-0.01*0.1)/12
ulsil_norm_1 <- -log(1-0.60)/12
ulsil_norm_2_3 <- -log(1-0.4)/12
ulsil_uhsil_1 <- -log(1-0.14)/12
ulsil_uhsil_2_3 <- -log(1-0.35)/12
uhsil_ulsil_1 <- -log(1-0.62)/12
uhsil_ulsil_2_3 <- -log(1-0.20)/12
uhsil_ucan <- -log(1-6.4/100000)/12

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
dcan_dcandeath <- -log(1-0.35)/12

# treatment
dhsil_norm <- -log(1-0.9*0.9)/12
dcan_norm <- -log(1-0.5*1)/12

# BRFSS rates!!
# screening
ulsil_dlsil <- -log(1-0.83/3)/12
uhsil_dhsil <- -log(1-0.83/3)/12
ucan_dcan <- -log(1-0.83/3)/12

# loss to follow up
dlsil_uhsil <- -log(1-0.17/3)/12
dhsil_ucan <- -log(1-0.17/3)/12

# hysterectomies
hyst_1 <- -log(1+log(0.99)/10)/12
hyst_2 <- -log(1+log(0.96)/10)/12

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
arr1[1,,] <- rmultinom(1, Pop_size_1+Pop_size_2+Pop_size_3, prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7,prev8)) 

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
result_2 <- apply(arr1, 2L, rowSums)
final_result_2 <- as.data.frame(result_2[nrow(result_2),])

# Log likelihood of result vs predicted
LL_2 <- sum(dpois(final_result[c(5,7,8),], c(196000*0.6,10510*0.6,3400*0.6), log=TRUE))

Results_2 <- final_result_2[c(5,7,8),]


### COMBINED RESULTS ###

Results_all = Results+Results_2

Results_all

