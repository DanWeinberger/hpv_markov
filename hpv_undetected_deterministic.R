### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Nov 28, 2021

### This model is for the undetected compartments ONLY

# 500 years
t = 500

# 1000 people
Pop_size = 1000

# 7 states: normal, infected, CIN 1, CIN 2, CIN 3, cancer, cancer deaths
N.states = 7


# Rough estimate based on proportion of US pop that is female, 10-14
# https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/

p_1_2 <- 0.1 # normal > infected

p_2_1 <- 0.5 # infected > normal
# Numbers vary widely for this parameter
p_2_3 = 0.15 # infected > CIN 1

p_3_2 <- 0.2 # CIN1 > infected
# These estimates vary widely
p_3_4 <- 0.18  # CIN1 > CIN2

p_4_5 <- 0.3 # CIN2 > CIN3
p_4_3 <- 0.13   # CIN2 > CIN1

p_5_4 <- 0.03 # CIN3 > CIN2
# Only one estimate for this, only for HPV 16/18
p_5_6 <- 0.2 # CIN 3 > cancer
# Two estimates for this, 10-fold difference between them

p_6_7 <- 0.5 # Cancer > Cancer death
# Don't have an estimate for this yet, but I think available on SEER

# Age into cohort
p_age_in <- 1/20
# Age out of cohort
p_age_out <- p_age_in

# Stay in compartment 
p_1_1 = 1-p_age_out-p_1_2 # Stay normal
p_2_2 = 1-p_age_out-p_2_1-p_2_3 # Stay infected
p_3_3 = 1-p_age_out-p_3_2-p_3_4 # Stay CIN1
p_4_4 = 1-p_age_out-p_4_5-p_4_3 # Stay CIN2
p_5_5 = 1-p_age_out-p_5_4-p_5_6 # Stay CIN3
p_6_6 = 1-p_age_out-p_6_7 # Stay cancer


# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:t)
years = paste(prefix, suffix, sep=" ")


# Create empty array
mat1 = matrix(NA, t, N.states, dimnames=list(years, c("Normal","Infected","CIN 1","CIN 2","CIN 3","Cancer","Death due to CC")))

# Assign starting prevalence of each state
prev1 = 0.29
prev2 = 0.5
prev3 = 0.1
prev4 = 0.05
prev5 = 0.05
prev6 = 0.01
prev7 = 0

# Assign starting states
set.seed(123)
mat1[1,] <- rmultinom(1, Pop_size, prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7)) 

# Run model
for(i in 2:t){
  # Normal
  mat1[i,1] <- (mat1[(i-1),1]) + #previous N people in normal
    Pop_size*p_age_in +      #aging in #0.05
    (mat1[(i-1),2])*p_2_1 -      #regress from inf #0.2
    (mat1[(i-1),1])*p_1_2 -      #progress to inf #0.25
    (mat1[(i-1),1])*p_age_out        #age out #0.1
  # Infected   
  mat1[i,2] <- (mat1[(i-1),2]) + #previous N people in infected
    (mat1[(i-1),1])*p_1_2 +      #progress from normal #0.15
    (mat1[(i-1),3])*p_3_2 -      #regress from CIN1 #0.2
    (mat1[(i-1),2])*p_2_1 -      #regress to normal #0.2
    (mat1[(i-1),2])*p_2_3 -      #progress to CIN1 #0.1
    (mat1[(i-1),2])*p_age_out        #age out #0.05
  # CIN1
  mat1[i,3] <- (mat1[(i-1),3]) + #previous N people in CIN1
    (mat1[(i-1),2])*p_2_3 +      #progress from infected #0.1
    (mat1[(i-1),4])*p_4_3 -      #regress from CIN2 #0.13
    (mat1[(i-1),3])*p_3_2 -      #regress to infected #0.2
    (mat1[(i-1),3])*p_3_4 -      #progress to CIN2 #0.1
    (mat1[(i-1),3])*p_age_out        #age out #1/21
  # CIN2
  mat1[i,4] <- (mat1[(i-1),4]) + #previous N people in CIN2
    (mat1[(i-1),3])*p_3_4 +      #progress from CIN1
    (mat1[(i-1),5])*p_5_4 -      #regress from CIN3
    (mat1[(i-1),4])*p_4_3 -      #regress to CIN1
    (mat1[(i-1),4])*p_4_5 -      #progress to CIN3
    (mat1[(i-1),4])*p_age_out        #age out
  # CIN3
  mat1[i,5] <- (mat1[(i-1),5]) + #previous N people in CIN3
    (mat1[(i-1),4])*p_4_5 -      #progress from CIN2
    (mat1[(i-1),5])*p_5_4 -      #regress to CIN2
    (mat1[(i-1),5])*p_5_6 -      #progress to cancer
    (mat1[(i-1),5])*p_age_out        #age out
  # Cancer
  mat1[i,6] <- (mat1[(i-1),6]) + #previous N people in cancer
    (mat1[(i-1),5])*p_5_6 -      #progress from CIN3
    (mat1[(i-1),6])*p_6_7 -      #die
    (mat1[(i-1),6])*p_age_out        #age out
  # Cumulative Cancer deaths
  mat1[i,7] <- (mat1[(i-1),7]) + 
    (mat1[(i-1),6])*p_6_7        #progress from cancer
  
}

mat1 <- round(mat1,0)


#Check Progression over time

States <- apply(mat1, c(1,2), sum)

plot(States[,'Normal'], type='l')

plot(States[,'Infected'], type='l')

plot(States[,'CIN 1'], type='l')

plot(States[,'CIN 2'], type='l')

plot(States[,'CIN 3'], type='l')

plot(States[,'Cancer'], type='l')

plot(States[,'Death due to CC'], type='l')

pop_size_t <- apply(States,1, sum)

plot(pop_size_t)

