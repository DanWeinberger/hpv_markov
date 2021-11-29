### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Nov 29, 2021

### This model is for the undetected compartments ONLY
### This model has not achieved an equilibirum state yet

# 500 years
t = 500

# 1000 people
Pop_size = 1000

# 5 states: 
# normal
# low-grade squamous intraepithelial lesion (LSIL, infection and CIN1)
# high-grade squamous intraepithelial lesion (HSIL, CIN2, CIN3, and AIS)
# cancer
# cancer deaths
N.states = 5


# Rough estimate based on proportion of US pop that is female, 10-14
# https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/

p_1_2 <- 0.1 # normal > LSIL

p_2_1 <- 0.5 # LSIL > normal
# Numbers vary widely for this parameter

p_2_3 <- 0.18  # LSIL > HSIL

p_3_2 <- 0.2   # HSIL > LSIL
# Numbers vary widely for this parameters

p_3_4 <- 0.04 # HSIL > cancer
# Two estimates for this, 10-fold difference between them

p_4_5 <- 0.5 # Cancer > Cancer death
# Don't have an estimate for this yet, but I think available on SEER

# Age into cohort
p_age_in <- 1/20
# Age out of cohort
p_age_out <- p_age_in

# Die of other cause
p_die <- 0.01


# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:t)
years = paste(prefix, suffix, sep=" ")


# Create empty array
mat1 = matrix(NA, t, N.states, dimnames=list(years, c("Normal","LSIL","HSIL","Cancer","Death due to CC")))

# Assign starting prevalence of each state
prev_norm = 0.35
prev_LSIL = 0.5
prev_HSIL = 0.1
prev_cancer = 0.05
prev_cancdeath = 0

# Assign starting states
set.seed(123)
mat1[1,] <- rmultinom(1, Pop_size, prob=c(prev_norm,prev_LSIL,prev_HSIL,prev_cancer,prev_cancdeath)) 

# Run model
for(i in 2:t){
  # Normal
  mat1[i,1] <- (mat1[(i-1),1]) + #previous N people in normal
    Pop_size*p_age_in +          #aging in
    (mat1[(i-1),2])*p_2_1 -      #regress from LSIL
    (mat1[(i-1),1])*p_1_2 -      #progress to LSIL
    (mat1[(i-1),1])*p_age_out -  #age out
    (mat1[(i-1),1])*p_die        #die
  # LSIL  
  mat1[i,2] <- (mat1[(i-1),2]) + #previous N people in LSIL
    (mat1[(i-1),1])*p_1_2 +      #progress from normal
    (mat1[(i-1),3])*p_3_2 -      #regress from LSIL
    (mat1[(i-1),2])*p_2_1 -      #regress to normal
    (mat1[(i-1),2])*p_2_3 -      #progress to HSIL
    (mat1[(i-1),2])*p_age_out -  #age out
    (mat1[(i-1),2])*p_die        #die
  # HSIL
  mat1[i,3] <- (mat1[(i-1),3]) + #previous N people in HSIL
    (mat1[(i-1),3])*p_2_3 -      #progress from LSIL
    (mat1[(i-1),4])*p_3_2 -      #regress to LSIL
    (mat1[(i-1),4])*p_3_4 -      #progress to cancer
    (mat1[(i-1),4])*p_age_out -  #age out
    (mat1[(i-1),3])*p_die        #die
  # Cancer
  mat1[i,4] <- (mat1[(i-1),4]) + #previous N people in cancer
    (mat1[(i-1),3])*p_3_4 -      #progress from HSIL
    (mat1[(i-1),4])*p_4_5 -      #die of cancer
    (mat1[(i-1),4])*p_age_out -  #age out
    (mat1[(i-1),4])*p_die        #die (other cause)
  # Cumulative Cancer deaths
  mat1[i,5] <- (mat1[(i-1),5]) + 
    (mat1[(i-1),4])*p_4_5        #progress from cancer
  
}

mat1 <- round(mat1,0)


#Check Progression over time

States <- apply(mat1, c(1,2), sum)

plot(States[,'Normal'], type='l')

plot(States[,'LSIL'], type='l')

plot(States[,'HSIL'], type='l')

plot(States[,'Cancer'], type='l')

plot(States[,'Death due to CC'], type='l')

pop_size_t <- apply(States,1, sum)

plot(pop_size_t)

