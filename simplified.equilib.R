### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Dec 5, 2021

### This model is for the undetected compartments ONLY

# 1000 years
t = 1000

N.states = 5


Pop_size =1000 

# aging into cohort
p_age_in <- 1/20
p_age_out <- p_age_in

p_die <- 0.01

# Rough estimate based on proportion of US pop that is female, 10-14
# https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/


p_1_2 <- 0.1 # normal > infected
p_2_1 <- 0.5 # infected > normal

p_2_3 <- 0.18  # LSIL > HSIL
p_3_2 <- 0.2   # HSIL > LSIL

p_3_4 <- 0.04 # HSIL > cancer
p_4_5 <- 0.5 # Cancer > Cancer death


# Starting number in each state

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

# Assign starting states
set.seed(123)
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
