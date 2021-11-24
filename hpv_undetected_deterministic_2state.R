### MARKOV MODEL HPV DRAFT
## Gwen Oliver
## Nov 12, 2021

### This model is for the undetected compartments ONLY

# 20 years
t = 1000

# 7 states: normal, infected, CIN 1, CIN 2, CIN 3, cancer, cancer deaths
N.states = 2


Pop_size =1000 
# Total prob leaving each state needs to add to 1?

# aging into cohort
p_age_in <- 1/120
p_age_out =p_age_in
# Rough estimate based on proportion of US pop that is female, 10-14
# https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/


p_1_2 <- 0.02 # normal > infected
p_2_1 <- 0.1 # infected > normal

p_1_1 = 1 - p_age_out - p_1_2  # Stay normal
p_2_2 = 1-  p_age_out -  p_2_1 # Stay infected

# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:t)
years = paste(prefix, suffix, sep=" ")

# Create empty array
mat1 = matrix(NA, t, N.states, dimnames=list(years, c("Normal","Infected")))

# Assign starting prevalence of each state
prev1 = 0.9
prev2 = 0.1

# Assign starting states
set.seed(123)
mat1[1,] <- rmultinom(1, Pop_size, prob=c(prev1,prev2)) 

# Run model
for(i in 2:t){
  # Normal
  mat1[i,1] <- (mat1[(i-1),1]) +      #N people in previous time step
    (Pop_size*p_age_in) +      #aging in #0.05  (A)
    (mat1[(i-1),2])*p_2_1 -      #regress from inf #0.2 (B)
    (mat1[(i-1),1])*p_1_2 -      #progress to inf #0.25 (B)
    (mat1[(i-1),1])*p_age_out        #age out #0.1 (A)
      # Infected   
  mat1[i,2] <-   (mat1[(i-1),2]) +      #N people in state 2 previous time step
    (mat1[(i-1),1])*p_1_2 -      #progress from normal state 1 #0.15 (B)
    (mat1[(i-1),2])*p_2_1 -      #regress to normal #0.2 (B)
    (mat1[(i-1),2])*p_age_out        #age out #0.05 (A)
 
}

# For steady state need (Pop_size*p_age_in) + (mat1[(i-1),1])*p_1_1 + (mat1[(i-1),2])*p_2_1 =
#                        (mat1[(i-1),2])*p_2_1 + (mat1[(i-1),1])*p_1_2 + (mat1[(i-1),1])*p_age_out 

mat1 <- round(mat1,0)

#Check Progression over time

States <- apply(mat1, c(1,2), sum)

plot(States[,'Normal'], type='l')

plot(States[,'Infected'], type='l')

pop_size_t <- apply(mat1,1, sum)

plot(pop_size_t)
