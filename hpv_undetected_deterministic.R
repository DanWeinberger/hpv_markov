### MARKOV MODEL HPV DRAFT 1
## Gwen Oliver
## Oct 22, 2021

### This model is for teh undetected compartments ONLY

# 20 years
t = 20

# 100 people
N.people = 100

# 6 states: normal, infected, CIN 1, CIN 2, CIN 3, cancer
N.states = 7


# Transition probabilities (made up)
# Total prob leaving each state needs to add to 1?
# Do I need to include p_1_1, etc or is this implied?

# births
p_0_1 = 0.3

# normal > normal
p_1_1 = 0.7

# normal > infected
p_1_2 <- 0.3

# infected > infected
p_2_2 = 0.5

# infected > normal
p_2_1 <- 0.2

# infected > CIN 1
p_2_3 = 0.3
p_3_3 = 0.5
p_3_2 = 0.2

# CIN 1 > CIN 2
p_3_4 = 0.3
p_4_4 = 0.5
p_4_3 = 0.2

# CIN 2 > CIN 3
p_4_5 = 0.3
p_5_5 = 0.5
p_5_4 = 0.2

# CIN 3 > cancer
p_5_6 = 0.3
p_6_6 = 1

# death
p_1_0
p_2_0
p_3_0
p_4_0
p_5_0
p_6_0 = 0.5


# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:20)
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
mat1[1,] <- rmultinom(1, 100, prob=c(prev1,prev2,prev3,prev4,prev5,prev6,prev7)) 

# Run model
for(i in 2:t){
                          #birth    #stay normal              #regress from inf     #progress to inf
      mat1[i,1] <- (mat1[(i-1),1]) + 10*p_0_1 + (mat1[(i-1),1])*p_1_1 + (mat1[(i-1),2])*p_2_1 - (mat1[(i-1),1])*p_1_2
    
      mat1[i,2] <- (mat1[(i-1),2]) + (mat1[(i-1),1])*p_1_2 + (mat1[(i-1),2])*p_2_2 + (mat1[(i-1),3])*p_3_2 - (mat1[(i-1),2])*p_2_1 - (mat1[(i-1),2])*p_2_3

      mat1[i,3] <- (mat1[(i-1),3]) + (mat1[(i-1),2])*p_2_3 + (mat1[(i-1),3])*p_3_3 + (mat1[(i-1),4])*p_4_3 - (mat1[(i-1),3])*p_3_2 - (mat1[(i-1),3])*p_3_4
    
      mat1[i,4] <- (mat1[(i-1),4]) + (mat1[(i-1),3])*p_3_4 + (mat1[(i-1),4])*p_4_4 + (mat1[(i-1),5])*p_5_4 - (mat1[(i-1),4])*p_4_3 - (mat1[(i-1),4])*p_4_5

      mat1[i,5] <- (mat1[(i-1),5]) + (mat1[(i-1),4])*p_4_5 + (mat1[(i-1),5])*p_5_5 - (mat1[(i-1),5])*p_5_4 - (mat1[(i-1),5])*p_5_6
      
      mat1[i,6] <- (mat1[(i-1),6]) + (mat1[(i-1),5])*p_5_6 + (mat1[(i-1),6])*p_6_6 - (mat1[(i-1),6])*p_6_0
      
      mat1[i,7] <- (mat1[(i-1),7]) + (mat1[(i-1),6])*p_6_0
      
}

mat1 <- round(mat1,0)


#Check Progression over time

States <- apply(mat1, c(1,2), sum)

plot(States[,'Normal'], type='l')

plot(States[,'Cancer'], type='l')





