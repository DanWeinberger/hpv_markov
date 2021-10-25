### MARKOV MODEL HPV DRAFT 1
## Gwen Oliver
## Oct 22, 2021

### This model is for teh undetected compartments ONLY

# 20 years
t = 20

# 100 people
N.people = 100

# 6 states: normal, infected, CIN 1, CIN 2, CIN 3, cancer
N.states = 6


# Transition probabilities (made up)
# Total prob leaving each state needs to add to 1?
# Do I need to include p_1_1, etc or is this implied?
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


# Starting number in each state

# Create years label
prefix = "Year"
suffix = seq(1:20)
years = paste(prefix, suffix, sep=" ")

# Create person label
prefix2 = "Person"
suffix2 = seq(1:100)
person = paste(prefix2, suffix2, sep=" ")

# Create empty array
mat1 = array(NA, dim=c(t, N.states, N.people), dimnames=list(years, c("Normal","Infected","CIN 1","CIN 2","CIN 3","Cancer"),person))
# Assign starting prevalence of each state
prev1 = 0.29
prev2 = 0.5
prev3 = 0.1
prev4 = 0.05
prev5 = 0.05
prev6 = 0.01
# Assign starting states
mat1[1,,] <- rmultinom(N.people, 1, prob=c(prev1,prev2,prev3,prev4,prev5,prev6)) 

# Run model
set.seed(123)

for(i in 2:t){
  for(j in 1:N.people){
    
    if(mat1[(i-1),1,j]==1){
      mat1[i,,j] <- (mat1[(i-1),1,j]==1)*rmultinom(1, 1, prob=c(p_1_1,p_1_2,0,0,0,0))
    }
    
    if(mat1[(i-1),2,j]==1){
      mat1[i,,j] <- (mat1[(i-1),2,j]==1)*rmultinom(1, 1, prob=c(p_2_1,p_1_2,p_2_3,0,0,0))
    }
    
    if(mat1[(i-1),3,j]==1){
      mat1[i,,j] <- (mat1[(i-1),3,j]==1)*rmultinom(1, 1, prob=c(0,p_3_2,p_3_3,p_3_4,0,0))
    }
    
    if(mat1[(i-1),4,j]==1){
      mat1[i,,j] <- (mat1[(i-1),4,j]==1)*rmultinom(1, 1, prob=c(0,0,p_4_3,p_4_4,p_4_5,0))
    }
    
    if(mat1[(i-1),5,j]==1){
      mat1[i,,j] <- (mat1[(i-1),5,j]==1)*rmultinom(1, 1, prob=c(0,0,0,p_5_4,p_5_5,p_5_6))
    }
    
    if(mat1[(i-1),6,j]==1){
      mat1[i,,j] <- (mat1[(i-1),6,j]==1)*rmultinom(1, 1, prob=c(0,0,0,0,0,p_6_6))
    }
  }
}


#Check Progression over time

States <- apply(mat1, c(1,2), sum)

plot(States[,'Normal'], type='l')







