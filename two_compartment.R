#test1

#One person
t <- 100
N.people <- 5
N.states <- 3

prev <- 0.7
p_0_1 <- 0.2
p_1_0 <- 0.5
p_1_2 <- 0.5
p_2_1 <- 0.5


mat1 <- array(NA, dim=c(t,N.states,N.people) )

set.seed(123)
mat1[1,1,] <- rbinom(N.people, 1, prob=prev) 
mat1[1,2,] <- 1 - mat1[1,1,]
mat1[1,3,] <- 0

for(i in 2:t){
  for(j in 1:N.people){
  mat1[i,,j] <- (mat1[(i-1),1,j]==1)*rbinom(1, 1, prob=p_0_1) +
    
               (mat1[(i-1),2,j]==1)*rmultinom(n, size=1, c(p_1_0, (1-(p_1_2+p_1_0)) , p_1_2) ) 

       # (mat1[(i-1),2,j]==1)*rmultinom(n, size=1, c(p_1_0, (1-(p_1_2+p_1_0)) , p_1_2) ) 
  
  }
}


