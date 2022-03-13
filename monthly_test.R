

### Simple model - birth and death, monthly


b <- 0.1 # birth

m <- 0.11 # death (mortality)

t <- 100 # years

df1 <- matrix(data=NA, nrow=100, ncol=1, dimnames=list(1:100,"population"))

df1 <- as.data.frame(df1)

df1[1,] <- 1000

for(i in 2:t){
  
  df1[i,] <- df1[(i-1),] + df1[(i-1),]*b - df1[(i-1),]*m
  
}

df1 <- round(df1,0)


### Monthly version

b_2 <- -log(1-0.1)/12 # birth

m_2 <- -log(1-0.11)/12 # death (mortality)

t_2 <- 12*100 # years

df2 <- matrix(data=NA, nrow=12*100, ncol=1, dimnames=list(1:(12*100),"population"))

df2 <- as.data.frame(df2)

df2[1,] <- 1000

for(i in 2:t_2){
  
  df2[i,] <- df2[(i-1),] + df2[(i-1),]*b_2 - df2[(i-1),]*m_2
  
}

df2 <- round(df1,2)



tail(df1)
tail(df2)

