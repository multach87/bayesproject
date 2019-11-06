#libraries
library(mvtnorm)
library(magrittr)
library(purrr)

#generate initial 108 data conditions
sim.structure1 <- as.data.frame(matrix(ncol = 4 , nrow = 108))
{
       colnames(sim.structure1) <- c("n" , "p" , "eta.x" , "eta.y")
       sim.structure1[ , "eta.x"] <- c(rep(c(0.0 , 0.1 , 0.2)))
       sim.structure1[ , "eta.y"] <- c(rep(0.0 , 3) , rep(0.1 , 3) , rep(0.2 , 3))
       sim.structure1[ , "p"] <- c(rep(8 , 9) , rep(30 , 9) , rep(50 , 9))
       sim.structure1[ , "n"] <- c(rep(25 , 27) , rep(50 , 27) , rep(100 , 27) , rep(200 , 27))
}
View(sim.structure1)
#generate repped conditions dataframe
sim.structure.repped <- as.data.frame(matrix(ncol = 5 , nrow = (108*1000)))
colnames(sim.structure.repped) <- c("n" , "p" , "eta.x" , "eta.y" , "seed")
for(i in 1:nrow(sim.structure1)) {
       sim.structure.repped[ ((1000*(i - 1)) + 1): (1000*i), (1:4)] <- 
              purrr::map_dfr(seq_len(1000) , ~sim.structure1[i , ])
}
#create checks indices just in case
checks.index <- numeric(216)
for(i in 1:108) {
       checks.index[((2*i) - 1)] <- (((i-1)*1000) + 1)
       checks.index[(2*i)] <- (i*1000)
}
sim.structure.repped[ , "seed"] <- rnorm((108*1000))
head(sim.structure.repped)
#get seed from sim.structure rather than internally-generating
data.gen <- function(n , p , eta.x , eta.y , seed) {      
       betas <- matrix(0 , nrow = p , ncol = 1)
       betas[1,1] <- 0.5
       betas[2,1] <- 1.0
       betas[3,1] <- 1.5
       betas[4,1] <- 2.0
       seed <- seed                       #set seed
       covar.X <- matrix(rep(0 , p^2) , ncol = p)  #generate covariance matrix
       diag(covar.X) <- 1                          #1's along cavariance diagonal
       X.UC <- rmvnorm(floor((1 - eta.x)*n) , mean = rep(0 , p) , sigma = covar.X)
       #generate uncontam. X values
       if(eta.x > 0) {                             #generate contam. X values
              X.C <- rmvnorm(ceiling(eta.x*n) , mean <- rep(10 , p) , sigma = covar.X)
              X <- rbind(X.UC , X.C)
       } else {
              X.C <- 0
              X <- X.UC
       }
       err.UC <- rnorm(floor((1-eta.y)*n) , mean = 0 , sd = 1)   #generate uncontom. residuals
       if(eta.y > 0) {                                           #generate contam. residuals
              err.C <- rnorm(ceiling(eta.y*n) , mean = 2 , sd = 5)
              err <- c(err.UC , err.C)
       } else {
              err.c <- 0
              err <- err.UC
       }
       Y <- X %*% betas[ , 1] + err                                    #generate Y values
       combine <- list(Y = Y , X = X , err = err)        #create combined list of all values
       return(combine)                       #save combined list of all values
}
#TESTING Again
data.test <- sim.structure.repped %>%   
       pmap(data.gen)








