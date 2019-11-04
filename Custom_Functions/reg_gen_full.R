reg_gen_full <- function(n , beta , eta.x , eta.y , seed) {      #Actual data generation
       beta <- beta  #store beta value
       seed <- seed  #store seed
       #generate uncontam. X values
       X.UC <- rnorm(floor((1 - eta.x)*n) , mean = 0 , sd = 1)
       if(eta.x > 0) {                             #generate contam. X values
              X.C <- rnorm(ceiling(eta.x*n) , mean <- 10 , sd = 1)
              X <- c(X.UC , X.C)
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
       Y <- X * beta + err                      #generate Y values
       reg.sum <- summary(lm(Y ~ X))            #generate regression results
       combine <- list(Y = Y , X = X , err = err , n = n ,
                       eta.x = eta.x , eta.y = eta.y ,
                       seed = seed , reg.summary = reg.sum)        #create combined list of all values
       return(combine)                       #save combined list of all values
}