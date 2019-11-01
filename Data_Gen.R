
reg.gen.simple <- function(n , beta) {
       X <- rnorm(n , mean = 0 , sd = 5)    #generate X values; randomly chose mean, sd
       res <- rnorm(n)                       #generate residuals
       Y <- X*beta + res                     #generate Y values
       combine <- list(Y = Y , X = X , res = res)        #create combined list of all values
       return(combine)                       #save combined list of all values
}

#reg.gen.complex
{
        #libraries
        library(magrittr)
        library(purrr)
        
        #generate initial 108 data conditions
        sim.structure1 <- as.data.frame(matrix(ncol = 4 , nrow = CONS)) #MAKE SURE TO REPLACE "CONS
        {
                colnames(sim.structure1) <- c("n" , "beta" , "eta.x" , "eta.y")
                sim.structure1[ , "eta.x"] <- c(rep(c(0.0 , 0.1 , 0.2)))
                sim.structure1[ , "eta.y"] <- c(rep(0.0 , 3) , rep(0.1 , 3) , rep(0.2 , 3))
                sim.structure1[ , "beta"] <- c(rep(1 , 9) , rep(3 , 9) , rep(10 , 9))
                sim.structure1[ , "n"] <- c(rep(25 , 27) , rep(50 , 27) , rep(100 , 27) , rep(200 , 27))
        }
        View(sim.structure1)
        
        #generate repped conditions dataframe
        sim.structure.repped <- as.data.frame(matrix(ncol = 5 , 
                                                     nrow = (CONS*N.ITERS))) #MAKE SURE TO REPLACE
                                                                             ##"CONS" AND "N.ITERS"
        colnames(sim.structure.repped) <- c("n", "eta.x" , "eta.y" , "seed")
        for(i in 1:nrow(sim.structure1)) {
                sim.structure.repped[ ((1000*(i - 1)) + 1): (1000*i), (1:4)] <- 
                        purrr::map_dfr(seq_len(1000) , ~sim.structure1[i , ])
        }
        sim.structure.repped[ , "seed"] <- rnorm((CONS*N.ITERS))  #MAKE SURE TO REPLACE "CONS"
                                                                  #AND "N.ITERS"
        
        data.gen <- function(n , beta , eta.x , eta.y , seed) {      #Actual data generation
                beta <- beta  #store beta value
                seed <- seed  #store seed
                #generate uncontam. X values
                X.UC <- rnorm(floor((1 - eta.x)*n) , mean = 0 , sd = 1)
                if(eta.x > 0) {                             #generate contam. X values
                        X.C <- rnorm(ceiling(eta.x*n) , mean <- 10 , sd = 1)
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
                Y <- X * beta + err                      #generate Y values
                combine <- list(Y = Y , X = X , err = err ,
                                eta.x = eta.x , eta.y = eta.y ,
                                seed = seed)        #create combined list of all values
                return(combine)                       #save combined list of all values
        }
}



reg.list <- list()                           #create empty list to fill with loop
for(i in 1:1000) {
       cat("i = " , i , "\n")                #tracker for how much of loop is done
       low.list <- reg.gen.simple(n = 25 , beta = 3.7)   #create sub-list of simulated regression 
                                                         ##data
                                                         ##Randomly chose n, beta
       reg.list[[i]] <- low.list             #fill in meta-list with simulation data
       if(i == 1000) {                       #save combiend list of all simulated data if 
              return(reg.list)               ##end of loop; if you change the number of iterations,
       }                                     ##make sure to change if(i == ???)
}

reg.list.summary <- function(data) {         #function to model data from sub-list elements
       summary(lm(data$Y ~ data$X))
}
summary.list <- lapply(reg.list , reg.list.summary)     #apply above function to each iteration 
                                                        ##of simulated data

names(summary.list[[1]])         ##lots of info in these summaries













