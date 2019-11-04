
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
        sim.structure1 <- as.data.frame(matrix(ncol = 4 , nrow = 135)) 
        {
                colnames(sim.structure1) <- c("n" , "beta" , "eta.x" , "eta.y")
                sim.structure1[ , "eta.x"] <- c(rep(c(0.0 , 0.1 , 0.2)))
                sim.structure1[ , "eta.y"] <- c(rep(0.0 , 3) , rep(0.1 , 3) , rep(0.2 , 3))
                sim.structure1[ , "beta"] <- c(rep(1 , 9) , rep(3 , 9) , rep(10 , 9))
                sim.structure1[ , "n"] <- c(rep(25 , 27) , rep(50 , 27) , rep(100 , 27) , 
                                            rep(200 , 27) , rep(1000 , 27))
        }
        View(sim.structure1)
        
        #generate repped conditions dataframe
        create.sim.repped <- function(sim.structure1 , num.conds = nrow(sim.structure1) , 
                                      num.pars ,  num.iters) {
                sim.structure.repped <- as.data.frame(matrix(ncol = (num.pars + 1) , 
                                                             nrow = (num.conds*num.iters))) 
                colnames(sim.structure.repped) <- c("n", "beta" , "eta.x" , "eta.y" , "seed")
                for(i in 1:nrow(sim.structure1)) {
                        sim.structure.repped[ ((num.iters*(i - 1)) + 1): (num.iters*i), (1:num.pars)] <- 
                                purrr::map_dfr(seq_len(num.iters) , ~sim.structure1[i , ])
                }
                sim.structure.repped[ , "seed"] <- rnorm((num.conds*num.iters))
                return(sim.structure.repped)
        }
        sim.repped <- create.sim.repped(sim.structure1 = sim.structure1 , num.pars = 4 , num.iters = 1000)
        which(is.na(sim.repped))                    #check that there are no missing values
        head(sim.structure.repped)                  #verify that everything populations appropriated
        
        reg.gen.full <- function(n , beta , eta.x , eta.y , seed) {      #Actual data generation
                paste("n = " , n , "beta = " , beta , "eta.x = " , eta.x , "eta.y = " , 
                    eta.y , "\n")
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

bayes.data.full <- sim.structure.repped %>%   
        pmap(reg.gen.full)

#save data to computer
saveRDS(data.full , "/Users/Matt Multach/Desktop/XXXXX/bayes_data_110419.RData")



reg.list <- list()                           #create empty list to fill with loop

reg.list.summary <- function(data) {         #function to model data from sub-list elements
       summary(lm(data$Y ~ data$X))
}
summary.list <- lapply(bayes.data.full , reg.list.summary)     #apply above function to each iteration 
                                                        ##of simulated data

names(summary.list[[1]])         ##lots of info in these summaries





for(i in 1:nrow(sim.structure.repped)) {
        cat("i = " , i , "\n")                #tracker for how much of loop is done
        low.list <- reg.gen.simple(n = 25 , beta = 3.7)   #create sub-list of simulated regression 
        ##data
        ##Randomly chose n, beta
        reg.list[[i]] <- low.list             #fill in meta-list with simulation data
        if(i == 1000) {                       #save combiend list of all simulated data if 
                return(reg.list)               ##end of loop; if you change the number of iterations,
        }                                     ##make sure to change if(i == ???)
}














