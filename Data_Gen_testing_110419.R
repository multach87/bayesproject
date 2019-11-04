#reg.gen.complex
{
        #libraries
        library(magrittr)
        library(purrr)
        
        #generate initial 108 data conditions
        sim.structure <- as.data.frame(matrix(ncol = 4 , nrow = 135)) 
        {
                colnames(sim.structure) <- c("n" , "beta" , "eta.x" , "eta.y")
                sim.structure[ , "eta.x"] <- c(rep(c(0.0 , 0.1 , 0.2)))
                sim.structure[ , "eta.y"] <- c(rep(0.0 , 3) , rep(0.1 , 3) , rep(0.2 , 3))
                sim.structure[ , "beta"] <- c(rep(1 , 9) , rep(3 , 9) , rep(10 , 9))
                sim.structure[ , "n"] <- c(rep(25 , 27) , rep(50 , 27) , rep(100 , 27) , 
                                            rep(200 , 27) , rep(1000 , 27))
        }
        
        #generate repped conditions dataframe
        create_sim_repped <- function(sim.structure , num.conds = nrow(sim.structure) , 
                                      num.pars ,  num.iters) {
                sim.structure.repped <- as.data.frame(matrix(ncol = (num.pars + 1) , 
                                                             nrow = (num.conds * num.iters))) 
                colnames(sim.structure.repped) <- c("n", "beta" , "eta.x" , "eta.y" , "seed")
                for(i in 1:nrow(sim.structure)) {
                        sim.structure.repped[ ((num.iters * (i - 1)) + 1): (num.iters * i), (1:num.pars)] <- 
                                purrr::map_dfr(seq_len(num.iters) , ~sim.structure[i , ])
                }
                sim.structure.repped[ , "seed"] <- rnorm((num.conds * num.iters))
                return(sim.structure.repped)
        }
        sim.repped <- create_sim_repped(sim.structure = sim.structure , num.pars = 4 , num.iters = 10)
        
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
}

bayes.data.full <- sim.repped %>%   
        pmap(reg_gen_full)


#one round
View(sim.repped)
n <- sim.repped[11 , "n"]
beta <- sim.repped[11 , "beta"]
eta.x <- sim.repped[11 , "eta.x"]
eta.y <- sim.repped[11 , "eta.y"]
seed <- sim.repped[11 , "seed"]

reg_gen_full(n = n , beta = beta , eta.x = eta.x , eta.y = eta.y , seed = seed)



#save data to computer
saveRDS(bayes.data.full , "/Users/Matt Multach/Desktop/XXXXX/bayes_data_110419.RData")



reg.list <- list()                           #create empty list to fill with loop

reg_list_summary <- function(data) {         #function to model data from sub-list elements
       summary(lm(data$Y ~ data$X))
}
summary.list <- lapply(bayes.data.full , reg_list_summary)     #apply above function to each iteration 
                                                        ##of simulated data

names(summary.list[[1]])         ##lots of info in these summaries




















