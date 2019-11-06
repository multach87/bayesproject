
reg.gen.simple <- function(n , beta) {
       X <- rnorm(n , mean = 0 , sd = 10)    #generate X values; randomly chose mean, sd
       res <- rnorm(n)                       #generate residuals
       Y <- X*beta + res                     #generate Y values
       combine <- list(Y = Y , X = X , res = res)        #create combined list of all values
       return(combine)                       #save combined list of all values
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













