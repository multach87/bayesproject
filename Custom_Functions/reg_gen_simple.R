reg_gen_simple <- function(n , beta) {
       X <- rnorm(n , mean = 0 , sd = 5)    #generate X values; randomly chose mean, sd
       res <- rnorm(n)                       #generate residuals
       Y <- X*beta + res                     #generate Y values
       combine <- list(Y = Y , X = X , res = res)        #create combined list of all values
       return(combine)                       #save combined list of all values
}