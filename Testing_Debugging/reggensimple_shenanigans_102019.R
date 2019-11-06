reg.list <- list()
for(i in 1:1000) {
       cat("i = " , i , "\n")
       #name <- paste(c("sim." , as.character(i)) , collapse = "")
       #cat("name = " , name , "\n")
       low.list <- reg.gen.simple(n = 25 , beta = 3.7)
       reg.list[[i]] <- low.list
       #cat("names(reg.list) = " , names(reg.list) , "\n")
       #names(reg.list[i]) <- name
       if(i == 1000) {
              return(reg.list)
       }
}

summary(lapply(reg.list , lm , Y~X))
summary(lm(reg.list[[1]]$Y ~ reg.list[[1]]$X))




names(reg.list[1]) <- name
names(reg.list[1])


high.list <- rep(low.list , 1000)



test.list <- rep(reg.gen.simple(n = 25 , beta = 2.5) , 1000)








test <- reg.gen.simple(n = 25 , beta = 2)
test

test.2 <- reg.gen.simple(n = 25 , beta = 2)
test.2

meta.test <- list(sim1 = test , sim2 = test.2)




beta <- 3.5
X <- rnorm(10 , mean = 0 , sd = 10)
Y <- X*beta + rnorm(1)
X
Y

fit.test <- lm(Y ~ X)
summary(fit.test)
