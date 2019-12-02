#Balanced kfold subsetting
kfold_subsetter <- function(data , k , seed = 7 , list = FALSE) {
       
       if(length(dim(data)) == 2) { ###For 2D data
              #determine number of larger subsets (when unequal subsets)
              nsams.large <- nrow(data) %% k
              
              #determine number of smaller subsets (total number when equal subsets)
              nsams.small <- k - nsams.large
              
              #determine sample size of larger subsets (when unequal subsets)
              samsize.large <- ceiling(nrow(data) / k) * (nsams.large != 0)
              
              #determine sample size of smaller subsets (all subset size when equal subsets)
              samsize.small <- floor(nrow(data) / k)
              
              #indicator for which subset
              subset.indicator <- c(rep((1 : k) , floor(nrow(data) / k)) ,
                                    rep((1 : (nsams.large) ) , (1 * (nsams.large != 0)) ))
              
              #fix random assignment process
              if(seed) {
                     set.seed(seed)
              }
              
              #combine subset indicator with original data  
              newdata <- cbind(data , subset = sample(subset.indicator))
              if(list) {
                     newdata <- return(split(newdata[ , -ncol(newdata)] ,
                                             f = newdata[ , ncol(newdata)]))
              } else {
                     newdata <- return(newdata)
              }
       } else if (length(dim(data)) == 0){   #for 1D data
              #determine number of larger subsets (when unequal subsets)
              nsams.large <- length(data) %% k
              
              #determine number of smaller subsets (total number when equal subsets)
              nsams.small <- k - nsams.large
              
              #determine sample size of larger subsets (when unequal subsets)
              samsize.large <- ceiling(length(data) / k) * (nsams.large != 0)
              
              #determine sample size of smaller subsets (all subset size when equal subsets)
              samsize.small <- floor(length(data) / k)
              
              #indicator for which subset
              subset.indicator <- c(rep((1 : k) , floor(length(data) / k)) ,
                                    rep((1 : (nsams.large) ) , (1 * (nsams.large != 0)) ))
              
              #fix random assignment process
              if(seed) {
                     set.seed(seed)
              }
              
              #combine subset indicator with original data
              #create split list if desired
              newdata <- matrix(cbind(data , 
                                      subset = sample(subset.indicator)) , 
                                ncol = 2)
              if(list) {
                     newdata <- return(split(newdata[ , -ncol(newdata)] ,
                                             f = newdata[ , ncol(newdata)]))
              } else {
                     newdata <- return(newdata)
              }
       }
}





#reg.gen.pres
{
        #libraries
        library(magrittr)
        library(purrr)
        
        #generate initial 45 data conditions
        sim.structure <- as.data.frame(matrix(ncol = 4 , nrow = 45)) 
        {
                colnames(sim.structure) <- c("beta" , "n" , "eta.x" , "eta.y")
                sim.structure[ , "eta.x"] <- c(rep(c(0.0 , 0.1 , 0.2)))
                sim.structure[ , "eta.y"] <- c(rep(0.0 , 3) , rep(0.1 , 3) , rep(0.2 , 3))
                sim.structure[ , "beta"] <- 1
                sim.structure[ , "n"] <- c(rep(25 , 9) , rep(50 , 9) , rep(100 , 9) , 
                                           rep(200 , 9) , rep(1000 , 9))
        }
        
        #generate repped conditions dataframe
        create_sim_repped <- function(sim.structure , num.conds = nrow(sim.structure) , 
                                      num.pars ,  num.iters) {
                sim.structure.repped <- as.data.frame(matrix(ncol = (num.pars + 1) , 
                                                             nrow = (num.conds * num.iters))) 
                colnames(sim.structure.repped) <- c("beta", "n" , "eta.x" , "eta.y" , "seed")
                for(i in 1:nrow(sim.structure)) {
                        sim.structure.repped[((num.iters * (i - 1)) + 1) : (num.iters * i), (1 : num.pars)] <- 
                                purrr::map_dfr(seq_len(num.iters) , ~sim.structure[i , ])
                }
                sim.structure.repped[ , "seed"] <- rnorm((num.conds * num.iters))
                return(sim.structure.repped)
        }
        
        reg_gen_pres <- function(n , beta , eta.x , eta.y , seed) {      #Actual data generation
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
                combine <- list(Y = Y , X = X , err = err , n = n , beta = beta , 
                                eta.x = eta.x , eta.y = eta.y ,
                                seed = seed)        #create combined list of all values
                return(combine)                       #save combined list of all values
        }
}

sim.repped <- create_sim_repped(sim.structure = sim.structure , num.pars = 4 , num.iters = 5)

bayes.data.pres <- sim.repped %>%   
        pmap(reg_gen_pres)






#Note 1: Simple LM vs. CV vs. Bootstrapping
#Note 2: Different freq vs. bayes k



data <- bayes.data.pres[[1]]

X <- data$X
Y <- data$Y
temp <- data.frame(X,Y)
bayes_K = 5
temp <- kfold_subsetter(temp , k = bayes_K)

for(i in 1:bayes_K) {
        
}

train.ind <- temp[ , "subset"]

train <- split(temp[ , -3], f = temp[ , "subset"])


kfold_cv <- function(data) {
       
}


FreqyBayes.CV <- function(data , freq_K , R , bayes_K ,
                         Diffusion , family_brms = c("gaussian","student","skew_normal")) {
       #####frequentist section
       
       ####training data
       #///this needs to be flexible to the data that goes in
       #///currently only does from object with "X" and "Y" columns
       {X <- data$X
       Y <- data$Y}
       temp <- data.frame(X,Y)
       
       temp <- kfold_subsetter(temp , k = bayes_K)
       
       train <- split(temp[ , -3] , f = temp[ , "subset"])
       

       ###Cross-validated regression
       library(caret)
       data_ctrl <- trainControl(method = "cv", number = freq_K)
       model_caret <- train(Y ~ X , data = train , trControl = data_ctrl ,
                            method = "lm" , na.action = na.pass)
       model_cv_final <- model_caret$finalModel
       model_cv_final <- summary(model_cv_final)
       sigma <- model_cv_final$sigma
       model_coeff <- model_cv_final$coefficients
       intercept_b <- model_coeff[1 , 1]
       intercept_SE <- model_coeff[1 , 2]
       X_b <- model_coeff[2 , 1]
       X_SE <- model_coeff[2 , 2]


       
       ###Diffusion parameter
       #\\\Diffusion = 0.1 #this should be changeable in the function (0 to 1)
       Diffusion <- Diffusion*100
       if (Diffusion == 0){
              Diffusion <- 1
       }
       SD_intercept <- intercept_SE*sqrt(Diffusion)
       SD_b <- intercept_b*sqrt(Diffusion)
       sigma <- sigma * sqrt(Diffusion)
       
       ###Bayesian section
       library(brms)
       #\\\family_brms = "gaussian" #this should be changable (gaussian,student,skew_normal)
       
       jj_intercept <- paste("normal(",intercept_b,',',SD_intercept, ")")
       jj_b <- paste("normal(",X_b,',',SD_b, ")")
       jj_sigma <- paste("student_t(", sigma, ',', 0, ',', 5, ")")
       
       
       m1priors <- c(
              prior_string(jj_intercept, class = "Intercept"),
              prior_string(jj_b, class = "b"),
              prior_string(jj_sigma, class = "sigma")
       )
       
       m1 <- brm(
              Y~X,
              data = test,
              prior = m1priors,
              family = family_brms,
              seed = 1,
              iter = 4000,
              chains = 4
       )
       message("the summary of the Baysian model on the test data appears below")
       return(summary(m1)) #return 
       message("Use 'prior_summary(m1)' to extract priors based on the training data")
       
}




#example
test <- MattMohammad(data = bayes.data.full[[1]], method = "cv", K = 10, train_percent = .50, Diffusion = 0.1, family_brms = "gaussian")


