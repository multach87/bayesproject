#FINAL
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
                newdata <- cbind(data , sample(subset.indicator))
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
                                        sample(subset.indicator)) , 
                                  ncol = 2)
                if(list) {
                        newdata <- return(split(newdata[ , -ncol(newdata)] ,
                                                f = newdata[ , ncol(newdata)]))
                } else {
                        newdata <- return(newdata)
                }
        }
}
