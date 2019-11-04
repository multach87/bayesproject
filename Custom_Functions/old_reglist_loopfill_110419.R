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