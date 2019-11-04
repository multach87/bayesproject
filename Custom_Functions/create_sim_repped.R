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