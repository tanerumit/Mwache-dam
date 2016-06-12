

#SAMPLE FROM THE NETWORK
BN_SAMPLE <- function(data, n, seed = FALSE) {

  #require(iterators)
  #require(doParallel)
  #require(foreach)

  m <- length(data)

  if(seed) set.seed(seed)
  LHS_mat <- sapply(1:m, function(x) sample(1:n)) %>% as.data.frame()
  colnames(LHS_mat) <- names(data)

  #Data frame to store samples
  LHS_set <- replace(x=LHS_mat, values=NA)

  #cl <- makeCluster(6)
  #registerDoParallel(cl)
  #Loop through each instatiation
  #OUT <- foreach(k=1:n, .combine = rbind) %dopar% {

  for(k in 1:n) {

    #Loop through each node
    for(s in 1:m) {

      #Current values in the iteration
      node  <- names(data)[s]
      state <- data[[node]]$states
      CPT_i <- data[[node]]$CPT

      # Set the CPT for root and non-root nodes
      if(length(names(CPT_i)) == 2) {

        CPT <- CPT_i

      } else {

        #parents <- c("n1", "n4")
        parents <- setdiff(names(CPT_i), c("s","p"))
        parent_vals <- sapply(parents, function(x) unlist(LHS_set[[x]][k]))

        dots <- sapply(1:length(parents),
                       function(i) lazyeval::interp( ~ variable == value,
                                                     variable = as.name(parents[i]), value = parent_vals[i]))

        CPT <- filter_(CPT_i, .dots = dots)

      }

      # Sample from the interval
      sintval <- c(0, cumsum(CPT$p)) * n + 0.1
      sindex  <- findInterval(LHS_mat[[node]][k], sintval, rightmost.closed = T)
      LHS_set[[node]][k] <- data[[node]]$state[sindex]

    }

    #as.vector(LHS_set[k,])
  }

  return(LHS_set %>% as_data_frame())
  #stopCluster(cl)
  #return(OUT %>% as_data_frame())

}
