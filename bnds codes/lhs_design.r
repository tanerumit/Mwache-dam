

# LHS Design -------------------------------------------------------------------

  #SAMPLE SIZE
  n <- 1000

  # Required libraries
  library(foreach, quietly = T, warn.conflicts = F)
  library(doParallel, quietly = T, warn.conflicts = F)

  ### INPUT PARAMETERS AND INTERVALS
  par <- list()
  par$nvar$intv  <- 1:n              # natural variability
  par$mcc$intv   <- 1:n               # climate change
  par$dem$intv   <- c(55,110)
  par$sed$intv   <- c(0.5, 1.0)
  par$dprc$intv  <- c(0.7, 1.3)
  par$iprc$intv  <- c(40, 80)
  par$disc$intv  <- c(0.02, 0.08)

  ### Define LHS matrix #(n rows x m columns) , uniform priors...
  seed <- 17
  m <- length(par)
  set.seed(seed)
  LHS <- sapply(1:m, function(x) sample(1:n))
  LHS %<>% as.data.frame() %>% as_data_frame()
  colnames(LHS) <- names(par)

# SAMPLE VALUES  ---------------------------------------------------------------

  #Data frame to store sample values
  SAMPLE <- data_frame(n = 1:n)

# (1) Natural variability samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  c <- "nvar"
  qtiles <- cumsum(rep(1/max(par[[c]]$intv),max(par[[c]]$intv)))
  sintval <- c(0, qtiles) * n + 0.1
  SAMPLE[[c]]  <- findInterval(LHS[[c]], sintval, rightmost.closed = T)

  # (2) Mean climate changes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate mean precip & temp change grid
  DelP_range <- c(0.4, 1.8)
  DelT_range <- c(0,4)
  DelP <- seq(DelP_range[1], DelP_range[2], length.out = 50)
  DelT <- seq(DelT_range[1], DelT_range[2], length.out = 20)
  mcc_grid <- expand.grid(Prec=DelP,Tavg=DelT) %>% mutate(Index = 1:n)
  t_step <- (DelT[2] - DelT[1])/2
  p_step <- (DelP[2] - DelP[1])/2
  c <- "mcc"
  out <- sapply(1:n, function(x) {
    sub <- filter(mcc_grid, Index == LHS[[c]][x])
    #set.seed(seed)
    t <- runif(1, min = sub$Tavg - t_step, max = sub$Tavg + t_step)
    #set.seed(seed)
    p <-runif(1, min = sub$Prec - p_step,  max = sub$Prec + p_step)
    c(t,p,as.integer(sub$Index))
  }) %>% t()

  SAMPLE$mcc <- out[,3]
  SAMPLE$temp <- out[,1]
  SAMPLE$prec <- out[,2]

  # (3) Demand samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Continuous variable, known interval
  c <- "dem"
  levels <- seq(min(par[[c]]$intv), max(par[[c]]$intv), length = n)
  center <- levels[LHS[[c]]]
  step <-  (levels[2] - levels[1])/2
  SAMPLE[[c]] <- sapply(1:n, function(x) {set.seed(seed)
      runif(1, min = center[x] - step, max = center[x] + step)})

  # (4) Sediment rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Continuous variable known interval
  c <- "sed"
  levels <- seq(min(par[[c]]$intv), max(par[[c]]$intv), length = n)
  center <- levels[LHS[[c]]]
  step <-  (levels[2] - levels[1])/2
  SAMPLE[[c]] <- sapply(1:n, {
    set.seed(seed)
    function(x) runif(1, min = center[x] - step, max = center[x] + step)})

  # (5) Change factor for unit price of waterv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # As dimensionless number, continuous variable
  c <- "dprc"
  levels <- seq(min(par[[c]]$intv), max(par[[c]]$intv), length = n)
  center <- levels[LHS[[c]]]
  step <-  (levels[2] - levels[1])/2
  SAMPLE[[c]] <- sapply(1:n, {set.seed(seed); function(x)
    runif(1, min = center[x] - step, max = center[x] + step)})

  # (5) Change factor for unit price of waterv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # As dimensionless number, continuous variable
  c <- "iprc"
  levels <- seq(min(par[[c]]$intv), max(par[[c]]$intv), length = n)
  center <- levels[LHS[[c]]]
  step <-  (levels[2] - levels[1])/2
  SAMPLE[[c]] <- sapply(1:n, {set.seed(seed); function(x)
    runif(1, min = center[x] - step, max = center[x] + step)})

  # (6) Economic discount rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # As percantage, continuous variable
  c <- "disc"
  levels <- seq(min(par[[c]]$intv), max(par[[c]]$intv), length = n)
  center <- levels[LHS[[c]]]
  step <-  (levels[2] - levels[1])/2
  SAMPLE[[c]] <- sapply(1:n, {set.seed(seed); function(x)
    runif(1, min = center[x] - step, max = center[x] + step)})

#-------------------------------------------------------------------------------

