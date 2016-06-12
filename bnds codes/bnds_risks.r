
  #Load vulnerability analysis results
  #load("./BN/BNDS_simulation.Rdata")

# POPSTERIOR LIKELIHOODS -------------------------------------------------------


  qtiles <- c(0.00001, seq(0.1,0.9,0.1), 0.99999)

  #SAMPLE <- rename(SAMPLE, cvar = var)
  #LHS <- rename(LHS, cvar = var)

  #Data frame to store marginal pmf of each variable
  MarjP <- data_frame(n=1:n, cvar=0, mcc=0, dem=0, sed=0, prc=0, dsc=0)

  CPT  <- list()
  NPT  <- list()
  PARM <- list()

  # CLIMATE VARIABILITY (var) ++++++++++++++++++++++++++++++++++++++++++++++++++
  var <- "cvar"
  NPT[[var]] <- data_frame(s = 1:n, p = 1/max(s))
  MarjP[[var]] <- NPT[[var]]$p[SAMPLE[[var]]]

  # ECONOMIC DEVELOPMENT (dev) +++++++++++++++++++++++++++++++++++++++++++++++++
  # Categorical, root node (low, med, high)
  var <- "dev"
  NPT[[var]] <- data_frame(s = as.numeric(1:3),  p = c(0.4,0.3,0.3))

  # Target population (pop) ++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Continuous, root node, reported in millions
  var <- "pop"
  PARM[[var]] <- data_frame(min = 1.6, max  = 2.6 , mean = 2.1, sd = 0.4)
  NPT[[var]]  <- data_frame(s = qunif(qtiles, PARM[[var]]$min, PARM[[var]]$max))
  NPT[[var]]$p <- dtruncnorm(x = NPT[[var]]$s, a = PARM[[var]]$min,
    b = PARM[[var]]$max, mean = PARM[[var]]$mean, sd = PARM[[var]]$sd)
  NPT[[var]]$p <- NPT[[var]]$p / sum(NPT[[var]]$p)

  # Per capita demand (pcd) ++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Continuous node, conditional on development, reported in l/p/d (pcd)
  var <- "pcd"
  PARM[[var]] <- data_frame(dev = 1:3, mean = c(100, 110, 120)) %>%
    mutate(min = mean - 10, max = mean + 10, sd = 10)
  states <- qunif(qtiles, 100, 120)

  CPT[[var]] <- lapply(1:3, function(x) {
    data_frame(s = states,
      p = dtruncnorm(x= s, a = PARM[[var]]$min[x], b = PARM[[var]]$max[x],
        mean = PARM[[var]]$mean[x], sd = PARM[[var]]$sd[x])) %>%
    mutate(p = p/sum(p))}) %>% bind_rows(.id = "dev") %>%
    mutate(dev = as.integer(dev))

  NPT[[var]] <- CPT[[var]]
  NPT[[var]]$jp <- with(CPT[[var]], p * NPT[["dev"]]$p[dev])
  NPT[[var]] %<>% group_by(s) %>% summarize(mp = sum(jp))

  # Mean climate change (mcc) ++++++++++++++++++++++++++++++++++++++++++++++++++
  #Bivariate, continuous variable
  var <- "mcc"
  gcm_mu <- c(1.023917, 1.200884)
  gcm_sd <- matrix(c(0.02236143, 0.01509313, 0.01509313, 0.60468366), nrow = 2)
  PARM[[var]] <- list(mean = gcm_mu, sd = gcm_sd)
  cc_pair <- as.matrix(SAMPLE[,c("dtemp","dprec")])
  p <- dmvc(x = cc_pair, mu = PARM[[var]]$mean, S = PARM[[var]]$sd)
  MarjP[[var]] <- p / sum(p)

  # Sediment flux (sed) ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Continuous node, conditional on development
  var <- "sed"
  PARM[[var]] <- data_frame(dev = 1:3, mean = c(0.7, 0.75, 0.8)) %>%
    mutate(sd = 0.1)

  CPT[[var]] <- lapply(1:3, function(x) {
    data_frame(s = SAMPLE[[var]],
      p = dnorm(x = s, mean = PARM[[var]]$mean[x], sd = PARM[[var]]$sd[x])) %>%
    mutate(p = p/sum(p))}) %>% bind_rows(.id = "dev")

  NPT[[var]] <- CPT[[var]]
  NPT[[var]]$jp <- CPT[[var]]$p * NPT[["dev"]]$p[as.numeric(CPT[[var]]$dev)]
  NPT[[var]] %<>% group_by(s) %>% summarize(mp = sum(jp))
  MarjP[[var]] <- NPT[[var]]$mp[LHS[[var]]]

  # Price of water (price) +++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Continuous, conditional on development
  var <- "prc"
  PARM[[var]] <- data_frame(dev = 1:3, mean = c(0.7, 1.0, 1.3)) %>%
    mutate(sd = 0.1)

  CPT[[var]] <- lapply(1:3, function(x) {
    data_frame(s = SAMPLE[[var]],
      p = dnorm(x = s, mean = PARM[[var]]$mean[x], sd = PARM[[var]]$sd[x])) %>%
    mutate(p = p/sum(p))}) %>% bind_rows(.id = "dev")

  NPT[[var]] <- CPT[[var]]
  NPT[[var]]$jp <- CPT[[var]]$p * NPT[["dev"]]$p[as.numeric(CPT[[var]]$dev)]
  NPT[[var]] %<>% group_by(s) %>% summarize(mp = sum(jp))
  MarjP[[var]] <- NPT[[var]]$mp[LHS[[var]]]

  # Economic discount rate +++++++++++++++++++++++++++++++++++++++++++++++++++++
  var <- "dsc"
  PARM[[var]] <- data_frame(dev = 1:3, mean = c(0.02, 0.05, 0.08)) %>%
    mutate(sd = 0.005)

  CPT[[var]] <- lapply(1:3, function(x) {
    data_frame(s = SAMPLE[[var]],
      p = dnorm(x = s, mean = PARM[[var]]$mean[x], sd = PARM[[var]]$sd[x])) %>%
    mutate(p = p/sum(p))}) %>% bind_rows(.id = "dev")

  NPT[[var]] <- CPT[[var]]
  NPT[[var]]$jp <- CPT[[var]]$p * NPT[["dev"]]$p[as.numeric(CPT[[var]]$dev)]
  NPT[[var]] %<>% group_by(s) %>% summarize(mp = sum(jp))
  MarjP[[var]] <- NPT[[var]]$mp[LHS[[var]]]

  # Total Demand +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Conditional on pop & pcd
  var <- "dem"

  PARM[[var]] <- expand_grid(pop = NPT[["pop"]]$s, pcd = NPT[["pcd"]]$s) %>%
    mutate(mean = pop * pcd * 0.365, sd = mean/5)

  CPT[[var]] <- expand_grid_df(PARM[[var]], data.frame(value = SAMPLE[[var]]))
  CPT[[var]]$p <- dnorm(CPT[[var]]$value, CPT[[var]]$mean, CPT[[var]]$sd)
  CPT[[var]] %<>% group_by(pop, pcd) %>% mutate(p = p / sum(p))

  NPT[[var]] <- CPT[[var]] %>% ungroup() %>%
    left_join(rename(NPT[["pop"]], p_pop =  p), by = c("pop" = "s")) %>%
    left_join(rename(NPT[["pcd"]], p_pcd =  mp), by = c("pcd" = "s"))
  NPT[[var]]$jp <- with(NPT[[var]], p * p_pop * p_pcd)

  NPT[[var]] <- group_by(NPT[[var]], value) %>% summarize(mp = sum(jp))
  MarjP[[var]] <- NPT[[var]]$mp

  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  MarjP %<>% mutate(P = cvar * mcc * dem * sed * prc * dsc, P = P/sum(P))
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  #write_csv(MarjP, './BN/MarjP_10000.csv')




#-------------------------------------------------------------------------------



  # PRIM ANALYSIS --------------------------------------------------------------

  LHS_df <- OUT %>%
    left_join(mutate(SAMPLE, n = 1:n()), by = "n") %>%
    select(-pop, -pcd, -mcc, -drel, -irel)

  df <- filter(LHS_df, size == 120) %>%
    select(nvar:dprec, NPV)

  library(prim)
  y.thold <- 2200
  test <- prim.box(x=as.matrix(df[,c("dprec","sed")]),
                   y=as.matrix(df[,"NPV"]),
                   threshold=y.thold, threshold.type= 1)

  summary(test)
  df$coloring <- ifelse(df$NPV >= y.thold, T, F)

  ggplot(df, aes(x=dprec, y = sed, color = coloring)) + geom_point()

  df <- results_upper_arun_new[1:3000,]
  y.thold <- 0
  test <- prim.box(peel.alpha=0.05,
    x=as.matrix(df[,2:6]),
                   y=as.matrix(df[,"NPV_BillUSD"]),
                   threshold=y.thold, threshold.type= 1)

  summary(test)
  df$coloring <- ifelse(df$NPV >= y.thold, T, F)


range(df$NPV)




# ANALYZE RESULTS --------------------------------------------------------------







  apply(LHS_df, 2, range)


  NPV_df <- select(LHS_df, size, n, NPV) %>%
    left_join(select(MarjP, n, P), by = "n") %>%
    group_by(n) %>%
    mutate(Regret = max(NPV) - NPV)

  df <- NPV_df %>%
    #filter(size %in% c(120,130,140)) %>%
    group_by(size) %>%
    arrange(Regret) %>%
    mutate(P2 = 1/1000, cumP2 = cumsum(P2)) %>%
    mutate(cumP = cumsum(P)) %>%
    ungroup() %>%
    mutate(size = as.factor(size))

  ggplot(filter(df, n %in% 90:100),
         aes(x = n, y = NPV, group = size, color = size)) +
    theme_bw() +
    geom_line()

  ggplot(df, aes(x = cumP, y = Regret, color = size)) +
    geom_smooth(se = FALSE)



  ggplot(df, aes(x = value, y = P, color = dem)) +
    theme_bw() +
    geom_point(alpha = 0.5) +
    #scale_x_continuous(limits = c(600,1800), breaks = seq(600,1800,300)) +
    #scale_y_log10() +
    labs(x= "NPV", y = "normalized probability") +
    geom_vline(xintercept = 1200, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed")

  perf_dat <- data_frame(value = OUT[,4]) %>%
    mutate(Prior = 1/n, Posterior = MarjP$P) %>%
    gather(key = variable, value = P, Prior:Posterior)

  df <- perf_dat %>%
    group_by(variable) %>%
    arrange(desc(value)) %>%
    mutate(ExP = cumsum(P)) %>%
    mutate(Pnorm = P/(1/n))

  df2 <- filter(df, variable == "Posterior")
  sum(df2$Pnorm)
  range(df2$Pnorm)

  library(scales)


################################################################################

  LHS_df <- OUT %>%
    left_join(mutate(SAMPLE, n = 1:n()), by = "n") %>%
    select(-pop, -pcd, -mcc, -drel, -irel)

  #write_csv(LHS_df, "LHS_output.csv")


  #  rowwise() %>%
  #  mutate(Kindex = which(size == dindex)) %>%
  #  ungroup() %>%



  SAMPLE

  #Vulnerability analysis
  vul_df <- read_csv("./BN/sim_50000.csv") %>%
    select(size = V1, n = V2, rel_dom = V3, rel_irr = V4, Q = V5, NPV = V6) %>%
    mutate(NPV = NPV + designs$PVCost[3] - designs$PVCost[Kindex],
     Q = Q * 12)

  sample_df <- read_csv("./BN/sample_10000.csv") %>% mutate(n = 1:n())
  MarjP_df <- read_csv("./BN/MarjP_10000.csv")

  parcoord_data <- vul_df %>%
    left_join(sample_df, by = "n") %>%
    select(n, size, var, dem, sed, dtemp, dprec, rel_dom, Q, NPV)

  write_csv(parcoord_data, "LHS_output.csv")




  p <- ggplot(df, aes(x = Q, y = NPV)) +
   theme_bw() +
   labs(x = "Q (MCM)", y = "NPV (MUSD)") +
   scale_x_continuous(limits = c(60,180), breaks = seq(60, 180, 20)) +
   geom_point(alpha = 0.5, size = 0.6)

  NPV <- vul_df %>%
    select(K, n, Q, NPV) %>%
    group_by(n) %>%
    mutate(Regret = max(NPV) - NPV) %>%
    ungroup()

  sum(MarjP_df$P)

  NPV_df <- NPV %>%
    left_join(sample_df, by = "n") %>%
    left_join(select(MarjP_df, n, P), by = "n") %>%
    mutate(Pnorm = P/(1/10000), Punif = 1/10000)

  ggplot(NPV_df, aes(x = Pnorm, y = Regret, color = dprec)) +
    theme_bw() +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 50, 100)) +
    geom_point(alpha = 0.4) +
    geom_vline(xintercept = 1) +
    geom_hline(yintercept = 50) +
    labs(x = "likelihood", y = "impact (regret)")

  RegCDF <- NPV_df %>%
    select(K, Regret, P) %>%
    group_by(K) %>%
    arrange(desc(Regret)) %>%
    mutate(CumP = cumsum(P)) %>%
    ungroup()

  ggplot(NPV_df, aes(x = Regret, y = P, group = as.factor(K), color = as.factor(K))) +
    geom_smooth(se = F)



  ggplot(RegCDF, aes(x = CumP, y = Regret, color = as.factor(K))) +
    theme_bw() +
    geom_line() +
    labs(x = "Exceedance P", y = "Regret (M.USD)")

    geom_point(alpha = 0.5) +
    scale_x_continuous(limits = c(600,1800), breaks = seq(600,1800,300)) +
    scale_y_log10() +
    labs(x= "NPV", y = "normalized probability") +
    geom_vline(xintercept = 1200, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed")





