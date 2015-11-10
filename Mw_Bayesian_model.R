
#------------------------------------------------------------------------------#
# Title: Bayesian Belief Network for the Mwache project design
# By:    M.Umit Taner
# Date:  September, 2015
#------------------------------------------------------------------------------#

  #install.packages("devtools")
  #library(devtools)
  #install_github("ecbrown/LaplacesDemon")

##############################################################################
#CLIMATE PROBABILITIES


  source("Mw_Input_pars.R")
  script_dir <- "/Users/umit/Dropbox/Research/Scripts/"
  source(paste0(script_dir, "BAYESNET.R"))


gcm_dat <- gcms_data %>%
  as_data_frame() %>%
  filter(Ensemble == "CMIP5") %>%
  select(Index, Temp_Abs, Prec_Abs) %>%
  mutate(Delta_Temp = Temp_Abs - hist_clim$Temp,
         Delta_Prec = Prec_Abs / hist_clim$Prec)

#Fit gcm data to a distribution, get mu and sigma
gcm_mu <- as.vector(apply(gcm_dat[,c("Delta_Temp","Delta_Prec")], 2, mean))
gcm_covv <- cov(cbind(gcm_dat$Delta_Temp, gcm_dat$Delta_Prec)) * 20

#climate probability matrix, ordered from driest
clim_priors <- as_data_frame(clim_matrix) %>%
  mutate(c = Clim, Delta_Prec = Prec/100, Delta_Temp = Tavg) %>%
  select(c, Delta_Prec, Delta_Temp, NVar) %>%
  mutate(P = dmvc(x = cbind(.$Delta_Temp, .$Delta_Prec), mu=gcm_mu, S=gcm_covv)) %>%
  mutate(P = P / length(unique(.$NVar))) %>%
  mutate(P = P / sum(P)) %>%
  arrange(Delta_Prec, desc(Delta_Temp)) %>%
  mutate(Climate = 1:n()) %>%
  select(Climate, c, Delta_Prec, Delta_Temp, NVar, P)

clim_prior_tbl <- clim_priors %>%
  select(Climate, Delta_Prec, Delta_Temp, P) %>%
  group_by(Delta_Prec, Delta_Temp) %>%
  summarize(P = sum(P))

p <- ggplot(clim_prior_tbl) +
  geom_tile(aes(x = Delta_Temp, y = Delta_Prec*100, fill = P), color = "black") +
  scale_x_continuous(expand=c(0,0), breaks = seq(0,2.5, 0.5)) +
  scale_y_continuous(expand=c(0,0), breaks = seq(70, 140, 10)) +
  scale_fill_gradient(
    low = "white", high = "black", name = "Probability", guide = "legend",
    limits = c(0, 0.06), breaks = seq(0, 0.06, 0.01)) +
  labs(x=alabels[1], y=alabels[2]) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_blank(),
    plot.title        = element_text(face="bold"),
    legend.background = element_rect(fill=NA),
    panel.border      = element_rect(color="black", fill=NA),
    strip.background  = element_rect(fill="white", color='black'),
    panel.margin      = unit(1, "lines"),
    plot.title = element_text(face="bold"))

ggsave(file = "climate_probs.png", plot = p, height =5, width = 6)

##############################################################################



# DATA -------------------------------------------------------------------------

  clim_response_dat <-  read.csv('./inputs/safeyield_Sep25.csv')
  bn <- read_excel(path="./inputs/MwBN_data.xlsx", sheet = "Nodes", skip =1)
  edges_dat <- read_excel(path="./inputs/MwBN_data.xlsx", sheet = "Edges", skip =1)
  signs_dat <- read_excel(path="./inputs/MwBN_data.xlsx", sheet = "Signs", skip =1)

  nodes <- bn$Variable

  edges <- setNames(sapply(1:length(nodes),
    function(x) nodes[which(edges_dat[x,-1] != 0)]), nodes)

  weights <- setNames(sapply(1:length(bn$Variable), function(x)
    unlist((edges_dat[x, which(edges_dat[x,] != 0)])[-1], use.names = FALSE)),
    bn$Variable)

  signs <- setNames(sapply(1:length(bn$Variable), function(x)
    unlist((signs_dat[x, which(signs_dat[x,] != 0)])[-1], use.names = FALSE)),
    bn$Variable)

  node_priors <- list()
  node_priors$Equake   <- c(0.001, 0.999)
  node_priors$Climate  <- clim_priors$P
  node_priors$Discount <- c(0.236, 0.378, 0.271, 0.106, 0.009)
  node_priors$Econdev  <- c(0.106, 0.247, 0.294, 0.247, 0.106)
  node_priors$Capacity <- c(0.166, 0.167, 0.167, 0.167, 0.167, 0.166)

  ############### NO-climate change!!!!
  #bn$min_bound[2] <- 1
  #bn$max_bound[2] <- 10
  #bn$state_num[2] <- 10
  #node_priors$Climate  <- rep(0.1, 10)
  #####################################

  #sapply(node_priors, sum)
  #sapply(node_priors, length)

  ##############################################################################

  sample_size <- 10000
  #variance_coef = 1

  #Arguments for BN sampling
  LHS_scaled <- list()
  LHS_actual <- list()

  #Probabilities are conditional here....
  BN_pars <- NPT_GENERATE(nodes = bn$Variable,
    min_bounds = bn$min_bound, max_bounds = bn$max_bound,
    state_nums = bn$state_num, node_priors  = node_priors,
    node_links = edges, edge_weights = weights,
    edge_signs = signs, variance_coef= 0.2)

  #sapply(BN_pars$NPT, function(x) sum(x[["P"]]))
  #sapply(node_priors, sum)
  #sapply(node_priors, length)

  ptm <- proc.time()

  #Apply sampling (rescaled interval 0-1)
  for (i in 1:6) {

    #Sampled nodes (scaled to 0-1 interval)
    LHS_scaled[[i]] <- BAYESNET_LHS(
      sample_size = sample_size,
      decision    = list(Capacity = i),
      nodes       = nodes,
      CPTs        = BN_pars$NPT,
      edges       = edges)

    LHS_index <- LHS_scaled[[i]][[2]]

    #Retrieve actual levels from the scaled intervals
    LHS_actual[[i]] <- sapply(1:length(BN_pars$Levels), function(x)
       BN_pars$Levels[[x]][LHS_index[[x]]]) %>%
       as.data.frame() %>% as_data_frame()

    colnames(LHS_actual[[i]]) <- nodes

  }

  proc.time() - ptm


################################################################################
################################################################################
################################################################################


#POSTERIOR DISTRIBUTIONS

  for (d in 1:6) {

    sample_dat <- LHS_actual[[d]] %>% gather(key = node, value = value)

    p <- list()

    for (i in 1:length(nodes)) {

      df <- filter(sample_dat, node == nodes[i])
      df$value <- as.factor(df$value)

      p[[i]] <-  ggplot(df, aes(x=value, fill = value)) +
        #theme_bw(base_size = 10) +
        ggtitle(paste0("(",i,") ", bn$Definition[i])) +
        geom_bar(aes(y = (..count..)/sum(..count..)),
                       fill ="gray55", color = "black", width = 0.8) +
        #scale_y_continuous(expand=c(0,0)) +
        scale_x_discrete(drop=FALSE) +
        labs(x = "levels", y = "frequency") +
        theme(plot.title = element_text(face="bold")) +
        guides(fill = FALSE)

    }

  fig_name  <- paste0("./posteriors/posteriors_design_",d,".png")
  fig_title <- paste0("Storage capacity:", designs$Labels[d])

  out <- arrangeGrob(p[[6]], p[[7]], p[[8]], p[[9]], ncol=2,
    top=textGrob(fig_title, gp=gpar(fontsize=20,font=3), just = "centre"))

  ggsave(filename = fig_name, out, width = 10, height = 10)

  }


################################################################################
################################################################################
################################################################################


  # PERFORMANCE METRIC PLOTS ---------------------------------------------------

  clim_subset <- clim_matrix %>% as_data_frame() %>%
    mutate(c = Clim) %>% select(-Clim)

  syield <- read.csv('./inputs/safeyield_Sep25.csv') %>%
    as_data_frame() %>%
    left_join(clim_priors, by = "c") %>%
    select(Climate, d, yield)
    #left_join(clim_subset, by = "c") %>%
    #filter(CC == 19) %>%
    #select(d, Climate = NVar, yield)


  NPV_tbl <- bind_rows(LHS_actual, .id ="d") %>%
    mutate(d = as.integer(d), Climate = as.integer(Climate)) %>%
    left_join(syield, by = c("d","Climate")) %>%
    mutate(AnnB = yield * Valuewater) %>%
    rowwise() %>%
    mutate(PVB  = PV_CALCULATE(rep(AnnB, Lifespan), r=Discount),
           NPV  = PVB - designs$PVCost[d]) %>%
    ungroup()


  #NPV_tbl <- vuly_df


  #Joint probability distribution of NPV
  NPV_posterior_cdf <- NPV_tbl %>%
    select(d, NPV) %>%
    mutate(d = ordered(.$d, labels = designs$Labels)) %>%
    group_by(d) %>%
    arrange(desc(NPV)) %>%
    mutate(Index = 1:n(), ExP = Index/(n()+1))

  options(scipen=3)

  ### Density function
  p1 <- ggplot(NPV_posterior_cdf, aes(x=NPV, group=d, color=d)) +
    theme_bw() +
    stat_density(position = "identity", fill=NA, trim = TRUE, size = 1) +
    scale_color_brewer(name="Designs Capacities", palette="Spectral") +
    labs(x="NPV (Million dollars)", y = "frequency"); p1

  p2 <- ggplot(NPV_posterior_cdf, aes(x=ExP, y=NPV, group=d, color=d)) +
    theme_bw() +
    geom_line(size=1) +
    scale_color_brewer(name="Designs Capacities", palette="Spectral") +
    scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.1), expand=c(0,0)) +
    geom_vline(xintercept=0.90, linetype='dashed') +
    labs(x="Probability of exceedance", y = "NPV (M$)") ; p2

  #ggsave(file="NPV_posterior_pdf.png", plot = p1, width=25, height=10, unit='cm')

  #save.image("Results_10000_nocc_posteriors.Rdata")

################################################################################
################################################################################
################################################################################



  # DECISION CRITERIA **********************************************************
  ##############################################################################

  #cVaR calculations ********************
  alpha_values <- c(0, 0.50, 0.75, 0.90, 0.95, 0.99)

  metrics <- data_frame(Metric = alpha_values,
    "40 MCM" = NA,  "60 MCM" = NA, "80 MCM" = NA,
    "100 MCM" = NA, "120 MCM" = NA, "140 MCM" = NA)

  for (i in 1:nrow(designs)) {

    data <- filter(NPV_posterior_cdf, d == designs$Labels[i])
    ind <- sapply(alpha_values, function(x) which.min(abs(data$ExP - x)):nrow(data))
    metrics[,i+1] <- sapply(ind, function(x) sum(data[x, "NPV"]) / length(x))
  }

  metrics_df <- metrics %>% gather(key = Designs, value = value, -Metric)

  metric_tbl <- metrics %>% slice(c(1,2,4))

  #write.csv(metric_tbl, file = "metric_10000_cc_uniform.csv")

  ggplot(metrics_df, aes(x=Metric, y=cVaR, color=Designs)) +
    theme_bw() +
    geom_line(size=1)

  filter(metrics_df, Metric %in% c(0,0.90))



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

  vuly_df <- expand.grid(Climate  = BN_pars$Levels$Climate,
    Discount = BN_pars$Levels$Discount,
    Valuewater = BN_pars$Levels$Valuewater,
    Lifespan  = BN_pars$Levels$Lifespan,
    d = 1:6) %>%
    as.data.frame() %>% as_data_frame() %>%
    left_join(syield, by = c("Climate","d")) %>%
    mutate(AnnB = yield * Valuewater) %>%
    rowwise() %>%
    mutate(PVB  = PV_CALCULATE(rep(AnnB, Lifespan), r=Discount),
           NPV  = PVB - designs$PVCost[d]) %>%
    ungroup()




################################################################################
################################################################################
################################################################################


#  save.image(file = "Sept29.Rdata")


  st_num  <- 6
  shape_1 <- 1
  shape_2 <- 1

  df <- data_frame(x = seq(1/(2*st_num), by = 1/(st_num), length.out = st_num)) %>%
    mutate(y = dbeta(x, shape_1, shape_2), y = round(y/sum(y),3))

  ggplot(df, aes(x=as.factor(x),y)) +
    theme_bw() +
    geom_bar(stat="identity")

