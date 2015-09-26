
#------------------------------------------------------------------------------#
# Title: Bayesian Belief (Decision) Network for Mwache Dam design
# By:    M.Umit Taner
# Date:  September, 2015                                                       
#------------------------------------------------------------------------------#

  #install.packages("devtools")
  #library(devtools)
  #install_github("ecbrown/LaplacesDemon")


  #Load associated data
  load("Inputs/Mb_StressTest.Rdata")
  load("Inputs/Mb_Climate.Rdata")
  source("/Users/umit/GDrive/Research/Scripts/SAMPLING.R")
  source("/Users/umit/GDrive/Research/Scripts/MULTIPLOT.R")
  require(devtools)
  require(LaplacesDemon)
  require(mvtnorm)
  require(truncnorm)
  require(lazyeval)
  require(readxl)

  
  
# BAYESIAN BELIEF NETWORK ******************************************************

# Step 1. |Define nodes, edges, and significance weights
# Step 2. |Define prior probabilities of the nodes
# Step 3. |Populate CPTs by weighted sums algrorithm
# Step 4. |Apply LHS sampling to approximate PMFs

  
################################################################################  
  
  
# STEP 1) NODES - EDGES - RANKING ----------------------------------------------

  
  #*****************************************************************************
  
  #Define nodes in the system
  nodes  <- c(
    "Equakes",  "ClimChange", "ClimVar",  "DiscRate",
    "Develop",  "Population", "Capacity", "Sediment",
    "PCDemand", "PriceWater", "LifeSpan")
  
  #*****************************************************************************
  
  #Define parent-child relationships between the nodes (NA means no parent)
  edges <- list(
    Equakes    = NA,
    ClimChange = NA,
    ClimVar    = NA,
    DiscRate   = NA,
    Develop    = NA,
    Population = NA,
    Capacity   = NA,
    Sediment   = c("ClimChange", "ClimVar"),
    PCDemand   = c("ClimChange", "ClimVar", "Develop", "Capacity"),
    PriceWater = c("Develop", "PCDemand", "Capacity"),
    LifeSpan   = c("Equakes", "Sediment", "Capacity"))

  #*****************************************************************************

  #Define the strengths of the parent-child relationshps
  edge_weights <- list(
    Equakes    = NA,
    ClimChange = NA,
    ClimVar    = NA,
    DiscRate   = NA,
    Develop    = NA,
    Population = NA,
    Capacity   = NA,
    Sediment   = c(2, 3),
    PCDemand   = c(1, 3, 2, 5),
    PriceWater = c(2, 4, 2),
    LifeSpan   = c(1, 3, 5))
  
  #*****************************************************************************
  
  

  
################################################################################
  
  
# STEP 2) PRIOR PROBABILITIES  -------------------------------------------------

  # First define a prior probability table for each node. The tables
  # should include Vlevel, Scale, Prior columns
  # All variable levels are scaled to the interval [0,1]

  ClimGCMs <- filter(gcms_data, Ensemble == "CMIP5") %>% 
    select(x=Temp_Abs, y= Prec_Abs)
  ClimGCMs_mean <- as.vector(apply(ClimGCMs, 2, mean))
  ClimGCMs_sdev <- cov(cbind(ClimGCMs$x, ClimGCMs$y))
  Prior = dmvc(x  = as.matrix(select(cc_space, Temp, Prec)), 
               mu = ClimGCMs_mean, S= ClimGCMs_sdev)
  Prior = Prior/sum(Prior)
  
  #Store node level information and prior probabilities
  node_priors <- list()

  #Node 1, Earthquake (catastrophic events) ************************************
  node_priors[["Equakes"]] <- data_frame(Level = c("True", "False")) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.001, 0.999))

  #Node 2, Mean Climate change (DeltaP, DeltaT change) *************************
  node_priors[["ClimChange"]] <- data_frame(
      Level = paste0("cc",1:48))%>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = Prior)

  #Node 3, Natural Climate variability *****************************************
  node_priors[["ClimVar"]] <- data_frame(
      Level = paste0("tr",1:10)) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = rep(1/n(), n()))

  #Node 4, economic discount rate **********************************************
  node_priors[["DiscRate"]] <- data_frame(
      Level = c("Low", "Med", "High"),
      Value = c(0.03, 0.05, 0.08)) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.35, 0.4, 0.25))

  #Node 5, future economic development *****************************************
  node_priors[["Develop"]] <- data_frame(
      Level = c("Low", "Med", "High"),
      Value = c("Low", "Med", "High")) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.3, 0.4, 0.3))

  #Node 6, future population ***************************************************
  node_priors[["Population"]] <- data_frame(
      Level = c("Low", "Med", "High"),
      Value = c(3e6, 3.5e6, 4e6)) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.3, 0.4, 0.3))

  #Node 7, Sedimentation and siltation *****************************************
  node_priors[["Sediment"]] <- data_frame(
    Level = c("Low", "Med", "High"),
    Value = c("NA","NA","NA")) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.5, 0.4, 0.1))

  #Node 8, Per Capita Water Demand *********************************************
  node_priors[["PCDemand"]] <- data_frame(
      Level = c("Low", "Med", "High"),
      Value = c("NA","NA","NA")) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.3, 0.5, 0.2))

  #Node 9, Per Capita Water Demand *********************************************
  node_priors[["PriceWater"]] <- data_frame(
      Level = c("Low", "Med", "High"),
      Value = c("NA","NA","NA")) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.3, 0.4, 0.3))

  #Node 10, Life span of the Reservoir *****************************************
  node_priors[["LifeSpan"]] <- data_frame(
      Level = c("Low", "Med", "High"),
      Value = c("NA","NA","NA")) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = c(0.6, 0.3, 0.1))

  #Node 11, Reservoir capacity (decision node) *********************************
  node_priors[["Capacity"]] <- data_frame(
      Level = seq(40,140,20),
      Value = c("NA","NA","NA","NA","NA","NA")) %>%
    mutate(Scale = seq(1/(2*n()), by = 1/(n()), length.out = n()),
           Prior = rep(1/length(Scale), length(Scale)))

  
################################################################################

  
# STEP 3) CONDITIONAL PROBABILITY TABLES (CPTs) --------------------------------
  mw_CPTs <- list()

  # Create CPT tables for each node in the system
  # Start with the first node, propogate CPT tables
  for (n in 1:length(nodes)) {

    #Current node
    node_name  <-  nodes[n]
    node_level <- node_priors[[node_name]]$Level
    node_scale <- node_priors[[node_name]]$Scale
    node_parents <- edges[[n]]

    #If node does not have any parents, use the prior table
    if(is.na(node_parents)[1]) {

      mw_CPTs[[node_name]] <- node_priors[[node_name]] %>%
        select(Level = Scale, Prob = Prior)

    #If node has parents, propogate the CPTs
    } else  {

      #Parent nodes
      parent_level <- lapply(node_parents, function(x) node_priors[[x]]$Level)
      parent_scale <- lapply(node_parents, function(x) node_priors[[x]]$Scale)
      parent_wghts <- edge_weights[[n]]

      #Template CPT table
      CPT_level <- expand.grid(parent_level)
      CPT_scale <- expand.grid(parent_scale)
      colnames(CPT_level) <- node_parents
      colnames(CPT_scale) <- node_parents

      CPT_template <- merge(CPT_scale, node_scale, by=NULL) %>%
        as_data_frame() %>% mutate(Prob = NA)
      colnames(CPT_template) <- c(node_parents, "Level", "Prob")

      #Loop through each level combination to propogate probabilities
      for (i in 1:nrow(CPT_template)) {

        output_vals <- CPT_template[[i,"Level"]]
        parent_vals <- as.vector(unlist(CPT_template[i,node_parents]))

        CPT_template[i,"Prob"] <- dtruncnorm(
          output_vals,
          mean = sum(parent_vals * parent_wghts)/sum(parent_wghts),
          sd = 1/sum(parent_wghts), a = 0, b = 1)

      }

    #Normalize CPT based on their conditionals
    dots <- lapply(node_parents, as.symbol)
    CPT_template2 <- suppressMessages(CPT_template %>%
      group_by_(.dots= dots) %>%
      summarize(normf = sum(Prob)) %>%
      left_join(CPT_template, .) %>%
      mutate(Prob = Prob/normf)) %>%
      select(-normf)

    #Final CPT for the node
    mw_CPTs[[node_name]] <- CPT_template2

    }
  }

  
################################################################################

  
# STEP 4) LATIN HYPERCUBE SAMPLING (LHS) ---------------------------------------

  edges <- edges
  CPTs  <- mw_CPTs
  ssize <- 200

  #Arguments for BN sampling
  LHS_scale <- list()
  LHS_org <- list()

  ptm <- proc.time()
  #Apply sampling (rescaled interval 0-1)
  for (i in 1:nrow(designs)) {

    #Sampled nodes (scaled to 0-1 interval)
    LHS_scale[[i]] <- BN_SAMPLING(sample_size = ssize,
      decision = list(Capacity = mw_CPTs$Capacity$Level[i]),
      nodes = nodes,
      CPTs = mw_CPTs,
      edges = edges)
  }
  proc.time() - ptm

  ptm <- proc.time()
  #Sampled nodes (Original internval)
  for (i in 1:nrow(designs)) {
    LHS_org[[i]] <- sapply(nodes, function(y)
      sapply(1:nrow(LHS_scale[[i]]), function(x)
        node_priors[[y]][[match(LHS_scale[[i]][[x,y]], node_priors[[y]]$Scale),"Level"]])) %>%
      as.data.frame() %>% as_data_frame()
  }
  proc.time() - ptm


  
  
  
  
    
  
################################################################################    
################################################################################    
################################################################################    
################################################################################  
################################################################################  
  
#STAKEHOLDER INFORMATION 

  #Plots to be saved    
  require(xtable)
  
  #Belief data 
  belief_df <- read_excel("inputs/Workshop_survey_v2.xlsx", sheet = 1) %>%
    gather(key = Particip, value = value, x1:x28) %>%
    mutate(value = as.factor(.$value))
  
  belief_sum <- belief_df %>%
    mutate(Value2 = as.numeric(value)) %>%
    group_by(Uncertainty, Short) %>%
    summarize(Mean = round(mean(Value2)))

  
  #Belief data plots
  p <- list()
  for (i in 1:12) {
    
    plot_name <- unique(belief_df$Uncertainty)[i]
    plot_data <- filter(belief_df, Uncertainty == plot_name)
    
    p[[i]] <- ggplot(plot_data, aes(x = value, fill = Short)) +
      ggtitle(plot_name) +
      geom_bar(width=.7) +
      facet_grid(.~ Short) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      labs(x="", y = "") +
      theme(plot.title = element_text(face="bold"),
            plot.margin = unit(c(0,0,0,0),"mm")) +
      guides(fill = FALSE)  

    ggsave(paste0("./Results_in_progress/factor_",i,".png"), p[[i]], width = 8, height = 2.5, units = "in")
    
  }
  
  
  
  
  
  
  
################################################################################
# #Profiling
# Rprof("BN_sample.out")
# 
# ptm1 <- proc.time()
# LHS1<- BN_SAMPLING(sample_size = 500,
#     decision = list(Capacity = mw_CPTs$Capacity$Level[i]),
#     nodes = nodes,
#     CPTs = mw_CPTs,
#     edges = edges)
# time1 <- proc.time() - ptm1
# 
# ptm2 <- proc.time()
# LHS_samples <- BN_SAMPLING2(sample_size = 2000,
#     decision = list(Capacity = mw_CPTs$Capacity$Level[i]),
#     nodes = nodes,
#     CPTs = mw_CPTs,
#     edges = edges)
# time2 <- proc.time() - ptm2

