

################################################################################
# PREPARATION OF DATA FOR THE BN MODEL
################################################################################

  # Required input data
  source("./Codes/Mw_input.R")

  cat3 <- data_frame(levels = c(1,2,3), labels = c("Low", "Med", "High"))

################################################################################

#-------------------------------------------------------------------------------
# Population projections -------------------------------------------------------

  # Consider the population projections of the three districts from Tahal 2012
  # Use historical pop information + future projections to define categorical
  # population levels: Low, Medium, High.

  # Historical and projected population (Tahal, 2012) ++++++++++++++++++++++++++
  pop_data <- read_csv("./Inputs/pop_projections.csv")

  # Population lumped
  pop_data_sum <- pop_data %>% group_by(year, scenario) %>%
    summarize(value = sum(value))

  pop_2015 <- pop_data_sum %>% filter(year==2015 & scenario=="Med") %$% value

  # Filter future pop projections & interpolate for every year
  pop_data_yr <- pop_data_sum %>%
    group_by(scenario) %>%
    filter(year > 2014) %>%
    do(Demand = as.data.frame(approx(x= .$year, y = .$value, xout = 2015:2042))) %>%
    unnest(Demand) %>%
    select(scenario, year = x, value = y) %>%
    arrange(scenario, year) %>%
    mutate(value = round(value) %>% as.integer(),
           scenario = factor(scenario, levels = c("Low","Med","High"),ordered = T)) %>%
    arrange(scenario, year)

  pop_change <- data_frame(
      scenario = c("Low","Med","High"),
      pop_2015 = rep(970000,3),
      pop_2034 =  pop_data_yr %>% filter(year == 2034) %$% value) %>%
    mutate(change = round((pop_2034 - pop_2015)/pop_2015,2))

  pop_2034 <- pop_data_yr %>%
    filter(year == 2034) %>%
    mutate(change = (value - pop_2015)/pop_2015)

  #Levels used in the analysis +++++++++++++++++++++++++++++++++++++++++++++++++

  pop_levels  <- data_frame(levels = cat3$levels, labels = cat3$labels,
    '2015' = as.numeric(pop_2015),'2034' = c(pop_2015 * c(1.50, 2.00, 2.50))) %>%
    gather(key = year, value = value, -labels, -levels)

  pop_levels_yr <- pop_levels %>%
    mutate(Year = as.integer(year)) %>%
    group_by(levels) %>%
    do(Demand = as.data.frame(approx(x = .$year, y = .$value, xout = 2015:2070))) %>%
    unnest(Demand) %>%
    mutate(levels = levels %>% as.factor()) %>%
    select(levels, year = x, value = y) %>%
    group_by(levels) %>%
    mutate(value = ifelse(year > 2034, value[year==2034], value))

  df <- pop_levels_yr

  p <- ggplot(df, aes(x=year, y=value/10^3, color = levels)) +
    theme_bw() +
    geom_line(size = 1) +
    labs(x = "", y = "thousand cap")

  #write_csv(x = df, path = "./pop_levels.csv")

# ------------------------------------------------------------------------------
# Demand projections (deterministic) -------------------------------------------

  # Baseline pc demand (2015)
  demand_2015 <- 106.7

  # Set demand level in 2034
  dem_levels  <- data_frame(levels = cat3$levels, labels = cat3$labels,
    '2015' = as.numeric(demand_2015),'2034' = c(demand_2015 * c(0.8,1,1.2))) %>%
    gather(key = year, value = value, -labels, -levels)

  dem_levels_yr <- dem_levels %>%
    mutate(Year = as.integer(year)) %>%
    group_by(levels) %>%
    do(Demand = as.data.frame(approx(x = .$year, y = .$value, xout = 2015:2070))) %>%
    unnest(Demand) %>%
    mutate(levels = levels %>% as.factor()) %>%
    select(levels, year = x, value = y) %>%
    group_by(levels) %>%
    mutate(value = ifelse(year > 2034, value[year==2034], value))

  df <- dem_levels_yr

  p <- ggplot(df, aes(x=year, y=value, color = levels)) +
    theme_bw() +
    geom_line(size = 1) +
    labs(x = "", y = "l per cpd")

  #write_csv(x = df, path = "./dempc_levels.csv")

# ------------------------------------------------------------------------------
 # Sedimentation flux ----------------------------------------------------------

  # Land use classes (sediment rate given m3/km2-year)
  land_use <- data_frame(
    Type  = c("Low_Forest", "Closed_Forest", "Woodland", "Closed_Grassland",
      "Open_Grassland", "Urban", "Cropland"),
    Area  = c(45, 434.25, 402.75, 1055.25, 2.25, 112.50, 225),
    Erosion  = c("Low", "Low", "Low", "Med", "Med", "High", "High"),
    SedRate = c(500, 500, 500, 1000, 1000, 1500, 1500))

  #Sediment flux per land-use class (in MCM/km2-yr)
  land_use %<>%
     mutate(SedFlux = (Area * SedRate) / 10^6)

  # Historical Mean annual sediment flux (in MCM/yr)
  check_dam_reduction <- 0.65
  sedRate_2015 <-  sum(land_use$SedFlux) * (1-0.65)

  # Set demand level in 2034
  sedRate_levels  <- data_frame(
    levels = levels_3lev,
    labels = labels_3lev,
    '2015' = as.numeric(sedRate_2015),
    '2070' = c(sedRate_2015 * c(1.00,1.25,1.50))) %>%
    gather(key = year, value = value, -labels, -levels)

  sedRate_levels_yr <- sedRate_levels %>%
    mutate(Year = as.integer(year)) %>%
    group_by(levels) %>%
    do(value = as.data.frame(approx(x = .$year, y = .$value, xout = 2015:2070))) %>%
    unnest(value) %>%
    mutate(levels = as.factor(levels),
           labels = labels_3lev[levels] %>% as.factor()) %>%
    select(levels, labels, year = x, value = y)

  df <- sedRate_levels_yr

  #p <- ggplot(df, aes(x=year, y=value, color = labels)) +
  #  theme_bw() +
  #  geom_line(size = 1) +
  #  labs(x = "", y = "l per cpd") +
  #  scale_x_continuous(limits = c(2015,2034), breaks = seq(2015,2034,2))
  #write_csv(x = df, path = "sedRate_levels.csv")

  #Historical mean_flow (in km3/year)
  Q_df <- read_csv("./Inputs/flow_gage_cms.csv", skip = 1) %>%
    filter(Year %in% 1977:1989) %>%
    mutate(Q_mcm = Flow_int * 3600 * 24 * 365 / (12 * 10^6)) %>%
    group_by(Year) %>%
    summarize(Q_mcm = sum(Q_mcm))

  #Trapping efficiency (Vorosmorthy et al., 2003)
  Rt <- designs$Volume/mean(Q_df$Q_mcm)
  TE_dam <- 1 - (0.05 / sqrt(Rt))
  sed_flux <- sedRate_2015 * TE_dam

  #Explore sediment flux across a wider range
  span_nClim <- 1:750
  span_K  <- c(80, 100, 120, 130, 140)

  out <- expand_grid(K=span_K, nClim = span_nClim, Q_mean=NA, TE=NA, Rt=NA, Flux=NA)
  out <- filter(out, K == 120)


  #abcd calibration parameters
  par_hydro <- hydro_pars$Princeton

  for (n in 1:nrow(out)) {

    c <- out$nClim[n]
    K <- out$K[n]

    #climate forcings
    clim <- forcing %>% filter((nClim == c) & (Dataset == "Princeton"))

    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3

   out$Q_mean[n] <- median(clim$Q_sim * 12)
   out$Rt[n] <- K/out$Q_mean[n]
   out$TE[n] <- 1 - (0.05 / sqrt(out$Rt[n]))
   out$Flux[n] <- sedRate_2015 * out$TE[n]

  }

  ggplot(out, aes(x = Q_mean, y = Flux, color = as.factor(K))) +
    geom_line()

# ------------------------------------------------------------------------------


################################################################################














################################################################################
# DAG - GRAPHICAL MODE----------------------------------------------------------



  #devtools::install_github('rich-iannone/DiagrammeR')
  library(DiagrammeR)

  #NET PRESENT VALUE
  grViz({"

    digraph DAG {

    # Intialization of graph attributes
    graph[overlap = false, ratio = 1, layout = dot, margin = 0, nodesep = 0.5]

    #Initial node options
    node[fontname = Helvetica,
    shape = plaintext,
    penwidth = 2,
    fontsize = 16,
    height = 1, width = 1.3,
    color = black,
    fixedsize = true]

    m1[label = 'weather\ngenerator', fontcolor = '#E80514']
    m2[label = 'hydrology\nmodel', fontcolor = '#E80514']
    m3[label = 'system\nmodel', fontcolor = '#E80514']

    #Decision nodes
    node[shape = box, peripheries = 1, style = bold];
    d1[label = 'reservoir\ncapacity'];

    #Initial edge attributes
    edge[fontcolor = Red,
    fontname = Helvetica,
    fontsize = 12,
    penwidth = 2,
    color = blue];

    # Chance nodes
    node[shape = ellipse, fillcolor = blue];

    n1[label = 'climate\ndata'];
    n2[label = 'climate\nrealiza-\ntions'];
    n3[label = 'climate\nchange'];
    n4[label = 'develop-\nment'];
    n5[label = 'sediment\nrate'];
    n6[label = 'price\nof water'];

    # Deterministic nodes
    node[shape = ellipse, peripheries = 2, style=filled, fillcolor = Gainsboro];
    o1[label = 'climate\ntraces'];
    o2[label = 'land-use'];
    o3[label = 'sediment\nflux'];
    o4[label = 'stream-\nflow'];
    o5[label = 'yield'];
    o6[label = 'discount\nrate'];
    o7[label = 'project\ncost'];

    #Utility nodes
    node[shape = diamond, height = 1.5, width = 1.5, peripheries = 1, style = bold];
    u1[label = 'NPV'];

    # Edge statements
    n1->m1; n1->m2;
    m1->n2;
    n2->o1;
    n3->o1;
    n4->o2; n4->n6;
    n5->o3;
    n6->u1;
    o1->o3; o1->m2;
    m2->o4;
    o2->n5;
    o3->o5;
    o4->m3;
    m3->o5;
    o5->u1;
    o6->u1;
    o7->u1;
    d1->m3; d1->o7; d1 -> n5;


    }
    "})

  #WATER SECURITY
  grViz({"

    digraph DAG {

    # Intialization of graph attributes
    graph[overlap = false, ratio = 1, layout = dot, margin = 0, nodesep = 0.5]

    #Initial node options
    node[fontname = Helvetica,
    shape = plaintext,
    penwidth = 2,
    fontsize = 16,
    height = 1, width = 1.3,
    color = black,
    fixedsize = true]

    #Initial edge attributes
    edge[fontcolor = Red,
    fontname = Helvetica,
    fontsize = 12,
    penwidth = 2,
    color = blue];

    # Chance nodes
    node [shape = ellipse, fillcolor = blue];
    CDat [label = 'climate\ndata'];
    CRel [label = 'climate\nrealiza-\ntions'];
    CC   [label = 'climate\nchange'];
    Dev  [label = 'develop-\nment'];
    Pop  [label = 'population']
    SedR [label = 'sediment\nrate'];

    # Deterministic nodes
    node [shape = ellipse, peripheries = 2, style=filled, fillcolor = Gainsboro];
    CTr  [label = 'climate\ntraces'];
    Luse [label = 'land-use'];
    SedF [label = 'sediment\nflux'];
    Flow [label = 'stream-\nflow'];
    SedR [label = 'sedimentation\nrate'];
    Inc  [label = 'income\nlevel'];
    PCD  [label = 'per capita\ndemand'];
    Dem  [label = 'domestic\ndemand'];

    #Utility nodes
    node[shape = diamond, height = 1.5, width = 1.5, peripheries = 1, style = bold];
    Rel[label = 'Reliability'];

    #Decision nodes
    node[shape = box, peripheries = 1, style = bold];
    Cap[label = 'reservoir\ncapacity'];

    # Edge statements
    CDat -> CRel;
    CRel -> CTr;
    CC -> CTr;
    CTr -> Flow;
    Dev -> Luse;
    Dev -> Inc;
    Luse -> SedR;
    SedR -> SedF;
    SedF -> ResOp;
    Inc -> PCD;
    PCD -> Dem;
    Pop -> Dem;
    Dem -> ResOp;
    Cap -> ResOp;
    Cap -> SedF;
    Flow -> SedF;
    Flow -> ResOp;
    ResOp-> Rel;

    }
    "})












  #Fit a linear model to annual population increase (%)
  #pop_incrate <- pop_yr %>%
  #  filter(year %in% c(2015,2016)) %>%
  #  group_by(scenario) %>%
  #  mutate(change = (value - value[which(year == 2015)])/ value[which(year == 2015)])

  #Rates 0.04, 0.05, 0.07

  #Apply coefficients from a truncated norm distribution (for the future replace this with a beta distribution)
  #library(truncnorm)
  #samples <- rtruncnorm(50, a = 0.03, b = 0.08, mean = 0.05, sd = 0.01)

  #Find demand time-series
  #years <- 2015:2040
  #pop_baseline <- 970538
  #sim_period <- 1:length(years)
  #demand <- sapply(samples, function(x) pop_baseline * (1 + x)^sim_period)
  #demand_df <- data.frame(Year = years, demand, check.names = FALSE) %>%
  #  as_data_frame() %>%
  #  gather(key = sample, value = value, -Year)

  #ggplot(demand_df, aes(x = Year, y = value, color = sample)) + geom_line()


# ------------------------------------------------------------------------------




