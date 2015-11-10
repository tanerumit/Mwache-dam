
#------------------------------------------------------------------------------#
#   Title   :    Decision Tree Project - MWACHE HYDROLOGY/SYSTEMS MODEL
#   By      :    Mehmet Umit Taner
#   Date    :    September 25, 2015
#   Version :    1.0
#   Notes   :    This version is only designed to implement baseline
#                simulation analysis
#
#------------------------------------------------------------------------------#

  source("Mw_Input_pars.R")
  require(readr)

  #Monthly precip forcings data global met datasets
  p_forcings <- read_csv("./inputs/grid_monthly_prec_long_mm.csv")

  #Abcd model paramters for each global met dataset
  cal_pars   <- read_csv("./Inputs/grid_met_calibration_pars.csv")

# HYDROLOGY SIMULATION ---------------------------------------------------------

  #Select met. dataset to use
  met_data <- "Princeton"

  #subset related parameters
  forcings_historic <- filter(forcings_gcm_monthly, ID == 1)
  forcings_metdata  <- filter(p_forcings, Source == met_data)

  prec_per <- forcings_metdata %>% select(value) %>% rowSums()
  date_per <- forcings_metdata$Date
  pars_per <- cal_pars[[met_data]]
  forcings_per <- forcings_historic[date_hist %in% date_per,]

  #Calculate PET using HARGREAVES (From Princeton data)
  Tdif <- forcings_per$Tmax - forcings_per$Tmin
  PE_per <- HARGREAVES(date_per, forcings_per$Tavg, Tdif, BasinLat = 3.5)

  #Vectors to store the simulated streamflow sequences
  Q_sim_mm  <- vector(mode = "numeric", length = length(PE_per))
  Q_sim_mcm <- vector(mode = "numeric", length = length(PE_per))

  Q_sim_mm <- ABCD_QEST_NO_SNOW( parm = pars_per, P  = prec_per, PE = PE_per,
    S_ini = 1, G_ini = 1)

  Q_sim_mcm <- Q_sim_mm * (catchment_area * 1e-3)

################################################################################
################################################################################
################################################################################

# YIELD-RELIABILITY RELATINSHIP ------------------------------------------------

  # Explore yield vs reliability relationship at different reservoir capacities
  # Plot & compare yield-reliabilty curves

  #Simulate reservoir performance over the given period
  yieldT <- c(seq(40, 80, 10), seq(85, 150, 5))
  envir_mon <- envir_monthly$demand
  inflows   <- Q_sim_mcm

  # Create the input-output matrix for the analysis.
  # Matrix size = # designs x # yield levels
  loopArray <- expand.grid(d = designs$Volume, y = yieldT) %>%
    mutate(iter = 1:n(), output = rep(list(NA), n()), rel = NA) %>%
    select(iter, d, y, rel, output)

  #Progress bar & timing
  pb <- txtProgressBar(min = 1, max = nrow(loopArray), style = 3)

  #Loop through each scenario
  for (i in 1:nrow(loopArray)) {

    setTxtProgressBar(pb, i)

    #Set current values for iteration
    capacity <- loopArray[[i,"d"]]
    yield_m  <- loopArray[[i,"y"]] * as.numeric(days_in_month(1:12)/365)

    loopArray[[i,"output"]] <- RESERVOIR_SIM(beg.y  = 1950, beg.m = 1,
      K = capacity, Q = inflows, K_d = 20,
      ratingc = rating_curve, evap.m = netevap,
      Tar_dom  = yield_m, Tar_irr  = rep(0,12), Tar_eco  = envir_mon)

    T_Dom <- loopArray[[i,"output"]]$Tar_dom
    R_Dom <- loopArray[[i,"output"]]$R_dom
    loopArray[[i,"rel"]] <- 100 * length(which((R_Dom - T_Dom) == 0))/ length(T_Dom)

    close(pb)
  }

  loopArray %<>% as.tbl() %>% mutate(d = as.factor(d))

  #Plot yield-reliability curve for all designs
  p1 <- ggplot(loopArray, aes(x = y, y = rel, group = d, color = d)) +
    theme_bw() +
    geom_line(size = 1) +
    geom_vline(xintercept=80.3, show.legend = F, linetype="dashed", size=1) +
    scale_y_continuous(limits = c(60,100), breaks = seq(60,100,10)) +
    scale_x_continuous(limits = c(50,120), breaks = seq(50,120,10)) +
    labs(x = "Yield (MCM/year)",  y= "Reliability (%)") +
    scale_color_brewer(name="Design\nCapacity\n(MCM)", palette="Spectral")

  # Compare yields for the de, in which the likelihood estimates of the
  # plausible future outcomes and their impacts are incomplete, sign options
  rel_levs <- c(0.85, 0.90, 0.95, 0.99)
  yield_df <- data_frame(Designs = designs$Volume,
    Yield_85 = NA, Yield_90 = NA, Yield_95 = NA, Yield_99 = NA)

  #Find yields at %95 reliability for each design
  for (i in 1:nrow(designs)) {

    df <- filter(loopArray, d == designs$Volume[i])
    yield_df[i,2:5] <- approx(df$rel, df$y, xout = rel_levs)$y

  }

  #Yields at given reliabities
  df <- gather(yield_df, key = Yield_level, value = value, -Designs) %>%
    mutate(Designs     = factor(Designs, levels = designs$Volume),
           Yield_level = factor(Yield_level, levels = colnames(yield_df)[-1]))

  p2 <- ggplot(df, aes(x=Designs, y = value, fill = Yield_level)) +
    theme_bw() +
    geom_bar(stat="identity", position = "dodge", width =0.8) +
    facet_wrap( ~ Yield_level, ncol = 1) +
    scale_y_continuous(limits= c(0, 80), breaks = seq(0,80,10)) +
    guides(fill = FALSE)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


# CROP CALCULATIONS ------------------------------------------------------------

  library(sirad)

  #Arguments needed to calculated ET0 (Based on PM- FAO 56)
  #  Tmax      = Vector of length n containing daily maximum temperature [C].
  #  Tmin      = Vector of length n containing daily minumum temperature [C].
  #  z         = Altitude above the sea level [m].
  #  lat       = Latitude in decimal degrees. Required only if extraT=NA.
  #  sol_rad   = Vector of length n of daily solar radiation [MJm-2d-1]
  #  uz        = Wind speed measured at heith ’meah’ [ms-1].
  #  meah      = The height (a.s.l) of the wind speed measurement [m].
  #  days      = Vector of class ’Date’ of length n.
  #  vap_pres  = Vector of length n of mean daily vapour pressure [kPa].
  #  tal       = Clear sky transmissivity [0-1].

  spr_eff <- 0.65     # efficiency for the sprinkler system
  dip_eff <- 0.90     # efficiency for the dip system

  #Crop coefficient values
  crop_data <- read_csv("./inputs/crop_data_irrig.csv", skip = 1) %>%
    gather(key= Month, value = Kc, Jun:May) %>%
    unite(Crop, Season, c) %>%
    mutate(Crop = as.factor(Crop))

  #Planned crop areas and irrigaiton efficiencies
  crop_area <- crop_data %>%
    select(Month, Crop, Area_ha, Method) %>%
    filter(Month == "Jun") %>% select(-Month) %>%
    mutate(Eff = ifelse(Method == "sprinkler", spr_eff, dip_eff))

  #Monthly parameters needed for the Penmann-Montheith calculation
  crop_pars <- read_csv("./inputs/crop_data_et0.csv", skip = 1) %>%
    select(-Unit) %>%
    gather(key   = Month, value, -Parameter) %>%
    rowwise() %>%
    mutate(Month = which(Month == month.abb)) %>%
    ungroup() %>%
    spread(key = Parameter, value = value) %>%
    select(Month, Humid = Humidity, Solar = Solar_rad, Suns = Sunshine, Wind)

  #Calculate ET0 and append to the met forcing table
  forcings <- forcings_per %>%
    mutate(Prec_e = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 -10))) %>%
    left_join(crop_pars, by = "Month") %>%
    mutate(
      Wind  = Wind * 0.01157407,
      v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
      vap_p = Humid/100 * v_sat,
      ET0   = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
      tal   = cst(Solar, date_per, radians(Latitude)),z = 150, uz = Wind,
      days  = date_per, lat = Latitude) * as.numeric(days_in_month(Month)))

  #Calculate CRW (mm/month)
  ETc <- crop_data %>%
    select(Crop, Month, Kc) %>%
    spread(key = Crop, value = Kc) %>%
    mutate(Month = c(6:12,1:5)) %>%
    left_join(select(forcings, Month), ., by = "Month") %>%
    select(-Month) %>%
    mutate_each(funs(. * forcings$ET0))

  #Net Irrigation and Gross Irrigation Requirements (mm/month)
  crop_irr <- ETc %>%
    mutate_each(funs(pmax(0,. - forcings$Prec_e))) %>%
    mutate(Index = 1:n()) %>%
    gather(key=Crop, value = NIR_mm, -Index) %>%
    left_join(crop_area, by = "Crop") %>%
    mutate(GIR_mm = NIR_mm / Eff) %>%
    mutate(GIR_mcm = GIR_mm * Area_ha / 10^5)

  #Monthly irrigation demand
  irr_monthly <- crop_irr %>%
    group_by(Index) %>%
    summarize(Irr_mcm = sum(GIR_mcm)) %>%
    bind_cols(select(forcings, Year, Month),.) %>%
    select(-Index)

  #Annual irrigation demand
  irr_annual <- irr_monthly %>%
    group_by(Year) %>%
    summarize(Irr_mcm = sum(Irr_mcm))

  #Plot monthly & annual irrigation demand
  p <- ggplot(data = irr_annual, mapping = aes(x = Year, y = Irr_mcm)) +
    geom_line() +
    geom_hline(aes(yintercept = 34), color = "red"); p


# RESERVOIR SIMULATION ---------------------------------------------------------

  inflows  <- Q_sim_mcm
  capacity <- 120
  domes_mon <- 70 * as.numeric(days_in_month(1:12) / 365)
  irrig_mon <- irr_monthly[["Irr_mcm"]]
  envir_mon <- envir_monthly$demand

  out <- RESERVOIR_SIM(beg.y  = 1950, beg.m = 1,
    K = 120, Q = inflows, K_d = 20, ratingc = rating_curve, evap.m = netevap,
    Tar_dom  = domes_mon, Tar_irr  = irrig_mon, Tar_eco  = envir_mon)


# RESULTS VISUALIZATION (SINGLE RUN) -------------------------------------------

  results <- results_sim %>%
    unite("Date", c(Year, Month) ,sep = "-") %>%
    mutate(Date = as.Date(paste0(Date,"-1"))) %>%
    mutate(Deficit_eco = Tar_eco - R_eco,
           Deficit_dom = Tar_dom - R_dom,
           Deficit_irr = Tar_irr - R_irr)

  p_storage <- ggplot(results, aes(x = Date, y = S)) +
    geom_line() + labs(x = "Date", y = "Storage (MCM)")

  #Deficit time-series
  df <- select(results, Date, Deficit_eco : Deficit_irr) %>%
    gather(key = Deficits, value = value, -Date)

  p_deficit <- ggplot(df, aes(x = Date)) +
    geom_bar(aes(y = value, fill = Deficits), stat = "identity") +
    #geom_line(aes(y = S), data = results) +
    labs(x = "Date", y = "MCM")


################################################################################







# Mass balance analysis --------------------------------------------------------

  #Mass-balance accounting for the projected demand & supplies
  year_set <- 2015:2065
  file_path <- "./inputs/BWSS_balance.xlsx"

  demand_dat <- read_excel(path = file_path, sheet = "Demand", skip = 1)

  demand_tbl <- data_frame(Year = year_set) %>%
    mutate(Pop = approxExtrap(
      x=demand_dat$Year, demand_dat$Population, xout = year_set)$y) %>%
    mutate(Demand =  approxExtrap(
      x=demand_dat$Year, demand_dat$Demand_mcm, xout = year_set)$y)

  supply_tbl <- read_excel(path = file_path, sheet = "Supply", skip = 1) %>%
    mutate_each(funs(. * (365/1e6)), -Year) %>%
    select(Year, Tiwi, Marere, Others, Mzima, Baricho, Msambweni, Mwache) %>%
    gather(key = Source, value = Yield, -Year)

  var_cols <- c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd",
                "#969696","#737373","#253494")

  #Plot cummulative yield from alternative sources
  p3 <- ggplot(supply_tbl, aes (x = Year)) +
    theme_bw() +
    geom_area(mapping=aes(y = Yield, fill = Source)) +
    scale_fill_brewer(palette = "Blues") +
    scale_x_continuous(expand = c(0,0),
      limits = c(2015, 2065), breaks = seq(2015,2065,10)) +
    scale_y_continuous(expand = c(0,0),
      limits = c(0,350), breaks = seq(0, 350, 50)) +
    geom_line(data = demand_tbl, mapping = aes(y=Demand),
              linetype = "longdash", size = 1) +
    labs( x = "Year", y = "Cummulative yield (MCM)")

  #Surplus storage from Mwache
  balance_tbl <- supply_tbl %>%
    spread(key = Source, value = Yield) %>%
    mutate(Other_supply = Tiwi + Marere + Mzima + Baricho + Others + Msambweni,
           Demand       = demand_tbl$Demand,
           Net          = Other_supply + Mwache - Demand,
           Balance      = as.factor(ifelse(Net > 0, 1, 0))) %>%
    select(Year, Demand, Other_supply, Mwache, Net, Balance)

  net.rle = rle(balance_tbl$Net < 0)
  balance_tbl$group = rep.int(1:length(net.rle$lengths), times=net.rle$lengths)

  p4 <- ggplot(
      data = balance_tbl,
      mapping = aes(x = Year, y = Net, group = group, fill = Balance)) +
    geom_area() +
    scale_fill_manual(values  = c('red', 'blue'), labels = c("Def.", "Surp")) +
    scale_y_continuous(limits = c(-150,50), breaks = seq(-150,50,50)) +
    labs(x = "Year", y = "Surplus supply (MCM)") +
    guides(fill = FALSE)




