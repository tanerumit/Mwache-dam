

# DATA  ------------------------------------------------------------------------

  #MODEL DATA
  source("./Codes/Mw_input.R")

  #LIBRARIES
  library(foreach, quietly = T, warn.conflicts = F)
  library(doParallel, quietly = T, warn.conflicts = F)

  ##### CONSTANTS
  Res_Kd <- 20
  par_calib <- hydro_pars[["Princeton"]]
  demand_bline <- 40
  price_dom_bline <- 0.7 # USD/m3
  price_irr_bline <- 40  # USD/ha (4000 KSH/ha)

  #Climate realizations considered
  clim_realizations <- read_csv(paste0("./BN/climvar_Princeton.csv")) %>%
    rename(nvar = nVar) %>% filter(nvar %in% 1:n)

  # ITERATE THROUGH SAMPLE -----------------------------------------------------

  Date = as.Date("2020-01-01") + months(0:599)
  sim_yrs <- length(unique(year(Date)))
  iterate <- expand.grid(K = c(3,4,5), n = 1:n)
  end <- nrow(iterate)

  #Parellell running
  cl <- makeCluster(6)
  registerDoParallel(cl)

  #SIMULTION (LOOP)
  OUT <- foreach(s=1:end, .combine = rbind,
    .packages = c("lubridate", "sirad")) %dopar% {

    #Current design and iteration
    i <- iterate$n[s]
    K <- designs$Volume[iterate$K[s]]

    #Current parameter values
    cli_var    <- SAMPLE[['nvar']][i]   # climate variability (ts index)
    del_temp   <- SAMPLE[['temp']][i]   # delta temperature (delta degC)
    del_prec   <- SAMPLE[['prec']][i]   # delpta precip (delta %)
    dem_dom    <- SAMPLE[['dem']][i]    # domestic demand (MCM/yr)
    sed_flux   <- SAMPLE[['sed']][i]    # sediment flux (MCM/yr)
    price_dom  <- SAMPLE[['dprc']][i]   # price of domestic use (USD/m3)
    price_irr  <- SAMPLE[['iprc']][i]   # price of irrigation (USD/ha
    disc_rate  <- SAMPLE[['disc']][i]   # econ discount rate (%)
    PVCost <- designs$PVCost[iterate$K[s]]

    # Climate forcing & inflow
    clim <- clim_realizations %>% filter(nvar == cli_var) %>%
      mutate(Prec = Prec * rep(seq(1, del_prec, length.out = sim_yrs), each = 12),
        Tavg = Tavg + rep(seq(0, del_temp, length.out = sim_yrs), each = 12),
        Tmin = Tmin + rep(seq(0, del_temp, length.out = sim_yrs), each = 12),
        Tmax = Tmax + rep(seq(0, del_temp, length.out = sim_yrs), each = 12),
        PE = HARGREAVES(Date, Tavg, Tmax - Tmin, lat = 3.5),
        Q = ABCD_QEST(par_calib, Prec, PE, S_ini=100, G_ini=2) %>%
      as.numeric() * area * 1e-3)

    ####### DOMESTIC DEMAND
    demand_fn <- approx(x=c(2020,2034), y=c(demand_bline,dem_dom), xout=2020:2069)
    demand_yr <- data_frame(year = demand_fn$x, value = demand_fn$y) %>%
      mutate(value = ifelse(year > 2034, value[year==2034], value))
    demand_dom <- demand_yr %>%
      expand_grid_df(., data.frame(month = 1:12)) %>%
      mutate(value = value * month_coef[month]) %>%
      select(year, month, value) %>%
      arrange(year, month) %$% value

    ###### IRRIGATION DEMAND
    #Calculate ET0 (Penmann - Monteith equation (FAO))
    irrig_cal <- clim %>%
      left_join(crop_pars_et0, by = "Month") %>%
      mutate(
        EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
        Wind  = Wind * 0.01157407,
        v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
        vap_p = Humid/100 * v_sat,
        ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
          tal = cst(Solar, Date, radians(Latitude)),
          z = 150, uz = Wind, days  = Date, lat = Latitude),
        ET0 = ET0 * as.numeric(days_in_month(Month))) %>%
      select(Year, Month, Prec, Tavg, Tmin, PE, EPrec, ET0)

    # #Calculate CRW (mm/month)
    demand_irr <- irrig_cal %>%
      left_join(crop_Kc, by = "Month") %>%
      left_join(crop_area, by = c("Season","Crop")) %>%
      mutate(
        Area_ha = ifelse((Year > 2033) & Crop == "Mango", 0, Area_ha)) %>%
      mutate(
        ETc = ET0 * Kc,
        NIR_mm = pmax(ETc - EPrec, 0),
        GIR_mm = NIR_mm / (Eff/100),
        GIR_mcm = GIR_mm * Area_ha / 10^5) %>%
      mutate(Month = as.integer(Month)) %>%
      group_by(Year, Month) %>%
      summarize(Irr_mcm = sum(GIR_mcm)) %$% Irr_mcm * 0.8

    #Sediment calculations
    #Trapping efficiency (Vorosmorthy et al., 2003)
    Rt <- K/(mean(clim$Q) * 12)
    TE <- 1 - (0.05 / sqrt(Rt))
    sed_trapped <- sed_flux * TE
    sedF <- c(0,rep(sed_trapped, sim_yrs-1))
    Kd_yr <- {sedF %>% cumsum() + Res_Kd} %>% rep(each = 12)

    # Simulate reservoir operations
    sim <- RESERVOIR_SIM(
       beg.y    = 2020,
       beg.m    = 1,
       Q        = clim$Q,
       K_dead   = Kd_yr,
       K        = K,
       K_cons   = K,
       K_buff   = K*0.9,
       buffer   = 1,
       Tar_dom  = demand_dom,
       Tar_irr  = demand_irr,
       Tar_eco  = demand_ces_mcm$env,
       evap.m   = (clim$PE - clim$Prec) * 1e-3,
       f_elev   = f_elev,
       f_vol    = f_vol,
       f_sarea  = f_sarea)

    # Time-based reliabilities
    nofail   <- length(which(sim$R_dom >= 1 * sim$Tar_dom))
    dom_rel  <- 100 * nofail / nrow(sim)
    nofail   <- length(which(sim$R_irr >= 1 * sim$Tar_irr))
    irr_rel  <- 100 * nofail / nrow(sim)

    #Net present value
    CBA <- sim %>% group_by(Year) %>%
      summarize(Dom = sum(R_dom), IrrR = sum(R_irr)/sum(Tar_irr)) %>%
      mutate(Area_ha = ifelse(Year < 2034, 2400, 1700),
        BDom = Dom * price_dom,
        BIrr = price_irr * Area_ha * IrrR * 10^-6,
        OMC  = PVCost * 0.005,
        NB   = BDom + BIrr - OMC)

     NPV <- PV_CALCULATE(CBA$NB, r = disc_rate) - PVCost
     Qmean <- mean(clim$Q)

     c(K,i, dom_rel, irr_rel, Qmean, NPV)
  }
  stopCluster(cl)

  OUT <- as_data_frame(as.data.frame(OUT)) %>%
    select(size = V1, n = V2, drel = V3, irel = V4, Q = V5, NPV = V6)

  save.image("./BN/simulation.Rdata")
  #write_csv(as_data_frame(as.data.frame(OUT)), "./BN/sim_50000.csv")
  save(OUT, file = "./BN/vulnerability.Rdata")
  #write_csv(SAMPLE, "./BN/sample_10000.csv")

################################################################################



