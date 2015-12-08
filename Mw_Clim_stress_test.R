
#------------------------------------------------------------------------------#
# Title: Climate stress test of the alternative Mwache system designs
# By:    M.Umit Taner
# Date:  August, 2015
#------------------------------------------------------------------------------#

  # Hydrologic simulations are based on the abcd model
  # Water system operations are based on WEAP

  source("Mw_Input_pars.R")

  #Read-in previously generated climate forcings
  forcings_st <- read_csv("./Inputs/forcings_st.csv", skip =1) %>% split(.$ID)

  #DATA structure
  results_tbl <- expand.grid(c=1:length(forcings_st), d=1:nrow(designs)) %>%
    as.tbl() %>%
    mutate(iter = 1:n(), rel = NA) %>%
    mutate(output = rep(list(NA), n())) %>%
    select(iter, d, c, rel, output)

  #Progress bar & timing
  pb <- txtProgressBar(min = 1, max = nrow(results_tbl), style = 3)

  end <- 100
  ptm <- proc.time()
  for (s in 1:end) {

    #Set climate (c) and design (d) counters
    c <- results_tbl[s,"c"] %>% as.numeric()
    d <- results_tbl[s,"d"] %>% as.numeric()

    #Calculate PET using HARGREAVES
    tdif <- forcings_st[[c]]$Tmax- forcings_st[[c]]$Tmin
    PE <- HARGREAVES(Dates = date_gcm,
      tavg=forcings_st[[c]]$Tavg, tdif = tdif, BasinLat = 3.5)

    #Simulate streamflow (mm)
    Q_sim_mm <- ABCD_QEST_NO_SNOW(
      par_calib,
      P = forcings_st[[c]]$Prec,
      PE = PE,
      S_ini = 1,
      G_ini = 1)

    #Simulate streamflow (MCM/mo)
    Q_sim <- Q_sim_mm * area * 1e-3

    OUT <- BIN_SEARCH(FUN = RESERVOIR_RELIABILITY, par = "Tar_dom",
      sign = "+", lower = 0, upper = 50, target = 95,
      beg.y = 2000, beg.m = 1, K = designs$Volume[d], K_d = 5, Q = Q_sim,
      Tar_irr = rep(0,12), Tar_eco = demand_ces_mcm$env,
      ratingc = rating_curve, evap.m  = netevap)


    #Find yield for the given reliability and inflow sequence
    result <- YIELD_ESTIMATE(
      dates         = sim_dates,
      Q             = Q_sim,
      capacity      = designs$Volume[d],
      top_conserv   = designs$Volume[d],
      top_buffer    = 35,
      top_inactive  = 5,
      buffer_coef   = 1,
      hVol_curve    = e_approx,
      Volh_curve    = s_approx,
      evap_m        = netevap,
      S_initial     = 1,
      reliability   = 0.95,
      max_iter      = 100,
      double_cycle  = TRUE)

    #Save results
    results_df[s, 1] <- c
    results_df[s, 2] <- d
    results_df[s, 3:5] <- result$Val

    setTxtProgressBar(pb, s)


  }






    close(pb)
  proc.time() - ptm

  head(results_df)

  results_out <- data.frame(results_df) %>%
    as_data_frame() %>%
    select(c = X1, d = X2, iter = X3, rel = X4, yield = X5)

  write.csv(x=results_out, file = "safeyield_Sep25.csv")





