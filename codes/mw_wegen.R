

## *****************************************************************************
## Mwache Weather Generator
## November, 2014
## M.Umit Taner
## *****************************************************************************
################################################################################

  scripts_dir <- "/Users/umit/Dropbox/Research/Scripts/"

  source("./Codes/Mw_input.R")
  source(paste0(scripts_dir,"PRIMARY.R"))
  source(paste0(scripts_dir,"PLOTTING.R"))
  library(car)

  dataset_list <- c("Princeton","GPCC","CRU")
  historic_period <- as.Date("1961-01-01") + c(0:467) * months(1)

  future_period <- as.Date("2020-01-1") + c(0:599) * months(1)
  future_years <- unique(year(future_period))
  forcing_gcm <- read_csv("./inputs/forcing_gcm.csv", skip = 1)

# ------------------------------------------------------------------------------
# WAVELET (AR) MODEL  ----------------------------------------------------------

  for (s in 2:length(dataset_list)) {

    dataset <- dataset_list[s]

    #Prepare annual and monthly climate data ***********************************
    prec_grid <- read_csv("./inputs/precip_grid_mm.csv")
    prec_data <- filter(prec_grid, Source == dataset) %>% .[["value"]]

    #Monthly precip with selected dataset
    clim_monthly <- filter(forcing_gcm, Index == 1) %>%
      mutate(Prec = prec_data) %>%
      filter(Year %in% c(1961:1999))

    #Annual precipitation (1950-1999)
    clim_annual  <- group_by(clim_monthly, Year) %>%
      summarise(Prec=sum(Prec),Tavg=mean(Tavg),Tmax=mean(Tmax),Tmin=mean(Tmin))
    ANNUAL_PRCP  <- clim_annual$Prec

    #### WAVELET DECOMPOSITION +++++++++++++++++++++++++++++++++++++++++++++++++

    #Apply box-cox transformation of the original data
    #lambda <- as.numeric(powerTransform(ANNUAL_PRCP)$lambda)
    #TRANSFORM_DATA <- bcPower(ANNUAL_PRCP,lambda)

    #Wavelet analysis on transformed data
    CLIMATE_VARIABLE <-  ANNUAL_PRCP
    source(paste0(scripts_dir, "/WAVELET_ANALYSIS.R"))
    POWER_SPECTRUM_PRCP_OBS <- GWS
    POWER_SPECTRUM_PRCP_SIGNIFICANCE <- signif_GWS
    POWER_SPECTRUM_PRCP_PERIOD <- period

    #Wavelet spectra
    df <- data_frame(Period = period, Sig85 = signif_GWS, GWS)
    #write_csv(x = df, path = paste0("gws_obs_",dataset,".csv"))

    #WARM Calibration parameters
    #Number of orthogonal component series that representing a low-freq signal
    NUM_FINAL_PERIODS <-  1
    #Length of the first component: This is an initial guess
    NUM_PERIODS_COMP1 <-  3
    #Length of all components
    NUM_PERIODS_ALL_COMPS <- c(NUM_PERIODS_COMP1) #, NUM_PERIODS_COMP2)
    #List indices of significant components (from all components)
    ALL_SIG_PERIODS <- c(5,6,7)

    #Wavelet decomposition
    source(paste0(scripts_dir,'WAVELET_DECOMP.R'))
    PRCP_LOW_FREQUENCY_COMPONENTS <- LOW_FREQUENCY_COMPONENTS
    PRCP_NOISE <- NOISE

    #GENERATE STOCHASTIC PRECIP DATA +++++++++++++++++++++++++++++++++++++++++++

    sim_length <- 50
    sim_num    <- 100000

    #sim_numate y times
    PRCP_FINAL_ANNUAL_SIM             <- array(NA,c(sim_length,sim_num))
    PRCP_FINAL_ANNUAL_SIM_TRANSFORMED <- array(NA,c(sim_length,sim_num))
    POWER_SPECTRUM_PRCP_ARIMA_SIM     <- array(NA,c(length(GWS)+1,sim_num))

    pb <- txtProgressBar(min = 1, max = sim_num, style = 3)

    ptm <- proc.time()
    for (y in 1:sim_num) {

      PRCP_NOISE_MEAN <- mean(PRCP_NOISE)
      PRCP_NOISE_CENTERED  <- PRCP_NOISE - PRCP_NOISE_MEAN
      MODEL1 <- auto.arima(PRCP_NOISE_CENTERED,max.p=5,max.q=0,max.P=0,max.Q=0,stationary=TRUE)
      AR1 <- as.vector(MODEL1$coef)[which(names(MODEL1$coef)!="intercept")]
      INTERCEPT <- 0
      if (length(which(names(MODEL1$coef)=="intercept"))>0) {
          INTERCEPT <- as.vector(MODEL1$coef)[which(names(MODEL1$coef)=="intercept")]}

      #Simulate noise component
      SIM1 <- arima.sim(n = (sim_length+1),
        list(ar = AR1),sd = sd(residuals(MODEL1)))[2:(sim_length+1)] + INTERCEPT +
        PRCP_NOISE_MEAN

      #Simulate low frequency component(s)
      SIM2 <- array(NA,c(sim_length,NUM_FINAL_PERIODS))
      for (k in 1:NUM_FINAL_PERIODS) {

        PRCP_LOW_FREQUENCY_COMPONENTS_MEAN <- mean(PRCP_LOW_FREQUENCY_COMPONENTS[,k])
        PRCP_LOW_FREQUENCY_COMPONENTS_CENTERED  <- PRCP_LOW_FREQUENCY_COMPONENTS[,k] -
          PRCP_LOW_FREQUENCY_COMPONENTS_MEAN

        MODEL2 <- auto.arima(PRCP_LOW_FREQUENCY_COMPONENTS_CENTERED,max.p=5,max.q=0,max.P=0,max.Q=0,stationary=TRUE)
        AR2 <- as.vector(MODEL2$coef)[which(names(MODEL2$coef)!="intercept")]
        INTERCEPT <- 0
        if (length(which(names(MODEL2$coef)=="intercept"))>0)
          {INTERCEPT <- as.vector(MODEL2$coef)[which(names(MODEL2$coef)=="intercept")]}
        SIM2[,k] <- arima.sim(n = sim_length, list(ar = AR2),
            sd = sd(residuals(MODEL2))) + INTERCEPT + PRCP_LOW_FREQUENCY_COMPONENTS_MEAN
      }

      #Aggregate all components and obtain the simulated time-series
      PRCP_FINAL_ANNUAL_SIM[,y] <- rowSums(cbind(SIM1,SIM2))

      #Inverse transform
      #PRCP_FINAL_ANNUAL_SIM[,y] <- InverseBoxCox(lambda,PRCP_FINAL_ANNUAL_SIM[,y])

      #Y <- PRCP_FINAL_ANNUAL_SIM[,y]
      #result <- (lambda*Y + 1)  ^ (1/lambda)


      #Wavelet spectra for the simulated time-series
      WAVELET_VARIABLE <- PRCP_FINAL_ANNUAL_SIM[,y]
      CLIMATE_VARIABLE <- WAVELET_VARIABLE
      source(paste0(scripts_dir,'/WAVELET_ANALYSIS.R'))
      POWER_SPECTRUM_PRCP_ARIMA_SIM[,y] <- GWS

      setTxtProgressBar(pb, y)

    }
    close(pb)
    proc.time() - ptm

    #write simulated time-series to file
    name1 <- paste0("./warm_",dataset,"_",sim_num,".csv")
    write.table(paste0(";",Sys.time()), file = name1,
                col.names=FALSE, row.names = FALSE, sep = ",")
    write.table(t(PRCP_FINAL_ANNUAL_SIM),
                file = name1, sep = ",", row.names = FALSE, append=TRUE)

    #Write simulated power-spectra to file
    name2 <- paste0("./gws_",dataset,"_",sim_num,".csv")
    write.table(paste0(";",Sys.time()), file = name2,
                col.names=FALSE, row.names = FALSE, sep = ",")
    write.table(t(POWER_SPECTRUM_PRCP_ARIMA_SIM),
                file = name2, sep = ",", row.names = FALSE, append=TRUE)

  }

  #### WAVELET OUTPUT PLOTS

  # df <- data.frame(x=1:length(CLIMATE_VARIABLE),
  #     Original=CLIMATE_VARIABLE,
  #     Component=LOW_FREQUENCY_COMPONENTS, Noise=NOISE) %>%
  #   gather(key = variable, value = value, -x)
  #
  # #Time_series components
  # pl_comp <- ggplot(df, aes(x,y=value)) +
  #   geom_line() +
  #   facet_grid(variable  ~ ., scales='free') +
  #   labs(x='Years',y='')
  #
  # #PowerSpectrum Plot
  # period_lower_limit <- 0
  # sig_periods <- which(GWS>signif_GWS & period > period_lower_limit)
  #
  # par(mfrow=c(1,2),font.axis=2,font.lab=2,cex.axis=1.1,cex.lab=1.1)
  # xx <- 1:n1
  # yy <- period
  # image(x=xx,y=yy,z=t(log(POWER,base=2)),xlab="Time (Years)",
  #   ylab="Fourier Period (Years)",ylim=rev(range(yy)),log="y",col=(heat.colors(12)))
  # lines(xx,coi,lty=2)
  # contour(x=xx,y=yy,z=t(sig95),add=TRUE,levels=c(-99,1),labels="")
  #
  # xmin <- (min(GWS,signif_GWS))
  # xmax <- (max(GWS,signif_GWS))
  # plot(GWS,yy,type="b",xlim=c(xmin,xmax),xlab="Global Wavelet Spectrum",
  #   ylab="Fourier Period (Years)",log="y",ylim=rev(range(yy)))
  # lines(signif_GWS,yy,col="red",lty=2)
  # if (background_noise=="white") {
  #   tt <- paste(siglvl*100,"% confidence \nlevel for white-noise\n spectrum",sep="")}
  # if (background_noise=="red") {
  #   tt <- paste(siglvl*100,"% confidence \nlevel for red-noise\n spectrum",sep="")}
  # legend(.2*max(GWS),.65*max(yy),tt,col="red",lty=2,box.lty=0,box.lwd=0,cex=.8)
  # mtext("SIGNIFICANT PERIODS:",line=2)
  # e <- paste(sig_periods[1])
  # if (length(sig_periods)==0) {
  #   mtext("NONE!",line=1)
  # } else if (length(sig_periods)==1) {
  # 	mtext(e,line=1)
  # } else {
  # 	for (h in 2:length(sig_periods)) {
  # 		e <- paste(e,",",paste(sig_periods[h]),sep="")
  # 	}
  # 	mtext(e,line=1)
  # }

  #NORMALITY CHECK **

  #Box-cox transformation of the original data
  # lambda <- powerTransform(ANNUAL_PRCP)$lambda
  # TRANSFORM_DATA <- bcPower(ANNUAL_PRCP,lambda)
  #
  # #Check normality of original ts
  # data <- TRANSFORM_DATA
  # model <- data.frame(x=1:length(data), y=data)
  # pl_1 <- ggplot(model, aes(x,y)) + geom_line() + labs(y='mm/year', x='') +
  #   ggtitle('time-series')
  # pl_2 <- ggplot(model, aes(x=y)) +
  #   geom_histogram(fill='gray90', color='blue', binwidth=0.05) +
  #   labs(x='mm/year', y='count') + ggtitle('histogram')
  # pl_3 <- ggplot(model, aes(sample = y)) + stat_qq()  +
  #   geom_abline(intercept = mean(model$y), slope = sd(model$y)) + ggtitle('qqplot')
  #
  # p <- arrangeGrob(pl_2, pl_3, ncol=1)
  # grid.draw(p)
# ------------------------------------------------------------------------------
# CLIMATE DATA GENERATION ------------------------------------------------------

  #Steps:
  #1. Subset annual warm outputs based on climate statistics
  #2. Generate monthly realizations from annual data
  #3. Add trends to realizations

  ########################
  dataset <- "Princeton"
  ########################

  #Prepare annual and monthly climate data ***********************************
  prec_grid <- read_csv("./inputs/precip_grid_mm.csv")
  prec_data <- filter(prec_grid, Source == dataset) %>% .[["value"]]

  #Monthly precip with selected dataset
  clim_monthly <- filter(forcing_gcm, Index == 1) %>%
    mutate(Prec = prec_data) %>%
    filter(Year %in% c(1961:1999))

  #Annual precipitation (1950-1999)
  clim_annual  <- group_by(clim_monthly, Year) %>%
    summarise(Prec=sum(Prec),Tavg=mean(Tavg),Tmax=mean(Tmax),Tmin=mean(Tmin))

  ANNUAL_PRCP  <- clim_annual$Prec
  prec_warm    <- read_csv("./Inputs/warm_Princeton_100k.csv", skip = 1)
  gws_sim      <- read_csv("./Inputs/gws_Princeton_100k.csv", skip = 1)
  gws_obs      <- read_csv('./Inputs/gws_obs_Princeton.csv')

  #prec_warm <- read_csv(paste0('./Inputs/warm_', dataset,'.csv'), skip=1)
  #gws_sim <- read_csv(paste0('./Inputs/gws_', dataset,'.csv'), skip=1)


  #ggplot(gws_obs, aes(x=Period)) +
  #  geom_line(aes(y = Sig85), color ="blue") +
  #  geom_line(aes(y = GWS), color ="red")

################################################################################
################################################################################
### 1) CLIMATE REALIZATIONS SUBSETTTING

  precip_neg <- sapply(1:100000, function(x) any(prec_warm[x,] < 0))
  index_neg <- which(precip_neg == TRUE)

  gws_sim <- gws_sim[-index_neg,]
  prec_warm <- prec_warm[-index_neg,]

  #Significant period(s) in the historic data
  gws_sig <- which(gws_obs$GWS > 1  * gws_obs$Sig85)
  gws_sig <- unique(c(sapply(gws_sig, function(x) c(x-1,x,x+1))))
  gws_nonsig <- setdiff((1:nrow(gws_obs)), gws_sig)

  #Identify significant periods in each realization
  gws_sig_sim <- sapply(1:nrow(gws_sim),
    function(x) which(gws_sim[x,] > gws_obs$Sig85))
  gws_sig_sim[sapply(gws_sig_sim, function(x) length(x)==0)] <- NA

  #Subset realizations with the near sig periods of the historic
  if (length(gws_sig) == 0) {
    gws_sub <- which(is.na(gws_sig_sim))
  } else {
    gws_sub1 <- which(sapply(gws_sig_sim, function(x) any(gws_sig %in% x)))
    gws_sub2 <- which(sapply(gws_sig_sim, function(x) all(!x %in% gws_nonsig)))
    gws_sub <- intersect(gws_sub1, gws_sub2)
  }
  #Exclude realizations with different sig periods

  #Plot subsetted realizations
  period <- gws_obs$Period
  signif <- gws_obs$Sig85
  obs <- gws_obs$GWS
  sim <- gws_sim[gws_sub,1:length(period)]

  #Plot wavelet spectra, means, sds and skew
  font_size <- 9
  p2 <- viz_wavelet(period, signif, obs, sim, Ribbon = TRUE); p2

  #Plot wavelet spectra & climate stats
  stats_obs <- data_frame(value = ANNUAL_PRCP) %>%
    summarize(mean = mean(value), st.dev = sd(value))

  stats_sim <- t(prec_warm) %>%
    as.data.frame() %>% as_data_frame() %>%
    gather(key = variable, value = value) %>%
    group_by(variable) %>%
    summarize(mean = mean(value), st.dev = sd(value))

  #Subsetting bounds
  mean_rng <- c(mean(ANNUAL_PRCP)*(1-0.05), mean(ANNUAL_PRCP)*(1+0.05))
  sdev_rng <- c(sd(ANNUAL_PRCP)*(1-0.10), sd(ANNUAL_PRCP)*(1+0.10))

  #Subset from simulated
  sub_mean  <- which(
    (stats_sim$mean > mean_rng[1]) & (stats_sim$mean < mean_rng[2]))
  sub_stdev <- which(
    (stats_sim$st.dev > sdev_rng[1]) & (stats_sim$st.dev< sdev_rng[2]))
  sub_stats <- Reduce(intersect, list(sub_mean, sub_stdev))

  #Intersect all
  sub_clim <- intersect(sub_stats, gws_sub)
  sim <- gws_sim[sub_clim,1:length(period)]
  p1 <- viz_wavelet(period, signif, obs, sim, Ribbon = T); p1

  #RANDOMLY SAMPLE from matching warm outputs
  sample_size <- 10
  sub_clim_sample <- sample(sub_clim, 10)

  p1 <- ggplot(mapping = aes(x=mean, y=st.dev)) +
    theme_bw(base_size = 9) +
    geom_point(data = stats_sim, color = "black", alpha = 0.1, size = 0.1) +
    geom_point(data = stats_sim[sub_clim_sample,], color = "blue", size = 1.5) +
    geom_point(data = stats_obs, color = "red", size = 1.5) +
    lims(x=c(600,1200), y = c(100,500))


  prec_warm_select <- prec_warm[sub_clim,]
  #write_csv(prec_warm_select, "warm_Princeton_select.csv")


  #ggsave(filename = paste0("warm_out_",dataset,".png"), plot = p1)
  #ggsave(filename = paste0("gws_plot_",dataset,".png"), plot = p2)

################################################################################
################################################################################
##### 2) KNN BOOTSTRAPPING (ANNUAL TO MONTHLY)

  #Lists for storing climate data set
  btstrap_annual  <- list()
  btstrap_monthly <- list()

  ##### selected realizations for each dataset
  #selection <-data_frame(
  #  Princeton = c(7655, 5769, 1233, 5253, 1727, 7029, 9825, 3478, 1872, 6470),
  #  GPCC = c(4564, 2688, 2074, 5976, 8408, 1786, 6470, 3537, 8233,  319),
  #  CRU  = c(7864, 1363,  272, 4457,  491, 5275, 3335, 4976, 8682, 4012))

  #WARM outputs
  prec_st <- prec_warm_select %>% t(.)
  #prec_st <- prec_warm[sub_clim_sample,] %>% t(.)
  num <- ncol(prec_st)

  #stochastic climate traces
  hist_years <- 1961:1999

  #KNN bootstrap historic years (annual time-series)
  KNN_yrs <- apply(prec_st, c(1,2), function(x) hist_years[KNN(x,7,ANNUAL_PRCP)])

  #Annual climate from the bootstrapped years **********************************
  btstrap_annual[['Prec']] <- prec_st
  btstrap_annual[['Tavg']] <- apply(KNN_yrs, c(1,2),
    function(x) filter(clim_annual, Year == x)$Tavg)
  btstrap_annual[['Tmax']] <- apply(KNN_yrs, c(1,2),
    function(x) filter(clim_annual, Year == x)$Tmax)
  btstrap_annual[['Tmin']] <- apply(KNN_yrs, c(1,2),
    function(x) filter(clim_annual, Year == x)$Tmin)

  #Monthly climate *************************************************************
  m_factors <- data.frame(Year = clim_monthly$Year, Month = clim_monthly$Month)
  m_factors$Prec <- apply(clim_monthly, 1,
    function(x) x["Prec"] / clim_annual$Prec[clim_annual$Year==x["Year"]])
  m_factors$Tavg <- apply(clim_monthly, 1,
    function(x) x["Tavg"] / clim_annual$Tavg[clim_annual$Year==x["Year"]])
  m_factors$Tmin <- apply(clim_monthly, 1,
    function(x) x["Tmin"] / clim_annual$Tmin[clim_annual$Year==x["Year"]])
  m_factors$Tmax <- apply(clim_monthly, 1,
    function(x) x["Tmax"] / clim_annual$Tmax[clim_annual$Year==x["Year"]])

  #APPLY monthly factors to annual climate data *******
  btstrap_monthly <- lapply(1:ncol(prec_st), function(x) {
    expand.grid(Month = 1:12, KNN_year = KNN_yrs[,x]) %>%
    as.data.frame() %>% as_data_frame() %>%
    mutate(Index = rep(1:50, each = 12)) %>%
    select(Index, KNN_year, Month) %>%
    left_join(m_factors, by = c("KNN_year" = "Year", "Month")) %>%
    mutate(Prec = Prec * btstrap_annual$Prec[Index,x],
          Tavg = Tavg * btstrap_annual$Tavg[Index,x],
          Tmin = Tmin * btstrap_annual$Tmin[Index,x],
          Tmax = Tmax * btstrap_annual$Tmax[Index,x]) %>%
    mutate(Year = year(future_period)) %>%
    select(Year, Month, Prec:Tmax)})

  #SAVE STOCHASTIC REALIZATIONS TO FILE
  dat <- bind_rows(btstrap_monthly, .id = "nVar")
  #write_csv(x = dat, path = paste0("clim_realizations_new",dataset,".csv"))

################################################################################
################################################################################

### 3) CLIMATE TRENDS ADDITION

  data <- read_csv("./inputs/WARM_realizations.csv")

  outmat <- expand.grid(dataset = dataset_lst, nClim = 1:nrow(clim_tbl)) %>%
    as.tbl() %>% left_join(clim_tbl, by ="nClim") %>%
    mutate(dataset = as.character(dataset)) %>%
    arrange(dataset, nClim)

  traject <- vector(mode = 'list', length = nrow(outmat))

  #Loop through each row of the outmat and modify the climate data
  end <- nrow(outmat)
  for (i in 1:end) {

    #Extract chance factors & variability index in each row
    dataset <- outmat$dataset[i]
    nClim <- outmat$nClim[i]
    Tavg <- outmat$Tavg[i]
    Prec <- outmat$Prec[i]
    nvar <- outmat$nVar[i]

    #Vector of transient climate change factors
    t_trend <- rep(seq(0, Tavg, length.out=length(future_years)),each=12)
    p_trend <- rep(seq(1, Prec, length.out=length(future_years)),each=12)

    #Monthly climate table
    cli <- data %>% filter((Dataset == dataset) & (nVar == nvar)) %>%
      select(-nVar)

    traject[[i]] <- cli %>%
      mutate(nClim = nClim) %>%
      mutate(Prec = Prec * p_trend, Tavg = Tavg + t_trend,
             Tmin = Tmin + t_trend, Tmax = Tmax + t_trend)
  }

  #Bind all trajectories
  #traject_all <- bind_rows(traject)
  #write_csv(x = traject_all, path = "./inputs/WARM_trajectories.csv")
# ------------------------------------------------------------------------------

