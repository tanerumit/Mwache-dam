

## *****************************************************************************
## Mwache Weather Generator
## November, 2014
## M.Umit Taner
## *****************************************************************************


################################################################################

# Data -------------------------------------------------------------------------

  #Historic and projected dates
  dates_hist <- seq.Date(as.Date("1950-1-1"),as.Date("1999-12-1"), by="month")
  dates_gcms <- seq.Date(as.Date("2001-1-1"),as.Date("2050-12-1"), by="month")
  years_hist <- unique(year(dates_hist))

  #Monthly and annual precipitation (1950-1999)
  clim_monthly <- read.csv(paste0(getwd(),'/inputs/forcings_historic.csv'))
  clim_annual  <- group_by(clim_monthly, Year) %>%
    summarise(Prec=sum(Prec), Tavg=mean(Tavg), Tmax=mean(Tmax), Tmin=mean(Tmin))
  ANNUAL_PRCP  <- as.numeric(as.matrix(clim_annual[,2]))


################################################################################

# Wavelet Analysis   -----------------------------------------------------------


  #WAVELET ANALYSIS AND DECOMPOSITION ******************************************

  #Apply box-cox transformation of the original data
  lambda <- powerTransform(ANNUAL_PRCP)$lambda
  TRANSFORM_DATA <- bcPower(ANNUAL_PRCP,lambda)

  #Wavelet analysis on transformed data
  CLIMATE_VARIABLE <-  TRANSFORM_DATA
  source(paste0(scripts_dir, "/WAVELET_ANALYSIS.R"))
  POWER_SPECTRUM_PRCP_OBS <- GWS
  POWER_SPECTRUM_PRCP_SIGNIFICANCE <- signif_GWS
  POWER_SPECTRUM_PRCP_PERIOD <- period

  #WARM Calibration parameters
  #Number of orthogonal component series that representing a low-frequency signal
  NUM_FINAL_PERIODS <-  1
  #Length of the first component: This is an initial guess
  NUM_PERIODS_COMP1 <-  3
  #Length of all components
  NUM_PERIODS_ALL_COMPS <- c(NUM_PERIODS_COMP1) #, NUM_PERIODS_COMP2)
  #List indices of significant components (from all components)
  ALL_SIG_PERIODS <- c(5,6,7)

  #Wavelet decomposition
  source(paste0(scripts_dir,'/WAVELET_DECOMPOSITION.R'))
  PRCP_LOW_FREQUENCY_COMPONENTS <- LOW_FREQUENCY_COMPONENTS
  PRCP_NOISE <- NOISE

  #WAVELET OUTPUT PLOTS ********************************************************

  #Time_series components
  pl_comp <- ggplot(data.frame(x=1:length(CLIMATE_VARIABLE),
      Original=CLIMATE_VARIABLE,
      Component=LOW_FREQUENCY_COMPONENTS, Noise=NOISE) %>%
      melt(id.vars=c('x')), aes(x,y=value)) +
    geom_line() +
    facet_grid(variable  ~ ., scales='free') +
    labs(x='Years',y='')

  #PowerSpectrum Plot
  period_lower_limit <- 0
  sig_periods <- which(GWS>signif_GWS & period > period_lower_limit)

  par(mfrow=c(1,2),font.axis=2,font.lab=2,cex.axis=1.1,cex.lab=1.1)
  xx <- 1:n1
  yy <- period
  image(x=xx,y=yy,z=t(log(POWER,base=2)),xlab="Time (Years)",
    ylab="Fourier Period (Years)",ylim=rev(range(yy)),log="y",col=(heat.colors(12)))
  lines(xx,coi,lty=2)
  contour(x=xx,y=yy,z=t(sig95),add=TRUE,levels=c(-99,1),labels="")

  xmin <- (min(GWS,signif_GWS))
  xmax <- (max(GWS,signif_GWS))
  plot(GWS,yy,type="b",xlim=c(xmin,xmax),xlab="Global Wavelet Spectrum",
    ylab="Fourier Period (Years)",log="y",ylim=rev(range(yy)))
  lines(signif_GWS,yy,col="red",lty=2)
  if (background_noise=="white") {
    tt <- paste(siglvl*100,"% confidence \nlevel for white-noise\n spectrum",sep="")}
  if (background_noise=="red") {
    tt <- paste(siglvl*100,"% confidence \nlevel for red-noise\n spectrum",sep="")}
  legend(.2*max(GWS),.65*max(yy),tt,col="red",lty=2,box.lty=0,box.lwd=0,cex=.8)
  mtext("SIGNIFICANT PERIODS:",line=2)
  e <- paste(sig_periods[1])
  if (length(sig_periods)==0) {
    mtext("NONE!",line=1)
  } else if (length(sig_periods)==1) {
  	mtext(e,line=1)
  } else {
  	for (h in 2:length(sig_periods)) {
  		e <- paste(e,",",paste(sig_periods[h]),sep="")
  	}
  	mtext(e,line=1)
  }



  #NORMALITY CHECK *************************************************************

  #Box-cox transformation of the original data
  lambda <- powerTransform(ANNUAL_PRCP)$lambda
  TRANSFORM_DATA <- bcPower(ANNUAL_PRCP,lambda)

  #Check normality of original ts
  data <- TRANSFORM_DATA
  model <- data.frame(x=1:length(data), y=data)
  pl_1 <- ggplot(model, aes(x,y)) + geom_line() + labs(y='mm/year', x='') +
    ggtitle('time-series')
  pl_2 <- ggplot(model, aes(x=y)) +
    geom_histogram(fill='gray90', color='blue', binwidth=0.05) +
    labs(x='mm/year', y='count') + ggtitle('histogram')
  pl_3 <- ggplot(model, aes(sample = y)) + stat_qq()  +
    geom_abline(intercept = mean(model$y), slope = sd(model$y)) + ggtitle('qqplot')

  p <- arrangeGrob(pl_2, pl_3, ncol=1)
  print(p)

  #DATA GENERATION *************************************************************

  #settings
  sim_length <- 50
  sim_num    <- 5000

  #sim_numate y times
  PRCP_FINAL_ANNUAL_SIM             <- array(NA,c(sim_length,sim_num))
  PRCP_FINAL_ANNUAL_SIM_TRANSFORMED <- array(NA,c(sim_length,sim_num))
  POWER_SPECTRUM_PRCP_ARIMA_SIM     <- array(NA,c(length(GWS),sim_num))

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
      if (length(which(names(MODEL2$coef)=="intercept"))>0) {INTERCEPT <- as.vector(MODEL2$coef)[which(names(MODEL2$coef)=="intercept")]}
      SIM2[,k] <- arima.sim(n = sim_length, list(ar = AR2),sd = sd(residuals(MODEL2))) + INTERCEPT + PRCP_LOW_FREQUENCY_COMPONENTS_MEAN
    }

    #Aggregate all components and obtain the simulated time-series
    PRCP_FINAL_ANNUAL_SIM[,y] <- rowSums(cbind(SIM1,SIM2))

    #Inverse transform
    PRCP_FINAL_ANNUAL_SIM[,y] <- InverseBoxCox(lambda,PRCP_FINAL_ANNUAL_SIM[,y])

    #Wavelet spectra for the simulated time-series
    WAVELET_VARIABLE <- PRCP_FINAL_ANNUAL_SIM[,y]
    CLIMATE_VARIABLE <- WAVELET_VARIABLE
    source(paste0(scripts_dir,'/WAVELET_ANALYSIS_NOPLOT.r'))
    POWER_SPECTRUM_PRCP_ARIMA_SIM[,y] <- GWS
  }
  proc.time() - ptm

  #Simulated WARM annual rainfall
  #name1 <- paste0(inputs_dir,"/WARM_SIM_",sim_length,"_",sim_num,".csv")
  #write.table(paste0(";",Sys.time()), file = name1, col.names=FALSE, row.names = FALSE, sep = ",")
  #write.table(t(PRCP_FINAL_ANNUAL_SIM), file = name1, sep = ",", row.names = FALSE, append=TRUE)

  #Simulated Power spectra
  #name2 <- paste0(inputs_dir,"/WARM_GWS_",sim_length,"_",sim_num,".csv")
  #write.table(paste0(";",Sys.time()), file = name2, col.names=FALSE, row.names = FALSE, sep = ",")
  #write.table(t(POWER_SPECTRUM_PRCP_ARIMA_SIM), file = name2, sep = ",", row.names = FALSE, append=TRUE)


source('/Users/umit/GDrive/Research/Scripts/WARM_PLOTTING.R')

################################################################################

# Realizations subsetting  -----------------------------------------------------

  #Load observed and simulated annual precipitation
  precip_obs <- ANNUAL_PRCP
  precip_sim <- data.frame(read.csv('./inputs/WARM_SIM_50_5000.csv', skip=1, header=TRUE))

  #Load Global Wavelet spectra(GWS) for observed and simulated annual precipitation
  gws_precip_obs <- data.frame(read.csv('./inputs/WARM_GWS_observed.csv', header=TRUE))
  gws_precip_sim <- data.frame(read.csv('./inputs/WARM_GWS_50_5000.csv', skip=1, header=TRUE))

  #Plot wavelet spectra & climate stats
  stats_obs <- data.frame(X='x', Type='Obs', Mean=mean(precip_obs),
    ST.Dev=sd(precip_obs), Skewness=skewness(precip_obs))
  stats_sim <- data.frame(X='x', Type='Sim',
    t(apply(precip_sim, 1, function(x) c(mean(x), sd(x), skewness(x)))))
  colnames(stats_sim) <- c('X','Type','Mean', 'ST.Dev', 'Skewness')

  gws_period <- gws_precip_obs$period
  gws_sig <- gws_precip_obs$signif85
  WaveletPlot (gws_obs = gws_precip_obs$GWS, gws_sim = gws_precip_sim, Ribbon = FALSE)
  ClimStatsPlot(var_obs = precip_obs, var_sim = precip_sim)

  #Trace subsetting based on GWS ***********************************************

  #Significant period(s) in historic data
  GWS_sig_index_obs <- which(gws_precip_obs$GWS > gws_precip_obs$signif85)
  GWS_sig_bound <- c(GWS_sig_index_obs-1, GWS_sig_index_obs, GWS_sig_index_obs+1)
  GWS_nonsig_index <- (1:length(gws_precip_obs$GWS))[!(1:length(gws_precip_obs$GWS)) %in% GWS_sig_bound]

  #Identify significant GWS periods
  Sigband_GWS  <- as.matrix(gws_precip_obs$signif85)
  GWS_data <- as.matrix(gws_precip_sim)
  GWS_sig_index  <- apply(GWS_data, 1, function(x) which(x > Sigband_GWS))
  GWS_sig_index[sapply(GWS_sig_index, function(x) length(x)==0)] <- NA

  #Select near bound of Sig GWS, but not the others
  GWS_subset1  <- which(sapply(GWS_sig_index, function(x) GWS_sig_index_obs %in% x))
  GWS_subset2  <- which(sapply(GWS_sig_index, function(x) all(!x %in% GWS_nonsig_index)))
  GWS_subset <- intersect(GWS_subset1, GWS_subset2)

  #Plot wavelet spectra, means, sds and skew
  WaveletPlot (gws_obs = gws_precip_obs$GWS, gws_sim = gws_precip_sim[GWS_subset,], Ribbon = TRUE)
  ClimStatsPlot(var_obs = precip_obs, var_sim = precip_sim[GWS_subset,])

  #Trace subsetting based on Climate statistics *********************************************************************

  #Subsetting bounds
  mean_bound <- c(mean(precip_obs)-30, mean(precip_obs)+30)
  sdev_bound <- c(sd(precip_obs)-30, sd(precip_obs)+30)
  skew_bound <- c(skewness(precip_obs)-0.2, skewness(precip_obs)+0.2)

  #Subset from simulated
  subset_mean  <- which((stats_sim$Mean > mean_bound[1]) & (stats_sim$Mean < mean_bound[2]))
  subset_sdev  <- which((stats_sim$ST.Dev > sdev_bound[1]) & (stats_sim$ST.Dev< sdev_bound[2]))
  subset_skew  <- which((stats_sim$Skewness > skew_bound[1]) & (stats_sim$Skewness < skew_bound[2]))
  stats_subset <- Reduce(intersect, list(subset_mean, subset_sdev, subset_skew))

  #Intersect all
  subset_clim <- intersect(stats_subset, GWS_subset)
  WaveletPlot (gws_obs = gws_precip_obs$GWS, gws_sim = gws_precip_sim[subset_clim,], Ribbon = TRUE)
  ClimStatsPlot(var_obs = precip_obs, var_sim = precip_sim[subset_clim,])

  #Plot means vs variance
  stats_sim_subset <- data.frame(X='x', Type='Sim', t(apply(precip_sim[stats_subset,], 1, function(x) c(mean(x), sd(x), skewness(x)))))
  colnames(stats_sim_subset) <- c('X','Type','Mean', 'ST.Dev', 'Skewness')
  pl_meanvar <- ggplot(mapping=aes(x=Mean, y=ST.Dev)) +
    geom_point(data=stats_sim, size=1, color="gray") +
    geom_point(data=stats_sim_subset, size=1, color="red") +
    geom_point(data = stats_obs, size=2, color="blue") +
    labs(x= 'mean (mm)', y='sd (mm)') +
    scale_x_continuous(expand=c(0,0), limits=c(750,950)) +
    scale_y_continuous(expand=c(0,0), limits=c(100,400)) +
    THEME_FONT_SIZE(x=font_size)


  #Dissagregate selected traces to monthly  *****************************************************************************

  #Lists for storing climate data set
  btstrap_clim_annual   <- list()     #annual climate set
  btstrap_clim_monthly  <- list()     #monthly climate set

  #stochastic climate traces
  selection <- subset_clim
  warm_prec_annual <- data.frame(read.csv('./inputs/WARM_SIM_50_5000.csv', skip=1, header=TRUE))
  warm_prec_annual <- t(warm_prec_annual[selection,])
  stochast_num <- ncol(warm_prec_annual)

  #climate change space
  #x=Temp. additive factor, y=Prec. multiplier, z=stochastic realization
  #KNN bootstrap historic years (annual time-series)
  KNN_years  <- apply(warm_prec_annual, c(1,2), function(x) years_hist[KNN(x,7,ANNUAL_PRCP)])
  KNN_order  <- apply(KNN_years, c(1,2), function(x) which(x == years_hist))

  #Prepare a annual climate data for all climate variables from the bootstrapped years
  btstrap_clim_annual[['Prec']] <- warm_prec_annual
  btstrap_clim_annual[['Tavg']] <- apply(KNN_years, c(1,2), function(x) filter(clim_annual, Year == x)$Tavg)
  btstrap_clim_annual[['Tmax']] <- apply(KNN_years, c(1,2), function(x) filter(clim_annual, Year == x)$Tmax)
  btstrap_clim_annual[['Tmin']] <- apply(KNN_years, c(1,2), function(x) filter(clim_annual, Year == x)$Tmin)

  #Monthly factors
  m_factors <- data.frame(Year = clim_monthly[,1], Month = clim_monthly[,2])
  m_factors$Prec <- apply(clim_monthly, 1, function(x) x[3] / clim_annual$Prec[clim_annual$Year==x[1]])
  m_factors$Tavg <- apply(clim_monthly, 1, function(x) x[4] / clim_annual$Tavg[clim_annual$Year==x[1]])
  m_factors$Tmax <- apply(clim_monthly, 1, function(x) x[5] / clim_annual$Tmax[clim_annual$Year==x[1]])
  m_factors$Tmin <- apply(clim_monthly, 1, function(x) x[6] / clim_annual$Tmin[clim_annual$Year==x[1]])

  #Monthly time-series
  btstrap_clim_monthly[['Prec']] <-  sapply(1:stochast_num, function(y)
    apply(m_factors, 1, function(x) btstrap_clim_annual[['Prec']][length(1950:x[1]),y] * x[3]))
  btstrap_clim_monthly[['Tavg']] <-  sapply(1:stochast_num, function(y)
    apply(m_factors, 1, function(x) btstrap_clim_annual[['Tavg']][length(1950:x[1]),y] * x[4]))
  btstrap_clim_monthly[['Tmax']] <-  sapply(1:stochast_num, function(y)
    apply(m_factors, 1, function(x) btstrap_clim_annual[['Tmax']][length(1950:x[1]),y] * x[5]))
  btstrap_clim_monthly[['Tmin']] <-  sapply(1:stochast_num, function(y)
    apply(m_factors, 1, function(x) btstrap_clim_annual[['Tmin']][length(1950:x[1]),y] * x[6]))

  #Organize data-set
  forcings <- lapply(1:stochast_num, function(x) data.frame(Year=year(dates_hist), Mon=month(dates_gcms),
    Prec = btstrap_clim_monthly[['Prec']][,x], Tavg = btstrap_clim_monthly[['Tavg']][,x],
    Tmax = btstrap_clim_monthly[['Tmax']][,x], Tmin = btstrap_clim_monthly[['Tmin']][,x]))

  #RUN TRACES IN WEAPSIM *************************************************************************************************************************

  #Simulation parameters
  netevap    <- c(132, 168, 135, 44, -23, 51, 45, 74, 104, 77, 45, 88)*1e-3 #m/month
  elev_curve <- data.frame(vol=c(0,1,5,15,30,50,81.5,110,160), elev=c(14, 24, 34, 44, 54, 64, 75.6, 81.1, 84))
  e_approx   <- approxfun(elev_curve[,1],elev_curve[,2])
  s_approx   <- approxfun(elev_curve[,2],elev_curve[,1])
  days_month = c(31,28,31,30,31,30,31,31,30,31,30,31)
  monthly_factors = days_month/sum(days_month)
  MonthlyDemand= 100 * monthly_factors  #78
  hydro_par = c(1.202155e-03, 2.517880e+02, 9.170601e-01, 1.089352e-04)

  #Calculate PET across scenarios
  YEARS <- as.numeric(format(dates_hist,"%Y"))
  MONTHS <- as.numeric(format(dates_hist,"%m"))
  trace_num <- length(forcings)
  metric <- vector()
  output <- list()

  #Run through all subsetted traces
  for (i in 1: trace_num) {

    PE <- Hargreaves(tavg=forcings[[i]]$Tavg, tdif=forcings[[i]]$Tmax-forcings[[i]]$Tmin, BasinLat=3.5) # mm/month

    output[[i]] <- WEAPSIM(YEARS, MONTHS, Area=2250, hydro_par=hydro_par, S_ini=800, GW_ini=1,
      Precip=forcings[[i]]$Prec, PET=PE, StorageCap=60, StorageBegin=50,
      TopConservation=60, TopBuffer=20, TopInactive=20,
      BufferCoef=1, MonthlyDemand=MonthlyDemand)[121:600,]

    metric[i] <- length(which((output[[i]]$Demand - output[[i]]$Release)==0))/length(output[[i]]$Demand)
  }

  #Run historic
  PE <- Hargreaves(tavg=clim_monthly$Tavg, tdif=clim_monthly$Tmax-clim_monthly$Tmin, BasinLat=3.5) # mm/month
  out_hist <- WEAPSIM(YEARS, MONTHS, Area=2250, hydro_par=hydro_par, S_ini=800, GW_ini=1,
    Precip=clim_monthly$Prec, PET=PE, StorageCap=60, StorageBegin=50,
    TopConservation=60, TopBuffer=20, TopInactive=20,
    BufferCoef=1, MonthlyDemand=MonthlyDemand)[121:600,]
  metric_hist <- data.frame(Trace='Historic', value=length(which((out_hist$Demand - out_hist$Release)==0))/length(out_hist$Demand))

  #Select quantiles
  sample_quantiles = c(0.05, seq(0.1,0.9,0.1), 0.95)
  qtile_select <- sapply(sample_quantiles, function(x) which(quantile(metric, x, names=FALSE, type = 3) == metric))
  qtile_select <- sapply(qtile_select, function(x) sample(x,1))


  #plot performance

  select_df <- data.frame(Trace=1:length(metric), value=metric, Type='ALL', stringsAsFactors=FALSE)
  select_df$Type[qtile_select] <- 'Select'

  #Plot selection
  #Boxplot
  pl_traceSelect <- ggplot(select_df, aes(x='', y=value)) + geom_boxplot() +
    geom_point(data=filter(select_df, Type=='Select'), color='red',size=3)+
    geom_point(data=metric_hist, color='blue', size=3)


  #Experimental violin plot
  pl_traceSelect <- ggplot(select_df, aes(x='', y=value)) + geom_violin() +
    geom_point(color='black',size=1)+
    geom_point(data=filter(select_df, Type=='Select'), color='red',size=1)+
    geom_point(data=metric_hist, color='blue', size=2)


  global_index <- subset_clim[qtile_select]
  #4553 1189  533 1450  311 3880  215  893 1659 4350 3555 (11 prior to Feb 7)
  #1384 1886  517 2161 4278  892  439 4151  938  524 (10 on Feb 7 2015)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

################################################################################

# KNN $ CLIMATE TRENDS  --------------------------------------------------------

  #Lists for storing climate data set
  btstrap_clim_annual   <- list()     #annual climate data set
  btstrap_clim_monthly  <- list()     #monthly climate data set
  stochast_clim_monthly <- list()     #final climate data set with trends

  #stochastic climate traces
  selection <- c(1384, 1886,  517, 2161, 4278,  892,  439, 4151,  938,  524)
  warm_prec_annual <- data.frame(
    read.csv('./inputs/WARM_SIM_50_5000.csv', skip=1, header=TRUE))
  warm_prec_annual <- t(warm_prec_annual[selection,])
  stochast_num <- ncol(warm_prec_annual)

  #Climate Space ***************************************************************

  #Mean changes
  cc_temp_mid <- seq(0,2.5,0.5)
  cc_prec_mid <- seq(0.7,1.4,0.1)

  #Final changes for linear trends
  cc_temp_end <- cc_temp_mid * 2
  cc_prec_end <- ifelse(cc_prec_mid <= 1,
    1-(1-cc_prec_mid)*2, 1 + (cc_prec_mid-1)*2)

  #Climate space
  cc_space_ind <- expand.grid(Temp=1:length(cc_temp_mid),
    Prec=1:length(cc_prec_mid))
  cc_space <- expand.grid(Temp=cc_temp_mid, Prec=cc_prec_mid)
  clim_space_ind <- expand.grid(Temp=1:length(cc_temp_mid),
    Prec=1:length(cc_prec_mid), TS=1:stochast_num)

  #KNN bootstrap historic years (annual time-series)
  KNN_years  <- apply(warm_prec_annual, c(1,2),
    function(x) years_hist[KNN(x,7,ANNUAL_PRCP)])
  KNN_order  <- apply(KNN_years, c(1,2),
    function(x) which(x == years_hist))

  #Climate data from bootstrapped years
  btstrap_clim_annual[['Prec']] <- warm_prec_annual
  btstrap_clim_annual[['Tavg']] <- apply(KNN_years, c(1,2),
    function(x) filter(clim_annual, Year == x)$Tavg)
  btstrap_clim_annual[['Tmax']] <- apply(KNN_years, c(1,2),
    function(x) filter(clim_annual, Year == x)$Tmax)
  btstrap_clim_annual[['Tmin']] <- apply(KNN_years, c(1,2),
    function(x) filter(clim_annual, Year == x)$Tmin)

  #Monthly factors for each climate variable
  m_factors <- data.frame(Year = clim_monthly[,1], Month = clim_monthly[,2])
  m_factors$Prec <- apply(clim_monthly, 1,
    function(x) x[3] / clim_annual$Prec[clim_annual$Year==x[1]])
  m_factors$Tavg <- apply(clim_monthly, 1,
    function(x) x[4] / clim_annual$Tavg[clim_annual$Year==x[1]])
  m_factors$Tmax <- apply(clim_monthly, 1,
    function(x) x[5] / clim_annual$Tmax[clim_annual$Year==x[1]])
  m_factors$Tmin <- apply(clim_monthly, 1,
    function(x) x[6] / clim_annual$Tmin[clim_annual$Year==x[1]])

  #Apply monthly factors to bootstrapped climate data
  btstrap_clim_monthly[['Prec']] <-  sapply(1:stochast_num,
    function(y) apply(m_factors, 1,
    function(x) btstrap_clim_annual[['Prec']][length(1950:x[1]),y] * x[3]))

  btstrap_clim_monthly[['Tavg']] <-  sapply(1:stochast_num,
    function(y) apply(m_factors, 1,
    function(x) btstrap_clim_annual[['Tavg']][length(1950:x[1]),y] * x[4]))

  btstrap_clim_monthly[['Tmax']] <-  sapply(1:stochast_num,
    function(y) apply(m_factors, 1,
    function(x) btstrap_clim_annual[['Tmax']][length(1950:x[1]),y] * x[5]))

  btstrap_clim_monthly[['Tmin']] <-  sapply(1:stochast_num,
    function(y) apply(m_factors, 1,
    function(x) btstrap_clim_annual[['Tmin']][length(1950:x[1]),y] * x[6]))

  #Monthly climate forcings (variability + trend)
  t_trend <- sapply(cc_temp_end, function(x)
    rep(seq(0, x, length.out=length(years_hist)),each=12))
  p_trend <- sapply(cc_prec_end, function(x)
    rep(seq(1, x, length.out=length(years_hist)),each=12))

  #Trends added to the stochastic realizations
  stochast_clim_monthly[['Prec']] <- apply(clim_space_ind, 1, function(x)
    btstrap_clim_monthly[['Prec']][,x[3]] * p_trend[,x[2]])
  stochast_clim_monthly[['Tavg']] <- apply(clim_space_ind, 1, function(x)
    btstrap_clim_monthly[['Tavg']][,x[3]] + t_trend[,x[1]])
  stochast_clim_monthly[['Tmax']] <- apply(clim_space_ind, 1, function(x)
    btstrap_clim_monthly[['Tmax']][,x[3]] + t_trend[,x[1]])
  stochast_clim_monthly[['Tmin']] <- apply(clim_space_ind, 1, function(x)
    btstrap_clim_monthly[['Tmin']][,x[3]] + t_trend[,x[1]])

  #Organize data-set
  forcings <- lapply(1:nrow(clim_space_ind),
    function(x) data.frame(Year=year(dates_gcms), Month=month(dates_gcms),
    Prec = stochast_clim_monthly[['Prec']][,x], Tavg = stochast_clim_monthly[['Tavg']][,x],
    Tmax = stochast_clim_monthly[['Tmax']][,x], Tmin = stochast_clim_monthly[['Tmin']][,x]))

  #Cross-check
  #check1 <- t(sapply(forcings, function(x) apply(x, 2, mean))) %>%
  #  cbind(clim_space,.)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  #Write results to directory
  for (i in 1:nrow(clim_space_ind)) {
    file_name   <- paste0('./inputs/ST_forcings_022015/clim_',i,'.csv')
    file_header <- rbind(paste0('Mombasa: ', Sys.time(), collapse=''),
      paste0('Precip/Temp/Realization: ', paste0(clim_space_ind[i,], collapse='-')))
    write.table(file_header, file_name, row.names=FALSE, col.names=FALSE)
    write.table(forcings[[i]], file_name, sep=',', row.names=FALSE, col.names=TRUE, append=TRUE)
  }



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

################################################################################



