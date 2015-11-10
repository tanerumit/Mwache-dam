
#-----------------------------------------------------------------------------#
# Title: Decision analysis model for Mwache Dam design
# By:    M.Umit Taner
# Date:  September, 2015                                                       
#-----------------------------------------------------------------------------# 
  

  source("Mw_Input_pars.R")

# Vulnerabilty Analysis -------------------------------------------------------

  #Socio-economic samples
  ValueWater = c(0.5, 0.7, 0.9)
  DiscRate   = c(0.02, 0.05, 0.08)
  ResLife    = c(40, 45, 50)
  
  NPV_tbl_sec <- expand.grid(ValueWater = ValueWater, 
      DiscRate = DiscRate, ResLife=ResLife) %>% as_data_frame() 
    
  
  #Simulated Yield: delta[cc] x delta[Variability]
  syield <- read.csv('./inputs/safeyield_Sep25.csv') %>% as_data_frame()
  
  #Look up table for climate parameters
  clim_index <- clim_matrix %>% as_data_frame() %>%
    dplyr::select(Clim, CC, NVar, Tavg, Prec)
  
  #Simulated yield for: delta[P] x delta[T] x delta[Variability]
  syield_sim <- syield %>%
    dplyr::select(Clim=c, d, yield) %>%
    arrange(d) %>%
    left_join(clim_index, "Clim") %>%
    dplyr::select(d, Clim, CC, Tavg, Prec, NVar, yield) %>%
    mutate(d = as.integer(d), CC = as.integer(CC), NVar = as.integer(NVar))
  
  #Sample_results table
  NPV_tbl <- expand.grid.df(syield_sim, NPV_tbl_sec) %>% 
    as_data_frame() %>% rowwise() %>%
    mutate(AnnB = yield * ValueWater) %>%
    mutate(PVB  = PV_CALCULATE(rep(AnnB, ResLife), r=DiscRate)) %>%
    mutate(NPV  = PVB - designs$PVCost[d])

  NPV_cdf <- NPV_tbl %>%
    ungroup() %>%
    mutate(d = ordered(.$d, labels = designs$Labels)) %>%
    group_by(d) %>%
    arrange(desc(NPV)) %>%
    mutate(Index =  row_number(), ExP = Index/(max(Index)+1))


  ##############################################################################
  
#   #PARALEL COORDINATE PLOTS | Interactive d3.js 
#   
#   
#   devtools::install_github("timelyportfolio/parcoords")
#   require(parcoords)
#   
#   df <- dplyr::select(NPV_tbl, d, 
#     Tavg, Prec, NVar, yield, ValueWater, DiscRate, ResLife, NPV)
#   
#   parcoords(df,
#     rownames = F # turn off rownames from the data.frame
#     , brushMode = "2D-strums"
#     , reorderable = T
#     , queue = T
#     , color = list(
#       colorBy = "d"
#       ,colorScale = htmlwidgets::JS("d3.scale.category10()")
#     )
#   )
#   
#   # Paralell Coordinate Plots ----------------------------------------------------
#   
#   pl_cols <- rev(c("#fc8d59","#fee08b","#ffffbf","#e6f598","#99d594","#3288bd"))
#   
#   for (i in 1:nrow(designs)) {
#     
#     dat <- dplyr::select(NPV_tbl, d, Tavg, Prec, NVar, ValueWater, DiscRate, ResLife, NPV) %>%
#       filter(d == i) %>% dplyr::select(-d)
#     
#     png(filename = paste0("plot_",i,".png"), width = 1400, height = 600)
#     parcoordlabel2(dat, result_var = "NPV", result_var_range = seq(-100,500, 100), 
#                    col = ifelse(dat[,"NPV"]> 400, pl_cols[1],
#                                 ifelse(dat[,"NPV"]> 300, pl_cols[2],
#                                        ifelse(dat[,"NPV"]> 200, pl_cols[3],
#                                               ifelse(dat[,"NPV"]> 100, pl_cols[4],     
#                                                      ifelse(dat[,"NPV"]> 0, pl_cols[5], pl_cols[6]))))))
#     
#     graphics.off()
#     
#   }

  
  ##############################################################################
  
  
  pl_cols <- c("#2b83ba","#abdda4","#ffffbf","#fdae61","#d7191c")
  pl_levs <- c(3000, 2500, 2000, 1000, 500)
  
  png(filename = paste0("par_coord_compare.png"), width = 800, height = 1000)  
  par(mfrow=c(2,1), mar=c(4,3,2,2)+0.1)
  
  
  #range(NPV_tbl$NPV)
  
  #min(NPV_tbl[which(NPV_tbl$NPV < 0),"NPV"])
  
  dat1 <- dplyr::select(NPV_tbl, d, Tavg, Prec, NVar, ValueWater, DiscRate, ResLife, NPV) %>%
    filter(d == 1) %>% dplyr::select(-d)
  parcoordlabel2(dat1, main = "Capacity: 40_MCM",
                 result_var = "NPV", result_var_range = seq(0, 4000, 1000), 
                 col = ifelse(dat1[,"NPV"]> pl_levs[1], pl_cols[1],
                              ifelse(dat1[,"NPV"]> pl_levs[2], pl_cols[1],
                                     ifelse(dat1[,"NPV"]> pl_levs[3], pl_cols[2],
                                            ifelse(dat1[,"NPV"]> pl_levs[4], pl_cols[3],     
                                                   ifelse(dat1[,"NPV"]> pl_levs[5], pl_cols[4], pl_cols[5]))))))
  
  
  dat2 <- dplyr::select(NPV_tbl, d, Tavg, Prec, NVar, ValueWater, DiscRate, ResLife, NPV) %>%
    filter(d == 6) %>% dplyr::select(-d)
  
  parcoordlabel2(dat2, main = "Capacity: 140_MCM",
                 result_var = "NPV", result_var_range = seq(0, 4000, 1000), 
                 col = ifelse(dat1[,"NPV"]> pl_levs[1], pl_cols[1],
                              ifelse(dat1[,"NPV"]> pl_levs[2], pl_cols[1],
                                     ifelse(dat1[,"NPV"]> pl_levs[3], pl_cols[2],
                                            ifelse(dat1[,"NPV"]> pl_levs[4], pl_cols[3],     
                                                   ifelse(dat1[,"NPV"]> pl_levs[5], pl_cols[4], pl_cols[5]))))))
  
  dev.off()
  
  ##############################################################################

  
  
# Climate Response Surfaces -----------------------------------------------
  
  
  metric_sim <- syield_sim %>%
    select(d, Tavg, Prec, NVar, yield)
    
  
  metric_avg <- metric_sim %>%
    group_by(d, Tavg, Prec) %>%
    summarize(yield = mean(yield)) %>%  
    spread(key = d, value=yield)


  #Intpolate for surface plots
  res <- 500
  metric_to_plot <- GRID_INTERPOLATE(
    x=metric_avg$Tavg,
    y=metric_avg$Prec,
    z=select(metric_avg, -Tavg, -Prec), resolution=res) %>%
    as_data_frame()
  
  #Set metric range & bins
  var_xy <- metric_to_plot %>% select(Temp = x, Prec = y)
  var_pl <- as_data_frame(metric_to_plot) %>% select(-x, -y)

  #Date-frame for plotting
  metric_df <- cbind(var_xy, var_pl) %>% as_data_frame() %>%
    gather(Designs, Value, -Temp, -Prec) %>%
    dplyr::filter(Designs %in% c(1,2,3,4,5,6))
  
  #Bins and colors
  var_bin <- pretty(c(max(var_pl), min(var_pl)), n=4)
  var_bin <- seq(30, 120, 10)
  var_col <- c(rev(brewer.pal(4, "Reds")), "#ffffff", brewer.pal(4, "Blues"))
    
    
  #var_col <- c("#e6550d", "#fdbe85", "#feedde", "#f0f9e8", "#7bccc4", "#0868ac")
  #var_col <- colorRampPalette(brewer.pal(length(var_bin), "Spectral"))(length(var_bin))
  metric_df$Bins <- cut(metric_df$Value, breaks=var_bin, dig.lab=5, include.lowest = TRUE)
  

  levels(metric_df$Designs) <- paste0("Capacity: ", designs$Labels)
  metric_df <- metric_df %>%
    group_by(Designs, Temp) %>%
    mutate(Thold = abs(Value - 75)) %>%
    mutate(Tline = ifelse(Thold == min(Thold), 1, NA)) %>%
    ungroup()
  
  
  filter(metric_df, Tline == 1)

  #Climate Surface Plot
  plot <- pl_surface +
    geom_tile(data=metric_df, aes(x=Temp, y=Prec, fill=Bins)) +
    geom_point(aes(x=Temp, y=Prec),
      data=dplyr::filter(gcms_data, Ensemble %in% 'CMIP5'), size=2) +
    geom_line(data=filter(metric_df, Tline == 1), aes(x=Temp, y=Prec), size = 1) +
    facet_wrap( ~  Designs, ncol = 2) +
    scale_fill_manual(values=var_col, name="Safe yield\n(MCM/year)")   +
    #scale_color_manual(values=c("red","black"), name="Climate\nProjections") +
    guides(fill = guide_legend(ncol = 1, order = 1, keyheight = 1.5), color = guide_legend(order=2))  

  ggsave(file="surface_syields_all2.png", plot = plot, width=17*1.2, height=22*1.2, unit='cm')
  

  #PLOT POSTERIOR PROBABILITIES.....
  clim_posteriors <- data.frame(cc_matrix, Value = Prior) %>%
    as_data_frame() %>%
    mutate(Tavg = as.factor(Tavg), Prec = as.factor(Prec))
  

  
  
  
  ##############################################################################


### Density function
p1 <- ggplot(NPV_cdf, aes(x=NPV, group=d, color=d)) +
  theme_bw() +
  stat_density(stats="identity", fill=NA) +
  scale_color_brewer(name="Designs Capacities", palette="Spectral") +
  labs(x="NPV", y = "frequency"); p1

### Cummulative probabilty
p2 <- ggplot(NPV_cdf, aes(x=ExP, y=NPV, group=d, color=d)) +
  theme_bw() +
  geom_line(size=1) +
  scale_color_brewer(name="Designs Capacities", palette="Spectral") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.1), expand=c(0,0)) +
  geom_vline(xintercept=0.90, linetype='dashed') +
  labs(x="Probability of exceedance", y = "NPV (M$)") ; p2

#PDF of variables
variables <- samples_sim %>%
  select(Dev, Pop, Price, DRate, PcD, Prec = round(Prec), Temp, NVar, d) %>%
  gather(key = x, value = Value, -d) %>%
  filter(d == 1) %>% select(-d)

p3 <- ggplot(variables, aes(x=Value)) +
  theme_bw() +
  geom_bar(position="dodge") +
  facet_wrap(~ x, scales = "free") +
  labs(x="Level", y = "frequency"); p3


################################################################################
# DECISION CRITERIA ************************************************************
################################################################################


#cVaR calculations ********************
alpha_values <- seq(0,1,0.01)

cVaRs <- data_frame(Alphas = alpha_values, "40 MCM" = NA,  "60 MCM" = NA,
                    "80 MCM" = NA, "100 MCM" = NA, "120 MCM" = NA, "140 MCM" = NA)

for (i in 1:nrow(designs)) {

  data <- filter(samples_cdf, d == Design_labels[i])
  ind <- sapply(alpha_values, function(x) which.min(abs(data$ExP - x)):nrow(data))
  cVaRs[,i+1] <- sapply(ind, function(x) sum(data[x, "NPV"]) / length(x))
}

cVaRs_df <- cVaRs %>%
  gather(key = Designs, value = cVaR, -Alphas)

ggplot(cVaRs_df, aes(x=Alphas, y=cVaR, color=Designs)) +
  theme_bw() +
  geom_line(size=1)


################################################################################
# DECISION-SCALING PLOTS *******************************************************
################################################################################



# SURFACE PLOTS ----------------------------------------------------------------



#Calculate NPV under given values of socio-economic parameters
Price     <- 0.50   # 0.35 - 0.81 USD per cubic m
Drate     <- 0.05   # discount rate
LifeSpan  <- 40     # project life-time (max. 50)

#Simulated Yield: delta[cc] x delta[Variability]
syield <- read.csv('./inputs/safeyield_Feb06.csv') %>% as_data_frame()

#Simulated yield for: delta[P] x delta[T] x delta[Variability]
metrics_sim <- select(syield, c, d, yield) %>%
  arrange(d) %>%
  left_join(clim_space, by=c("c" ="Clim")) %>%
  select(Temp, Prec, NVar, c, d, yield)

#Calculate NPV under given socio-economic conditions
metrics_sim %<>% rowwise() %>%
  mutate(AnnB = yield * 0.5) %>%
  mutate(PVB  = PV_CALCULATE(rep(AnnB, 40), r=0.05)) %>%
  mutate(NPV  = PVB- designs$Cost[d])



# CHOOSE METRIC TO ANALYZE:
metric <- "NPV"   #yield

#Metrics as list of variability runs
metric_lst <- select(metrics_sim, c, d, NVar, metric = get(metric)) %>%
  spread(d, metric) %>%
  select(-c) %>%
  split(., .$NVar)

#Metrics averaged-over variability runs
metric_avg <- (Reduce(`+`, metric_lst) / length(metric_lst)) %>%
  as.data.frame() %>% as_data_frame() %>% select(-NVar)

#SURFACE PLOTS +++++++++++++


# BEST & WORST DESIGNS PLOT ++++++++++++++++++++++++++++++++++++++++++++++++++

criteria <- 'max'
Value <- if(criteria == 'min') {
  as.factor(apply(var_pl, 1, which.min))
} else {as.factor(apply(var_pl, 1, which.max))}

optim_df <- data.frame(var_xy, Value=Value)
levels(optim_df$Value) <- Design_labels

plot <- pl + geom_tile(data=optim_df, aes(x=Temp, y=Prec, fill=Value)) +
  theme_bw(base_size = 12) +
  geom_point(aes(x=Temp, y=Prec, color=Ensemble),
             data=dplyr::filter(gcms_data, Ensemble %in% c('CMIP3', 'CMIP5')), size=2) +
  scale_color_manual(values=c("red","black"), name="Climate\nProjections") +
  geom_hline(yintercept = 0, linetype ="dashed", size = 0.5) +
  scale_fill_manual(values = brewer.pal(6, "Set1")[3:6] , name="Designs") +
  guides(fill = guide_legend(ncol = 1))

ggsave(file="surface_2.png", path="C:\\Users\\umit\\GDrive",
       plot = plot, width=20, height=15, unit='cm')


# PARALLEL COORDINATES ---------------------------------------------------------

#INTERACTIVE PARALELL COORDINATES PLOTS **************************************

##Clean-up the data table for parcoords plot
output_var <- "NPV"

ParCoord_df <- NPV_sim %>%
  select (d, Dev, Price, DRate, NVar, Temp, Prec, get(output_var)) %>%
  distinct()

rChartsLib <- "C:/Users/Umit/Documents/R/win-library/3.1/rCharts/libraries/parcoords"

######## Prepare for each design

for (i in 1:nrow(designs))  {

  Data <- ParCoord_df %>%
    filter(d == i) %>%
    select(-d) %>%
    rbind(data_frame(Dev="Low", Price=NA, DRate=NA, NVar=NA,
                     Temp ="NA", Prec = "NA", NPV=max(ParCoord_df$NPV)))

  p1 <- rCharts$new()
  p1$setLib(rChartsLib)
  options(RCHART_WIDTH = 1400, RCHART_HEIGHT = 450)

  p1$set(padding = list(top = 24, left = 100, bottom = 12, right = 100),
         data    = toJSONArray(Data, json = F),
         colorby = output_var,
         range   = c(0,100),
         colors  = c('gold', 'darkred'))

  chart_name <- paste0("ParCoords_", i,".html")
  p1$save(chart_name, standalone=TRUE)

}


################################################################################









