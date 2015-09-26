
#-----------------------------------------------------------------------------#
# Title: Decision analysis model for Mwache Dam design
# By:    M.Umit Taner
# Date:  September, 2015                                                       
#-----------------------------------------------------------------------------# 
  
 










# Decision Analysis ------------------------------------------------------------


################################################################################


#Create Sample table for all designs
samples <- as.data.frame(LHS_samples) %>%
  as_data_frame() %>%
  replicate(6, ., simplify=FALSE) %>%
  bind_rows() %>%
  mutate(d = rep(1:6, each = nrow(.)/6))

#Frequency analysis
#dat <- gather(LHS_samples, key = Par, value = value, - Sample)
#p <- ggplot(data = dat, aes(x=value)) +
#  geom_histogram() +
#  facet_wrap( ~ Par, scales = "free")

# STEP 3: DERIVE FULL JOINT-DISTRIBITIONS (CDFs) -------------------------------

#Simulated Yield: delta[cc] x delta[Variability]
syield <- read.csv('./inputs/safeyield_Feb06.csv') %>% as_data_frame()

#Look up table for climate parameters
clim_index <- clim_space %>% as_data_frame() %>%
  select(Clim, CC, NVar, Temp, Prec)

#Simulated yield for: delta[P] x delta[T] x delta[Variability]
syield_sim <- select(syield, Clim=c, d, yield) %>%
  arrange(d) %>%
  left_join(clim_index, "Clim") %>%
  select(d, Clim, CC, Temp, Prec, NVar, yield) %>%
  mutate(d = as.integer(d), CC = as.integer(CC), NVar = as.integer(NVar))

#Sample_results table
samples_df <- samples %>%
  left_join(syield_sim, by = c("CC", "NVar", "d"))

#Calculate NPV for each sample:
samples_sim <- samples_df %>%
  rowwise() %>%
  mutate(
    AnnB = yield * Price,
    PVB  = PV_CALCULATE(rep(AnnB, 40), r=DRate),
    NPV  = PVB - designs$Cost[d],
    Demand = Pop * PcD * (10**-9) *365) %>%
  select(-AnnB, -PVB)

samples_cdf <- ungroup(samples_sim) %>%
  as_data_frame() %>%
  mutate(d = ordered(.$d, labels = Design_labels)) %>%
  group_by(d) %>%
  arrange(desc(NPV)) %>%
  mutate(Index =  row_number(),
         ExP = Index/(max(Index)+1))


################################################################################




### Density function
p1 <- ggplot(samples_cdf, aes(x=NPV, group=d, color=d)) +
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


# PLOT SETTINGS ----------------------------------------------------------------

#Plotting variables
alabels <- c(expression(paste("Mean Temp. Increase (", ~degree~C, ")")),
             paste("Mean Precip. Change (%)"))
aticks <- list(x = seq(0,2.5,0.5), y = seq(-30,40,10))

#Color & shape schemes
color_2 <- c("#0072B2", "#D55E00")
color_4 <- c('#000000', '#666666', "#D55E00","#0072B2")
ColorQual6 <- brewer.pal(6, "Spectral")

shape_4 <- c(1,0,3,8)
shape_2 <- c(1,3,8)
shape_2 <- c(3,16)

#Template
pl <- ggplot(mapping=aes(x,y)) +
  scale_x_continuous(expand=c(0,0), limits=range(aticks$x), breaks=aticks$x) +
  scale_y_continuous(expand=c(0,0), limits=range(aticks$y), breaks=aticks$y) +
  labs(x=alabels[1], y=alabels[2]) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(face="bold"),
    legend.background = element_rect(fill=NA),
    panel.margin = unit(2, "lines"),
    panel.border = element_rect(color="black", fill=NA),
    strip.background = element_rect(fill="gray90", color='black'))


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

#Intpolate for surface plots
res <- 100
metric_to_plot <- GRID_INTERPOLATE(
  x=cc_space$Temp_Del,
  y=cc_space$Prec_Del,
  z=metric_avg, resolution=res) %>%
  as_data_frame()

#Set metric range & bins
var_xy <- metric_to_plot %>% select(Temp = x, Prec = y)
var_pl <- as_data_frame(metric_to_plot) %>% select(-x, -y)
colnames(var_pl) <- 1:nrow(designs)

#Date-frame for plotting
metric_df <- cbind(var_xy, var_pl) %>% as_data_frame() %>%
  gather(Designs, Value, -Temp, -Prec) %>%
  dplyr::filter(Designs %in% c(1,3,5,6))

#Bins and colors
var_bin <- pretty(c(max(var_pl), min(var_pl)), n=6)
var_col <- colorRampPalette(brewer.pal(6, "Blues"))(length(var_bin))
metric_df$Bins <- cut(metric_df$Value, breaks=var_bin,
                      dig.lab=5, include.lowest = TRUE)
levels(metric_df$Designs) <- paste0("Design Capacity: ", Design_labels)

#Climate Surface Plot
plot <- pl + geom_tile(data=metric_df, aes(x=Temp, y=Prec, fill=Bins)) +
  theme_bw(base_size = 12) +
  facet_wrap( ~  Designs) +
  geom_point(aes(x=Temp, y=Prec, color=Ensemble),
             data=dplyr::filter(gcms_data, Ensemble %in% c('CMIP3', 'CMIP5')), size=2) +
  scale_fill_manual(values=var_col, name='NPV\n(Million USD)') +
  #scale_fill_manual(values=var_col, name='Guaranteed supply\n(MCM/year)') +
  scale_color_manual(values=c("red","black"), name="Climate\nProjections") +
  geom_hline(yintercept = 0, linetype ="dashed", size = 0.5) +
  guides(fill = guide_legend(ncol = 1, order = 1), color = guide_legend(order=2)) +
  theme(panel.margin = unit(1.5, "lines"))

ggsave(file="surface_1.png", path="C:\\Users\\umit\\GDrive",
       plot = plot, width=18, height=15, unit='cm')

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









