


# CLIMATE DATA CLEANING --------------------------------------------------------


#OBSERVED CLIMATE DATA (Monthly, from 10 rain gauge stations)

obs_prec_long <- read.csv("./Inputs/Obs_monthly_precip_mm.csv") %>%
  gather(key = Month, value = value, -Year, -Station_ID, -Station_Name) %>%
  mutate(Month = as.vector(sapply(as.character(Month), function(x) which(x == month.abb)))) %>%
  transform(Date = as.Date(paste(Year, Month, 1, sep="-"))) %>%
  select(Station_ID, Station_Name, Date, value) %>%
  arrange(Station_ID, Date) %>%
  as_data_frame()

#GRIDDED CLIMATE DATA
GPCC_3.25 <- read.table("./Repository/Climate/GPCC/GPCC_-3.25_39.25") %>%
  as_data_frame() %>%
  select(Year = V1, Month = V2, Value = V3) %>%
  filter(Year %in% 1950:1999) %>%
  mutate(Source = "GPCC_3.25")

GPCC_3.75 <- read.table("./Repository/Climate/GPCC/GPCC_-3.25_39.75") %>%
  as_data_frame() %>%
  select(Year = V1, Month = V2, Value = V3) %>%
  filter(Year %in% 1950:1999) %>%
  mutate(Source = "GPCC_3.75")

CRU_3.25 <- read.table("./Repository/Climate/CRU/CRU_-3.25_39.25") %>%
  as_data_frame() %>%
  select(Year = V1, Month = V2, Value = V3) %>%
  filter(Year %in% 1961:1999) %>%
  mutate(Source = "CRU_3.25")

CRU_3.75 <- read.table("./Repository/Climate/CRU/CRU_-3.25_39.75") %>%
  as_data_frame() %>%
  select(Year = V1, Month = V2, Value = V3) %>%
  filter(Year %in% 1961:1999) %>%
  mutate(Source = "CRU_3.75")

grid_prec_long <- bind_rows(GPCC_3.25, GPCC_3.75, CRU_3.25, CRU_3.75) %>%
  spread(key = Source, value = Value) %>%
  mutate(Princeton = hist_forcings$Prec, Date = hist_dates) %>%
  transmute(Date, CRU = CRU_3.75, GPCC = GPCC_3.75, Princeton) %>%
  gather(key = Source, value = value, -Date)


#write.csv(x = grid_prec_long, file = "./inputs/grid_prec_long_mm.csv", row.names=FALSE)
#write.csv(x = obs_prec_long,  file = "./inputs/obs_prec_long_mm.csv", row.names=FALSE)


################################################################################
################################################################################
################################################################################



# STAKEHOLDER DATA -------------------------------------------------------------


  #Packages to load
  require(readxl)
  require(dplyr)
  require(tidyr)
  require(ggplot2)

  p <- list()

  belief_df <- read_excel("inputs/Workshop_survey_v2.xlsx", sheet = 1) %>%
    gather(key = Particip, value = value, x1:x28) %>%
    mutate(value = as.factor(.$value))

  for (i in 1:12) {

    plot_name <- unique(belief_df$Uncertainty)[i]
    plot_data <- filter(belief_df, Uncertainty == plot_name)

  p[[i]] <- ggplot(plot_data, aes(x = value, fill = Factor)) +
    ggtitle(plot_name) +
    geom_bar() +
    facet_grid(.~ Factor) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    labs(x="Significance score", y = "count (out of 28)") +
    theme(plot.title = element_text(face="bold"),
          strip.text = element_text(size = 6)) +
    guides(fill = FALSE)

    ggsave(paste0("factor_",i,".png"), p[[i]], width = 8, height = 4, units = "in")

  }




################################################################################
################################################################################
################################################################################
