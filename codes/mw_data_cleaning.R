

# GRIDDED DATA CLEANING --------------------------------------------------------

#OBSERVED CLIMATE DATA (Monthly, from 10 rain gauge stations)
obs_prec_long <- read.csv("./Inputs/Obs_monthly_precip_mm.csv") %>%
  gather(key = Month, value = value, -Year, -Station_ID, -Station_Name) %>%
  mutate(Month = as.vector(sapply(as.character(Month),
    function(x) which(x == month.abb)))) %>%
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

precip_obs_grid <- bind_rows(GPCC_3.25, GPCC_3.75, CRU_3.25, CRU_3.75) %>%
  spread(key = Source, value = Value) %>%
  mutate(Princeton = hist_forcings$Prec, Date = hist_dates) %>%
  transmute(Date, CRU = CRU_3.75, GPCC = GPCC_3.75, Princeton) %>%
  gather(key = Source, value = value, -Date)


#write.csv(x = precip_obs_grid,
#  file = "./inputs/precip_obs_grid_mm.csv", row.names=FALSE)
#write.csv(x = obs_prec_long,
#file = "./inputs/obs_prec_long_mm.csv", row.names=FALSE)




# # GCM DATA PREPARATION --------------------------------------------------


#Historical period: 1950 - 1999
#Projection period: 2055 - 2085

prec_hist <- read_excel("Inputs/GCM_deltaclim.xlsx", "pr_historical")
temp_hist <- read_excel("Inputs/GCM_deltaclim.xlsx", "tas_historical")
prec_rcp26 <- read_excel("Inputs/GCM_deltaclim.xlsx", "pr_rcp26")
temp_rcp26 <- read_excel("Inputs/GCM_deltaclim.xlsx", "tas_rcp26")
prec_rcp45 <- read_excel("Inputs/GCM_deltaclim.xlsx", "pr_rcp45")
temp_rcp45 <- read_excel("Inputs/GCM_deltaclim.xlsx", "tas_rcp45")
prec_rcp60 <- read_excel("Inputs/GCM_deltaclim.xlsx", "pr_rcp60")
temp_rcp60 <- read_excel("Inputs/GCM_deltaclim.xlsx", "tas_rcp60")
prec_rcp85 <- read_excel("Inputs/GCM_deltaclim.xlsx", "pr_rcp85")
temp_rcp85 <- read_excel("Inputs/GCM_deltaclim.xlsx", "tas_rcp85")

prec_hist_mean  <- prec_hist %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(hist = mean(value))

temp_hist_mean  <-  temp_hist %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(hist = mean(value))

prec_rcp26_mean <- prec_rcp26 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp26 = mean(value) * 0.01)

temp_rcp26_mean <- temp_rcp26 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp26 = mean(value))

prec_rcp45_mean <- prec_rcp45 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp45 = mean(value) * 0.01)

temp_rcp45_mean <- temp_rcp45 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp45 = mean(value))

prec_rcp60_mean <- prec_rcp60 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp60 = mean(value) * 0.01)

temp_rcp60_mean <- temp_rcp60 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp60 = mean(value))

prec_rcp85_mean <- prec_rcp85 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp85 = mean(value) * 0.01)

temp_rcp85_mean <- temp_rcp85 %>% gather(key = Mon, value = value, -GCM) %>%
  group_by(GCM) %>% summarize(rcp85 = mean(value))

GCM_temp <- temp_hist_mean %>%
  left_join(temp_rcp26_mean, by = "GCM") %>%
  left_join(temp_rcp45_mean, by = "GCM") %>%
  left_join(temp_rcp60_mean, by = "GCM") %>%
  left_join(temp_rcp85_mean, by = "GCM") %>%
  gather(key = scenario, value = temp, -GCM)

GCM_prec <- prec_hist_mean %>%
  left_join(prec_rcp26_mean, by = "GCM") %>%
  left_join(prec_rcp45_mean, by = "GCM") %>%
  left_join(prec_rcp60_mean, by = "GCM") %>%
  left_join(prec_rcp85_mean, by = "GCM") %>%
  gather(key = scenario, value = prec, -GCM)

GCM_means <- GCM_temp %>%
  left_join(GCM_prec, by = c("GCM","scenario"))


#write_csv(x = GCM_means, path = "CMIP5_.csv")

