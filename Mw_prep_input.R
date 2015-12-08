

### PREPARE IRRIGATION DEMAND --------------------------------------------------

  require(sirad)

  #Crop efficiency and Kc values for each crop type
  crop_data <- read_csv("./inputs/crop_data_irrig.csv", skip = 1)

  #Monthly parameters for Penmann - Montheith calculation
  crop_pars_et0 <- read_csv("./inputs/crop_pars_et0.csv", skip = 1)

  #Crop coefficients
  crop_Kc <- crop_data %>%select(Season, Crop = c, Jun:May) %>%
    gather(key = Month, value = Kc, -Season, -Crop) %>%
    rowwise() %>%
    mutate(Month = which(Month == month.abb)) %>%
    mutate(Month = as.numeric(Month))

 #Crop area & efficiency
  crop_area <- crop_data %>%
    select(Season, Crop = c, Area_ha, Method) %>%
    mutate(Eff = ifelse(Method == "sprinkler", par$irr_spr_eff, par$irr_dip_eff)) %>%
    select(-Method)

  #Calculate ET0 (Penmann - Monteith equation (FAO))
  forcing <- hydroclim_i %>%
    mutate(Year = as.numeric(year(Date)), Month = as.numeric(month(Date))) %>%
    left_join(crop_pars_et0, by = "Month") %>%
    mutate(
      EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
      Wind  = Wind * 0.01157407,
      v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
      vap_p = Humid/100 * v_sat,
      ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
        tal = cst(Solar, Date, radians(Latitude)),z = 150, uz = Wind,
        days  = Date, lat = Latitude) * as.numeric(days_in_month(Month)))

  #Calculate CRW (mm/month)
  irr_per_crop <- forcing %>%
    left_join(crop_Kc, by = "Month") %>%
    left_join(crop_area, by = c("Season","Crop")) %>%
    mutate(
      ETc = ET0 * Kc,
      NIR_mm = pmax(ETc - EPrec, 0),
      GIR_mm = NIR_mm / Eff,
      GIR_mcm = GIR_mm * Area_ha * par$irr_scale / 10^5)

  # Monthly total irrigation demand
  demand.m <- irr_per_crop %>%
    group_by(Year, Month) %>%
    summarize(Irr_mcm = sum(GIR_mcm)) %>%
    ungroup() %>%
    mutate(Year = as.numeric(Year), Month = as.numeric(Month))

  #DOMESTIC
  dom_y <- par$dom_baseline * par$dom_scale * (1 - par$dom_recycle)
  demand.m %<>% mutate(Dom_mcm = dom_y * coeff.m[Month])

  #ENVIRONMENTAL
  demand.m %<>% mutate(Eco_mcm = par$eco_baseline[Month] * par$eco_scale)

# FINAL DATA --------------------------------------------------------------

  input <- forcing %>%
    left_join(demand.m, by=c("Year","Month")) %>%
    select(key, Year, Month, Prec:Tmax, PE, Q_mcm, ET0, Irr_mcm, Dom_mcm, Eco_mcm)





