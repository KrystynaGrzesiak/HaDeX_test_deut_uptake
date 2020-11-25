library(HaDeX)
library(dplyr)
library(tidyr)

dat_raw <- read_hdx(system.file(package = "HaDeX", "HaDeX/data/KD_180110_CD160_HVEM.csv"))


u <- function(x) {
  if(length(x) == 1) {
    0
  }else {
    sqrt(sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)/((length(x) - sum(is.na(x)))*(length(x) - 1 - sum(is.na(x)))))
  }
}



data = dat_raw %>% 
  filter(Sequence == "ARSQKSGIRLQGHF", State == "CD160")


calculate_err_deut_uptake <- function(data) {
  
  proton_mass <- 1.0072764668799
  
  df_pepMass = data %>% 
    complete(Exposure, nesting(Protein, Sequence, State, Start, End, MaxUptake, z)) %>% 
    mutate(pepMass = z*(Center - 1.0072764668799)) %>% 
    select(- RT, -Modification, -Fragment, -z)
  
  df_pepMass %>% 
    group_by(Protein, Sequence, State, Exposure) %>% 
    mutate(aggMass = weighted.mean(pepMass, Inten, na.rm = TRUE)) %>% 
    mutate(avgMass = mean(pepMass)) %>% 
    mutate(err_avgMass = u(pepMass)) %>% 
    select(-File, -pepMass, -Center, -Inten) %>% 
    unique() %>% 
    ungroup() %>% 
    group_by(Protein, Sequence, State) %>% 
    mutate(theo_deut_uptake = avgMass - MHP) %>% 
    mutate(theo_frac_deut_uptake = 100* (avgMass - MHP)/(MaxUptake * proton_mass)) %>% 
    mutate(err_theo_deut_uptake = abs(err_avgMass/(MaxUptake * proton_mass))) %>% 
    select(-MHP) %>% 
    unique() %>% 
    mutate(deut_uptake = avgMass - avgMass[Exposure == 0.001]) %>%
    mutate(frac_deut_uptake = 100*(avgMass - avgMass[Exposure == 0.001])/
             (avgMass[Exposure == max(Exposure, na.rm = TRUE)] - avgMass[Exposure == 0.001])) %>% 
    mutate(err_frac_deut_uptake_mt  = err_avgMass/
                                         (avgMass[Exposure == max(Exposure, na.rm = TRUE)] - avgMass[Exposure == 0.001])) %>%
    
    mutate(err_frac_deut_uptake_m0 = err_avgMass[Exposure == 0.001]*(avgMass - avgMass[Exposure == max(Exposure, na.rm = TRUE)])/
                                        ((avgMass[Exposure == max(Exposure, na.rm = TRUE)] - avgMass[Exposure == 0.001])^2)) %>% 
    
    mutate(err_frac_deut_uptake_m100 = err_avgMass[Exposure == max(Exposure, na.rm = TRUE)]*(avgMass[Exposure == 0.001] - avgMass)/
             ((avgMass[Exposure == max(Exposure, na.rm = TRUE)] - avgMass[Exposure == 0.001])^2)) %>% 
    
    mutate(err_frac_deut_uptake = sqrt(err_frac_deut_uptake_mt^2 + err_frac_deut_uptake_m0^2 + err_frac_deut_uptake_m100^2)) %>% 
    select(Protein, Sequence, State, frac_deut_uptake, err_frac_deut_uptake, deut_uptake, err_avgMass, Exposure, theo_deut_uptake, theo_frac_deut_uptake, err_theo_deut_uptake)
}




states_from_file <- unique(dat_raw[["State"]]) #reactive
chosen_protein <- unique(dat_raw[["Protein"]]) #input
in_time <- 0.001 #input
chosen_time <- 5 #input
out_time <- 1440 #input
deut_concentration <- 100 #input



prepared_dat <- bind_rows(lapply(states_from_file,
                                 function(i) calculate_state_deuteration(dat_raw,
                                                                         protein = chosen_protein,
                                                                         state = i,
                                                                         time_0 = in_time,
                                                                         time_t = chosen_time,
                                                                         time_100 = out_time,
                                                                         deut_part = 0.01*deut_concentration))) 

calculate_err_deut_uptake(dat_raw) %>% 
  filter(Exposure == 5, Sequence == "ARSQKSGIRLQGHF")


prepared_dat %>% 
  filter(Sequence == "ARSQKSGIRLQGHF")







