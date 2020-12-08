# Adjusting the current code to fit the test result

library(dplyr)
library(tidyr) # spread

reference_table

##
dat_raw <- read_hdx(system.file(package = "HaDeX", "HaDeX/data/KD_180110_CD160_HVEM.csv"))
proton_mass <- 1.00727647
time_0 <- 0.001
time_t <- 5
time_100 <- 1440
sequence <-  "FTISQVTPLHSGT"
protein <- "db_CD160"
state <- "CD160"
deut_part <- 1
##

dat <- dat_raw %>%
  filter(Sequence == sequence,
         Protein == protein,
         State == state,
         Exposure %in% c(time_0, time_t, time_100))

##

dat %>%
  mutate(exp_mass = Center*z - z*proton_mass) %>%
  select(-Center, -z, -Modification, -Fragment) %>%
  # mutate(source = "tmp") %>%
  # as.data.frame() %>%
  # results(reference_table_exp_mass)
  group_by(Sequence, Start, End, MHP, MaxUptake, State, Exposure, Protein, File) %>%
  summarize(avg_exp_mass = weighted.mean(exp_mass, Inten, na.rm = TRUE)) %>%
  ungroup(.) %>% 
  mutate(Exposure = case_when(Exposure == time_0 ~ "time_0",
                              Exposure == time_t ~ "time_t",
                              Exposure == time_100 ~ "time_100")) %>%
  spread(key = Exposure, value = avg_exp_mass) %>%
  group_by(Sequence, Start, End, MaxUptake, MHP, Protein, State) %>%
  summarize(m_0 = mean(time_0, na.rm = TRUE),
            err_m_0 = coalesce(sd(time_0, na.rm = TRUE)/sqrt(sum(!is.na(time_0))), 0),
            m_t = mean(time_t, na.rm = TRUE),
            err_m_t = coalesce(sd(time_t, na.rm = TRUE)/sqrt(sum(!is.na(time_t))), 0),
            m_100 = mean(time_100, na.rm = TRUE),
            err_m_100 = coalesce(sd(time_100, na.rm = TRUE)/sqrt(sum(!is.na(time_100))), 0)) %>%
  ##
  # mutate(source = "tmp") %>%
  # as.data.frame() %>%
  # select(Sequence, m_0, err_m_0, m_t, err_m_t, m_100, err_m_100, source) %>%
  # results(reference_table_mass)
  # 
  ## dotad ok
  # odtad nieok
  mutate(# experimental calculations below - fractional
    frac_deut_uptake = 100*(m_t - m_0)/(m_100 - m_0),
    err_frac_deut_uptake = 100*sqrt((err_m_t/(m_100 - m_0))^2 + ((m_t - m_100)*err_m_0/(m_100 - m_0^2))^2 + ((m_0 - m_t)*err_m_100/(m_100 - m_0^2))^2),
    # experimental calculations below - absolute
    deut_uptake = (m_t - m_0),
    err_deut_uptake = sqrt(err_m_t^2 + err_m_0^2),
    # theoretical calculations below - fractional
    theo_frac_deut_uptake  = 100*(m_t - MHP)/(MaxUptake * proton_mass * deut_part),
    err_theo_frac_deut_uptake  = 100*err_m_t*(1/(MaxUptake * proton_mass * deut_part)),
    # theoeretical calculations below - absolute
    theo_deut_uptake = (m_t - MHP),
    err_theo_deut_uptake = err_m_t,
    # helper values
    Med_Sequence = Start + (End - Start)/2) %>%
  ungroup(.) %>%
  arrange(Start, End) %>%
  mutate(source = "tmp") %>%
  as.data.frame() %>%
  select(Sequence, frac_deut_uptake, err_frac_deut_uptake, 
         deut_uptake, err_deut_uptake, 
         theo_frac_deut_uptake, err_theo_frac_deut_uptake,
         theo_deut_uptake, err_theo_deut_uptake, source) %>%
  results(reference_table_deut_uptake)


###
mutate(source = "tmp") %>%
  as.data.frame() %>%
  results(reference_table_mass)

reference_table
tmp

results <- function(...){
  bind_rows(...) %>%
    gather(type, value, -Sequence, -source)  %>%
    spread(source, value) %>%
    mutate(delta_tmp = test - tmp) %>%
    arrange(delta_tmp)
}

results(reference_table, tmp)

####


  select(Sequence, m_0, err_m_0, m_t, err_m_t, m_100, err_m_100, source) 
  
  
  mutate(# experimental calculations below - fractional
    frac_deut_uptake = 100*(m_t - m_0)/(m_100 - m_0),
    err_frac_deut_uptake = 100*sqrt((err_m_t*(1/(m_100 - m_0)))^2 + (err_m_0*((m_t - m_100 )/((m_100 - m_0)^2)))^2 + (err_m_100*((m_0 - m_t)/((m_100 - m_0)^2)))^2),
    # experimental calculations below - absolute
    deut_uptake = (m_t - m_0),
    err_deut_uptake = sqrt(err_m_t^2 + err_m_0^2),
    # theoretical calculations below - fractional
    theo_frac_deut_uptake  = 100*(m_t - MHP)/(MaxUptake * proton_mass * deut_part),
    err_theo_frac_deut_uptake  = 100*err_m_t*(1/(MaxUptake * proton_mass * deut_part)),
    # theoeretical calculations below - absolute
    theo_deut_uptake = (m_t - MHP),
    err_theo_deut_uptake = err_m_t,
    # helper values
    Med_Sequence = Start + (End - Start)/2)

