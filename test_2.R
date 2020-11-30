# Test code from meeting, 30/11/2020

library(HaDeX)
library(dplyr)

options(digits=10)

dat_raw <- read_hdx(system.file(package = "HaDeX", "HaDeX/data/KD_180110_CD160_HVEM.csv"))

unique(dat_raw[["Sequence"]])
unique(dat_raw[["State"]])
unique(dat_raw[["Exposure"]])

dat <- dat_raw %>%
  filter(Sequence == "FTISQVTPLHSGT", 
         State == "CD160",
         Exposure %in% c(0.001, 5, 1440))

dat_tmp <- dat %>% select(-Protein, -Start, -End, -Modification, -Fragment, -RT, -State) 

time_0 <- 0.001
time_t <- 5
time_100 <- 1440
Sequence <-  "FTISQVTPLHSGT"

proton_mass <- 1.00727647

tmp_mass <- 
  dat_tmp %>%
  mutate(expMass = z*(Center-proton_mass)) %>%
  group_by(Sequence, MaxUptake, MHP, Exposure, File) %>%
  mutate(aggMass = weighted.mean(expMass, Inten, na.rm = TRUE)) %>%
  ungroup(.) %>%
  select(-File, -z, -Inten, -Center, -expMass) %>%
  unique(.) %>%
  mutate(aggMass = round(aggMass, 5)) %>% 
  as.data.frame()

agg_tmp_mass <- tmp_mass %>%
  group_by(Sequence, MaxUptake, MHP, Exposure) %>%
  summarise(pepMass = mean(aggMass, na.rm = TRUE),
            # err_pepMass = coalesce(sd(aggMass, na.rm = TRUE), 0)) %>% ## do poprawy
            err_pepMass = coalesce(sd(aggMass, na.rm = TRUE), 0)/4) %>% ## do poprawy
  as.data.frame()

agg_tmp_mass_old <- agg_tmp_mass

m_0 <- agg_tmp_mass[1, 5]
m_t <- agg_tmp_mass[2, 5]
m_100 <- agg_tmp_mass[3, 5]

err_m_0 <- agg_tmp_mass[1, 6]
err_m_t <- agg_tmp_mass[2, 6]
err_m_100 <- agg_tmp_mass[3, 6]

MaxUptake <- 11
MHP <- unique(dat_tmp[["MHP"]])

deut_uptake <- m_t - m_0 # Da
err_deut_uptake <- sqrt(err_m_t^2 + err_m_0^2)

frac_deut_uptake <- 100*(m_t - m_0)/(m_100 - m_0) # %
max_up <- m_100 - m_0
x_1 <- err_m_t/(max_up)
x_2 <- (m_t - m_100)*err_m_0/(max_up^2)
x_3 <- (m_0 - m_t)*err_m_100/(max_up^2)
err_frac_deut_uptake <- 100*sqrt(x_1^2 + x_2^2 + x_3^2)

theo_deut_uptake <- m_t - MHP
err_theo_deut_uptake <- err_m_t

theo_frac_deut_uptake <- 100*(m_t - MHP)/(MaxUptake*proton_mass)
err_theo_frac_deut_uptake <- 100*err_m_t/(MaxUptake*proton_mass)

result <- data.frame(Sequence, deut_uptake, err_deut_uptake, frac_deut_uptake, err_frac_deut_uptake, theo_deut_uptake, err_theo_deut_uptake, theo_frac_deut_uptake, err_theo_frac_deut_uptake) %>% mutate(source = "test") %>% as.data.frame()

res_fun <- calculate_state_deuteration(dat,
                                       protein = "db_CD160",
                                       state = "CD160",
                                       time_0 = time_0,
                                       time_t = time_t,
                                       time_100 = time_100,
                                       deut_part = 1) %>%
  select(-Protein, -Start, -End, -State, -Med_Sequence) %>%
  as.data.frame() %>%
  mutate(source = "function")


str(result)
str(res_fun)
bind_rows(result, res_fun)

bind_rows(x_1, result, res_fun)
