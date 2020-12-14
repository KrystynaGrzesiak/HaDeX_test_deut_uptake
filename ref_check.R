library(readr)
library(HaDeX)
library(dplyr)
library(tidyr)

dat <- read_hdx(system.file(package = "HaDeX", "HaDeX/data/KD_180110_CD160_HVEM.csv"))
ref_dat <- read_csv("C:/Users/pucha/Desktop/ref_table.csv")

head(dat)

chosen_peptide <- "INITSSASQEGTRLN"
chosen_protein <- "db_CD160"
chosen_state <- "CD160"
as.data.frame(ref_dat)

times <- unique(dat[ dat[["Exposure"]] > 0.001 & dat[["Exposure"]] != 1440, ][["Exposure"]])


fun_dat <- lapply(times, function(t){calculate_state_deuteration(dat = dat,
                                                                 protein = chosen_protein, 
                                                                 state = chosen_state, 
                                                                 time_0 <- 0.001,
                                                                 time_t <- t,
                                                                 time_100 <- 1440) %>%
    mutate(Exposure = t)
}) %>% bind_rows() %>%
  filter(Sequence == chosen_peptide) %>%
  mutate(source = "func") %>%
  select(Exposure, everything(), -Protein, -Start, -End, -State) %>%
  select(Exposure, deut_uptake, err_deut_uptake, frac_deut_uptake, err_frac_deut_uptake, source)

ref_dat <- ref_dat %>%
  filter(Exposure > 0.001, Exposure != 1440) %>%
  select(-Mass, -sd) %>%
  mutate(source = "ref")

options(digits=10)

bind_rows(fun_dat, ref_dat) %>%
  gather(type, value, -source, -Exposure)  %>%
  spread(source, value) %>%
  as.data.frame() %>%
  mutate(func = round(func, 7),
         ref = round(ref, 7)) %>%
  mutate(is_diff = !(ref == func)) %>%
  filter(is_diff)

# so it is ok

