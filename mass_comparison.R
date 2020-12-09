library(HaDeX)
library(dplyr)
library(tidyr)
library(data.table)

dat <- read_hdx(system.file(package = "HaDeX", "HaDeX/data/KD_180110_CD160_HVEM.csv"))

dat_raw <- read.csv("C:/Users/pucha/Desktop/CD160_wszystko_statedata.csv") %>%
  mutate(Center.SD = as.numeric(Center.SD),
         Uptake.SD = as.numeric(Uptake.SD)) %>%
  as.data.frame()
        
str(dat_raw)
tail(dat_raw)

calculate_mass_no_inten(dat) %>% head()
calculate_mass_new(dat) %>% head()
calculate_mass_old(dat) %>% head()


##

dat_dyn_all <- dat_raw %>%
  select(Protein, Sequence, Start, End, MHP, MaxUptake, State, Exposure, Center, Center.SD) %>%
  mutate(source = "dyn") %>%
  rename(mass = Center, 
         err_mass = Center.SD) %>%
  select(Sequence, Start, End, State, Exposure, mass, err_mass, source)

dat_new_all <- calculate_mass_new(dat) %>%
  arrange(Start, End, State, Exposure) %>%
  mutate(source = "new") %>%
  select(Sequence, Start, End, State, Exposure, mass, err_mass, source)

dat_old_all <- calculate_mass_old(dat) %>%
  arrange(Start, End, State, Exposure) %>%
  mutate(source = "old") %>%
  select(Sequence, Start, End, State, Exposure, mass, err_mass, source)

dat_no_in_all <- calculate_mass_no_inten(dat) %>%
  arrange(Start, End, State, Exposure) %>%
  mutate(source = "no_in") %>%
  select(Sequence, Start, End, State, Exposure, mass, err_mass, source)

##

chosen_peptide_start <- 34
chosen_peptide_end <- 48

dat_dyn <- filter(dat_dyn_all, Start == chosen_peptide_start, End == chosen_peptide_end)
dat_new <- filter(dat_new_all, Start == chosen_peptide_start, End == chosen_peptide_end)
dat_old <- filter(dat_old_all, Start == chosen_peptide_start, End == chosen_peptide_end)
dat_no_in <- filter(dat_no_in_all, Start == chosen_peptide_start, End == chosen_peptide_end)

## old and new mass calculations
bind_rows(dat_new, dat_old) %>%
  gather(type, value, -Sequence, -Start, -End, -source, -mass, -State, -Exposure) %>%
  spread(source, value) %>%
  mutate(no_diff = new == old) 

bind_rows(dat_new, dat_dyn) %>%
  gather(type, value, -Sequence, -Start, -End, -source, -State, -Exposure) %>%
  select(State, Exposure, source, type, value) %>%
  spread(source, value) %>%
  mutate(diff = (dyn - new)) %>%
  arrange(type)


###################################################
# intensity and no intensity
res <- bind_rows(dat_new_all, dat_no_in_all) %>%
  gather(type, value, -Sequence, -Start, -End, -source, -State, -Exposure) %>%
  spread(source, value) %>%
  mutate(diff = (new - no_in)) %>%
  arrange(desc(type)) 

# sredni roznica w masie
mean(filter(res, type == "mass")[["diff"]]) # 0.004992344
# srednia roznica w odchyleniu
mean(filter(res, type == "err_mass")[["diff"]]) # -0.009250437
###################################################

str(dat_raw)
str(tmp)
str(dat_new)

tmp <- dat_raw %>% 
  select(Sequence, Start, End, State, Exposure, Center, Center.SD, Uptake, Uptake.SD, RT, RT.SD)

# porownanie wartsci z drugiego liku dynamxowego i z naszymi policzonymi masami
merge(x = tmp, y = dat_new_all, by = c("Sequence", "Start", "End", "State", "Exposure"), all = T) %>%
  select(-source) %>%
  arrange(Start, End, State) %>%
  filter(!is.na(mass))
  