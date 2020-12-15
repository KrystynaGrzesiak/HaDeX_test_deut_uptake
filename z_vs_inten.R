## analysing the relation between z and Inten

library(HaDeX)
library(dplyr)
library(ggplot2)

dat_raw <- read_hdx(system.file(package = "HaDeX", "HaDeX/data/KD_180110_CD160_HVEM.csv"))

proton_mass <- 1.00727647

# state CD160
dat <- dat_raw %>%
  filter(State == "CD160") %>%
  select(-Protein, -Modification, -Fragment, -MHP, -RT) %>%
  mutate(Mass = z*(Center - proton_mass),
         Length = nchar(Sequence)) %>%
  select(-Start, -End) %>%
  select(Length, everything())

cor(dat[["z"]], dat[["Inten"]])

cor(dat_raw[["z"]], dat_raw[["Inten"]])

ggplot(dat, aes(x = z, y = Inten)) +
  geom_point()


dat %>% 
  group_by(Sequence, Exposure, z) %>%
  summarise(length = n(),
            mean_Inten = mean(Inten),
            err_mean_Inten = sd(Inten)) %>%
  filter(Exposure>0.001) %>%
  arrange(desc(nchar(Sequence)))
  