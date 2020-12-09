calculate_mass_new <- function(dat){
  
  proton_mass <- 1.00727647
  
  dat %>%
    mutate(exp_mass = Center*z - z*proton_mass) %>%
    select(-Center, -z, -Modification, -Fragment) %>%
    group_by(Sequence, Start, End, MHP, MaxUptake, State, Exposure, Protein, File) %>%
    summarize(avg_exp_mass = weighted.mean(exp_mass, Inten, na.rm = TRUE)) %>%
    ungroup(.) %>% 
    group_by(Sequence, Start, End, MHP, MaxUptake, State, Exposure, Protein) %>%
    summarize(mass = mean(avg_exp_mass, na.rm = TRUE),
              err_mass = coalesce(sd(avg_exp_mass, na.rm = TRUE)/sqrt(sum(!is.na(avg_exp_mass))), 0),
              num = (sum(!is.na(avg_exp_mass)))) %>%
    ungroup(.) %>%
    arrange(Start, End, Start, Exposure) %>%
    as.data.frame()
  
}

calculate_mass_old <- function(dat){
  
  proton_mass <- 1.00727647
  
  dat %>%
    mutate(exp_mass = Center*z - z*proton_mass) %>%
    select(-Center, -z, -Modification, -Fragment) %>%
    group_by(Sequence, Start, End, MHP, MaxUptake, State, Exposure, Protein, File) %>%
    summarize(avg_exp_mass = weighted.mean(exp_mass, Inten, na.rm = TRUE)) %>%
    ungroup(.) %>% 
    group_by(Sequence, Start, End, MHP, MaxUptake, State, Exposure, Protein) %>%
    summarize(mass = mean(avg_exp_mass, na.rm = TRUE),
              err_mass = coalesce(sd(avg_exp_mass, na.rm = TRUE)/sqrt(length(avg_exp_mass)), 0),
              num = length(avg_exp_mass)) %>%
    ungroup(.) %>%
    arrange(Start, End, Start, Exposure) %>%
    as.data.frame()
  
}

calculate_mass_no_inten <- function(dat){
  
  proton_mass <- 1.00727647
  
  dat %>%
    mutate(exp_mass = Center*z - z*proton_mass) %>%
    select(-Center, -z, -Modification, -Fragment) %>%
    group_by(Sequence, Start, End, MHP, MaxUptake, State, Exposure, Protein, File) %>%
    summarize(avg_exp_mass = mean(exp_mass, na.rm = TRUE)) %>%
    ungroup(.) %>% 
    group_by(Sequence, Start, End, MHP, MaxUptake, State, Exposure, Protein) %>%
    summarize(mass = mean(avg_exp_mass, na.rm = TRUE),
              err_mass = coalesce(sd(avg_exp_mass, na.rm = TRUE)/sqrt(sum(!is.na(avg_exp_mass))), 0)) %>%
    ungroup(.) %>%
    arrange(Start, End, Start, Exposure) %>%
    as.data.frame()
  
}

