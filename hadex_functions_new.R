## new versions of calculate_state_deuteration, with small
## variations. this is used for comparison of methods - the 
## results are in show_differences.R

calculate_state_deuteration_new <- function(dat,
                                        protein, 
                                        state, 
                                        time_0,
                                        time_t, 
                                        time_100,
                                        deut_part = 1){
  proton_mass <- 1.00727647
  dat <- dat[dat[["Protein"]] == protein & dat[["State"]] == state & dat[["Exposure"]] %in% c(time_0, time_t, time_100), ]
  
  dat %>%
    mutate(exp_mass = Center*z - z*proton_mass) %>%
    select(-Center, -z, -Modification, -Fragment) %>%
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
    mutate(# experimental calculations below - fractional
      frac_deut_uptake = 100*(m_t - m_0)/(m_100 - m_0),
      err_frac_deut_uptake = 100*sqrt((err_m_t/(m_100-m_0))^2 + (err_m_0*(m_t-m_100)/((m_100-m_0)^2))^2 + (err_m_100*(m_0-m_t)/((m_100-m_0)^2))^2),
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
    select(Protein, Sequence, Start, End, State, 
           frac_deut_uptake, err_frac_deut_uptake, 
           deut_uptake, err_deut_uptake, 
           theo_frac_deut_uptake, err_theo_frac_deut_uptake,
           theo_deut_uptake, err_theo_deut_uptake, 
           Med_Sequence)
  
}


generate_differential_data_set_new <- function(dat,
                                               states,
                                               protein,
                                               time_0,
                                               time_t,
                                               time_100,
                                               deut_part = 1){
  
  bind_rows(lapply(states, function(i) calculate_state_deuteration_new(dat, 
                                                                       protein = protein, 
                                                                       state = i, 
                                                                       time_0 = time_0,
                                                                       time_t = time_t, 
                                                                       time_100 = time_100,
                                                                       deut_part = deut_part))) %>%
    droplevels() %>% 
    mutate(State = factor(State, levels = states, labels = c("1", "2"))) %>%
    gather(variable, value, -c(Protein:End, State, Med_Sequence)) %>%
    unite(tmp, variable, State) %>%
    spread(tmp, value)  %>%
    mutate(diff_frac_deut_uptake = frac_deut_uptake_1 - frac_deut_uptake_2,
           err_diff_frac_deut_uptake = sqrt(err_frac_deut_uptake_1^2 + err_frac_deut_uptake_2^2),
           diff_deut_uptake = deut_uptake_1 - deut_uptake_2,
           err_diff_deut_uptake = sqrt(err_deut_uptake_1^2 + err_deut_uptake_2^2),
           diff_theo_frac_deut_uptake = theo_frac_deut_uptake_1 - theo_frac_deut_uptake_2, 
           err_diff_theo_frac_deut_uptake = sqrt(err_theo_frac_deut_uptake_1^2 + err_theo_frac_deut_uptake_2^2),
           diff_theo_deut_uptake = theo_deut_uptake_1 - theo_deut_uptake_2,
           err_diff_theo_deut_uptake = sqrt(err_theo_deut_uptake_1^2 + err_theo_deut_uptake_2^2)) %>%
    arrange(Start, End) %>%
    select(Protein, Start, End, Med_Sequence, everything(), -contains("1"), -contains("2"))
}

