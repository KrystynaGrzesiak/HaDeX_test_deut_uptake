show_differences(dat = dat,  time_t = 5)

show_differences <- function(dat, 
                             states = c("CD160", "CD160_HVEM"),
                             protein = "db_CD160",
                             time_0 = 0.001,
                             time_t,
                             time_100 = 1440,
                             confidence_limit = 0.98,
                             confidence_limit_2 = 0.99,
                             deut_part = 1){
  
  ## new
  
  woods_data_set_new <- generate_differential_data_set_new(dat = dat,
                                     states = states, 
                                     protein = protein,
                                     time_0 = time_0,
                                     time_t = time_t,
                                     time_100 = time_100,
                                     deut_part = deut_part)
  avg_err_diff_frac_deut_uptake_new <- mean(woods_data_set_new[["err_diff_frac_deut_uptake"]])
  avg_err_diff_deut_uptake_new <- mean(woods_data_set_new[["err_diff_deut_uptake"]])
  avg_err_diff_theo_frac_deut_uptake_new <- mean(woods_data_set_new[["err_diff_theo_frac_deut_uptake"]])
  avg_err_diff_theo_deut_uptake_new <- mean(woods_data_set_new[["err_diff_theo_deut_uptake"]])
  
  cl_1_diff_frac_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit, fractional = T, theoretical = F)[2]
  cl_2_diff_frac_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit_2, fractional = T, theoretical = F)[2]
  
  cl_1_diff_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit, fractional = F, theoretical = F)[2]
  cl_2_diff_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit_2, fractional = F, theoretical = F)[2]
  
  cl_1_diff_theo_frac_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit, fractional = T, theoretical = T)[2]
  cl_2_diff_theo_frac_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit_2, fractional = T, theoretical = T)[2]
  
  cl_1_diff_theo_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit, fractional = F, theoretical = T)[2]
  cl_2_diff_theo_deut_uptake_new <- calculate_confidence_limit_values(woods_data_set_new, confidence_limit = confidence_limit_2, fractional = F, theoretical = T)[2]
  
  num_cl_1_diff_frac_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_frac_deut_uptake"]]) > cl_1_diff_frac_deut_uptake_new)
  num_cl_2_diff_frac_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_frac_deut_uptake"]]) > cl_2_diff_frac_deut_uptake_new)
  num_cl_1_diff_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_deut_uptake"]] > cl_1_diff_deut_uptake_new))
  num_cl_2_diff_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_deut_uptake"]] > cl_2_diff_deut_uptake_new))
  num_cl_1_diff_theo_frac_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_theo_frac_deut_uptake"]] > cl_1_diff_theo_frac_deut_uptake_new))
  num_cl_2_diff_theo_frac_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_theo_frac_deut_uptake"]] > cl_2_diff_theo_frac_deut_uptake_new))
  num_cl_1_diff_theo_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_theo_deut_uptake"]] > cl_1_diff_theo_deut_uptake_new))
  num_cl_2_diff_theo_deut_uptake_new <- sum(abs(woods_data_set_new[["diff_theo_deut_uptake"]] > cl_2_diff_theo_deut_uptake_new))
  
  ## old
  
  woods_data_set_old <- generate_differential_data_set(dat = dat,
                                                           states = states, 
                                                           protein = protein,
                                                           time_0 = time_0,
                                                           time_t = time_t,
                                                           time_100 = time_100,
                                                           deut_part = deut_part)
  avg_err_diff_frac_deut_uptake_old <- mean(woods_data_set_old[["err_diff_frac_deut_uptake"]])
  avg_err_diff_deut_uptake_old <- mean(woods_data_set_old[["err_diff_deut_uptake"]])
  avg_err_diff_theo_frac_deut_uptake_old <- mean(woods_data_set_old[["err_diff_theo_frac_deut_uptake"]])
  avg_err_diff_theo_deut_uptake_old <- mean(woods_data_set_old[["err_diff_theo_deut_uptake"]])
  
  cl_1_diff_frac_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit, fractional = T, theoretical = F)[2]
  cl_2_diff_frac_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit_2, fractional = T, theoretical = F)[2]
  
  cl_1_diff_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit, fractional = F, theoretical = F)[2]
  cl_2_diff_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit_2, fractional = F, theoretical = F)[2]
  
  cl_1_diff_theo_frac_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit, fractional = T, theoretical = T)[2]
  cl_2_diff_theo_frac_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit_2, fractional = T, theoretical = T)[2]
  
  cl_1_diff_theo_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit, fractional = F, theoretical = T)[2]
  cl_2_diff_theo_deut_uptake_old <- calculate_confidence_limit_values(woods_data_set_old, confidence_limit = confidence_limit_2, fractional = F, theoretical = T)[2]
  
  num_cl_1_diff_frac_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_frac_deut_uptake"]]) > cl_1_diff_frac_deut_uptake_old)
  num_cl_2_diff_frac_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_frac_deut_uptake"]]) > cl_2_diff_frac_deut_uptake_old)
  num_cl_1_diff_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_deut_uptake"]] > cl_1_diff_deut_uptake_old))
  num_cl_2_diff_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_deut_uptake"]] > cl_2_diff_deut_uptake_old))
  num_cl_1_diff_theo_frac_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_theo_frac_deut_uptake"]] > cl_1_diff_theo_frac_deut_uptake_old))
  num_cl_2_diff_theo_frac_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_theo_frac_deut_uptake"]] > cl_2_diff_theo_frac_deut_uptake_old))
  num_cl_1_diff_theo_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_theo_deut_uptake"]] > cl_1_diff_theo_deut_uptake_old))
  num_cl_2_diff_theo_deut_uptake_old <- sum(abs(woods_data_set_old[["diff_theo_deut_uptake"]] > cl_2_diff_theo_deut_uptake_old))
  
  data.frame(Exposure = time_t, CL = 0.98, cl_1_diff_frac_deut_uptake_old, cl_1_diff_frac_deut_uptake_new, cl_1_diff_deut_uptake_old, cl_1_diff_deut_uptake_new, cl_1_diff_theo_frac_deut_uptake_old, cl_1_diff_theo_frac_deut_uptake_new, cl_1_diff_theo_deut_uptake_old, cl_1_diff_theo_deut_uptake_new, num_cl_1_diff_frac_deut_uptake_old, num_cl_1_diff_frac_deut_uptake_new,  num_cl_1_diff_deut_uptake_old, num_cl_1_diff_deut_uptake_new, num_cl_1_diff_theo_frac_deut_uptake_old, num_cl_1_diff_theo_frac_deut_uptake_new, num_cl_1_diff_theo_deut_uptake_old, num_cl_1_diff_theo_deut_uptake_new)
  
}

times <- unique(dat[["Exposure"]][dat[["Exposure"]] > time_0 & dat[["Exposure"]] < time_100])

bind_rows(lapply(times, function(t){
  show_differences(dat, time_t = t)
}))
