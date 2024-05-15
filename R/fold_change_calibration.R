uncertainty_fold_change_group <- function(df_prop_with_preliminary) {
  df_for_model <- df_prop_with_preliminary %>%
    mutate(PatientDesc = paste(Patient_ID, Condition, sep = "___"),
                         isFinal = as.numeric(!isPreliminary)) %>%
    group_by(PatientDesc) %>%
    filter(any(isPreliminary)) %>%
    ungroup()

  fit_calibration <- MASS::glm.nb(n ~ offset(-log(sum_all / 1000)) + PatientDesc + isFinal, data = df_for_model)

  fold_change_with_uncertainty <- df_for_model %>%
    select(n, Condition, sum_all, isPreliminary, PatientDesc, Sample_ID) %>%
    mutate(fitted = log(fitted(fit_calibration)),
           isPreliminary = if_else(isPreliminary, "preliminary", "final")) %>%
    pivot_wider(names_from = isPreliminary, values_from = c("Sample_ID", "n", "sum_all", "fitted")) %>%
    mutate(fold_change = ((n_final + 1) / sum_all_final) /   ((n_preliminary + 1) / sum_all_preliminary),
           fold_change_avg = exp(fit_calibration$coefficients["isFinal"])) %>%
    rowwise() %>%
    mutate(qfc = uncertainty_fold_change_single(fitted_final, sum_all_final, fitted_preliminary, sum_all_preliminary, fit_calibration$theta, fit_calibration$SE.theta)) %>%
    unnest("qfc")

  return(fold_change_with_uncertainty)
}


uncertainty_fold_change_single <- function(fitted1, sum_all1, fitted2, sum_all2, theta, SE.theta, N = 6000) {
  samp_theta <- rlnorm(N, log(theta), log(theta + SE.theta) - log(theta))
  samp1 <- rnbinom(N, mu = exp(fitted1), size = samp_theta)
  samp2 <- rnbinom(N, mu = exp(fitted2), size = samp_theta)
  fc_ratio = ((samp1 + 1) / sum_all1) /   ((samp2 + 1) / sum_all2)

  fc_q <- quantile(fc_ratio, prob = c(0.025,0.975), names = FALSE)
  return(
    data.frame(low = fc_q[1], high = fc_q[2])
  )
}
