fold_change_calibration_fit <- function(df_prop_with_preliminary, cache_dir, include_age = TRUE) {

  df_for_model <- df_prop_with_preliminary %>%
    mutate(PatientDesc = paste(Patient_ID, Condition, sep = "___"),
                         isFinal = as.numeric(!isPreliminary),
           AgeStd = (10 - Age)/5) %>%
    group_by(PatientDesc) %>%
    filter(any(isPreliminary)) %>%
    ungroup()

  if(include_age) {
    calibration_formula <- n ~ offset(-log(sum_all / 1000)) + PatientDesc + isFinal + isFinal:Condition + isFinal:AgeStd
  } else {
    calibration_formula <- n ~ offset(-log(sum_all / 1000)) + PatientDesc + isFinal + isFinal:Condition
  }

  model_label <- paste0(unique(df_prop_with_preliminary$Name), collapse = "-")

  cache_basename <- paste0(model_label, "_", rlang::hash(df_for_model), "_", rlang::hash(as.character(calibration_formula)), ".rds")
  cache_file <- file.path(cache_dir, cache_basename)

  if(mean(df_prop_with_preliminary$n == 0) > 0.5) {
    cat(model_label, " - more than 50% zeroes, skipping\n")
    return(NULL)
  }


  if(file.exists(cache_file)) {
    cat(model_label, " - calibration read from cache\n")
    return(readRDS(cache_file))
  }

  cat(model_label, " - running calibration fit\n")

  fit_calibration <- rstanarm::stan_glm.nb(calibration_formula, data = df_for_model, iter = 6000, warmup = 1000)

  saveRDS(fit_calibration, cache_file)

  return(fit_calibration)
}


calibration_fit_to_df <- function(fit_calibration) {
  if(is.null(fit_calibration)) {
    return(NULL)
  }
  df_for_model <- fit_calibration$data
  df_result <- tidybayes::add_predicted_rvars(df_for_model, fit_calibration, offset = -log(df_for_model$sum_all / 1000)) %>%
    mutate(pred_ratio = (.prediction +  1) / sum_all,
           pred_ratio_renorm = (.prediction +  1) / (sum_all - n + .prediction)) %>%
    select(n, Condition, sum_all, isPreliminary, PatientDesc, Sample_ID, Patient_ID, pred_ratio, pred_ratio_renorm, L1, L2,L3,Name)  %>%
    mutate(isPreliminary = if_else(isPreliminary, "preliminary", "final")) %>%
    pivot_wider(names_from = isPreliminary, values_from = c("Sample_ID", "n", "sum_all", "pred_ratio", "pred_ratio_renorm")) %>%
    mutate(fold_change = ((n_final + 1) / sum_all_final) /   ((n_preliminary + 1) / sum_all_preliminary),
           fc_rvar = pred_ratio_final / pred_ratio_preliminary,
           fc_rvar_renorm = pred_ratio_renorm_final / pred_ratio_renorm_preliminary,
           low = as.numeric(quantile(fc_rvar, 0.025)),
           high = as.numeric(quantile(fc_rvar, 0.975)),
           low_renorm = as.numeric(quantile(fc_rvar_renorm, 0.025)),
           high_renorm = as.numeric(quantile(fc_rvar_renorm, 0.975))
           ) %>%
    select(-starts_with("fc_rvar"),-starts_with("pred_ratio"))

  stopifnot(nrow(df_for_model) == 2 * nrow(df_result))

  return(df_result)
}

calibration_fit_to_fc_df <- function(fit_calibration) {
  if(is.null(fit_calibration)) {
    return(NULL)
  }
  df_for_model <- fit_calibration$data
  df_meta <- df_for_model %>% group_by(L1, L2,L3, Level, Name) %>%
    summarise(median_n_prelim = median(n[isPreliminary]),
              median_n_final = median(n[!isPreliminary]),
              median_sum_all_final = median(sum_all[!isPreliminary]))

  stopifnot(nrow(df_meta) == 1)
  df_meta$median_n <- median(df_for_model$n)

  fit_rvars <- posterior::as_draws_rvars(fit_calibration)
  overall_control <- exp(fit_rvars$isFinal)
  overall_dia <- exp(fit_rvars$isFinal + fit_rvars$`isFinal:ConditionDia T0`)
  overall_ratio <- exp(fit_rvars$`isFinal:ConditionDia T0`)

  df_control <- df_meta %>% mutate(Condition = "Ctrl T0", fc_rvar = overall_control)
  df_dia <- df_meta %>% mutate(Condition = "Dia T0", fc_rvar = overall_dia)
  df_ratio <- df_meta %>% mutate(Condition = "Ratio Dia T0/Ctrl T0", fc_rvar = overall_ratio)

  if(!is.null(fit_rvars$`isFinal:AgeStd`)) {
    overall_age <- exp(fit_rvars$`isFinal:AgeStd`)
    df_age <- df_meta %>% mutate(Condition = "Ratio Age (5 years)", fc_rvar = overall_age)
  } else {
    df_age <- NULL
  }

  df_result <- rbind(df_control, df_dia, df_ratio, df_age) %>%
    mutate(
      fc_rvar_renorm = ((median_n_prelim * fc_rvar + 1) / (median_sum_all_final - median_n_final + median_n_prelim * fc_rvar)) / (median_n_prelim / median_sum_all_final),
      fold_change = as.numeric(median(fc_rvar)),
      fold_change_renorm = as.numeric(median(fc_rvar_renorm)),
      low = as.numeric(quantile(fc_rvar, 0.025)),
      high = as.numeric(quantile(fc_rvar, 0.975)),
      low_renorm = as.numeric(quantile(fc_rvar_renorm, 0.025)),
      high_renorm = as.numeric(quantile(fc_rvar_renorm, 0.975))
  ) %>%
    select(-starts_with("fc_rvar"))

  return(df_result)
}

calibration_fit_all <- function(df_prop_with_preliminary, Level, cache_dir, include_age = TRUE) {
  df_prop_with_preliminary %>%
    filter(Level == !!Level) %>%
    group_by(Name) %>%
    group_map(\(df, group) fold_change_calibration_fit(cbind(group,df), cache_dir = cache_dir, include_age = include_age))
}

calibration_df_all <- function(df_prop_with_preliminary, Level, cache_dir, include_age = TRUE) {

  df_filtered <- df_prop_with_preliminary %>%
    filter(Level == !!Level)

  cache_basename <- paste0("df_all_", if_else(include_age, "_age", ""), rlang::hash(df_filtered), ".rds")
  cache_file <- file.path(cache_dir, cache_basename)

  if(file.exists(cache_file)) {
    cat("Reading from cache ", cache_basename)
    return(readRDS(cache_file))
  }


  fit_all <- calibration_fit_all(df_filtered, Level, cache_dir, include_age = include_age)
  df_all <- fit_all %>% purrr::map_df(calibration_fit_to_df)
  fc_df_all <- fit_all %>% purrr::map_df(calibration_fit_to_fc_df)

  result <- list(df = df_all, fc_df = fc_df_all)

  saveRDS(result, cache_file)
  return(result)
}

calibration_plot <- function(calibration_df, fc_df, renorm = FALSE) {
  fc_df  <- fc_df %>% filter(Condition %in% unique(calibration_df$Condition))

  if(renorm) {
    range_aes <- aes(ymin = low_renorm, ymax = high_renorm)
    hline_aes <- aes(yintercept = fold_change_renorm, color = Condition)
    rect_aes <- aes(ymin = low_renorm, ymax = high_renorm, fill = Condition)
  } else {
    range_aes <- aes(ymin = low, ymax = high)
    hline_aes <- aes(yintercept = fold_change, color = Condition)
    rect_aes <- aes(ymin = low, ymax = high, fill = Condition)
  }

  calibration_df %>%
    ggplot(aes(x = PatientDesc, y = fold_change, color = Condition)) +
    geom_hline(data = fc_df, hline_aes, linewidth = 1) +
    geom_rect(data = fc_df, rect_aes, inherit.aes = FALSE, xmin = -Inf, xmax = Inf, alpha = 0.3) +
    no_change_line() +
    geom_linerange(range_aes, alpha = 0.8) +
    geom_point(alpha = 0.8) +
    scale_y_fold_change() +
    scale_x_discrete("") +
    scale_color_condition() +
    scale_fill_condition() +
    coord_flip() +
    theme(axis.text.y = element_blank()) +
    facet_wrap(~Name, scales = "free_x")
}

fc_plot <- function(fc_df, renorm = FALSE) {
  if(renorm) {
    main_aes <- aes(x = Name, ymin = low_renorm, y = fold_change_renorm, ymax = high_renorm, color = Condition)
  } else {
    main_aes <- aes(x = Name, ymin = low, y = fold_change, ymax = high, color = Condition)
  }
  fc_df %>%
    mutate(group = if_else(grepl("Ratio", Condition), "Ratios", "Base")) %>%
    ggplot(main_aes) +
    no_change_line() +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    scale_y_fold_change() +
    scale_color_condition() +
    coord_flip() +
    facet_wrap(~group, scales = "free_x")
}
