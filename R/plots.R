scale_color_condition <- function(...) {
  scale_color_brewer(palette = "Dark2", ...)
}

scale_fill_condition <- function(...) {
  scale_fill_brewer(palette = "Dark2", ...)
}

fold_change_breaks_major <- function(range) {
  min_log <- ceiling(log10(range[1]))
  max_log <- floor(log10(range[2]))
  if(max_log - min_log < 2) {
    min_log <- ceiling(log2(range[1]))
    max_log <- floor(log2(range[2]))
    if(max_log - min_log < 2) {
      if(range[2] - range[1] <= 0.5) {
        lin_step <- 0.1
        start_exp <- 0
      } else {
        lin_step <- 0.2
        start_exp <- 1
      }
      if(range[2] >= 1.1) {
        high <- 1 + 2^seq(start_exp, log2((range[2] - 1) * 10)) / 10
      } else {
        high <- c()
      }
      low <- seq(1, round(range[1],2) - lin_step / 2, by = -0.2)
      return(sort(unique(c(low,high))))
    } else if(max_log - min_log <= 4) {
      return(2^(min_log:max_log))
    } else {
      return(sort(unique(c(1, 2^(seq(min_log,max_log, by = 2))))))
    }
  } else if(max_log - min_log <= 4) {
    return(10^(min_log:max_log))
  } else {
    return(sort(unique(c(1, 10^(seq(min_log,max_log, by = 2))))))
  }
  stop("Problem")
}

fold_change_breaks_minor <- function(range, major) {
  steps_linear <- unique(round(diff(major),2))
  if(length(steps_linear) == 1) {
    step <- steps_linear
    brks <- seq(min(major) - step, max(major) + step, by = 0.5 * step)
  } else {
    #step <- max(diff(log(major)))
    #brks <- exp(seq(min(log(major)) - step, max(log(major)) + step, by = 0.5 * step))
    log_diffs <- diff(log(major))
    log_brks_within <- log(major)[-length(major)] + 0.5 * log_diffs
    log_brk_low <- log(major[1]) - 0.5 * log_diffs[1]
    log_brk_high <- log(major[length(major)]) + 0.5 * log_diffs[length(major) - 1]
    brks <- exp(c(log_brk_low, log_brks_within, log_brk_high))
  }
  brks[brks > range[1] & brks < range[2]]
}

labels_log10_simple <- function(breaks) {
  labs <- format(breaks, scientific = FALSE, nsmall = 0,drop0trailing = TRUE, trim = TRUE)
  steps_linear <- unique(round(diff(breaks),2))
  if(length(steps_linear) > 1) {
    step <- max(diff(log(breaks)))
    if(round(step, 1) %in% c(2,4,8,16)) {
      base <- 2
    } else {
      base <- 10
    }
    extreme <- str_length(labs) > 4
    labs[extreme] <- paste0(base,"^{", format(log(breaks, base = base),scientific = FALSE, nsmall = 0,drop0trailing = TRUE, trim = TRUE),"}")[extreme]
    parse(text = labs)
  } else {
    labs
  }
}

scale_x_fold_change <- function(name = "Fold change", breaks =  fold_change_breaks_major, ...) {
  scale_x_continuous(name = name, breaks = breaks, minor_breaks = fold_change_breaks_minor,
                labels = labels_log10_simple,
                transform = "log10", ...)
}

scale_y_fold_change <- function(name = "Fold change", breaks =  fold_change_breaks_major, ...) {
  scale_y_continuous(name = name, breaks = breaks, minor_breaks = fold_change_breaks_minor,
                     labels = labels_log10_simple,
                     transform = "log10", ...)
}

no_change_line <- function(yintercept = 1, ...) {
  geom_hline(yintercept = yintercept, color = "darkviolet")
}

get_l2_draws <- function(fit) {
  tidybayes::spread_draws(fit, r_L2[population, effect]) %>%
    mutate(level = "L2") %>%
    rename(value = r_L2)
}

get_l3_draws <- function(fit) {
  tidybayes::spread_draws(fit, r_L3[population, effect])  %>%
    mutate(level = "L3") %>%
    rename(value = r_L3)
}

get_combined_draws <- function(fit, dataset) {
  l2_map_df <- dataset %>% ungroup() %>% select(L2,L3) %>% distinct()
  l2_map <- l2_map_df$L2
  names(l2_map) <- l2_map_df$L3

  l2 <- tidybayes::spread_draws(fit, r_L2[l2_pop, effect])
  l3 <- tidybayes::spread_draws(fit, r_L3[population, effect]) %>%
    mutate(l2_pop = l2_map[population])

  l2$l2_pop <- gsub(".", " ", l2$l2_pop, fixed = TRUE)

  l2 %>% inner_join(l3, by = c("l2_pop", "effect", ".chain", ".iteration", ".draw"),
                    relationship = "one-to-many",
                    unmatched = "error") %>%
    mutate(level = "Combined", value = r_L3 + r_L2) %>%
    select(-r_L2, -r_L3)
}

add_diaT1_draws <- function(effects_data) {
  effects_data_wide <- effects_data %>% pivot_wider(names_from = "effect", values_from = "value")

  control_t1_effect <- effects_data_wide %>%
    mutate(effect = "isDia + isT1", value = isDia + isT1) %>%
    select(-isDia, -isT1)

  stopifnot(nrow(control_t1_effect) == nrow(effects_data) / 2)

  combined_draws <- rbind(effects_data, control_t1_effect) %>%
    mutate(effect = factor(effect, levels = c("isDia", "isT1", "isDia + isT1")))

  stopifnot(all(!is.na(combined_draws$effect)))
  return(combined_draws)
}

plot_population_effects <- function(effects_data, title, order_populations = TRUE, fc_breaks = fold_change_breaks_major) {
  data_for_plot <- effects_data %>%
    mutate(value = exp(value)) %>%
    group_by(level, population, effect) %>%
    summarise(Estimate = median(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975),
              lower50 = quantile(value, 0.25),
              upper50 = quantile(value, 0.75),
              one_loc = mean(value < 1),
              .groups = "drop") %>%
    ungroup() %>%
    mutate(CI_excl_one = abs(one_loc - 0.5) * 2)

  if(order_populations) {
    # Hack to allow different ordering in different subpanels
    # Addint spaces to make population names different between panels
    data_for_plot <- data_for_plot %>% mutate(
      population = paste0(population, strrep(" ", as.integer(as.factor(effect)) - 1)),
      population = fct_reorder(population, Estimate, .desc = TRUE)
    )
    facet <- facet_wrap(~effect, scales = "free", nrow = 1)
  } else {
    facet <- facet_grid(~effect, scales = "free_x")
  }

  data_for_plot %>%
    ggplot(aes(color = CI_excl_one, x = population)) +
    geom_vline(xintercept = 1, color = "darkred") +
    geom_segment(aes(x = lower, xend = upper, y = population, yend = population)) +
    geom_segment(aes(x = lower50, xend = upper50, y = population, yend = population), linewidth = 2) +
    facet +
    scale_color_gradientn("CI excluding 1" , limits = c(0,1), colours = c("#303030","#303030","#0571b0","#ca0020","#ca0020"),
                          values = c(0,0.4,0.6,0.95, 1), breaks = c(0,0.5,0.95),
                          labels = c("0%","50%","95%")) +
    scale_x_fold_change(breaks = fc_breaks) +
    ggtitle(title)
}

compute_preliminary_comparison <- function(all_prop_with_preliminary){
  res <- all_prop_with_preliminary %>%
    mutate(isPreliminary = if_else(isPreliminary, "preliminary", "final")) %>%
    pivot_wider(names_from = isPreliminary, values_from = c("Sample_ID", "Experiment_ID", "n", "sum_parent", "sum_all", "prop_parent", "prop_all")) %>%
    filter(!is.na(Sample_ID_preliminary))

  stopifnot(!any(is.na(res$Sample_ID_final)))
  return(res)
}

plot_preliminary_comparison <- function(preliminary_comparison, Level = "L2") {
  preliminary_comparison <- preliminary_comparison %>%
    filter(Level == !!Level) %>%
    group_by(Patient_ID, Condition) %>%
    mutate(fold_change = ((n_final + 1) / sum_all_final) /   ((n_preliminary + 1) / sum_all_preliminary))

  preliminary_comparison %>%
    ggplot(aes(x = Name, y = fold_change, color = Condition)) +
    no_change_line() +
    stat_summary(fun.data = "mean_se", position = position_dodge(width = 0.3), linewidth = 1, fatten = 6) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.4) +
    scale_y_fold_change()  +
    scale_color_condition() +
    coord_flip()
}
