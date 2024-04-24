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
    mutate(level = "Combined", value = r_L3 + r_L2)
}

plot_population_effects <- function(effects_data, title) {
  effects_data %>%
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
    mutate(CI_excl_one = abs(one_loc - 0.5) * 2) %>%
    ggplot(aes(color = CI_excl_one, x = population)) +
    geom_vline(xintercept = 1, color = "darkred") +
    geom_segment(aes(x = lower, xend = upper, y = population, yend = population)) +
    geom_segment(aes(x = lower50, xend = upper50, y = population, yend = population), linewidth = 2) +
    facet_grid(~effect, scales = "free_y") +
    scale_color_gradientn("CI excluding 1" , limits = c(0,1), colours = c("#303030","#303030","#0571b0","#ca0020","#ca0020"),
                          values = c(0,0.4,0.6,0.95, 1), breaks = c(0,0.5,0.95),
                          labels = c("0%","50%","95%")) +
    scale_x_log10("Fold change") +
    ggtitle(title)
}
