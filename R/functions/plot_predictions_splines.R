## Plot model output functions
# Functions to plot prediction from combined models (i.e. fit_combined.R)
# July 21, 2020
# Modified from plot_predictions to include splines

# Vector of month labels
month_labs <- c("Feb", "Apr", "Jun", "Aug", "Oct", "Dec")

# Plot aggregate abundance data
plot_abund <- function(dat, ylab = NULL) {
  ggplot(data = dat, aes(x = month_n)) +
    geom_line(aes(y = pred_est / 1000, colour = region_c)) +
    geom_ribbon(aes(ymin = pred_low / 1000, ymax = pred_up / 1000, 
                    fill = region_c), 
                alpha = 0.5) +
    labs(x = "", y = ylab) +
    ggsidekick::theme_sleek() +
    theme(axis.text=element_text(size=9),
          # axis.text.y = element_text(angle = 90),
          plot.margin=unit(c(5.5, 10.5, 5.5, 5.5), "points")) +
    scale_fill_manual(name = "Region", values = pal, guide = FALSE) +
    scale_colour_manual(name = "Region", values = pal, guide = FALSE) +
    scale_x_continuous(breaks = seq(2, 12, by = 2), limits = c(1, 12),
                       labels = month_labs, expand = FALSE) +
    coord_cartesian(expand = FALSE, ylim = c(0, NA))
}

# Plot area CPUE data with observations
plot_cpue_area <- function(dat, obs_dat) {
  ggplot(data = dat, aes(x = month_n)) +
    geom_line(aes(y = pred_est_logcpue, colour = region_c)) +
    geom_ribbon(aes(ymin = pred_low_logcpue, ymax = pred_up_logcpue, 
                    fill = region_c), 
                alpha = 0.75) +
    geom_point(data = obs_dat, aes(x = month_n, y = log_area_cpue,
                                   fill = region_c),
               shape = 21, alpha = 0.5) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    scale_x_continuous(breaks = seq(2, 12, by = 2), limits = c(1, 12),
                       labels = month_labs, expand = c(0, 0)) +
    facet_wrap(~fct_reorder(area, as.numeric(region))) +
    labs(x = "Month", y = "ln(catch) / ln(effort)") +
    ggsidekick::theme_sleek()
}

# Plot composition data 
plot_comp <- function(comp_pred, raw_prop, raw = TRUE, 
                      ncol = NULL, facet_scales = "free_y") {
  p <- ggplot(data = comp_pred, aes(x = month_n)) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    labs(y = "Predicted Stock Composition", x = "Month") +
    facet_wrap(~ stock, ncol = ncol, scales = facet_scales) +
    ggsidekick::theme_sleek() +
    theme(legend.position = "top",
          axis.text=element_text(size=9),
          # axis.text.y = element_text(angle = 90),
          plot.margin = unit(c(5.5, 10.5, 5.5, 5.5), "points"),
          panel.spacing = unit(0.75, "lines")) +
    scale_x_continuous(breaks = seq(2, 12, by = 2), limits = c(1, 12),
                       labels = month_labs, expand = FALSE)
    
  if (raw == TRUE) {
    p2 <- p +
      geom_line(aes(y = pred_prob_est, colour = region_c)) +
      geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = region_c), 
                  alpha = 0.5) +
      geom_point(data = raw_prop,
                 aes(x = month_n, y = samp_g_ppn, fill = region_c),
                 shape = 21, alpha = 0.4, position = position_dodge(0.6))
  } else {
    p2 <- p +
      geom_line(aes(y = pred_prob_est, colour = region_c)) +
      geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = region_c), 
                  alpha = 0.5)
  }
  
  p2 + 
    coord_cartesian(expand = FALSE, ylim = c(0, NA))
}

# Plot stacked composition data 
plot_comp_stacked <- function(comp_pred, grouping_col) {
  palette_name <- ifelse(grouping_col == "pst_agg", "sunset", "midnight")
  stock_pal <- disco::disco(palette_name, n = length(unique(comp_pred$stock)))
  p <- ggplot(data = comp_pred, aes(x = month_n)) +
    geom_area(aes(y = pred_prob_est, colour = stock, fill = stock), 
              stat = "identity") +
    scale_fill_manual(name = "Stock", values = stock_pal) +
    scale_colour_manual(name = "Stock", values = stock_pal) +
    labs(y = "Predicted Stock Composition", x = "Month") +
    ggsidekick::theme_sleek() +
    theme(legend.position = "top",
          axis.text=element_text(size=9),
          plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points")
          ) +
    scale_x_continuous(breaks = seq(2, 12, by = 2), 
                       limits = c(1, 12),
                       labels = month_labs) +
    facet_wrap(~region_c) +
    coord_cartesian(expand = FALSE, ylim = c(0, NA))
}


plot_ss_abund <- function(comp_pred, raw_abund, raw = FALSE,
                          ncol = NULL, facet_scales = "free_y") {
  #adjustment to account for massive CIs on certain stocks in NSoG
  if ("NSoG" %in% comp_pred$region) {
    comp_pred <- comp_pred %>%
      group_by(stock) %>%
      mutate(temp_id = paste(stock, region, sep = "_")) %>% 
      filter(!temp_id %in% c("Fraser_Spring_4.2_NSoG", "CA/OR-coast_NSoG"))
  }
  
  p <- ggplot(data = comp_pred, aes(x = month_n)) +
    geom_line(aes(y = comp_abund_est, colour = region_c)) +
    geom_ribbon(aes(ymin = comp_abund_low, ymax = comp_abund_up, 
                    fill = region_c), 
                alpha = 0.5) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    labs(y = "Predicted Catch", x = "Month") +
    facet_wrap(~ stock, ncol = ncol, scales = facet_scales) +
    ggsidekick::theme_sleek() +
    theme(legend.position = "top",
          axis.text=element_text(size=9),
          plot.margin=unit(c(5.5, 10.5, 5.5, 5.5), "points")) +
    scale_x_continuous(breaks = seq(2, 12, by = 2), limits = c(1, 12),
                       labels = month_labs, expand = c(0, 0))

  if (raw == TRUE) {
    p <- p +
      geom_point(data = raw_abund,
                 aes(x = as.numeric(as.character(month)), y = cpue_g / 10,
                     fill = region),
                 shape = 21, alpha = 0.3, position = position_dodge(0.6))
  }
  
  p +
    coord_cartesian(expand = FALSE, ylim = c(0, NA))
    
}

