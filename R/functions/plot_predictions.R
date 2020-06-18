## Plot model output functions
# Functions to plot prediction from combined models (i.e. fit_combined.R)
# June 17, 2020

# Plot aggregate abundance data
plot_abund <- function(dat, ylab) {
  ggplot() +
    geom_pointrange(data = dat, aes(x = as.numeric(as.character(month)), 
                                    y = pred_est,
                                    ymin = pred_low, ymax = pred_up)) +
    labs(x = "", y = ylab) +
    ggsidekick::theme_sleek() +
    facet_wrap(~region) +
    scale_x_continuous(breaks = seq(1, 12, by = 1))
}

# Plot composition data 
plot_comp <- function(comp_pred, raw_prop, raw = TRUE, ncol = NULL) {
  p <- ggplot() +
    geom_pointrange(data = comp_pred,
                    #adjust to num-char to add spaces on x labels
                    aes(x = as.numeric(as.character(month)), y = pred_prob_est,
                        ymin = pred_prob_low, ymax = pred_prob_up,
                        fill = region),
                    shape = 21, size = 0.4, position = position_dodge(0.6)) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_x_continuous(breaks = seq(1, 12, by = 1)) +
    labs(y = "Probability", x = "Month") +
    facet_wrap(~ stock, ncol = ncol) +
    ggsidekick::theme_sleek()
  
  if (raw == TRUE) {
    p + 
      geom_point(data = raw_prop,
               aes(x = as.numeric(as.character(month)), y = samp_g_ppn, 
                   fill = region),
               shape = 21, alpha = 0.4, position = position_dodge(0.6))
  } else {
    p
  }
}

# Plot stock specific abundance data
plot_ss_abund <- function(comp_pred, raw_abund, raw = TRUE) {
  p <- ggplot() +
    geom_pointrange(data = comp_pred,
                    aes(x = as.numeric(as.character(month)), y = comp_abund_est,
                        ymin = comp_abund_low, ymax = comp_abund_up, 
                        fill = region),
                    shape = 21, position = position_dodge(0.6)) +
    facet_wrap(~stock, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = seq(1, 12, by = 1)) +
    scale_fill_manual(name = "Region", values = pal) +
    labs(x = "Month", y = "Predicted Catch") +
    ggsidekick::theme_sleek() +
    theme(legend.position="top")
  
  if (raw == TRUE) {
    p + 
      geom_point(data = raw_abund, 
                 aes(x = as.numeric(as.character(month)), y = catch_g, 
                     fill = region),
                 shape = 21, alpha = 0.3, position = position_dodge(0.6))
  } else {
    p
  }
}

