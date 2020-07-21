## Plot model output functions
# Functions to plot prediction from combined models (i.e. fit_combined.R)
# July 21, 2020
# Modified from plot_predictions to include splines

# Plot aggregate abundance data
plot_abund <- function(dat, ylab) {
  months = range(dat$month_n)
  ggplot(data = dat, aes(x = month_n)) +
    geom_line(aes(y = pred_est, colour = region)) +
    geom_ribbon(aes(ymin = pred_low, ymax = pred_up, fill = region), 
                alpha = 0.5) +
    labs(x = "", y = ylab) +
    ggsidekick::theme_sleek() +
    # facet_wrap(~region) +
    scale_x_continuous(breaks = seq(months[1], months[2], by = 1)) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal)
}

# Plot composition data 
plot_comp <- function(comp_pred, raw_prop, raw = TRUE, 
                      ncol = NULL, facet_scales = "free_y") {
  months = range(comp_pred$month_n)
  p <- ggplot(data = comp_pred, aes(x = month_n)) +
    geom_line(aes(y = pred_prob_est, colour = region)) +
    geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = region), 
                alpha = 0.5) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    scale_x_continuous(breaks = seq(months[1], months[2], by = 1)) +
    labs(y = "Predicted Encounter Probability", x = "Month") +
    facet_wrap(~ stock, ncol = ncol, scales = facet_scales) +
    ggsidekick::theme_sleek() +
    theme(legend.position="top")
  
 if (raw == TRUE) {
   p +
     geom_point(data = raw_prop,
              aes(x = month_n, y = samp_g_ppn, fill = region),
              shape = 21, alpha = 0.4, position = position_dodge(0.6))
 } else {
    p
  }
}

# Plot stock specific abundance data
plot_ss_abund <- function(comp_pred, raw_abund, raw = FALSE,
                          ncol = NULL, facet_scales = "free_y") {
  months = range(comp_pred$month_n)
  p <- ggplot(data = comp_pred, aes(x = month_n)) +
    geom_line(aes(y = comp_abund_est, colour = region)) +
    geom_ribbon(aes(ymin = comp_abund_low, ymax = comp_abund_up, fill = region), 
                alpha = 0.5) +
    scale_fill_manual(name = "Region", values = pal) +
    scale_colour_manual(name = "Region", values = pal) +
    scale_x_continuous(breaks = seq(months[1], months[2], by = 1)) +
    labs(y = "Predicted Catch", x = "Month") +
    facet_wrap(~ stock, ncol = ncol, scales = facet_scales) +
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

