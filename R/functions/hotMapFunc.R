## Function to create heat maps from HS database
hotMap <- function(fishDat, longRange = c(-129.5, -123), 
                   latRange = c(48, 52), facet = NULL) {
  nAm <- map_data("world") %>% 
    filter(region %in% c("Canada", "USA"))
  
  p <- ggplot(fishDat) +
    stat_density2d(aes(x = long, y = lat, fill = ..level..), 
                   geom = "polygon",
                   alpha = 0.5) +
    scale_fill_gradient(low = "green", high = "red") +
    scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
    geom_map(data = nAm, map = nAm, aes(long, lat, map_id = region), 
             color = "black", fill = "gray80") + 
    coord_fixed(xlim = longRange, ylim = latRange, ratio = 1.3) + 
    # lims(x = longRange, y = latRange) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    samSim::theme_sleekX()
  if (!is.null(facet)) {
    if (facet == "month") {
      nSize <- fishDat %>% 
        group_by(month) %>% 
        tally()
      p <- p +
        facet_wrap(~month) 
    }
    if (facet == "year") {
      nSize <- fishDat %>% 
        group_by(year) %>% 
        tally()
      p <- p +
        facet_wrap(~year) 
    }
    p <- p +
      geom_text(data = nSize, 
                aes(x = -Inf, y = min(latRange), label = n), 
                vjust = 0, hjust = -1)
  }
  return(p)
}
