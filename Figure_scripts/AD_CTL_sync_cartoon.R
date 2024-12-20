library(ggplot2)
library(gridExtra)
library(grid)

# Generate data
x <- seq(0, 2 * pi, length.out = 100) # Shared x-axis for sine waves
colors <- c("red", "blue", "green", "orange", "purple") # Colors for the waves
phase_shifts <- c(pi, pi/2, 0, 2*pi/3-0.4, pi/5-.4)# Different phase shifts for out-of-sync waves

# Calculate the dot positions
dot_positions <- lapply(seq_along(colors), function(i) {
  dot_time <-  (i - 1) * (2*pi / 4) # Increment pi/4 for each row
  list(
    sync = data.frame(x = dot_time, y = sin(dot_time), color = colors[i]),
    async = data.frame(x = dot_time, y = sin(dot_time + phase_shifts[i]), color = colors[i])
  )
})

# Function to generate sine wave plots with a dot
create_wave_plot <- function(x, y, y_label, color, dot_data) {
  ggplot(data.frame(x, y), aes(x, y)) +
    geom_line(color = color, linewidth = 1.2) +
    geom_point(data = dot_data, aes(x, y), color = color, size = 4) +
    theme_gray() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      plot.title = element_text(size = 20),      # Title font size
      legend.text = element_text(size = 12),    # Legend text font size
      legend.title = element_text(size = 14)    # Legend title font size
    ) +
    ylab(y_label)+
    ylim(c(-1.2,1.2))
  
}

# Function to create a plot for all dots
create_all_dots_plot <- function(y_label, dot_data) {
  
  ggplot(dot_data, aes(x, y, color = color)) +
    geom_point(size = 4) +
    scale_color_identity() + # Use colors as-is
    theme_gray() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      # plot.margin = margin(10, 10, 30, 30), # Add space around the plot
      plot.title = element_text(size = 20),      # Title font size
      axis.title = element_text(size = 14),     # Axis title font size
      legend.text = element_text(size = 12),    # Legend text font size
      legend.title = element_text(size = 14)    # Legend title font size
    ) +
    ylab(y_label)+
    ylim(c(-1.2,1.2))
  
  
}

# Create the plots
plots <- list()
for (i in 1:6) {
  # In-sync sine wave with dot
  if (i < 6) {
    y_in_sync <- sin(x)
    dot_time <-  (i - 1) * (2*pi / 4) # Increment pi/4 for each row
    plots[[i * 2 - 1]] <- create_wave_plot(
      x, y_in_sync, paste("CTL", i), colors[i], dot_positions[[i]]$sync
    )
    
    # Out-of-sync sine wave with dot
    y_out_of_sync <- sin(x + phase_shifts[i])
    plots[[i * 2]] <- create_wave_plot(
      x, y_out_of_sync, paste("AD", i), colors[i], dot_positions[[i]]$async
    )
  } else {
    # Row 6: All dots
    dot_data_ctl <- do.call(rbind, lapply(dot_positions, function(pos) pos$sync))
    plots[[i * 2 - 1]] <- create_all_dots_plot("CTL Measured", dot_data_ctl)
    
    dot_data_ad <- do.call(rbind, lapply(dot_positions, function(pos) pos$async))
    plots[[i * 2]] <- create_all_dots_plot("AD Measured", dot_data_ad)
  }
}

# Arrange the plots in a 5x2 grid
out = grid.arrange(grobs = plots, nrow = 6, ncol = 2) # X-axis title for the entire plot

ggsave("~/Desktop/AD_CTL_sync_cartoon_fig1.pdf", plot = out, width = 4.5, height = 5.4, units = "in")