library(tidyverse)
per =read_csv("~/Desktop/per_differences.csv", show_col_types = F)

# Summarize data: compute mean and standard error
per_summary <- per %>%
  group_by(group, time) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            Q1 = quantile(value, 0.25, na.rm = TRUE),
            Q3 = quantile(value, 0.75, na.rm = TRUE),
            .groups = "drop")

# Plot
per_plot = ggplot(per_summary, aes(x = (time), y = mean_value, fill = group)) +
  # geom_line(aes(group = group), position = position_dodge(width = .5), linewidth = .7) +  # Connect mean points
  # geom_point(position = position_dodge(width = .5), size = 2) +
  geom_col(position = position_dodge())+
  geom_errorbar(aes(ymin = mean_value - sd, ymax = mean_value + sd),
                width = .5, position = position_dodge(width = .9), linewidth = .5) +
  labs(x = "Relative Phase", y = "Expression", title = "Per mRNA") +
  scale_fill_brewer(palette = "Pastel1") +  
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20),      # Title font size
    axis.title = element_text(size = 14),     # Axis title font size
    axis.text = element_text(size = 14),      # Axis text font size
    legend.text = element_text(size = 12),    # Legend text font size
    legend.title = element_text(size = 14)    # Legend title font size
  )
  

# ggplot(per, aes(x = jitter(time, amount = .1), y = jitter(value), color = group) )+
#   geom_point(alpha = 0.5)

ggsave("~/Desktop/oscillator_perCryMrna_riboSlow1_riboNoise03to06_bar.pdf", per_plot, device = "pdf", width = 4, height = 5, units = "in")


# Two-way ANOVA with interaction
anova_model_per <- aov(value ~ group * as.factor(time), data = per)

# View ANOVA table
summary(anova_model_per)

bmal =read_csv("~/Desktop/bmal_differences.csv", show_col_types = F)
# Summarize data: compute mean and standard error
bmal_summary <- bmal %>%
  group_by(group, time) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            Q1 = quantile(value, 0.25, na.rm = TRUE),
            Q3 = quantile(value, 0.75, na.rm = TRUE),
            .groups = "drop")

# Plot
bmal_plot = ggplot(bmal_summary, aes(x = (time), y = mean_value, fill = group)) +
  # geom_line(aes(group = group), position = position_dodge(width = .5), linewidth = .7) +  # Connect mean points
  geom_col(position = position_dodge())+
  # geom_point(position = position_dodge(width = .5), size = 2) +
  geom_errorbar(aes(ymin = mean_value - sd, ymax = mean_value + sd),
                width = .5, position = position_dodge(width = .9), linewidth = .5) +
  labs(x = "Relative Phase", y = "Expression", title = "Bmal1 mRNA") +
  scale_fill_brewer(palette = "Pastel1") +  
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20),      # Title font size
    axis.title = element_text(size = 14),     # Axis title font size
    axis.text = element_text(size = 14),      # Axis text font size
    legend.text = element_text(size = 12),    # Legend text font size
    legend.title = element_text(size = 14)    # Legend title font size
  )

ggsave("~/Desktop/oscillator_bmalMrna_riboSlow1_riboNoise03to06_bar.pdf", bmal_plot, device = "pdf", width = 4, height = 5, units = "in")

# Two-way ANOVA with interaction
anova_model <- aov(value ~ group * as.factor(time), data = bmal)

# View ANOVA table
summary(anova_model)
