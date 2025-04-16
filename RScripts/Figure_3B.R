## -- Figure S6D

# Determine exclusion criteria
model <- lm(n ~ log2(area_corrected), data = d_stats_n)
coefficients <- coef(model)

d_stats_n %>%
  mutate(residual = (coefficients[1] + coefficients[2] * log2(area_corrected))/n) -> d_stats_residuals

Q1 = quantile(d_stats_residuals$residual, 0.25)
Q3 = quantile(d_stats_residuals$residual, 0.75)
IQR = Q3 - Q1
lower_fence = Q1 - 2.5 * IQR
upper_fence = Q3 + 2.5 * IQR

d_stats_residuals %>%
  mutate(included = !(residual < lower_fence | residual > upper_fence)) -> d_stats_residuals

exclusion_region <- data.frame(area_corrected = seq(100, 1400, by = 1)) %>%
  mutate(n = (coefficients[1] + coefficients[2] * log2(area_corrected))/upper_fence)

# Plot
ggplot() +
  geom_point(data = d_stats_residuals, aes(x = area_corrected, y = n, shape = included), fill = "darkblue")+
  geom_line(data = exclusion_region, aes(x = area_corrected, y = n), lty = "dotted")+
  theme_classic()+
  geom_smooth(data = d_stats_residuals, aes(x = area_corrected, y = n),
              formula = y ~ log2(x), se = F, method = "lm", lty = "dotted", color = "darkgrey")+
  labs(x = "Area (um2)", y = "Number of proteins")+
  scale_shape_manual(values = c(4, 21)) +
  scale_x_continuous(limits = c(123,1500), breaks = seq(100,1400, by = 100))

ggsave(file = "../output/Figures/Figure_3B.pdf", width = 5, height = 5)
