## -- Figure S1H

main_dir <- "../data/AATD_BIAS-meta/"

subfolders <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE)

combined_data <- list()

# Process each subfolder
for (subfolder in subfolders) {
  
  features_file <- file.path(subfolder, "ALL_FEATURES.csv")
  names_file <- file.path(subfolder, "ALL_CLASS_NAME.csv")
  
  if (file.exists(features_file) && file.exists(names_file)) {
    
    features_data <- read.csv(features_file, header = FALSE)
    print(features_file)
    names_data <- read.csv(names_file, header = FALSE)
    
    colnames(features_data) <- names_data$V1
    
    features_data$Folder <- basename(subfolder)
    
    combined_data[[length(combined_data) + 1]] <- features_data
  }
}

hist_data <- bind_rows(combined_data)

combine_columns <- function(data, pattern, new_column_name) {
  cols_to_combine <- grep(pattern, colnames(data), value = TRUE, ignore.case = TRUE)
  
  if (length(cols_to_combine) > 0) {
    data[[new_column_name]] <- rowSums(data[cols_to_combine], na.rm = TRUE)
  }
  return(data)
}

hist_data <- combine_columns(final_data, "low|Low", "Combined_Low")

hist_data <- combine_columns(final_data, "medium-|Medium-", "Combined_Medium_Minus")

hist_data <- combine_columns(final_data, "medium\\+|Medium\\+", "Combined_Medium_Plus")

hist_data <- combine_columns(final_data, "high-|High-", "Combined_High_Minus")

hist_data <- combine_columns(final_data, "high\\+|High\\+", "Combined_High_Plus")

hist_data %>%
  dplyr::select(Folder, Combined_Low, Combined_Medium_Minus, Combined_Medium_Plus, Combined_High_Minus, Combined_High_Plus) %>%
  dplyr::rename(sample = Folder) %>%
  mutate(Combined_Low = ifelse(Combined_Low == 0, Combined_Low, 1)) %>%
  group_by(sample) %>%
  summarise(
    Low = sum(Combined_Low, na.rm = TRUE),
    Medium_Minus = sum(Combined_Medium_Minus, na.rm = TRUE),
    Medium_Plus = sum(Combined_Medium_Plus, na.rm = TRUE),
    High_Minus = sum(Combined_High_Minus, na.rm = TRUE),
    High_Plus = sum(Combined_High_Plus, na.rm = TRUE),
    Total = Low + Medium_Minus + Medium_Plus + High_Minus + High_Plus
  ) %>%
  mutate(
    Percent_Low = Low / Total * 100,
    Percent_Medium_Minus = Medium_Minus / Total * 100,
    Percent_Medium_Plus = Medium_Plus / Total * 100,
    Percent_High_Minus = High_Minus / Total * 100,
    Percent_High_Plus = High_Plus / Total * 100
  ) %>%
  left_join(meta_biopsies %>% distinct(sample, kleiner_score)) %>%
  gather(class, fraction, contains("Percent")) %>%
  drop_na(kleiner_score) -> hist_data_summary

normalize_per_group <- function(x) {
  q5 <- quantile(x, 0.01, na.rm = TRUE)
  q95 <- quantile(x, 0.99, na.rm = TRUE)
  scaled <- (x - q5) / (q95 - q5)
  scaled <- pmax(0, pmin(1, scaled)) # Clamp values to [0, 1]
  return(scaled)
}

hist_data_norm <- hist_data %>%
  group_by(Folder) %>%
  mutate(
    Normalized_Intensity = normalize_per_group(`CELLMASK INTENSITY-MEAN AF568`)
  ) %>%
  ungroup() %>%
  dplyr::select(sample = Folder, Object_ID = `Object ID`, File_ID = `File ID`, Normalized_Intensity)

hist_data_norm <- hist_data_norm %>%
  left_join(meta_biopsies %>% distinct(sample, kleiner_score)) %>%
  drop_na(kleiner_score)

ggplot(hist_data_norm , aes(x = Normalized_Intensity, y = as.factor(kleiner_score))) +
  ggridges::geom_density_ridges(fill = "#00AFBB")+
  labs(
    title = "Stacked Density Plot by Kleiner Score",
    x = "Normalized Intensity",
    y = "Density",
    fill = "Kleiner Score"
  ) +
  theme_minimal()

ggsave(file = "../output/Figures/Figure_S1H.pdf", width = 5, height = 5)

ggplot(hist_data_norm , aes(x = Normalized_Intensity, y = as.factor(kleiner_score))) +
  geom_boxplot()
