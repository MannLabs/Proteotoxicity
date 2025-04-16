############################
## Single Shapes (Fig. 3) ##
############################

## -- Prepare Workspace
cat("\014")
rm(list=ls())
set.seed(123)

setwd(dirname(getActiveDocumentContext()$path))

## -- Package setup
packages <- c("rstudioapi",
              "tidyverse",
              "limma",
              "viridis",
              "ggrepel",
              "pheatmap",
              "XML",
              "sf",
              "Biostrings",
              "stringr")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

## -- Additional functions

transform_string <- function(s) {
  # Separate the letter and number parts of the string
  letter <- substr(s, 1, 1)
  number <- as.numeric(substr(s, 2, nchar(s)))
  
  # Convert letter to number, subtract 1, and convert back to letter
  new_letter <- intToUtf8(utf8ToInt(letter) - 1)
  
  # Subtract 1 from the number
  new_number <- number - 1
  
  # Combine new letter and new number into a new string
  paste0(new_letter, new_number)
}

xml_to_polygon <- function(x){
  data <- xmlParse(x)
  
  xml_data <- xmlToList(data) # Read xml data
  
  n_elements <- as.numeric(xml_data$ShapeCount) # Extract number of shapes
  n_shapes <- n_elements
  
  ## XML format correct
  paste("Shape names correct:", n_elements == sum(grepl("Shape_", names(xml_data))))
  paste("XML length correct:", n_elements == length(names(xml_data)) - 8)
  
  ##  Put shapes into st_sf format one by one
  
  for(i in seq(from = 9, to = length(names(xml_data)), by = 1)){
    
    if(i == 9){polygons_list <- c()} # Empty variable toappend data as it comes
    
    # Wrangling
    tibble::enframe(unlist(xml_data[[i]])) %>%
      mutate(dim = str_replace(name, "_.*", "")) %>%
      dplyr::filter(dim %in% c("X", "Y")) %>%
      mutate(value = as.numeric(value)) %>%
      mutate(name = as.numeric(str_replace(name, ".*_", ""))) %>%
      spread(dim, value) %>%
      column_to_rownames("name") -> pol_i
    
    
    pol_i_mat <- as.matrix((pol_i))
    pol_i_list <- st_polygon(list(rbind(pol_i_mat, pol_i_mat[1,]))) # Close the polygon, and put into correct data format
    
    polygons_list <- append(polygons_list, pol_i_list)
  }
  
  # Put into geom_sf readable format
  polygons_sfc = lapply(polygons_list, FUN = st_polygon)
  
  return(polygons_sfc)
}

polygon_to_map <- function(xml, data, protein, width = 5, height = 5, metadata = meta, subset = NULL, filepath = NULL){
  polygon <- xml_to_polygon(xml)
  
  sample_tmp <- str_replace(str_replace(str_replace(xml, ".*\\/", ""), "__.*", ""),"\\..*","")
  
  meta <- metadata %>%
    filter(sample_short %in% sample_tmp)
  
  level_of <- data.frame(rank = c(1:length(polygon))) %>%
    left_join(meta, by = "rank") %>%
    left_join(data %>%
                filter(sample_short %in% sample_tmp) %>%
                filter(Protein.Group == protein) %>%
                dplyr::select(well, int_norm))
  
  st_sf(polygon, ID = as.factor(c(1:length(polygon))), level_of) -> sf_object
  
  # Conditionally apply filtering if subset is not NULL
  if (!is.null(subset)) {
    sf_object<- sf_object %>% filter(biology == subset)
  }
  
  values <- data %>%
    filter(Protein.Group == protein) %>%
    pull(int_norm)
  
  min <- quantile(values, probs = c(0.05), na.rm = TRUE)
  max <- quantile(values, probs = c(0.95), na.rm = TRUE)
  
  ggplot(data = sf_object)+
    geom_sf(aes(fill = int_norm))+
    theme(legend.position="none") +
    theme_classic() +
    labs(title = protein) +
    scale_fill_viridis_c(limits = c(min, max), oob = scales::squish)+
    scale_y_continuous(limits = c(max(unlist(sf_object$polygon)) - 25000, c(max(unlist(sf_object$polygon)) + 3000)))+
    scale_x_continuous(limits = c(min(unlist(sf_object$polygon)) + 25000, c(min(unlist(sf_object$polygon)) - 3000)))+
    theme_void() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()) -> plot
  
  if (!is.null(filepath)) {
    ggsave(plot, file = filepath, width = width, height = height)
  }
  
  return(plot)
}

normalize_core <- function(data, SA_incl){
  
  data %>%
    filter(sample %in% SA_incl) -> data
  
  data %>%
    select(Protein.Group, sample, int) %>%
    spread(sample, int) %>%
    column_to_rownames("Protein.Group") %>%
    filter(complete.cases(.)) %>%
    rownames_to_column("protein") %>%
    pull(protein) -> core_proteome
  
  data %>%
    filter(Protein.Group %in% core_proteome) %>%
    group_by(sample) %>%
    summarise(median = median(int, na.rm = T)) %>%
    ungroup() %>%
    summarise(median_total = median(median)) %>%
    pull(median_total) -> proteome_core_median
  
  data  %>%
    filter(Protein.Group %in% core_proteome) %>%
    group_by(sample) %>%
    summarise(median_sample = median(int)) %>%
    ungroup() %>%
    mutate(median_total = proteome_core_median) %>%
    mutate(norm_factor = median_total / median_sample) %>%
    dplyr::select(sample, norm_factor) -> proteome_core_norms
  
  data %>%
    left_join(proteome_core_norms) %>%
    mutate(int_norm = int * norm_factor) %>%
    dplyr::rename(int_nonorm = int) %>%
    mutate(int_norm = log2(int_norm)) -> d_norm
  
  return(d_norm)
}

## -- Load data
read_tsv("../data/protein_quantification/scDVP_report.tsv") %>%
  select(contains(c( "Protein.Group", "OA"))) %>%
  gather(sample, int, contains("scDVP")) %>%
  drop_na(int) %>%
  filter(int != 0) %>%
  mutate(batch = str_extract(sample, "(SoSt|FAR)")) %>%
  mutate(well = ifelse(batch == "SoSt",
                       sapply(str_extract(sample, "(?<=cW-)[^_]+(?=_Box)"), transform_string),
                       str_replace_all(sample, ".*_", ""))) -> d_long_noexcl

meta <- read_tsv("../data/meta/meta_scDVP.txt") %>%
  mutate(
    well_letter = str_extract(well, "[A-Z]+"),
    well_number = as.numeric(str_extract(well, "\\d+")),
    well_id = paste(batch, well, sep = "_")) %>%
  group_by(sample_short) %>%
  arrange(well_letter, well_number) %>%
  mutate(rank = row_number()) 

read_tsv("../data/protein_quantification/scDVP_report.tsv") %>%
  select(Protein.Group, Genes) -> meta_pg

## -- Exclusion criteria

d_long_noexcl %>%
  mutate(well_id = paste(batch, well, sep = "_")) %>%
  left_join(meta) %>%
  mutate(area_corrected = area_gpt*(0.345*0.345)/1000) %>%
  group_by(sample_short, sample, area_corrected, batch) %>%
  summarise(n = n()) -> d_stats_n

# Determine criteria
model <- lm(n ~ log2(area_corrected), data = d_stats_n)
coefficients <- coef(model)

d_stats_n %>%
  mutate(residual = (coefficients[1] + coefficients[2] * log2(area_corrected))/n) -> d_stats_residuals

Q1 = quantile(d_stats_residuals$residual, 0.25)
Q3 = quantile(d_stats_residuals$residual, 0.75)
IQR = Q3 - Q1
lower_fence = Q1 - 1.5 * IQR
upper_fence = Q3 + 1.5 * IQR

d_stats_residuals %>%
  mutate(include = !(residual < lower_fence | residual > upper_fence)) -> d_stats_residuals

d_stats_residuals %>%
  filter(include == TRUE) %>%
  pull(sample) -> SA_incl

d_stats_n %>%
  ungroup() %>%
  filter(sample %in% SA_incl) %>%
  summarise(N = median(n))

d_stats_n %>%
  ungroup() %>%
  filter(sample %in% SA_incl) %>%
  group_by(batch) %>%
  summarise(N = median(n))

d_long_noexcl %>%
  mutate(well_id = paste(batch, well, sep = "_")) %>%
  left_join(meta) %>%
  mutate(area_corrected = area_gpt*(0.345*0.345)/1000) %>%
  filter(sample %in% SA_incl) -> d_long

## -- Data normalization

d_norm <- normalize_core(data = d_long, SA_incl = SA_incl)

## -- Call figures

source("Figure_3A.R")
source("Figure_3B.R")
source("Figure_S6B.R")
source("Figure_S6C.R")
source("Figure_S6-truncations.R")
source("Figure_S6D.R")
source("Figure_S6E.R")
source("Figure_S6G.R")
source("Figure_S6H.R")
source("Figure_3C.R")
source("Figure_S7B.R")
source("Figure_S7C.R")
source("Figure_3D.R")