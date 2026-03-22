library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(writexl)

# Read the LSUs in each subbasin
subbasin_text <- readLines("LSUs in each subbasin.txt")

subbasin_list <- list()
current_subbasin <- NULL

for (line in subbasin_text) {
  line <- str_trim(line)
  if (line == "") next  # skip blank lines
  if (grepl("^Subbasin", line)) {
    current_subbasin <- line
    subbasin_list[[current_subbasin]] <- c()
  } else {
    subbasin_list[[current_subbasin]] <- c(subbasin_list[[current_subbasin]], line)
  }
}

# Read the sediment and runoff CSV
sediment_data <- read_csv("Sediment_LSU.csv")  
# Assumes columns: name, mon, day, yr, runoff (mm), sed (t)

# Read the LSU areas
lsu_areas <- read_excel("Lsu areas.xlsx")  # expects columns: name, area (ha)

# Join area information into sediment_data
sediment_data <- sediment_data %>%
  left_join(lsu_areas, by = "name")  # adds the area column (ha)

# Assign subbasin names to each row
sediment_data$Subbasin <- NA

for (subbasin in names(subbasin_list)) {
  lsu_names <- subbasin_list[[subbasin]]
  sediment_data$Subbasin[sediment_data$name %in% lsu_names] <- subbasin
}

sediment_data <- sediment_data %>% filter(!is.na(Subbasin))

output_folder <- "subbasin_daily_summary"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

mapping <- do.call(rbind, lapply(names(subbasin_list), function(sub) {
  data.frame(name = subbasin_list[[sub]], Subbasin = sub, stringsAsFactors = FALSE)
}))

area_by_subbasin <- mapping %>%
  left_join(lsu_areas, by = "name") %>%
  group_by(Subbasin) %>%
  summarise(Total_Area_m2 = sum(area * 10000, na.rm = TRUE))  # convert ha to m²

# Compute daily runoff_volume, flow, and SSC for each subbasin
for (subbasin in unique(sediment_data$Subbasin)) {
  
  sub_data <- sediment_data %>% filter(Subbasin == subbasin)
  
  # Convert area from ha to m²
  sub_data <- sub_data %>% mutate(area_m2 = area * 10000)
  
  summary <- sub_data %>%
    group_by(yr, mon, day) %>%
    summarise(
      weighted_runoff = sum(runoff * area_m2, na.rm = TRUE) / sum(area_m2, na.rm = TRUE),
      runoff_volume = sum((runoff / 1000) * area_m2, na.rm = TRUE),  # runoff in mm → m, then volume in m³
      total_sed = sum(sed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      
      weighted_runoff = ifelse(weighted_runoff < 0.1, 0, weighted_runoff),
      
      runoff_volume = ifelse(weighted_runoff < 0.1, 0, runoff_volume),
      
      total_sed = ifelse(total_sed < 0.01, 0, total_sed),
      
      flow_m3s = runoff_volume / 86400,   
      
      SSC = ifelse(flow_m3s == 0, 0, (total_sed * 1e6) / (flow_m3s * 86400))
    )
  
  clean_name <- gsub(" ", "_", subbasin)
  file_name <- file.path(output_folder, paste0(clean_name, "_daily.xlsx"))
  write_xlsx(summary, file_name)
}


library(ggplot2)
library(readxl)
library(dplyr)
library(scales)

# Generate scatter plots
files <- list.files("subbasin_daily_summary", pattern = "\\.xlsx$", full.names = TRUE)

for (file in files) {
  summary <- read_excel(file)
  summary <- summary %>%
    filter(weighted_runoff >= 0.1, total_sed >= 0.01)
  if (nrow(summary) < 2) next
  
  base_name <- basename(file)
  subbasin_name <- gsub("_daily\\.xlsx", "", base_name)
  subbasin_title <- gsub("_", " ", subbasin_name)
  
  log_model <- lm(log10(total_sed) ~ log10(weighted_runoff), data = summary)
  intercept <- coef(log_model)[1]
  slope <- coef(log_model)[2]
  r_squared <- summary(log_model)$r.squared
  
  eq_label <- paste0(
    "y = ", round(10^intercept, 2), " x^", round(slope, 2),
    "\nR² = ", round(r_squared, 2)
  )
  
  p <- ggplot(summary, aes(x = weighted_runoff, y = total_sed)) +
    geom_point() +
    stat_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      aes(x = weighted_runoff, y = 10^(predict(log_model))),
      color = "blue"
    ) +
    scale_x_log10(labels = label_number()) +
    scale_y_log10(labels = label_number()) +
    annotate("text", 
             x = Inf, y = Inf, 
             label = eq_label,
             hjust = 1.1, vjust = 1.2,
             size = 4, color = "black") +
    labs(
      title = subbasin_title,
      x = "Runoff (mm)",
      y = "Sediment (tons)"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(dirname(file), paste0(subbasin_name, "_loglog.png")),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
}

