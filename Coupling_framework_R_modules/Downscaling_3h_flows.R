library(readr)
library(dplyr)
library(lubridate)
library(stringr)

# Process the 3-hourly reference file
reference <- read_csv("North Opuha 3h.csv", show_col_types = FALSE) %>%
  rename(datetime = Date, flow_3h = Flow) %>%
  mutate(
    datetime = parse_date_time(datetime, orders = "dmy HM"),
    date = as.Date(datetime)
  )

# Process the daily flows file with multiple channels
daily_flows <- read_csv("North sub basins flow.csv", show_col_types = FALSE) %>%
  rename(date = Date) %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y"))

# Build mean-normalized factors (mean factor per day = 1)
factors <- reference %>%
  group_by(date) %>%
  mutate(
    day_mean = mean(flow_3h, na.rm = TRUE),
    factor = ifelse(is.finite(day_mean) & day_mean > 0, flow_3h / day_mean, 1)
  ) %>%
  ungroup() %>%
  select(datetime, date, factor)

# Disaggregate each channel flow
channel_names <- setdiff(names(daily_flows), "date")

for (channel in channel_names) {
  message("Processing channel: ", channel)

  disagg <- factors %>%
    left_join(daily_flows %>% select(date, all_of(channel)), by = "date") %>%
    rename(daily_flow = !!channel) %>%
    mutate(flow3h = daily_flow * factor) %>%
    select(datetime, flow3h)

  # Verification: daily average of 3-hour values should equal the original daily flow
  check <- disagg %>%
    mutate(date = as.Date(datetime)) %>%
    group_by(date) %>%
    summarise(daily_avg = mean(flow3h, na.rm = TRUE), .groups = "drop") %>%
    left_join(daily_flows %>% select(date, !!channel), by = "date") %>%
    mutate(diff = daily_avg - .data[[channel]])

  output_name <- paste0(channel, "_3h.csv")
  write_csv(disagg, output_name)
}
