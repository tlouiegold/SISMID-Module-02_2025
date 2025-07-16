library("nfidd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("purrr")
library("lubridate")
library("tidybayes")

set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)

df <- add_delays(
  infection_times,
  hosp_params = list(meanlog = 1.0, sdlog = 0.5)
)

head(df)
