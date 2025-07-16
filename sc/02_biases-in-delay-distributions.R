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

# Use the floor() function to round down to integers
df_dates <- df |>
  mutate(
    infection_time = floor(infection_time),
    onset_time = floor(onset_time),
    hosp_time = floor(hosp_time)
  )
head(df_dates)

df_dates <- df_dates |>
  mutate(
    incubation_period = onset_time - infection_time,
    onset_to_hosp = hosp_time - onset_time
  )

summary(df$onset_time - df$infection_time)

summary(df_dates$incubation_period)

mod <- nfidd_cmdstan_model("lognormal")
res <- nfidd_sample(mod,
                    data = list(
                      n = nrow(na.omit(df_dates)),
                      y = na.omit(df_dates)$onset_to_hosp + 0.01
                    )
)

res

res %>%
  summarise_lognormal()

cmod <- nfidd_cmdstan_model("censored-delay-model")
cmod

cres <- nfidd_sample(cmod,
                     data = list(
                       n = nrow(na.omit(df_dates)),
                       onset_to_hosp = na.omit(df_dates)$onset_to_hosp
                     )
)

cres %>%
  summarise_lognormal()

set.seed(123)
df <- add_delays(infection_times,
                 hosp_params = list(meanlog = 1.75, sdlog = 0.5))

library(ggplot2)
library(dplyr)

# Set parameters for the truncation demo
set.seed(890)
n <- 5000
meanlog <- 1.75
sdlog <- 0.5

# Generate delays
true_delays <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)

# Simulate truncation at different time points
truncation_times <- c(6, 10, 15)
truncated_data <- map_dfr(truncation_times, function(t) {
  truncated_delays <- true_delays[true_delays <= t]
  data.frame(
    delay = truncated_delays,
    truncation = paste("Truncated at", t, "days"),
    mean_delay = mean(truncated_delays)
  )
})

# Create comparison plot
ggplot() +
  # True distribution (black line)
  geom_function(
    fun = dlnorm,
    args = list(meanlog = meanlog, sdlog = sdlog),
    color = "black",
    linewidth = 1.2,
    xlim = c(0, 20)
  ) +
  # Empirical densities for truncated data
  geom_density(
    data = truncated_data,
    aes(x = delay, fill = truncation),
    alpha = 0.6
  ) +
  # Mean lines
  geom_vline(
    data = truncated_data %>% distinct(truncation, mean_delay),
    aes(xintercept = mean_delay, color = truncation),
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_fill_manual(values = c("Truncated at 6 days" = "#E31A1C",
                               "Truncated at 10 days" = "#FF7F00",
                               "Truncated at 15 days" = "#1F78B4")) +
  scale_color_manual(values = c("Truncated at 6 days" = "#E31A1C",
                                "Truncated at 10 days" = "#FF7F00",
                                "Truncated at 15 days" = "#1F78B4")) +
  labs(
    title = "Effect of Right Truncation on Delay Distributions",
    x = "Delay (days)",
    y = "Density",
    fill = "Truncated data",
    color = "Mean delay",
    caption = "Black line: True log-normal distribution\nColored densities: Truncated data\nDashed lines: Mean delays"
  ) +
  xlim(0, 20) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  guides(colour = "none")

df_realtime <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  filter(hosp_time <= 70)

# truncated mean delay
mean(df_realtime$onset_to_hosp)

# compare with the mean delay over the full outbreak
mean(df$hosp_time - df$onset_time, na.rm=TRUE)

res <- nfidd_sample(mod,
                    data = list(
                      n = nrow(na.omit(df_realtime)),
                      y = na.omit(df_realtime)$onset_to_hosp
                    )
)

res

res %>%
  summarise_lognormal()

tmod <- nfidd_cmdstan_model("truncated-delay-model")
tmod

tres <- nfidd_sample(tmod,
                     data = list(
                       n = nrow(df_realtime),
                       onset_to_hosp = df_realtime$onset_to_hosp,
                       time_since_onset = 70 - df_realtime$onset_time
                     )
)

tres
tres %>% summarise_lognormal()
