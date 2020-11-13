
# pop counts plot
# stacked percentage bar plot


library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)


plot_drug <-
  melt(pop[2, , ]) %>%
  mutate(state = factor(state,
                        levels = c("Dead", "Progressive_disease", "Asymptomatic_disease")))

ggplot(plot_drug, aes(fill = state, y = value, x = cycle)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = 1:44) +
  ylab("Percentage of cohort")

plot_free <-
  melt(pop[1, , ]) %>%
  mutate(state = factor(state,
                        levels = c("Dead", "Progressive_disease", "Asymptomatic_disease")))

ggplot(plot_free, aes(fill = state, y = value, x = cycle)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = 1:44) +
  ylab("Percentage of cohort")


# combine
total_dat <- rbind(cbind(strat = "drug", plot_drug),
                   cbind(strat = "free", plot_free))

facet_plot <-
  ggplot(total_dat, aes(fill = state, y = value, x = cycle)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  scale_x_continuous(breaks = seq(0, 44, by = 2)) +
  ylab("Percentage of cohort") +
  facet_grid(.~strat) +
  theme_bw()

facet_plot

# ggsave(facet_plot, file = "images/state_pop_over_time.png")


