library(ggplot2)

# df1 <- read.csv('~/GitHub/CommunityImmunity/output/test_100sims_t1000/summary.csv')
# df1 <- read.csv('~/GitHub/CommunityImmunity/output/test_100sims_t10000/summary.csv')
df1 <- read.csv('~/GitHub/CommunityImmunity/output/20211013_0856_1000sims_t5000/summary.csv')
df2 <- read.csv("~/GitHub/CommunityImmunity/output/test_1sims_t10000/param_evolve.csv")

final_rate <- df1$m_final

hist(log10(final_rate))
boxplot(log10(final_rate))
plot(ecdf(log10(final_rate))) # CPD

p1 <- ggplot(NULL, aes(x=log10(final_rate))) +
  geom_histogram() +
  xlab('Final Average Mutation rate') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(5,"pt"),)
# ggsave(plot = p1, filename = '~/20210324_labmeetingp1.png', width = 7, height = 5)

p2 <- ggplot(df2, aes(x=X, y=m)) +
  geom_line() +
  scale_y_continuous('average mutation rate',limits=c(1e-7,10^-6)) +
  xlab('time') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(5,"pt"),)

# ggsave(plot = p2, filename = '~/20210324_labmeetingp2.png', width = 7, height = 5)
