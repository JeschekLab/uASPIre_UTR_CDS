# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(RColorBrewer)
library(scales)

# read data
print('Reading data ...')
data <- read.table('./data/data_combined.txt',
  header = T,
  sep = '\t',
  colClasses = c(
    rep('character', 5),
    rep('integer', 2),
    rep('numeric', 21),
    'integer'),
  check.names = F)

# set constants
TPS <- c('0', '95', '225', '290', '360', '480')
N <- 1000 # lines for line plot



############################# Fig. 1 E #############################
### lineplot ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 1 E: Lineplot of library ', current_library))

# filter data and make subset
data_lib <- data %>%
  filter(lib == current_library) %>%
  arrange(desc(sum)) %>%
  select('seq', all_of(TPS))

# make subset for line plot
data_lib <- data_lib %>%
  head(., N)

# make longer format and fraction flipped as %
data_plot <- data_lib %>%
  pivot_longer(
    cols = -seq,
    names_to = 'tp',
    values_to = 'fraction') %>%
  mutate(tp = as.numeric(tp)) %>%
  mutate(fraction = fraction * 100)

# generate plot
p <- ggplot(data_plot, aes(x = tp, y = fraction, group = seq)) + 
  geom_line(
    col = 'black',
    alpha = 0.5,
    linewidth = 0.1)+
  scale_x_discrete(
    name = 'time after induction (h)',
    breaks = seq(0, 480, 120),
    labels = c(0, 2, 4, 6, 8),
    expand = c(0,0)) +
  scale_y_continuous(
    name ='fraction flipped (%)',
    breaks = seq(0, 100, 25),
    limits = c(0, 100),
    expand = c(0,0)) +
  theme_SH() +
  theme(legend.position = 'none') +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(plot = p, file = 'Fig_01_E_lines.png', width = 3, height = 3)



############################# Fig. 1 F #############################
### heatmap ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 1 F: Heatmap of library ', current_library))

# filter data and make subset
data_lib <- data %>%
  filter(lib == current_library) %>%
  arrange(desc(sum)) %>%
  select('seq', all_of(TPS), rTR)

# make subset for heatmap and arrange by rTR
data_lib <- data_lib %>%
  head(., N) %>%
  arrange(rTR)

# make matrix
M <- as.matrix(data_lib[, TPS])
rownames(M) <- 1:N

# adapt matrix to add more fields for better color resolution
n_points <- length(seq(0, 480, 10))
plot_data <- as.data.frame(matrix(nrow = N, ncol = n_points),
  stringsAsFactors = F)

# change colnames and rownames
colnames(plot_data) <- seq(0, 480, 10)
rownames(plot_data) <- 1:N

# vectors
y_known <- as.numeric(TPS)
x_known <- (y_known / 10) + 1

# loop through time points and calculate intermediate values
for (i in 1:5) {
  print(i)
  x1 <- x_known[i]
  x2 <- x_known[i + 1]
  data_x1 <- M[, i]
  data_x2 <- M[, i + 1]
  len <- x2 - x1 + 1
  for (j in 1:N) {
    plot_data[j,x1:x2] <- seq(data_x1[j], data_x2[j], length = len)
  }
}

plot_data$ID <- 1:nrow(plot_data)
plot_data$rTR <- data_lib$rTR

# make long format
plot_data <- plot_data %>%
  select(-rTR) %>%
  pivot_longer(cols = -ID,
    names_to = 'tp',
    values_to = 'fraction') %>%
  mutate(tp = as.numeric(tp)) %>%
  mutate(fraction = fraction * 100)

# add logfraction
plot_data <- plot_data %>%
  mutate(logfraction = log2(fraction + 1)) # add 1 for log scale

# set colors
cols <- rev(brewer.pal(n = 11, name = 'Spectral'))

# generate plot
p <- ggplot(plot_data,
  aes(x = as.factor(tp), y = as.factor(ID), fill = fraction)) + 
  geom_tile() + 
  scale_x_discrete(
    name = 'time after induction (h)',
    breaks = seq(0, 480, 120),
    labels = c(0, 2, 4, 6, 8),
    expand = c(0,0)) +
  scale_y_discrete(
    name = 'ranked variants',
    breaks = c(1, seq(N/4, N, N/4)),
    labels = c(1, seq(
      floor(nrow(data_lib)/4),
      floor(nrow(data_lib)),
      floor(nrow(data_lib)/4))),
    expand = c(0,0)) +
  scale_fill_gradientn(
    name = 'fraction\nflipped (%)',
    colours = cols, 
    values = rescale(c(0, 2.5, 5.3, 8.7, 12.7, 17.6, 23.8, 32.4, 46.8, 100)),
    guide = 'colorbar',
    limits = c(0, 100)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(plot = p, file = 'Fig_01_F_heatmap.png', width = 10, height = 10)



############################# Fig. 1 H #############################
### IFP to GFP ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 1 H: IFP to rTR of library ', current_library))

# read calibration data
df_GFP <- read.table(
  file = './data/data_calibration_GFP.txt',
  header = T,
  sep = '\t',
  colClasses = c('character', 'numeric', 'numeric'))

# copy the lowest GFP value and add AUC 0 for loess fit
df_GFP <- rbind(df_GFP, df_GFP[df_GFP$GFP == min(df_GFP$GFP), ])
df_GFP$AUC[32] <- 0
rownames(df_GFP) <- 1:nrow(df_GFP)

# compute pseudo R2
GFP_loess <- loess(data = df_GFP, formula = GFP ~ AUC)
ss.dist <- sum(scale(df_GFP$GFP, scale = FALSE)^2)
ss.resid <- sum(resid(GFP_loess)^2)
paste0('Pseudo R2 = ', round((1 - ss.resid / ss.dist), 4))

# generate plot
p <- ggplot(df_GFP, aes(x = AUC/480, y = GFP)) +
  geom_point(
    size = 2,
    color = 'black',
    shape = 1) +
  scale_x_continuous(
    name = 'IFP (-)',
    limits = c(0, 1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = 'slope sfGFP0-290min (a.u.)',
    limits = c(0, 70),
    breaks = seq(0, 70, 10),
    expand = c(0, 0)) +
  geom_smooth(
    method = 'loess',
    se = F,
    fullrange = T,
    linewidth = 0.5,
    color = 'black') +
  theme_SH() +
  theme(legend.position = 'none') +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave('Fig_01_H_IFP_GFP.png', plot = p,
  width = 3, height = 3, units = c('in'), scale = 1)



############################# Fig. 1 I #############################
### histogram ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 1 I: Histogram of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == 1)

# generate plot
p <- ggplot(data_lib, aes(x = rTR)) +
  geom_histogram(
    binwidth = 0.04,
    fill = 'grey90',
    color = 'black') +
  scale_y_continuous(
    name = 'variant count',
    trans = 'log10',
    expand = c(0,0),
    limits = c(1, 10^5),
    breaks = c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5),
    labels = c(100, 101, 102, 103, 104, 105)) +
  scale_x_continuous(
    name = 'rTR (-)',
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.25),
    expand = c(0,0)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave('Fig_01_I_hist.png', plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)
