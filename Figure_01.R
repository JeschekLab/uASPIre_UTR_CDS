# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(RColorBrewer)
library(scales)

# read data
data <- read.table('./data/data_combined.txt',
  header = T,
  sep = '\t',
  colClasses = c(
    rep('character', 5),
    rep('integer', 2),
    rep('numeric', 21),
    'integer'),
  check.names = F)



############################# Fig. 1 E #############################
### lineplot ###
current_library <- 1
print(paste0('Lineplot of library ', current_library))

# set temp data
data_lib <- data %>%
  select('0', '95', '225', '290', '360', '480', 'rTR', 'lib', 'sum', 'seq') %>%
  filter(lib == current_library) %>%
  arrange(desc(sum))

# make fraction flipped as %
data_lib[, c('0', '95', '225', '290', '360', '480')] <- 
  data_lib[, c('0', '95', '225', '290', '360', '480')] * 100

# subset for line plot
n <- 1000
subsample <- data_lib[1:n, ] %>%
  arrange(rTR)

# make matrix for diverging color
M <- as.matrix(subsample[, 1:6])
rownames(M) <- 1:n

# adapt matrix for diverging color
n_points <- length(seq(0, 480, 5))
plot_data <- as.data.frame(matrix(nrow = n, ncol = n_points),
  stringsAsFactors = F)

# change col and rownames
colnames(plot_data) <- seq(0, 480, 5)
rownames(plot_data) <- 1:n

# vectors
y_known <- c(0, 95, 225, 290, 360, 480)
x_known <- (y_known / 5) + 1

for (i in 1:5) {
	x1 <- x_known[i]
	x2 <- x_known[i + 1]
	data_x1 <- M[, i]
	data_x2 <- M[, i + 1]
	len <- x2 - x1 + 1
	for (j in 1:n) {
		plot_data[j,x1:x2] <- seq(data_x1[j], data_x2[j], length = len)
	}
}

plot_data$ID <- 1:nrow(plot_data)

plot_melted <- reshape2::melt(plot_data, id = 'ID')
names(plot_melted) <- c('ID', 'tp', 'fraction')

# plot data
p <- ggplot(plot_melted, aes(x = tp, y = fraction, group = ID)) + 
  geom_line(col = 'black', alpha = 0.5, linewidth = 0.1)+
  scale_x_discrete('time after induction (h)',
    breaks = seq(0, 480, 120), labels = c(0, 2, 4, 6, 8), expand = c(0,0)) +
  scale_y_continuous('fraction flipped (%)',
    breaks = seq(0, 100, 25), limits = c(0, 100), expand = c(0,0)) +
  theme_SH() +
  theme(legend.position = 'none') +
  coord_cartesian(clip = 'off')

ggsave('Fig_01_E_lines.png',
  plot = p, width = 3, height = 3, units = c('in'), scale = 1)



############################# Fig. 1 F #############################
### heatmap ###
current_library <- 1
print(paste0('Heatmap of library ', current_library))

# set temp data
data_lib <- data %>%
  filter(lib == current_library) %>%
  arrange(desc(sum)) %>%
  select('0', '95', '225', '290', '360', '480', 'rTR')

# make fraction flipped as %
data_lib[, c('0', '95', '225', '290', '360', '480')] <- 
  data_lib[, c('0', '95', '225', '290', '360', '480')] * 100

# subset for line plot
n <- 1000 # nrow(data_lib) # all data, can be adjusted for faster plotting
subsample <- data_lib[1:n, ] %>% arrange(rTR)

# make matrix
M <- as.matrix(subsample[, 1:6])
rownames(M) <- 1:n

# adapt matrix for diverging color
n_points <- length(seq(0, 480, 10))
plot_data <- as.data.frame(matrix(nrow = n, ncol = n_points),
  stringsAsFactors = F)

# change col and rownames
colnames(plot_data) <- seq(0, 480, 10)
rownames(plot_data) <- 1:n

# vectors
y_known <- c(0, 95, 225, 290, 360, 480)
x_known <- (y_known / 10) + 1

for (i in 1:5) {
  print(i)
  x1 <- x_known[i]
  x2 <- x_known[i + 1]
  data_x1 <- M[, i]
  data_x2 <- M[, i + 1]
  len <- x2 - x1 + 1
  for (j in 1:n) {
    plot_data[j,x1:x2] <- seq(data_x1[j], data_x2[j], length = len)
  }
}

plot_data$ID <- 1:nrow(plot_data)
plot_data$rTR <- subsample$rTR

# split into low and high
melted <- reshape2::melt(plot_data %>% select(-rTR), id = 'ID')

names(melted) <- c('ID', 'tp', 'fraction')

# add logfraction
melted$logfraction <- log2(melted$fraction + 1)

# plot data
cols <- rev(brewer.pal(n = 11, name = 'Spectral'))

# plot
p <- ggplot(melted, aes(
    x = as.factor(tp), y = as.factor(ID), fill = fraction)) + 
  geom_tile() + 
  scale_fill_gradientn(
    name = 'fraction\nflipped (%)',
    colours = cols, 
    values = rescale(c(0, 2.5, 5.3, 8.7, 12.7, 17.6, 23.8, 32.4, 46.8, 100)),
    guide = 'colorbar',
    limits=c(0, 100)) +
  scale_x_discrete('time after induction (h)',
    breaks = seq(0, 480, 120), labels = c(0, 2, 4, 6, 8), expand = c(0,0)) +
  scale_y_discrete('ranked variants',
    breaks = c(1, seq(n/4, n, n/4)),
    labels = c(1, seq(floor(nrow(data_lib)/4), floor(nrow(data_lib)), floor(nrow(data_lib)/4))),
    expand = c(0,0)) +
  theme_SH() +
  theme(panel.grid.major = element_blank()) +
  coord_cartesian(clip = 'off')

ggsave('Fig_01_F_heatmap.png',
  plot = p, width = 10, height = 10)



############################# Fig. 1 H #############################
### IFP to GFP ###
df_GFP <- read.table(
  './data/data_calibration_GFP.txt',
  header = T, sep = '\t', colClasses = c('character', 'numeric', 'numeric'))

# copy the lowest GFP value and add IFP 0 for loess fit
df_GFP_0 <- rbind(df_GFP, df_GFP[df_GFP$GFP == min(df_GFP$GFP), ])
df_GFP_0$AUC[32] <- 0
rownames(df_GFP_0) <- 1:nrow(df_GFP_0)

# compute pseudo R2
GFP_loess <- loess(data = df_GFP_0, formula = GFP ~ AUC)
ss.dist <- sum(scale(df_GFP_0$GFP, scale = FALSE)^2)
ss.resid <- sum(resid(GFP_loess)^2)
paste0('Pseudo R2 = ', round((1 - ss.resid / ss.dist), 4))

p <- ggplot(df_GFP_0, aes(x = AUC/480, y = GFP)) +
  geom_point(size = 2, color = 'black', shape = 1) +
  scale_x_continuous('IFP (-)',
    limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous('slope sfGFP0-290min (a.u.)',
    limits = c(0, 70), breaks = seq(0, 70, 10), expand = c(0, 0)) +
  geom_smooth(method = 'loess', se = F, fullrange = T,
    linewidth = 0.5, color = 'black') +
  theme_SH() +
  theme(legend.position = 'none') +
  coord_cartesian(clip = 'off')

ggsave('Fig_01_H_IFP_GFP.png', plot = p,
  width = 3, height = 3, units = c('in'), scale = 1)



############################# Fig. 1 I #############################
### histogram ###
data_lib <- data %>% filter(lib == 1)

# histogram
p <- ggplot(data_lib, aes(x = rTR)) +
  geom_histogram(binwidth = 0.04, fill = 'grey90', color = 'black') +
  scale_y_continuous('variant count',
    trans = 'log10', expand = c(0,0), limits = c(1, 10^5),
    breaks = c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5),
    labels = c(100, 101, 102, 103, 104, 105)) +
  scale_x_continuous('rTR (-)',
    limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = c(0,0)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

ggsave('Fig_01_I_hist.png', plot = p, width = 3,
  height = 3, units = c('in'), scale = 1)
