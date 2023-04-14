# set wd
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(car)
library(pals)
library(ggseqlogo)

# read data
print('Reading data ...')
data <- read.table(
  file = './data/data_combined.txt',
  header = T,
  sep = '\t',
  colClasses = c(
    rep('character', 5),
    rep('integer', 2),
    rep('numeric', 21),
    'integer'),
  check.names = F)

# set constants
BASES <- c('A', 'C', 'G', 'U')
LABELS <- c(
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'U', 'G', 'C', 'G',
  'N', 'G', 'C', 'N', 'C', 'U', 'N', 'G', 'U', 'N', 'G', 'U', 'N', 'A', 'U',
  'H', 'C', 'G', 'N', 'C', 'U', 'N', 'U', 'C', 'N', 'C', 'G', 'N', 'G', 'U',
  'N', 'A', 'C', 'N', 'G', 'A', 'Y', 'G', 'C', 'N', 'A', 'C', 'N')



############################# Fig. 2 A #############################
### ANOVA ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 2 A: ANOVA of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library) %>%
  select('seq', 'rTR')

# split into N40
data_lib_split <- matrix(nrow = nrow(data_lib),
  ncol = 40, unlist(strsplit(data_lib$seq, '')), byrow = T) %>%
  as.data.frame(., stringsAsFactors = FALSE)

# add rTR
data_lib_split <- data_lib_split %>%
  mutate(rTR = data_lib$rTR)

# set formula for ANOVA
full_formula <- as.formula(
  paste('rTR ~', paste(paste0('V', 1:39), '+', collapse = ' '), 'V40'))

# run ANOVA
res.aov_full <- Anova(
  lm(
    data = data_lib_split,
    formula = full_formula),
  type = 'II')

# extract sum of squares
SOS <- res.aov_full[1] %>%
  as.data.frame(., stringsAsFactors = FALSE) %>%
  rename(SOS = "Sum Sq") %>%
  mutate(SS_rel = SOS / sum(SOS) * 100) %>%
  filter(rownames(.) != 'Residuals') %>%
  mutate(position = c(3:27, seq(33, 75, 3))) %>%
  mutate(id = c(paste0('0', 1:9), paste0(10:40)))

# get max und sum
print(paste0('maximum SOS: ', round(max(SOS$SS_rel), 2)))
print(paste0('total SOS: ', round(sum(SOS$SS_rel), 2)))

# generate plot
p <- ggplot(data = SOS, aes(x = position, y = SS_rel)) +
  geom_bar(
    stat = 'identity',
    fill = 'grey90',
    colour = 'black',
    width = 0.80) +
  scale_x_continuous(
    name = '',
    breaks = 3:75,
    labels = LABELS,
    expand = c(0, 0)) +
  scale_y_continuous(
    name = 'contribution to rTR\nvariance (%)',
    limits = c(0, 2.5),
    expand = c(0, 0)) +
  theme_SH()  +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_02_A_anova.pdf',
  width = 10,
  height = 3,
  scale = 0.8)



############################# Fig. 2 B #############################
### single base effect ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 2 B: single base effect of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library) %>%
  select(full, rTR)

# initialize empty data frame
df_bases <- as.data.frame(matrix(nrow = 4, ncol = 75))
names(df_bases) <- 1:75
rownames(df_bases) <- BASES

# loop through data and bases and calculate absolute base effect
for (pos in 1:75) {
  for (base in 1:4) {
    b <- BASES[base]
    rows <- substring(data_lib$full, pos, pos) == b
    df_bases[b, pos] <- log2(
      mean(data_lib$rTR[rows]) / mean(data_lib$rTR[!rows])
      )
    print(paste0('Base: ', b, ', position: ', pos))
  }
}

# add base column
df_bases$base <- BASES

# make long format
data_plot <- df_bases %>%
  pivot_longer(
    cols = -base,
    names_to = 'position',
    values_to = 'rTR') %>%
  mutate(position = as.integer(position)) %>%
  filter(position >= 3)

# generate plot
p <- ggplot(data_plot, aes(x = position, y = base, fill = rTR)) + 
  geom_tile() +
  scale_fill_gradientn(
    name = 'effect\non rTR (log2 FC)', 
    colours = rev(ocean.curl(256)),
    limits = c(-0.5, +0.5),
    na.value = NA) +
  theme_SH() +
  scale_x_continuous(
    name = '',
    breaks = 3:75,
    labels = LABELS,
    expand = c(0, 0)) +
  scale_y_discrete(
    name = 'base',
    limits = rev,
    expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_02_B_base_effects.pdf',
  width = 10,
  height = 3,
  scale = 0.8)



############################# Fig. 2 C #############################
### enrichment ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 2 C: enrichment of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == 1) %>%
  select(full, rTR)

# split into low and high rTR data sets
df_low <- data_lib %>% filter(rTR < 0.5)
df_high <- data_lib %>% filter(rTR >= 0.5)

# find length of library
N <- nchar(data_lib$full[1])

# count base occurence in df_low
df_count_low <- as.data.frame(matrix(nrow = 4, ncol = N+1))
names(df_count_low) <- c('Base', paste0('p', 1:N))
rownames(df_count_low) <- BASES
df_count_low$Base <- BASES

# loop through data and bases
for (pos in 1:N) {
  for (k in 1:4) {N
    b <- BASES[k]
    df_count_low[b, pos+1] <- sum(substring(df_low$full, pos, pos) == b) /
      nrow(df_low)
    print(paste0('Base: ', b, ', position: ', pos))
  }
}

# count base occurence in df_high
df_count_high <- as.data.frame(matrix(nrow = 4, ncol = N+1))
names(df_count_high) <- c('Base', paste0('p', 1:N))
rownames(df_count_high) <- BASES
df_count_high$Base <- BASES

# loop through data and bases
for (pos in 1:N) {
  for (k in 1:4) {
    b <- BASES[k]
    df_count_high[b, pos+1] <- sum(substring(df_high$full, pos, pos) == b) /
      nrow(df_high)
    print(paste0('Base: ', b, ', position: ', pos))
  }
}

# get difference
df_count_diff <- log2(df_count_high[, -(1:3)] / df_count_low[, -(1:3)])

# generate plot
p <- ggplot() +
  geom_logo(
    data = as.matrix(df_count_diff),
    method = 'custom',
    seq_type = 'rna') +
  theme_SH() +
  scale_y_continuous(
    name = 'enrichment\namongst strong\nvariants (log2 FC)',
    limits = c(-1.6, 1.2),
    breaks = seq(-1.6, 1.2, 0.4),
    expand = c(0, 0)) +
  scale_x_continuous(
    name = '',
    breaks = 1:73,
    labels = LABELS,
    expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_02_C_enrichment.pdf',
  width = 10,
  height = 3,
  scale = 0.8)

# done!
