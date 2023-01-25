# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(car)
library(pals)
library(ggseqlogo)

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

# set constant
BASES <- c('A', 'C', 'G', 'U')



############################# Fig. 2 A #############################
data_lib <- data %>% filter(lib == 1) %>% select('seq', 'rTR')

# split into N40
data_lib_split <- as.data.frame(matrix(nrow = nrow(data_lib),
  ncol = 40, unlist(strsplit(data_lib$seq, '')), byrow = T),
  stringsAsFactors = FALSE)

data_lib_split$rTR <- data_lib$rTR

# set formula
full_formula <- as.formula(
  paste('rTR ~', paste(paste0('V', 1:39), '+', collapse = ' '), 'V40'))

# test with all positions
res.aov_full <- Anova(lm(data = data_lib_split, formula = full_formula), type = 'II')

# extract Sum Sq
SOS <- res.aov_full[1]
SOS$SS_rel <- SOS$'Sum Sq' / sum(SOS$'Sum Sq')
SOS$position <- 1:41
SOS$id <- c(paste0('0', 1:9), paste0(10:41))

data_plot_all <- data.frame(position = 1:75, SS_rel = NA, id = c(paste0('0', 1:9), paste0(10:75)))
part_of_library <- c(3:27, seq(33, 75, 3))
data_plot_all$SS_rel[part_of_library] <- SOS$SS_rel[1:40]

# get max und sum
print(paste0('maximum SOS: ', max(SOS$SS_rel)*100))
print(paste0('total SOS: ', sum(SOS$SS_rel)*100))

# remove first 2 positions
data_plot_cut <- data_plot_all[-c(1:2), ]

p <- ggplot(data = data_plot_cut, aes(x = id, y = SS_rel * 100)) +
  geom_bar(stat = 'identity', fill = 'grey90', colour = 'black', width = 0.80) +
  scale_x_discrete('', expand = c(0, 0)) +
  scale_y_continuous('contribution to rTR\nvariance (%)', limits = c(0, 2.5),
    expand = c(0, 0)) +
  theme_SH()  +
  coord_cartesian(clip = 'off')

ggsave('Fig_02_A_anova.png', plot = p,
  width = 10, height = 3, units = c('in'), scale = 0.8)



############################# Fig. 2 B #############################
df_bases_expression_abs <- as.data.frame(matrix(nrow = 4, ncol = 75))
names(df_bases_expression_abs) <- 1:75
rownames(df_bases_expression_abs) <- BASES

base_effects <- list()

# loop through data and bases
data_lib <- data %>% filter(lib == 1)
for (i in 1:75) {
  for (k in 1:4) {
    b <- BASES[k]
    rows <- substring(data_lib$full, i, i) == b
    df_bases_expression_abs[b, i] <- log2(mean(data_lib$rTR[rows]) / mean(data_lib$rTR[!rows]))
    print(paste0(b,i))
  }
}

# fill NA with 0
df_bases_expression_abs[is.na(df_bases_expression_abs)] <- NA
df_bases_expression_abs$Base <- BASES
base_effects <- df_bases_expression_abs

df_long_abs <- reshape2::melt(base_effects, id = 'Base')
names(df_long_abs) <- c('Base', 'Position', 'rTR')

df_long_abs <- df_long_abs[!(df_long_abs$Position %in% c(1:2)), ]

p <- ggplot(df_long_abs, aes(x = Position, y = Base, fill = rTR)) + 
  geom_tile() +
  labs(x = '', y = 'base') +
  scale_fill_gradientn('effect\non rTR (log2FC)', 
    colours = rev(ocean.curl(256)), limits = c(-0.5, +0.5), na.value = NA) +
  theme_SH() +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

ggsave('Fig_02_B_base_effects.png', plot = p,
  width = 10, height = 3, units = c('in'), scale = 0.8)



############################# Fig. 2 C #############################
data_lib <- data %>%
  filter(lib == 1) %>%
  select(full, rTR)

# add expr and bin
df_low <- data_lib %>% filter(rTR < 0.5)
df_high <- data_lib %>% filter(rTR >= 0.5) 

n <- 75

# variants below 50% activity
# count base occurence in df_bin
df_count_bases_low <- as.data.frame(matrix(nrow = 4, ncol = n+1))
names(df_count_bases_low) <- c('Base', paste0('p', 1:n))
rownames(df_count_bases_low) <- BASES
df_count_bases_low$Base <- BASES

# loop through data and bases
for (i in 1:n) {
  for (k in 1:4) {
    b <- BASES[k]
    df_count_bases_low[b,i+1] <- sum(substring(df_low$full, i, i) == b)
    print(paste0(b,i))
  }
}

# make columns relative
for (i in 1:n) {
  df_count_bases_low[,i+1] <- round(df_count_bases_low[,i+1] / sum(df_count_bases_low[,i+1]), 3)
  print(i)
}

# variants above 50% activity
# count base occurence in df_bin
df_count_bases_high <- as.data.frame(matrix(nrow = 4, ncol = n+1))
names(df_count_bases_high) <- c('Base', paste0('p', 1:n))
rownames(df_count_bases_high) <- BASES
df_count_bases_high$Base <- BASES

# loop through data and bases
for (i in 1:n) {
  for (k in 1:4) {
    b <- BASES[k]
    df_count_bases_high[b,i+1] <- sum(substring(df_high$full, i, i) == b)
    print(paste0(b,i))
  }
}

# make columns relative
for (i in 1:n) {
  df_count_bases_high[,i+1] <- round(df_count_bases_high[,i+1] / sum(df_count_bases_high[,i+1]), 3)
  print(i)
}

# get difference
df_count_bases_diff <- as.data.frame(matrix(nrow = 4, ncol = n+1))
names(df_count_bases_diff) <- c('Base', paste0('p', 1:n))
rownames(df_count_bases_diff) <- BASES
df_count_bases_diff$Base <- BASES
df_count_bases_diff[, 2:(n+1)] <- log2(df_count_bases_high[, 2:(n+1)] / df_count_bases_low[, 2:(n+1)])
df_count_bases_diff$Base <- NULL
df_count_bases_diff <- as.matrix(df_count_bases_diff)

df_count_bases_diff <- df_count_bases_diff[, -c(1, 2)]

p <- ggplot() +
  geom_logo(data = df_count_bases_diff, method = 'custom', seq_type = 'rna') +
  theme_SH() +
  scale_y_continuous('enrichment\namongst strong\nvariants (log2 FC)',
    limits = c(-1.6, 1.2), breaks = round(seq(-1.6, 1.2, 0.4), 1), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = 1:73) +
  coord_cartesian(clip = 'off')

ggsave('Fig_02_C_enrichment.png', plot = p,
  width = 10, height = 3, units = c('in'), scale = 0.8)
