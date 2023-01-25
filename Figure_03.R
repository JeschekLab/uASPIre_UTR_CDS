# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(caret)
library(h2o)
library(viridis)

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



############################# Fig. 3 A #############################
### effect of folding ###
data_lib <- data %>%
  filter(lib == 1)

df_corr_structure <- data.frame()

rTR <- data_lib$rTR

for (current_model in c('mfeC', 'mfeT', 'efeC', 'efeT', 'accC', 'accT')) {
  current_energy <- data_lib[, current_model]
  p <- cor(rTR, current_energy, method = 'spearman')
  R <- cor(rTR, current_energy, method = 'pearson')
  df_temp <- data.frame(
    model = current_model,
    cor = c(p, R),
    stat = c('p', 'R'))
  df_corr_structure <- rbind(df_corr_structure, df_temp)
  print(current_model)
}

# calculate GC
p <- cor(data_lib$rTR, data_lib$GC_all, method = 'spearman')
R <- cor(data_lib$rTR, data_lib$GC_all, method = 'pearson')
df_temp <- data.frame(
  model = 'GC',
  cor = c(p, R),
  stat = c('p', 'R'))

df_corr_structure <- rbind(df_temp, df_corr_structure)
df_corr_structure$explainability <- (df_corr_structure$cor^2)*100

# plot
p <- ggplot(data = df_corr_structure,
    aes(y = explainability, x = model, fill = stat)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  scale_fill_discrete(name = '') +
  scale_x_discrete('') +
  scale_y_continuous('p2 / R2 (%)',
    limits = c(0, 40), expand = c(0, 0)) +
  theme_SH() +
  coord_flip(clip = 'off')

ggsave('Fig_03_A_folding.png', plot = p,
  width = 3, height = 3, units = c('in'), scale = 1)



############################# Fig. 3 B #############################
### scatter efeC ###
data_lib <- data %>%
  filter(lib == 1)

# plot data
p <- ggplot(data_lib, aes(x = efeC, y = rTR)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_x_continuous('efeC (kcal x mol-1)',
    limits = c(-25, 0), expand = c(0, 0), breaks = seq(-25, 0, 5)) +
  scale_y_continuous('rTR (-)', 
    limits = c(0, 1), expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

ggsave('Fig_03_B_efeC.png', plot = p,
  width = 5, height = 5, units = c('in'), scale = 1)



############################# Fig. 3 C #############################
### accC/accT scanning ###
df_corr_acc <- data.frame()

for (current_mode in c('accT', 'accC')) {
  df_wider <- read.table(paste0('./data/data_', current_mode, '_10.txt'),
      sep = '\t', header = T)

  # add library
  data_lib <- cbind(data, df_wider) %>% filter(lib == 1)

  # calculate correlation
  rTR <- data_lib$rTR

  for (current_position in 1:71) {
    current_name <- paste0('X', current_position)
    current_acc <- data_lib[, current_name]
    p <- cor(rTR, current_acc, method = 'spearman')
    df_temp <- data.frame(
      pos = current_position,
      cor = p,
      stat = 'p',
      mode = current_mode
      )
    df_corr_acc <- rbind(df_corr_acc, df_temp)
    print(current_name)
  }
}

# correct position and add center
df_corr_acc$pos_corrected <- df_corr_acc$pos - 28 + 4.5

# make points
p <- ggplot(df_corr_acc, aes(x = pos_corrected, y = cor, color = mode)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  theme_SH() +
  scale_color_discrete(name = '') +
  scale_x_continuous('center position of accessibility window',
    limits = c(-25, 50), 
    breaks = seq(-25, 50, 10),
    expand = c(0, 0)) +
  scale_y_continuous('p (-)',
    expand = c(0, 0), limits = c(0, -0.5), trans = 'reverse') +
  coord_cartesian(clip = 'off')

ggsave('Fig_03_C_acc.png', plot = p,
  width = 5, height = 3, units = c('in'), scale = 1)



############################# Fig. 3 D #############################
### hybridization scanning ###
data_lib <- data %>% filter(lib == 1)

# read hybridization
scanning_energy <- read.table(paste0('./data/data_hybridization.txt'),
  sep = '\t', header = T, colClasses = 'numeric')

n <- ncol(scanning_energy)

scanning_energy$rTR <- data_lib$rTR

scanning_correlation <- data.frame()

for (i in 1:n) {
  print(paste0('Current substring: ', i))
  R <- cor(scanning_energy[, i], scanning_energy$rTR, method = 'pearson')
  p <- cor(scanning_energy[, i], scanning_energy$rTR, method = 'spearman')
  df_temp <- data.frame(pos = i - 28, cor = c(R, p), stat = c('R', 'p'))
  scanning_correlation <- rbind(scanning_correlation, df_temp)
}

p <- ggplot(scanning_correlation, aes(x = pos+4, y = cor, fill = stat)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  theme_SH() +
  scale_x_continuous('center position of sliding window',
    limits = c(-25, 5),
    breaks = seq(-25, 5, by = 5)) +
  scale_y_reverse('p / R (-)',
    limits = c(0.2, -0.3),
    breaks = seq(0.2, -0.3, -0.1),
    expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

ggsave('Fig_03_D_hyb.png', plot = p,
  width = 3, height = 3, units = c('in'), scale = 1)



############################# Fig. 3 E #############################
### random forest ###
df_wider <- read.table('./data/data_accC_1.txt',
  sep = '\t', header = T)

data_lib <- cbind(data, df_wider) %>%
  filter(lib == 1) %>%
  select(seq, rTR, paste0('X', 1:80), hyb_opt,
    GC_all, mfeC, mfeT, efeC, efeT, accC, accT)

# make one hot
data_lib$onehot <- data_lib$seq %>% 
  gsub('A', '1000', .) %>%
  gsub('C', '0100', .) %>%
  gsub('G', '0010', .) %>%
  gsub('U', '0001', .)

df_one_hot <- as.data.frame(matrix(nrow = nrow(data_lib),
  ncol = nchar(data_lib$onehot[1]), unlist(strsplit(data_lib$onehot, '')), byrow = T),
  stringsAsFactors = FALSE)

for (i in 1:ncol(df_one_hot)) {
  class(df_one_hot[, i]) <- 'integer'
}

# acombine
df_RF <- cbind(data_lib, df_one_hot) %>%
  select(-seq, -onehot)

# split
set.seed(12345)
index <- createDataPartition(df_RF$rTR, p = 0.90, list = FALSE)

train <- df_RF[index, ]
test <- df_RF[-index, ]

# initialize h20 cluster
h2o.init()

# copy to cluster
train.h2o <- as.h2o(train)
test.h2o <- as.h2o(test)

y.dep <- 'rTR'

# overall model
rforest.model_cv <- h2o.randomForest(
  y = y.dep,
  training_frame = train.h2o,
  seed = 12345,
  keep_cross_validation_predictions = TRUE,
  nfolds = 10
  )

test$predict <- as.data.frame(h2o.predict(rforest.model_cv, test.h2o))$predict

# calculate correlation
cor(test$rTR, test$predict, method = 'pearson')^2

p <- ggplot(test, aes(x = rTR, y = predict)) +
  geom_point(size = 0.3, alpha = 0.3) +
  ggtitle(paste0(R_model)) +
  scale_x_continuous('measured rTR',
    limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous('predicted rTR', 
    limits = c(0, 1), expand = c(0, 0)) +
  geom_smooth(method = 'lm', se = T, fullrange = F) +
  theme_SH()

ggsave('03_RF_scatter.png', plot = p,
  width = 5, height = 5, units = c('in'), scale = 1)

### plot feature importance ###
varimp <- as.data.frame(h2o.varimp(rforest.model_cv))

varimp$name <- NA

# correct variable names for acc
varimp_acc <- varimp[varimp$variable %in% paste0('X', 1:80), ]
varimp <- varimp[!(varimp$variable %in% paste0('X', 1:80)), ]

varimp_acc$name <- paste0(
  'accC at position ',
  as.integer(gsub('X', '', varimp_acc$variable)) - 28)

# correct variable names for position
varimp_pos <- varimp[varimp$variable %in% paste0('V', 1:160), ]
varimp <- varimp[!(varimp$variable %in% paste0('V', 1:160)), ]

varimp_pos$ID <- as.integer(gsub('V', '', varimp_pos$variable))

df_temp <- data.frame(base = BASES, ID = 1:160, pos = rep(c(-25:-1, seq(6, 48, 3)), each = 4))

varimp_pos <- merge(varimp_pos, df_temp, by = 'ID', all = F)

varimp_pos$name <- paste0(
  varimp_pos$base,
  ' at position ',
  varimp_pos$pos)

varimp$name <- varimp$variable %>%
  gsub('hyb_opt', 'optimised hybridisation', .) %>%
  gsub('GC_all', 'GC-content', .)

varimp_all <- rbind(
  varimp %>% select(name, relative_importance, scaled_importance, percentage),
  varimp_pos %>% select(name, relative_importance, scaled_importance, percentage),
  varimp_acc %>% select(name, relative_importance, scaled_importance, percentage)
  ) %>%
  arrange(desc(percentage)) %>%
  mutate(relative_importance = round(relative_importance, 0)) %>%
  mutate(scaled_importance = round(scaled_importance, 3)) %>%
  mutate(percentage = round(percentage, 4)) %>%
  mutate(Rank = 1:nrow(.)) %>%
  select(Rank, name, relative_importance, scaled_importance, percentage) %>%
  rename('Name' = 'name') %>%
  rename('Relative importance' = 'relative_importance') %>%
  rename('Scaled importance' = 'scaled_importance') %>%
  rename('Percentage' = 'percentage')

# make barplots
n <- 10
labels <- head(varimp_all$Name, n)

p <- ggplot(data = head(varimp_all, n), aes(x = rev(Rank), y = Percentage * 100)) +
  geom_bar(stat = 'identity', fill = 'black') +
  scale_x_continuous('Ranked ', breaks =1:n, labels = rev(labels),
    expand = c(0, 0)) +
  scale_y_continuous('RF feature importance (%)',
    limits = c(0, 25), expand = c(0, 0)) +
  theme_SH() +
  coord_flip()

ggsave('03_RF_features.png', plot = p,
  width = 3, height = 3, units = c('in'), scale = 1)



############################# Fig. 3 F #############################
### heatmap efeC and hyb_opt ###
data_lib <- data %>%
  filter(lib == 1) %>%
  select(rTR, efeC, hyb_opt)

# bin efeC
data_lib$bin_efeC <- NA

data_lib$bin_efeC[data_lib$efeC <= -15] <- 1
data_lib$bin_efeC[data_lib$efeC > -15] <- 2
data_lib$bin_efeC[data_lib$efeC > -12.5] <- 3
data_lib$bin_efeC[data_lib$efeC > -10] <- 4
data_lib$bin_efeC[data_lib$efeC > -7.5] <- 5
data_lib$bin_efeC[data_lib$efeC > -5] <- 6

# bin hyb_opt
data_lib$bin_hyb <- NA

data_lib$bin_hyb[data_lib$hyb_opt <= -10] <- 1
data_lib$bin_hyb[data_lib$hyb_opt > -10] <- 2
data_lib$bin_hyb[data_lib$hyb_opt > -7.5] <- 3
data_lib$bin_hyb[data_lib$hyb_opt > -5] <- 4
data_lib$bin_hyb[data_lib$hyb_opt > -2.5] <- 5
data_lib$bin_hyb[data_lib$hyb_opt > 0] <- 6

# calculate mean
plot_data <- data_lib %>%
  group_by(bin_efeC, bin_hyb) %>%
  summarize(mean = mean(rTR))

# plot
p <- ggplot(plot_data, aes(x = bin_efeC, y = bin_hyb, fill = mean)) + 
  geom_tile() +
  theme_SH() +
  scale_fill_viridis(name = 'rTR (-)',
    option = 'turbo', limits = c(0, 1)) +
  scale_x_continuous('efeC (kcal x mol-1)',
    expand = c(0, 0), breaks = c(1:6),
    labels = c('<=-15', '>-15', '>-12.5', '>-10', '>-7.5', '>-5')) +
  scale_y_continuous('hybopt (kcal x mol-1)',
    expand = c(0, 0), breaks = c(1:6),
    labels = c('<=-10', '>-10', '>-7.5', '>-5', '>-2.5', '>0')) +
  coord_cartesian(clip = 'off')

ggsave('Fig_03_F_efeC_hyb.png', plot = p,
  width = 4, height = 3, units = c('in'), scale = 1)
