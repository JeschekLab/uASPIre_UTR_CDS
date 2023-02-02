# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(caret)
library(h2o)
library(viridis)

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



############################# Fig. 3 A #############################
### folding ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 3 A: folding of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library)

# initialize data frame
df_corr <- data.frame()

# loop through energies and calculate both p and R
for (current_model in
  c('mfeC', 'mfeT', 'efeC', 'efeT', 'accC', 'accT', 'GC_all')) {
  p <- cor(data_lib$rTR, data_lib[, current_model], method = 'spearman')
  R <- cor(data_lib$rTR, data_lib[, current_model], method = 'pearson')
  df_temp <- data.frame(
    model = current_model,
    stat = c('p', 'R'),
    cor = c(p, R))
  df_corr <- rbind(df_corr, df_temp)
  print(current_model)
}

# calculate explainability as R/p^2 in %
df_corr <- df_corr %>%
  mutate(explainability = (cor^2)*100) %>%
  mutate(stat2 = paste0(stat, 2))

# generate plot
p <- ggplot(data = df_corr,
  aes(y = explainability, x = model, fill = stat2)) +
  geom_bar(
    position = 'dodge',
    stat = 'identity',
    color = 'black') +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90')) +
  scale_x_discrete(
    name = '') +
  scale_y_continuous(
    name = 'p2 / R2 (%)',
    limits = c(0, 40),
    expand = c(0, 0)) +
  theme_SH() +
  coord_flip(
    clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_03_A_folding.pdf',
  width = 3,
  height = 3)



############################# Fig. 3 B #############################
### scatter efeC ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 3 B: scatter efeC of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library) %>%
  select(seq, rTR, efeC)

# print statement
print(paste0(
  'p = ', round(cor(data_lib$efeC, data_lib$rTR, method = 'spearman'), 3),
  '; R = ', round(cor(data_lib$efeC, data_lib$rTR, method = 'pearson'), 3)))

# plot data
p <- ggplot(data_lib, aes(x = efeC, y = rTR)) +
  geom_point(
    size = 0.3,
    alpha = 0.3) +
  scale_x_continuous(
    name = 'efeC (kcal x mol-1)',
    limits = c(-25, 0),
    breaks = seq(-25, 0, 5),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = 'rTR (-)', 
    limits = c(0, 1),
    expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(
    clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_03_B_efeC.png',
  width = 5,
  height = 5)



############################# Fig. 3 C #############################
### accC/accT scanning ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 3 C: scanning accC/accT of library ', current_library))

# initialize data frame
df_corr_acc <- data.frame()

# loop through data frames and calculate correlation
for (current_mode in c('accT', 'accC')) {
  # read additional data
  df_wider <- read.table(paste0('./data/data_', current_mode, '_10.txt'),
      sep = '\t', header = T)

  # add rTR
  data_lib <- cbind(data, df_wider) %>%
    filter(lib == current_library)

  # calculate correlation
  for (current_position in 1:71) {
    current_name <- paste0('X', current_position)
    current_acc <- data_lib[, current_name]
    p <- cor(data_lib$rTR, current_acc, method = 'spearman')
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
df_corr_acc <- df_corr_acc %>%
  mutate(pos_corrected = pos - 28 + 4.5)

# generate plot
p <- ggplot(df_corr_acc, aes(x = pos_corrected, y = -cor, color = mode)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  theme_SH() +
  scale_x_continuous(
    name = 'center position of accessibility window',
    limits = c(-25, 50), 
    breaks = seq(-25, 50, 10),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = 'p (-)',
    expand = c(0, 0),
    limits = c(0, 0.5)) +
  scale_color_manual(
    name = '',
    values = c('grey90', 'grey30')) +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_03_C_acc.pdf',
  width = 4,
  height = 3)



############################# Fig. 3 D #############################
### enrichment ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 3 D: hybridisation of library ', current_library))

# select data
data_lib <- data %>% filter(lib == current_library)

# read hybridization
scanning_energy <- read.table(file = './data/data_hybridization.txt',
  sep = '\t', header = T, colClasses = 'numeric')

# add rTR
scanning_energy$rTR <- data_lib$rTR

# initialize data frame
scanning_correlation <- data.frame()

# loop through all positions and calculate correlation
for (i in 1:ncol(scanning_energy)) {
  print(paste0('Current substring: ', i))
  R <- cor(scanning_energy[, i], scanning_energy$rTR, method = 'pearson')
  p <- cor(scanning_energy[, i], scanning_energy$rTR, method = 'spearman')
  df_temp <- data.frame(
    pos = i - 28,
    cor = c(R, p),
    stat = c('R', 'p'))
  scanning_correlation <- rbind(scanning_correlation, df_temp)
}

# generate plot
p <- ggplot(scanning_correlation, aes(x = pos+4, y = cor, fill = stat)) +
  geom_bar(
    position = 'dodge',
    stat = 'identity') +
  theme_SH() +
  scale_x_continuous(
    name = 'center position of sliding window',
    limits = c(-25, 5),
    breaks = seq(-25, 5, by = 5)) +
  scale_y_reverse(
    name = 'p / R (-)',
    limits = c(0.2, -0.3),
    breaks = seq(0.2, -0.3, -0.1),
    expand = c(0, 0)) +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90')) +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_03_D_hyb.pdf',
  width = 4,
  height = 3)



############################# Fig. 3 E #############################
### enrichment ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 3 E: feature importance of library ', current_library))

# read accC data
df_wider <- read.table(
  file = './data/data_accC_1.txt',
  sep = '\t', header = T)

# select lib random
data_lib <- cbind(data, df_wider) %>%
  filter(lib == current_library) %>%
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

# combine into one data frame
df_RF <- cbind(data_lib, df_one_hot) %>%
  select(-seq, -onehot)

# split randomly in test and train
set.seed(12345)
index <- createDataPartition(df_RF$rTR, p = 0.90, list = FALSE)

train <- df_RF[index, ]
test <- df_RF[-index, ]

# initialize h20 cluster
h2o.init()

# copy to cluster
train.h2o <- as.h2o(train)
test.h2o <- as.h2o(test)

# set variable that will be predicted
y.dep <- 'rTR'

# overall model
rforest.model_cv <- h2o.randomForest(
  y = y.dep,
  training_frame = train.h2o,
  seed = 12345,
  keep_cross_validation_predictions = TRUE,
  nfolds = 10
  )

# predict
test$predict <- as.data.frame(h2o.predict(rforest.model_cv, test.h2o))$predict

# calculate correlation
cor(test$rTR, test$predict, method = 'pearson')^2

# generate plot
p <- ggplot(test, aes(x = rTR, y = predict)) +
  geom_point(
    size = 0.3,
    alpha = 0.3) +
  ggtitle(paste0(R_model)) +
  scale_x_continuous(
    name = 'measured rTR',
    limits = c(0, 1),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = 'predicted rTR', 
    limits = c(0, 1),
    expand = c(0, 0)) +
  geom_smooth(
    method = 'lm',
    se = F,
    fullrange = F) +
  theme_SH()

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_03_E_RF_scatter.pdf',
  width = 5,
  height = 5)


### plot feature importance ###
# extract feature importance
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

df_temp <- data.frame(
  base = BASES,
  ID = 1:160,
  pos = rep(c(-25:-1, seq(6, 48, 3)), each = 4))

varimp_pos <- merge(varimp_pos, df_temp, by = 'ID', all = F)

varimp_pos$name <- paste0(
  varimp_pos$base,
  ' at position ',
  varimp_pos$pos)

varimp$name <- varimp$variable %>%
  gsub('hyb_opt', 'optimised hybridisation', .) %>%
  gsub('GC_all', 'GC-content', .)

varimp_all <- rbind(
  varimp %>%
    select(name, relative_importance, scaled_importance, percentage),
  varimp_pos %>%
    select(name, relative_importance, scaled_importance, percentage),
  varimp_acc %>%
    select(name, relative_importance, scaled_importance, percentage)
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

# generate plot
p <- ggplot(
  data = head(varimp_all, n), aes(x = rev(Rank), y = Percentage * 100)) +
  geom_bar(
    stat = 'identity',
    fill = 'grey30') +
  scale_x_continuous(
    name = 'Ranked ',
    breaks =1:n,
    labels = rev(labels),
    expand = c(0, 0)) +
  scale_y_continuous(
    name = 'rel. feature importance (%)',
    limits = c(0, 25),
    expand = c(0, 0)) +
  theme_SH() +
  coord_flip()

# save plot to file
ggsave(
  plot = p,
  file = '03_RF_features.pdf',
  width = 3,
  height = 3)



############################# Fig. 3 F #############################
### heatmap efeC and hyb_opt ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 3 F: efeC-hybopt of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library) %>%
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
  summarize(mean = mean(rTR)) %>%
  ungroup()

# generate plot
p <- ggplot(plot_data, aes(x = bin_efeC, y = bin_hyb, fill = mean)) + 
  geom_tile() +
  theme_SH() +
  scale_fill_viridis(
    name = 'rTR (-)',
    option = 'turbo',
    limits = c(0, 1)) +
  scale_x_continuous(
    name = 'efeC (kcal x mol-1)',
    expand = c(0, 0),
    breaks = c(1:6),
    labels = c('<=-15', '>-15', '>-12.5', '>-10', '>-7.5', '>-5')) +
  scale_y_continuous(
    name = 'hybopt (kcal x mol-1)',
    expand = c(0, 0),
    breaks = c(1:6),
    labels = c('<=-10', '>-10', '>-7.5', '>-5', '>-2.5', '>0')) +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_03_F_efeC_hyb.pdf',
  width = 4,
  height = 3)

# done!
