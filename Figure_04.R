# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)
library(data.table)
library(car)
library(caret)
library(h2o)

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



############################# Fig. 4 B #############################
### UTR and CDS ###
### enrichment ###

### calculate for lib_comb1 ###
# set library and make print statement
current_library <- 2
print(paste0('Figure 4 B: UTR-CDS of library ', current_library))

# select data
data_lib <- data %>% filter(lib == current_library)

# table UTRs and CDSs and give them an ID
UTRs <- data_lib %>%
  group_by(UTR) %>%
  summarize(Freq_U = n(), mean_UTR = mean(rTR)) %>%
  arrange(desc(Freq_U)) %>%
  mutate(ID_U = 1:nrow(.))

# table UTRs and CDSs and give them an ID
CDSs <- data_lib %>%
  group_by(CDS) %>%
  summarize(Freq_C = n(), mean_CDS = mean(rTR)) %>%
  arrange(desc(Freq_C)) %>%
  mutate(ID_C = 1:nrow(.))

# merge with data
data_lib <- data_lib %>%
  full_join(UTRs, by = 'UTR') %>%
  full_join(CDSs, by = 'CDS')

# add to data_comb
data_comb <- data_lib


### calculate for lib_comb2 ###
# set library and make print statement
current_library <- 3
print(paste0('Figure 4 B: UTR-CDS of library ', current_library))

# select data
data_lib <- data %>% filter(lib == current_library)

# table UTRs and CDSs and give them an ID
UTRs <- data_lib %>%
  group_by(UTR) %>%
  summarize(Freq_U = n(), mean_UTR = mean(rTR)) %>%
  arrange(desc(Freq_U)) %>%
  mutate(ID_U = 1:nrow(.))

# table UTRs and CDSs and give them an ID
CDSs <- data_lib %>%
  group_by(CDS) %>%
  summarize(Freq_C = n(), mean_CDS = mean(rTR)) %>%
  arrange(desc(Freq_C)) %>%
  mutate(ID_C = 1:nrow(.))

# merge with data
data_lib <- data_lib %>%
  full_join(UTRs, by = 'UTR') %>%
  full_join(CDSs, by = 'CDS')

# add to data_comb
data_comb <- rbind(data_comb, data_lib)

# combine
data_comb$comb <- data_comb$mean_CDS + data_comb$mean_UTR

### calculate for lib_fact ###
# set library and make print statement
current_library <- 4
print(paste0('Figure 4 B: UTR-CDS of library ', current_library))

# select data
data_lib <- data %>% filter(lib == current_library)

# read sanger sequencing
pools <- read.table(
  file = './data/data_pools_lib_fact.txt',
  header = T,
  sep = '\t',
  colClasses = c('character', 'integer'))

# extract variable positions
pools$seq <- paste0(
  substring(pools$full, 1, 25),
  substring(pools$full, 31, 31),
  substring(pools$full, 34, 34),
  substring(pools$full, 37, 37),
  substring(pools$full, 40, 40),
  substring(pools$full, 43, 43),
  substring(pools$full, 46, 46),
  substring(pools$full, 49, 49),
  substring(pools$full, 52, 52),
  substring(pools$full, 55, 55),
  substring(pools$full, 58, 58),
  substring(pools$full, 61, 61),
  substring(pools$full, 64, 64),
  substring(pools$full, 67, 67),
  substring(pools$full, 70, 70),
  substring(pools$full, 73, 73))

# check if Sanger reads are in NGS data
pools$correct <- pools$seq %in% data$seq

# table UTRs and CDSs and give them an ID
UTRs <- data_lib %>% group_by(UTR) %>%
  summarize(Freq_U = n()) %>% arrange(desc(Freq_U)) %>%
  mutate(ID_U = 1:nrow(.))

# table UTRs and CDSs and give them an ID
CDSs <- data_lib %>% group_by(CDS) %>%
  summarize(Freq_C = n()) %>% arrange(desc(Freq_C)) %>%
  mutate(ID_C = 1:nrow(.))

# merge with data
data_lib <- data_lib %>%
  full_join(UTRs, by = 'UTR') %>%
  full_join(CDSs, by = 'CDS')

# calculate sum of RBS and CDS
data_lib$appearance <- data_lib$Freq_U + data_lib$Freq_C
data_lib <- data_lib %>% arrange(desc(appearance))

# split into pools
pool_list <- list()

for (current_pool in 1:10) {
  print(paste0('Current pool: ', current_pool))
  pool_x <- data_lib[data_lib$seq %in% pools$seq[pools$pool == current_pool], ]

  # select CDSs and UTRs
  sanger_UTRs <- unique(pool_x$UTR)
  sanger_CDSs <- unique(pool_x$CDS)

  # select all combinations with those UTRs and CDSs
  sanger_variants <- data_lib[(data_lib$UTR %in% sanger_UTRs | data_lib$CDS %in% sanger_CDSs) & data_lib$Freq_U <= 125 & data_lib$Freq_C <= 150, ]
  pool_UTRs <- unique(sanger_variants$UTR)
  pool_CDSs <- unique(sanger_variants$CDS)

  # select all combinations with all UTRs and CDSs from the pool
  pool_variants <- data_lib[(data_lib$UTR %in% pool_UTRs | data_lib$CDS %in% pool_CDSs) & data_lib$appearance <= 250, ] # 23,003

  # make table
  pool_table_UTR <- as.data.frame(table(pool_variants$UTR), stringsAsFactors = F) %>%
    arrange(desc(Freq))
  pool_table_CDS <- as.data.frame(table(pool_variants$CDS), stringsAsFactors = F) %>%
    arrange(desc(Freq))

  # give ID
  pool_table_UTR$ID_R_pool <- 1:nrow(pool_table_UTR)
  pool_table_CDS$ID_C_pool <- 1:nrow(pool_table_CDS)
  names(pool_table_UTR)[1:2] <- c('UTR', 'Freq_pool_R')
  names(pool_table_CDS)[1:2] <- c('CDS', 'Freq_pool_C')

  # merge with data
  pool_temp <- merge(data_lib, pool_table_UTR, by = 'UTR', all = F)
  pool <- merge(pool_temp, pool_table_CDS, by = 'CDS', all = F)

  # calculate pool appearance
  pool$appearance_pool <- pool$Freq_pool_R + pool$Freq_pool_C

  # remove all that appear less than 100 times (50 each)
  cutoff <- 50
  pool_cutoff <- pool %>%
    filter(appearance_pool >= (cutoff*2)) %>%
    filter(Freq_pool_R >= cutoff) %>%
    filter(Freq_pool_C >= cutoff)

  pool_final <- pool_cutoff
  pool_final$pool <- current_pool

  # calculate average UTR and CDS expression
  data_aggregated_UTR <- pool_final %>%
    group_by(UTR) %>%
    summarize(mean_UTR = mean(rTR))
  data_aggregated_CDS <- pool_final %>%
    group_by(CDS) %>%
    summarize(mean_CDS = mean(rTR))

  data_merge_UTR <- merge(pool_final, data_aggregated_UTR, by = 'UTR')
  data_merge_both <- merge(data_merge_UTR, data_aggregated_CDS, by = 'CDS')

  pool_list[[current_pool]] <- data_merge_both
}

data_fact <- as.data.frame(rbindlist(pool_list))
data_fact$lib <- 4

### relative rTR change (UTR and CDS) ###
temp1 <- data_comb %>% select(seq, lib, rTR, mean_UTR, mean_CDS)
temp2 <- data_fact %>% select(seq, lib, rTR, mean_UTR, mean_CDS)

data_lib <- rbind(temp1, temp2)

# make log2 FC
data_lib$FC_CDS <- log2(data_lib$rTR / data_lib$mean_UTR)
data_lib$FC_UTR <- log2(data_lib$rTR / data_lib$mean_CDS)

# make long format
data_plot <- data_lib %>%
  select(seq, lib, FC_UTR, FC_CDS) %>%
    pivot_longer(
      cols = c(-seq, -lib),
      names_to = 'variable',
      values_to = 'value') %>%
    mutate(lib = as.character(lib))

# generate plot
p <- ggplot(data_plot,
  aes(x = lib, y = abs(value), fill = variable)) +
  geom_violin(
    scale = 'width') +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 2) +
  theme_SH() +
  scale_x_discrete(
    name = '') +
  scale_y_continuous(
    name = 'rel. rTR change (abs. log2 FC)',
    limits = c(0, 4),
    expand = c(0, 0)) +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90')) +
  coord_cartesian(
    clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_04_B_violins.pdf',
  width = 3,
  height = 3)



############################# Fig. 4 C #############################
### ANOVA of effects ###

# make print statement
print(paste0('Figure 4 C: ANOVA of library ', current_library))

# select data
aov_pools <- list()
for (current_pool in 1:10) {
  pool_x <- data_fact %>% filter(pool == current_pool)

  # make linear model
  mod.lm <- lm(data = pool_x, formula = rTR~mean_UTR+mean_CDS)
  res.aov <- Anova(mod.lm, type = 'II')

  res <- as.data.frame(res.aov[1])
  res$effect <- c('UTR', 'CDS', 'unexplained')
  res$pool <- current_pool
  res$rel <- res$'Sum Sq' / sum(res$'Sum Sq')
  aov_pools[[current_pool]] <- res
}

aov_pools_combined <- as.data.frame(rbindlist(aov_pools))

# make bar plots
aov_combined <- aov_pools_combined %>%
  group_by(effect) %>%
  summarize(mean = mean(rel), sd = sd(rel))

# generate plot
p <- ggplot(aov_combined, aes(x = effect, y = mean, fill = effect)) + 
  geom_bar(
    stat = 'identity',
    colour = 'black') +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .2) +
  theme_SH() +
  scale_x_discrete(
    name = '') +
  scale_y_continuous(
    name = 'contribution to rTR variance (%)',
    limits = c(0, 0.65),
    expand = c(0, 0)) +
  theme(
    legend.position = 'none') +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90', 'white')) +
  coord_cartesian(
    clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_04_C_anova.pdf',
  width = 3,
  height = 3)



############################# Fig. 4 D #############################
### effect size CAI and tAI ###

# set library and make print statement
print(paste0('Figure 4 D: CAI/tAI'))

# select data
df_corr_CAI <- data.frame()

# loop through library and calculate effect sizes of CAI and tAI
for (current_lib in 1:4) {
  data_lib <- data %>%
    filter(lib == current_lib) %>%
    select(rTR, tAI, CAI)

  rTR <- data_lib$rTR

  for (current_model in c('tAI', 'CAI')) {
    current_energy <- data_lib[, current_model]
    p <- cor(rTR, current_energy, method = 'spearman')
    R <- cor(rTR, current_energy, method = 'pearson')
    df_temp <- data.frame(
      model = current_model,
      cor = c(p, R),
      stat = c('p2', 'R2'),
      lib = current_lib)
    df_corr_CAI <- rbind(df_corr_CAI, df_temp)
    print(current_model)
  }
}

df_corr_CAI$explainability <- (df_corr_CAI$cor^2)*100
class(df_corr_CAI$lib) <- 'integer'

# generate plot
p <- ggplot(data = df_corr_CAI,
    aes(y = explainability, x = lib, fill = stat)) +
  geom_bar(
    position = 'dodge',
    stat = 'identity',
    color = 'black') +
  scale_x_reverse(
    name = '') +
  scale_y_continuous(
    name = 'p2 / R2 (%)',
    limits = c(0, 1.2),
    breaks = seq(0, 1.2, 0.4),
    expand = c(0, 0),
    ) +
  theme_SH() +
  coord_flip(clip = 'off') +
  facet_wrap(~model) +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90'))

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_04_D_CAI_tAI.pdf',
  width = 3,
  height = 3)



############################# Fig. 4 E #############################
### bin CAI and tAI ###

# set library and make print statement
print(paste0('Figure 4 E: low CAI/tAI'))

# select data
data_lib <- data %>%
  select(rTR, CAI, tAI, efeC)

cutoff <- 0.10

# bin CAI
CAI_10 <- data_lib %>%
  filter(CAI <= cutoff) %>%
  mutate(label = 'CAI_10')
CAI_90 <- data_lib %>%
  filter(CAI > cutoff) %>%
  mutate(label = 'CAI_90')

# bin tAI
tAI_10 <- data_lib %>%
  filter(tAI <= cutoff) %>%
  mutate(label = 'tAI_10') 
tAI_90 <- data_lib %>%
  filter(tAI > cutoff) %>%
  mutate(label = 'tAI_90')

# combine
data_plot <- rbind(CAI_10, CAI_90, tAI_10, tAI_90)
data_plot$efeC_norm <- data_plot$efeC / -25

# make long format
data_plot <- data_plot %>%
  select(rTR, label, efeC_norm) %>%
    pivot_longer(
      cols = -label,
      names_to = 'group',
      values_to = 'value')

# generate plot
p <- ggplot(data_plot, aes(x = label, y = value, fill = group, color = group)) +
  geom_violin(
    scale = 'width',
    color = 'black') +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 2) +
  scale_x_discrete(
    name = '') +
  scale_y_continuous(
    name = 'rTR (-)',
    limits = c(0, 1.1),
    expand = c(0,0),
    breaks = seq(0, 1, 0.2),
    sec.axis = sec_axis(~.*-25, name = 'efeC')) +
  theme_SH() +
  coord_cartesian(clip = 'off') +
  scale_color_manual(
    name = '',
    values = c('grey30', 'grey90')) +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90'))

# save plot to file
ggsave(
  plot = p,
  'Fig_04_E_CAI_efeC.pdf',
  width = 4,
  height = 3)

# make test statistic
t.test(CAI_10$rTR, CAI_90$rTR) # 0.213 vs. 0.217, 1.00e-12
t.test(tAI_10$rTR, tAI_90$rTR) # 0.296 vs. 0.216, 6.67e-226

t.test(tAI_10$efeC, tAI_90$efeC, alternative = 'greater') # -6.17 vs. -7.23, 0
t.test(CAI_10$efeC, CAI_90$efeC, alternative = 'greater') # -7.29 vs. -7.21, 1.38e-46



############################# Fig. 4 F #############################
### ANOVA and linear model ###

# set library and make print statement
print(paste0('Figure 4 F: ANOVA of CAI/tAI'))

# select data
data_lib <- data

# make linear model
model <- lm(data = data_lib, formula = rTR ~ efeC + CAI + tAI)

model_aov <- as.data.frame(Anova(model, type = 'II')[1])
model_aov$name <- rownames(model_aov)
model_aov$SOS <- model_aov$'Sum Sq'
model_aov$SOS_rel <- round(model_aov$SOS / sum(model_aov$SOS) * 100, 5)

temp_plot <- model_aov[-4, ]
temp_plot$rel <- (temp_plot$SOS_rel / sum(temp_plot$SOS_rel)) * 100

# generate plot
p <- ggplot(data = temp_plot, aes(x = name, y = SOS_rel)) +
  geom_bar(
    position = 'dodge',
    stat = 'identity',
    fill = 'grey90',
    color = 'black') +
  scale_x_discrete(
    name = '') +
  scale_y_continuous(
    name = 'contribution to rTR variance (%)',
    expand = c(0, 0),
    limits = c(0, 15)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_04_F_efeC_CAI_tAI.pdf',
  width = 3,
  height = 3)



############################# Fig. 4 G #############################
### random forest ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 4 G: Feature importance of library ', current_library))

# read accC data
df_wider <- read.table(paste0('./data/data_accC_1.txt'),
  sep = '\t', header = T)

# combine with data and select only library 1
data_lib <- cbind(data, df_wider) %>%
  filter(lib == current_library) %>%
  select(seq, rTR, paste0('X', 1:80), hyb_opt,
    GC_all, mfeT, mfeC, efeT, efeC, accT, accC, CAI, tAI)

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

# set seeds
seeds <- c(12345, 275635, 535361, 843675, 617942)

y.dep <- 'rTR'

df_cor_RF <- data.frame()

# initialize h20 cluster
h2o.init()

# loop through seeds and train random forest
for (current_replicate in 1:5) {
  current_seed <- seeds[current_replicate]

  set.seed(current_seed)

  index <- createDataPartition(df_RF$rTR, p = 0.90, list = FALSE)

  train <- df_RF[index, ]
  test <- df_RF[-index, ]

  # copy to cluster
  train_RF.h2o <- train %>%
    select(-CAI, -tAI) %>%
    as.h2o(.)
  train_RF_CAI.h2o <- train %>%
    as.h2o(.)
  train_RF_noFold.h2o <- train %>%
    select(rTR, paste0('V', 1:160), hyb_opt, GC_all) %>%
    as.h2o(.)
  train_RF_noFold_CAI.h2o <- train %>%
    select(rTR, paste0('V', 1:160), hyb_opt, GC_all, CAI, tAI) %>%
    as.h2o(.)

  # train
  RF <- h2o.randomForest(
    y = y.dep,
    training_frame = train_RF.h2o,
    nfolds = 10,
    seed = current_seed)
  RF_CAI <- h2o.randomForest(
    y = y.dep,
    training_frame = train_RF_CAI.h2o,
    nfolds = 10,
    seed = current_seed)
  RF_noFold <- h2o.randomForest(
    y = y.dep,
    training_frame = train_RF_noFold.h2o,
    nfolds = 10,
    seed = current_seed)
  RF_noFold_CAI <- h2o.randomForest(
    y = y.dep,
    training_frame = train_RF_noFold_CAI.h2o,
    nfolds = 10,
    seed = current_seed)

  # predict
  test.h2o <- test %>% as.h2o(.)

  test$RF <- as.data.frame(h2o.predict(RF, test.h2o))$predict
  test$RF_CAI <- as.data.frame(h2o.predict(RF_CAI, test.h2o))$predict
  test$RF_noFold <- as.data.frame(h2o.predict(RF_noFold, test.h2o))$predict
  test$RF_noFold_CAI <- as.data.frame(h2o.predict(RF_noFold_CAI, test.h2o))$predict

  # calculate correlations
  rTR <- test$rTR

  for (current_name in c('RF', 'RF_CAI', 'RF_noFold', 'RF_noFold_CAI')) {
    R <- cor(rTR, test[, current_name], method = 'pearson')
    df_temp <- data.frame(name = current_name, rep = current_replicate, R = R)
    df_cor_RF <- rbind(df_cor_RF, df_temp)
  }
}

# make R2
df_cor_RF$R2 <- df_cor_RF$R^2

# make t-tests
t.test(df_cor_RF %>% filter(name == 'RF') %>% .$R2,
  df_cor_RF %>% filter(name == 'RF_CAI') %>% .$R2)
t.test(df_cor_RF %>% filter(name == 'RF_noFold') %>% .$R2,
  df_cor_RF %>% filter(name == 'RF_noFold_CAI') %>% .$R2)

# plot
df_cor_RF_combined <- df_cor_RF %>% group_by(name) %>%
  summarize(mean = mean(R2), sd = sd(R2))

# generate plot
p <- ggplot(df_cor_RF_combined, aes(x = name, y = mean, fill = name)) + 
  geom_bar(
    stat = 'identity',
    colour = 'grey30') +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .2) +
  theme_SH() +
  scale_x_discrete(
    name = 'training features RF') +
  scale_y_continuous(
    name = 'R2 on test set (%)',
    limits = c(0, 0.7),
    breaks = seq(0, 0.6, 0.2),
    expand = c(0, 0)) +
  theme(
    legend.position = 'none') +
  coord_cartesian(
    clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_04_G_RF.pdf',
  width = 3,
  height = 3)

# done!
