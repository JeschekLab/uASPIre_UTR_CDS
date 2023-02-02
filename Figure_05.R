# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)

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

# set constant
BASES <- c('A', 'C', 'G', 'U')



############################# Fig. 5 A #############################
### POI ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 5 A.1: rTR of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library) %>%
  select(seq, rTR, efeC)

data_lib$U1 <- substring(data_lib$seq, 25, 25)
data_lib$G1 <- substring(data_lib$seq, 26, 26)

data_lib$UU <- data_lib$U1 == 'U'
data_lib$GG <- data_lib$G1 == 'G'

# generate plot
p <- ggplot(data_lib, aes(x = GG, y = rTR, fill = GG)) +
  geom_violin(
    scale = 'width') +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 2) +
  scale_x_discrete(
    name = 'G at CDS position +6') +
  scale_y_continuous(
    name = 'rTR (-)',
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(
    clip = 'off') +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90'))

ggsave(
  plot = p,
  file = 'Fig_05_A1_rTR.pdf',
  width = 3,
  height = 3)

# generate plot
p <- ggplot(data_lib, aes(x = GG, y = efeC, fill = GG)) +
  geom_violin(
    scale = 'width') +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 2) +
  scale_x_discrete(
    name = 'G at CDS position +6') +
  scale_y_continuous(
    name = 'efeC (kcal x mol-1)',
    limits = c(-25, 0),
    expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(
    clip = 'off') +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90'))

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_05_A2_efeC.pdf',
  width = 3,
  height = 3)

# t-test for rTR and efeC
t.test(data_lib$rTR[data_lib$GG], data_lib$rTR[!data_lib$GG])
t.test(data_lib$efeC[data_lib$GG], data_lib$efeC[!data_lib$GG])



############################# Fig. 5 B #############################
### relative triplet frequency CDS+6 ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 5 B: triplet frequency of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library) %>%
  select(seq, rTR)

data_lib$CDS6 <- substring(data_lib$seq, 26, 26)

data_plot <- data_lib %>%
  group_by(CDS6) %>%
  summarize(mean = mean(rTR)) %>%
  as.data.frame()

# manually add codon usage
RSCU_arg <- data.frame(
  CDS6 = c('A', 'C', 'G', 'U'),
  RSCU = c(0.09, 0.26, 0.15, 0.30))

data_plot <- merge(data_plot, RSCU_arg, by = 'CDS6')

# generate plot
p <- ggplot(data_plot, aes(x = RSCU, y = mean, color = CDS6)) +
  geom_point(
    size = 3) +
  theme_SH() +
  scale_y_continuous(
    name = 'mean rTR (-)',
    limits = c(0, 0.3),
    breaks = seq(0, 0.3, 0.1),
    expand = c(0, 0)) +
  scale_x_continuous(
    name = 'E. coli triplet frequency (-)',
    limits = c(0, 0.4),
    expand = c(0, 0)) +
  coord_cartesian(clip = 'off') +
  scale_color_manual(
    name = '',
    values = c('#55934E', '#3D5D96', '#E5B139', '#B72C3B'))

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_05_B_triplet.pdf',
  width = 3,
  height = 3)



############################# Fig. 5 C #############################
### POI ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 5 C.1: rTR of library ', current_library))

# select data
data_lib <- data %>%
  filter(lib == current_library) %>%
  select(seq, rTR, efeC)

data_lib$U1 <- substring(data_lib$seq, 25, 25)
data_lib$G1 <- substring(data_lib$seq, 26, 26)

data_lib$UU <- data_lib$U1 == 'U'
data_lib$GG <- data_lib$G1 == 'G'

# generate plot
p <- ggplot(data_lib, aes(x = UU, y = rTR, fill = UU)) +
  geom_violin(
    scale = 'width') +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 2) +
  scale_x_discrete(
    name = "U at 5'-UTR position -1") +
  scale_y_continuous(
    name = 'rTR (-)',
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(
    clip = 'off') +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90'))

ggsave(
  plot = p,
  file = 'Fig_05_C1_rTR.pdf',
  width = 3,
  height = 3)

# generate plot
p <- ggplot(data_lib, aes(x = UU, y = efeC, fill = UU)) +
  geom_violin(
    scale = 'width') +
  stat_summary(
    fun = mean,
    geom = 'point',
    size = 2) +
  scale_x_discrete(
    name = "U at 5'-UTR position -1") +
  scale_y_continuous(
    name = 'efeC (kcal x mol-1)',
    limits = c(-25, 0),
    expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(
    clip = 'off') +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90'))

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_05_C2_efeC.pdf',
  width = 3,
  height = 3)

# t-test for rTR and efeC
t.test(data_lib$rTR[data_lib$UU], data_lib$rTR[!data_lib$UU])
t.test(data_lib$efeC[data_lib$UU], data_lib$efeC[!data_lib$UU])



############################# Fig. 5 E #############################
### growth of E. coli ###

# set library and make print statement
print(paste0('Figure 5 E: tRNA growth'))

# read growth data
df_growth <- read.table(
  file = './data/data_OD_tRNA.txt',
  header = T,
  sep = '\t',
  fill = T)

# calculate mean and sd
data_plot <- df_growth %>%
  group_by(Name) %>%
  mutate(
    mean = mean(c(DT1, DT2, DT3, DT4), na.rm = T),
    sd = sd(c(DT1, DT2, DT3, DT4), na.rm = T))

# generate plot
p <- ggplot(data_plot, aes(x = Plasmid, y = mean, fill = Plasmid)) + 
  geom_bar(
    position = 'dodge',
    stat = 'identity',
    color = 'black') +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .2) +
  theme_SH() +
  facet_wrap(~GT) +
  scale_x_discrete('tRNAfMet pos. 37') +
  scale_y_continuous("doubling time (min)",
    limits = c(0, 60),
    breaks = seq(0, 60, by = 10),
    expand = c(0, 0)) +
  scale_fill_manual(
    name = '',
    values = c('#55934E', '#E5B139', '#B72C3B', '#555555')) +
  coord_cartesian(clip = 'off') +
  theme(
    legend.position = 'none')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_05_E_growth.pdf',
  width = 3,
  height = 3)



############################# Fig. 5 G #############################
### single base effects of tRNA ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 5 G: Singel base effects of library ', current_library))

# read tRNA data
data_tRNA <- read.table(
  file = './data/data_tRNA_combined.txt',
  sep = '\t',
  header = T)

# filter variants that occur in original data set
temp <- data %>%
  filter(lib == current_library) %>%
  select(-lib) 

# combine with tRNA data
data_tRNA <- data_tRNA %>%
  inner_join(., temp, by = 'seq')

# select important rows only
data_tRNA <- data_tRNA %>%
  select(seq, UTR, CDS, rTR, rTR_tRNA, p1, genotype, tRNA)

# initialize
df_tRNA <- data.frame()

# loop through genotypes and bases and calculate -1 base effect
for (current_genotype in c('WT', 'KO')) {
  for (current_tRNA in c('A', 'G', 'U')) {
    for (current_UTR in BASES) {
      temp <- data_tRNA %>%
        filter(genotype == current_genotype) %>%
        filter(tRNA == current_tRNA)
      has <- temp %>%
        filter(p1 == current_UTR) %>%
        pull(rTR_tRNA) %>%
        mean()
      not <- temp %>%
        filter(p1 != current_UTR) %>%
        pull(rTR_tRNA) %>%
        mean()
      df_temp <- data.frame(
        genotype = current_genotype,
        tRNA = current_tRNA,
        UTR = current_UTR,
        has = has,
        not = not)
      df_tRNA <- rbind(df_tRNA, df_temp)
    }
  }
}

# calculate difference
df_tRNA <- df_tRNA %>%
  mutate(FC = log2(has/not))

# generate plot
p <- ggplot(df_tRNA, aes(x = tRNA, y = FC, fill = UTR)) + 
  geom_bar(
    position = 'dodge',
    stat = 'identity',
    color = 'black') +
  scale_x_discrete(
    name = 'tRNAfMet pos. 37') +
  scale_y_continuous(
    name = 'effect on rTR (log2 FC)',
    expand = c(0, 0),
    limits = c(-0.5, 0.5)) +
  theme_SH() +
  facet_wrap(~genotype, nrow = 1) +
  scale_fill_manual(
    name = '',
    values = c('#55934E', '#3D5D96', '#E5B139', '#B72C3B')) +
  coord_cartesian(clip = 'off')

# save plot to file
ggsave('Fig_05_G_positions.pdf', plot = p,
  width = 6, height = 3, units = c('in'), scale = 1)



############################# Fig. 5 H #############################
### complementarity ###

# set library and make print statement
current_library <- 1
print(paste0('Figure 5 H: Complementary of library ', current_library))

# read tRNA data
data_tRNA <- read.table(
  file = './data/data_tRNA_combined.txt',
  sep = '\t',
  header = T)

# filter variants that occur in original data set
temp <- data %>%
  filter(lib == current_library) %>%
  select(-lib)

# add tRNA data
data_tRNA <- data_tRNA %>%
  inner_join(., temp, by = 'seq')

# select important rows only
data_tRNA <- data_tRNA %>%
  select(seq, UTR, CDS, rTR, rTR_tRNA, p1, genotype, tRNA)

# find complementary
data_tRNA$comp <- '0'
data_tRNA$comp[data_tRNA$p1 == 'A' & data_tRNA$tRNA == 'U'] <- 1
data_tRNA$comp[data_tRNA$p1 == 'C' & data_tRNA$tRNA == 'G'] <- 1
data_tRNA$comp[data_tRNA$p1 == 'G' & data_tRNA$tRNA == 'C'] <- 1
data_tRNA$comp[data_tRNA$p1 == 'U' & data_tRNA$tRNA == 'A'] <- 1

# calculate fold changes
data_tRNA <- data_tRNA %>%
  mutate(FC = rTR_tRNA / rTR) %>%
  mutate(logFC = log2(FC)) # %>%
  # mutate(log_rTR_tRNA = log2(rTR_tRNA)) %>%
  # mutate(log_rTR = log2(rTR)) %>%
  # mutate(FC_log = log_rTR_tRNA / log_rTR) %>%
  # mutate(logFC_log = log2(FC_log))

# calculate fold change
df_FC <- data_tRNA %>%
  group_by(genotype, comp, tRNA) %>%
  summarise(mean_comp = mean(rTR_tRNA)) %>%
  ungroup()

df_rTR <- data_tRNA %>%
  group_by(genotype, tRNA) %>%
  summarise(mean_all = mean(rTR_tRNA)) %>%
  ungroup()

df_combined <- inner_join(df_FC, df_rTR,
  by = c('genotype', 'tRNA')) %>%
  mutate(FC = mean_comp / mean_all) %>%
  mutate(logFC = log2(FC))

df_mean <- df_combined %>%
  group_by(comp, genotype) %>%
  summarize(mean = mean(logFC), sd = sd(logFC))

# generate plot
p <- ggplot(df_mean, aes(x = comp, y = mean, fill = comp)) +
  geom_bar(
    position = 'dodge',
    stat = 'identity',
    color = 'black') +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),
    width = 0.2,
    position = position_dodge(.9)) +
  geom_point(
    data = df_combined, aes(x = comp, y = logFC, color = tRNA), size = 3) +
  scale_x_discrete(
    name = 'complementarity') +
  scale_y_continuous(
    name = 'effect on rTR (log2 FC)',
    limits = c(-0.2, 0.4),
    breaks = seq(-0.2, 0.4, 0.2),
    expand = c(0,0)) +
  theme_SH() +
  scale_fill_manual(
    name = '',
    values = c('grey30', 'grey90')) +
  scale_color_manual(
    name = 'tRNAfMet\npos. 37',
    values = c('#55934E', '#E5B139', '#B72C3B')) +
  facet_wrap(~genotype) +
  coord_cartesian(
    clip = 'off')

# save plot to file
ggsave(
  plot = p,
  file = 'Fig_05_H_comp.pdf',
  width = 4,
  height = 3)

# done!
