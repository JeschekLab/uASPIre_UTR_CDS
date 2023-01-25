# set wd and load libraries
ROOT_DIR <- '~/gitlab/uASPIre_UTR_CDS'
setwd(ROOT_DIR)

# load libraries
library(tidyverse)

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

# read tRNA data
data_tRNA <- read.table('./data/data_tRNA_combined.txt',
  sep = '\t', header = T)

# set constant
BASES <- c('A', 'C', 'G', 'U')



############################# Fig. 5 A #############################
### POI ###
data_lib <- data %>% filter(lib == 1) %>% select(seq, rTR, efeC)

data_lib$U1 <- substring(data_lib$seq, 25, 25)
data_lib$G1 <- substring(data_lib$seq, 26, 26)

data_lib$UU <- data_lib$U1 == 'U'
data_lib$GG <- data_lib$G1 == 'G'

p <- ggplot(data_lib, aes(x = GG, y = rTR)) +
  geom_violin(scale = 'width') +
  stat_summary(fun = mean, geom = 'point', size = 3) +
  scale_x_discrete('G at CDS position +6') +
  scale_y_continuous('rTR (-)',
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(ylim = c(0.00, 1.1), clip = 'off')

ggsave('Fig_05_A1_rTR.png', plot = p, width = 3, height = 6,
  units = c('in'), scale = 1)

p <- ggplot(data_lib, aes(x = GG, y = efeC)) +
  geom_violin(scale = 'width') +
  stat_summary(fun = mean, geom = 'point', size = 3) +
  scale_x_discrete('G at CDS position +6') +
  scale_y_continuous('efeC (kcal x mol-1)',
    limits = c(-25, 0), expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(clip = 'off')

ggsave('Fig_05_A2_efeC.png', plot = p, width = 3, height = 3,
  units = c('in'), scale = 1)

# t-test for rTR and efeC
t.test(data_lib$rTR[data_lib$GG], data_lib$rTR[!data_lib$GG])
t.test(data_lib$efeC[data_lib$GG], data_lib$efeC[!data_lib$GG])



############################# Fig. 5 B #############################
### relative triplet frequency CDS+6 ###
data_lib <- data %>%
  filter(lib == 1) %>%
  select(seq, rTR)

data_lib$CDS6 <- substring(data_lib$seq, 26, 26)

data_plot <- data_lib %>%
  group_by(CDS6) %>%
  summarize(mean = mean(rTR)) %>%
  as.data.frame()

RSCU_arg <- data.frame(CDS6 = c('A', 'C', 'G', 'U'),
  RSCU = c(0.09, 0.26, 0.15, 0.30))

data_plot <- merge(data_plot, RSCU_arg, by = 'CDS6')

p <- ggplot(data_plot, aes(x = RSCU, y = mean, color = CDS6)) +
  geom_point(size = 3) +
  theme_SH() +
  scale_y_continuous('mean rTR (-)',
    limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1),
    expand = c(0, 0)) +
  scale_x_continuous('E. coli triplet frequency (-)',
    limits = c(0, 0.4), expand = c(0, 0)) +
  coord_cartesian(clip = 'off')

ggsave('Fig_05_B_triplet.png', plot = p, width = 3, height = 3,
  units = c('in'), scale = 1)



############################# Fig. 5 C #############################
### POI ###
data_lib <- data %>% filter(lib == 1) %>% select(seq, rTR, efeC)

data_lib$U1 <- substring(data_lib$seq, 25, 25)
data_lib$G1 <- substring(data_lib$seq, 26, 26)

data_lib$UU <- data_lib$U1 == 'U'
data_lib$GG <- data_lib$G1 == 'G'

# make boxplots for rTR
p <- ggplot(data_lib, aes(x = UU, y = rTR)) +
  geom_violin(scale = 'width') +
  stat_summary(fun = mean, geom = 'point', size = 3) +
  scale_x_discrete("U at 5'-UTR position -1") +
  scale_y_continuous('rTR (-)',
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)) +
  theme_SH() +
  coord_cartesian(ylim = c(0.00, 1.1))

ggsave('Fig_05_C1_rTR.png', plot = p, width = 3, height = 6,
  units = c('in'), scale = 1)

# make boxplots for efeC
p <- ggplot(data_lib, aes(x = UU, y = efeC)) +
  geom_violin(scale = 'width') +
  stat_summary(fun = mean, geom = 'point', size = 3) +
  scale_x_discrete("U at 5'-UTR position -1") +
  scale_y_continuous('efeC (kcal x mol-1)',
    limits = c(-25, 0), expand = c(0, 0)) +
  theme_SH()

ggsave('Fig_05_C2_efeC.png', plot = p, width = 3, height = 3,
  units = c('in'), scale = 1)

# t-test for rTR and efeC
t.test(data_lib$efeC[data_lib$UU], data_lib$efeC[!data_lib$UU])
t.test(data_lib$rTR[data_lib$UU], data_lib$rTR[!data_lib$UU])



############################# Fig. 5 E #############################
### growth of E. coli ###
df_growth <- read.table('./data/data_OD_tRNA.txt',
  header = T, sep = "\t", fill = T)

data_plot <- df_growth %>%
  select(-Name) %>% 
  reshape2::melt(., id = c('GT', 'Plasmid')) %>%
  group_by(Plasmid, GT) %>%
  summarize(mean = mean(value, na.rm = T), sd = sd(value, na.rm = T)) %>%
  ungroup()

p <- ggplot(data_plot, aes(x = Plasmid, y = mean)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .2) +
  theme_SH() +
  facet_wrap(~GT) +
  scale_x_discrete('tRNAfMet pos. 37') +
  scale_y_continuous("doubling time (min)",
    limits = c(0, 60),
    breaks = seq(0, 60, by = 10),
    expand = c(0, 0))

ggsave('Fig_05_E_growth.png', plot = p,
  width = 3, height = 3, units = c("in"), scale = 1)



############################# Fig. 5 G #############################
### single base effects of tRNA ###
# filter variants that occur in original data set
temp <- data %>% filter(lib == 1) %>% select(-lib) 
data_tRNA <- data_tRNA %>%
  inner_join(., temp, by = 'seq')

# select important rows only
data_tRNA <- data_tRNA %>%
  select(seq, UTR, CDS, rTR, rTR_tRNA, p1, genotype, tRNA)

# initialize
df_tRNA <- data.frame()

# loop
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

p <- ggplot(df_tRNA, aes(x = tRNA, y = FC, fill = UTR)) + 
  geom_bar(position = 'dodge', stat = 'identity') +
  labs(x = 'tRNAfMet pos. 37', y = 'effect on rTR (log2 FC)') +
  theme_SH() +
  scale_y_continuous(expand = c(0, 0),
    limits = c(-0.5, 0.5)) +
  facet_wrap(~genotype, nrow = 1) +
  coord_cartesian(clip = 'off')

ggsave('Fig_05_G_positions.png', plot = p,
  width = 6, height = 3, units = c('in'), scale = 1)



############################# Fig. 5 H #############################
### complementarity ###
# filter variants that occur in original data set
temp <- data %>% filter(lib == 1) %>% select(-lib) 
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
  mutate(logFC = log2(FC)) %>%
  mutate(log_rTR_tRNA = log2(rTR_tRNA)) %>%
  mutate(log_rTR = log2(rTR)) %>%
  mutate(FC_log = log_rTR_tRNA / log_rTR) %>%
  mutate(logFC_log = log2(FC_log))

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

p <- ggplot(df_mean, aes(x = comp, y = mean)) +
  geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2, position = position_dodge(.9)) +
  geom_point(data = df_combined, aes(x = comp, y = logFC, color = tRNA), size = 3) +
  facet_wrap(~genotype) +
  scale_fill_discrete(name = 'tRNAfMet\npos. 37') +
  scale_x_discrete(
    name = '') +
  scale_y_continuous(
    name = 'effect on rTR (log2 FC)',
    limits = c(-0.2, 0.4),
    breaks = seq(-0.2, 0.4, 0.2),
    expand = c(0,0)) +
  theme_SH() +
  coord_cartesian(
    clip = 'off')

ggsave('Fig_05_H_comp.pdf', plot = p,
  width = 4, height = 3, units = c('in'), scale = 1)
