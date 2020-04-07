# ##############################################################################
#
##  Meta-Analysis Confounder Analysis
#
# ##############################################################################

# Packages
library("tidyverse")
library("coin")

alpha.meta <- 0.05

CRC_relative<- as.matrix(read.csv(file = 'species_L7.csv', header=TRUE,
                                  row.names =1,stringsAsFactors = FALSE,check.names = FALSE))

metadata <- read.csv(file = 'metadata_4study.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)


cn.pval.adj <- read.csv(file = 'taxa_result/AN/p.value_an_block_L7.csv',row.names = 1,header=TRUE,
                        stringsAsFactors = FALSE,check.names = FALSE)

metadata_4study <- metadata %>% filter(!is.na(Age))

  
metadata_4study <- metadata_4study %>%
  # age
  mutate(age_factor=as.factor(
    cut(metadata_4study$Age, breaks = quantile(metadata_4study$Age), labels=c(1,2,3,4)))) %>%
  # bmi
  mutate(bmi_factor=as.factor(
    cut(metadata_4study$BMI, breaks = c(0, 25, 30, 100),
        labels=c('lean', 'overweight', 'obese')))) 

metadata_cn <- rbind(metadata_4study %>% filter(Group=='Adenoma'),metadata_4study %>% 
                       filter(Group=='Normal'))
CRC_relative_cn<- CRC_relative[,metadata_cn$Sample_ID]

# 
# write.csv(CRC_relative_cn, file = 'CRC_relative_cn_L7.csv')
# write.csv(CRC_relative_an, file = 'CRC_relative_an_L7.csv')
# write.csv(CRC_relative_ca, file = 'CRC_relative_ca_L7.csv')

#write.csv(metadata_4study,file = '../run/metadata_4study.csv')
#metadata_new 是把metadata_4study 的最后两列加入metadata中的meta信息

#metadata_new <- read.csv("metadata_new.csv",header = TRUE, stringsAsFactors = FALSE)


# ##############################################################################
#  variance explained by disease status
CRC_relative_cn <- read.csv('CRC_relative_an_L7.csv', header = TRUE, row.names = 1,
                            stringsAsFactors = FALSE, check.names = FALSE)

ss.disease <- apply(CRC_relative_cn, 1, FUN=function(x, label){
  rank.x <- rank(x)/length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l){
    sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)
  return(1-ss.o.i/ss.tot)
}, label=metadata_cn %>% pull(Group))

# calculate trimmed mean abundance
t.mean <- apply(CRC_relative_cn, 1, mean, trim=0.1)


CRC_cn_plot.all <- tibble(
  species=rownames(CRC_relative),
  disease=ss.disease,
  t.mean=t.mean,
  adj.p.val=cn.pval.adj[rownames(CRC_relative), 'all'],
  meta.significance=cn.pval.adj[rownames(CRC_relative), 'all'] < 0.05)

# ##############################################################################
# Test all possible confounder variables


df.list <- list()
for (meta.var in c('age_factor','bmi_factor','Gender', 'Study','Race',
                   'Platform','Diabetes','NSAID'
                   )){
  
  cat('###############################\n', meta.var, '\n')
  metadata.c <- metadata_cn %>%
    filter(!is.na(eval(parse(text=meta.var))))
  
  cat('After filtering, the distribution of variables is:\n')
  print(table(metadata.c$Group, metadata.c %>% pull(meta.var)))
  print(table(metadata.c$Study))
  crc.red <- CRC_relative[,metadata.c$Sample_ID]
  
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(crc.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
  }, label=metadata.c %>% pull(meta.var))
  CRC_cn_plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (metadata.c %>% pull(meta.var) %>% unique %>% length > 2){
    meta.significance <- apply(crc.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=metadata.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(crc.red, 1, FUN=function(x, var){
      wilcox.test(x~as.factor(var))$p.value
    }, var=metadata.c %>% pull(meta.var))
  }
  CRC_cn_plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}
  write.csv(CRC_cn_plot.all,file = 'taxa_result/AN/confounder_an_block_L7.csv')



# plot only study
CRC_cn_plot.all <- read.csv('taxa_result/AN/confounder_an_block.csv',row.names = 1,header=TRUE,
                            stringsAsFactors = FALSE,check.names = FALSE)
df.plot.study <- CRC_cn_plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  filter(type=='Study')

g2 <- df.plot.study %>%
  ggplot(aes(x=disease, y=meta)) +
  geom_point(aes(size=t.mean, fill=meta.significance), shape=21,
             col=alpha(c('black'), alpha=0.4)) +
  xlab(paste0('Variance explained by Disease\n','species',' average: ',
              formatC(mean(df.plot.study$disease)*100, digits=2), '%')) +
  ylab(paste0('Variance explained by Study\n','species',' average: ',
              formatC(mean(df.plot.study$meta)*100, digits=2), '%')) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from=0, to=0.05, by=0.01)) +
  scale_y_continuous(breaks=seq(from=0, to=0.4, by=0.02)) +
  scale_fill_manual(values = alpha(c('grey', '#CC071E'),
                                   alpha=c(0.4, .8)),
                    name=paste0('Significance\n(', alpha.meta,')')) +
  scale_size_area(name='Trimmed mean abundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')

ggsave(g2, filename = 'taxa_result/figures/confounder_an_plot_study_L7.pdf',
       width = 6, height = 6)

cat('Successfully computed confounder effects in',
    proc.time()[1]-start.time, 'second...\n')

# ##############################################################################
# plot
CRC_an_plot_new <- read.csv('taxa_result/AN/confounder_an_block_L7.csv',row.names = 1,header=TRUE,
                            stringsAsFactors = FALSE,check.names = FALSE)
g <- CRC_an_plot_new %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  mutate(facet=case_when(type=='Gender' ~ 'Sex',
                         type=='age_factor' ~ 'Age',
                         type=='bmi_factor' ~ 'BMI',
                      
                         TRUE ~ type)) %>%
  ggplot(aes(x=disease, y=meta, size=t.mean+1e-08, col=meta.significance)) +
  geom_point(shape=19) +
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from=0, to=0.2, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
  scale_colour_manual(values = alpha(c('black', '#CC071E'),
                                     alpha=c(0.1, .75)),
                      name=paste0('Significance\n(', alpha.meta, ')')) +
  scale_size_area(name='Trimmed mean\nabundance',
                  breaks=c(1e-05, 1e-03, 1e-02)) +
  guides( size = "legend", colour='legend')

ggsave(g, filename = 'taxa_result/figures/confounder_an_plot_pvalue_L7.pdf',
       width = 6, height = 6)

