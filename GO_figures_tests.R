library(tidyverse)
library(topGO)
library(ggrepel)
library(funfuns)
library(vegan)
library(GO.db)


########## GO and reference STUFF ###########

GO_terms = select(GO.db, keys(GO.db, "GOID"), c("TERM", "ONTOLOGY"))


# Here I'll get a df of my detected protein accessions and their corresponding GO terms

# used interpro to assign GO terms to all proteins #

interpro1 <- read_delim('NCBI_O157.tsv.txt', delim = '\t')
interpro2 <- read_delim('Uniprot_ecoli_pan.tsv', delim = '\t', col_names = colnames(interpro1))

interpro <- rbind(interpro1, interpro2)

GO <- interpro[!is.na(interpro$GO),] # this uses the middle of the Uniprot seq.id (might not work well for itraq?)

GO_all <- GO %>%
  dplyr::select(accno, GO) %>% separate_rows(GO, sep = '\\|') %>%
  unique() %>% 
  group_by(accno) %>%
  summarise(GO=paste(GO, sep = '', collapse = ',')) %>%
  unique()


# protein annotations
all_descriptions <- read_delim('./reference_files/prot_descriptions2.csv', delim = '\t', col_names = c('accno', 'description'))

# psort locations
psort <- read_delim('reference_files/psort_classifications_final.txt', delim = '\t', col_names = c('accno', 'psort_loc')) %>% 
  mutate(psort_loc=stringr::str_remove_all(psort_loc,' ')) %>% 
  mutate(psort_loc=stringr::str_remove_all(psort_loc,'\\(.*\\)')) %>% 
  mutate(psort_loc=stringr::str_remove_all(psort_loc,'[0-9]?[0-9]\\.[0-9][0-9]')) 

# to fix inconsistencies in accnos

ID_MAPPER <- read_delim('IDMAPPER.tsv', delim = '\t', col_names = c('ID1', 'ID2'))
ID_dict <- ID_MAPPER$ID1
names(ID_dict) <- ID_MAPPER$ID2

## fixing psort accnos

psort$accno2 <- ID_dict[psort$accno]

psort <- psort %>%
  mutate(accno=ifelse(is.na(accno2), accno, accno2)) %>% 
  dplyr::select(-accno2) %>%
  unique()

################################

### iTRAQ stuff ###
# Scaffold normalized iTRAQ intensities as well as Scaffold statistical test

lact <- read_csv('iTRAQ_lact.csv')
maint <- read_csv('iTRAQ_maint.csv')


maint <- maint[-grep('\\.', as.character(maint$num)),]
# lact <- lact[-grep('\\.', as.character(lact$num)),]

# lact$accno %in% MQ_iBAQ$accno



maint$accno <- sub('(.*) (.*)','\\1',maint$accno)
lact$accno <- sub('(.*) (.*)','\\1',lact$accno)

### FIX ACCNOS HERE ###


# ID_dict['A0A2H9ENY2_ECOLX']

maint$accno2 <- ID_dict[maint$accno]
lact$accno2 <- ID_dict[lact$accno]


maint <- maint %>% 
  mutate(accno=ifelse(is.na(accno2), accno, accno2)) %>%
  dplyr::select(-accno2)

lact <- lact %>%
  mutate(accno=ifelse(is.na(accno2), accno, accno2)) %>%
  dplyr::select(-accno2)


####



colnames(lact)[5:10] <- c('86-24_vitro', 'EDL933_vitro', 'SS17_vitro', '86-24_vivo', 'EDL933_vivo', 'SS17_vivo')
colnames(maint)[5:10] <- c('86-24_vitro', 'EDL933_vitro', 'SS17_vitro', '86-24_vivo', 'EDL933_vivo', 'SS17_vivo')

lact_gather <- lact %>% 
  dplyr::select(ends_with('o')) %>% 
  gather(key='sample', value='value', -accno) %>% 
  mutate(value=as.numeric(value)) %>%
  separate(sample, into = c('strain', 'condition'), sep = '_')


lact_gather$value[is.na(lact_gather$value)]


lact_gather %>% filter(is.na(value))


lact_gather %>%
  filter(accno == c('AAG59207.1') & condition == 'vivo') %>%
  dplyr::select(value) %>%
  summarise(mean=mean(value,na.rm=TRUE))

lact_gather[lact_gather$accno == 'AAG59207.1' &
              lact_gather$condition == 'vivo' &
              lact_gather$strain == 'EDL933',]$value <- 17.6


lact_gather %>%
  filter(accno == c('A0A2X2HSZ0') & condition == 'vivo') %>%
  dplyr::select(value) %>% summarise(mean=mean(value,na.rm=TRUE))

lact_gather[lact_gather$accno == 'A0A2X2HSZ0' &
              lact_gather$condition == 'vivo' &
              lact_gather$strain == 'EDL933',]$value <- 14.8


lact_mat <- lact_gather %>% 
  mutate(sample=paste(strain, condition, sep = '_')) %>%
  dplyr::select(-strain, -condition) %>% 
  spread(key=sample, value = value) %>%
  column_to_rownames(var='accno')

lact_matt <- t(lact_mat)

maint_mat <- maint %>%
  dplyr::select(ends_with('o')) %>%
  column_to_rownames(var='accno')

maint_matt <- t(maint_mat)



lact_meta <- data_frame(sample=row.names(lact_matt), 
                        strain = sub('(.*)_(.*)','\\1',rownames(lact_matt)), 
                        condition=sub('(.*)_(.*)','\\2',rownames(lact_matt)))

maint_meta <- data_frame(sample=row.names(maint_matt), 
                         strain = sub('(.*)_(.*)','\\1',rownames(maint_matt)), 
                         condition=sub('(.*)_(.*)','\\2',rownames(maint_matt)))




###

# rownames(lact_matt)
rownames(lact_meta) <- lact_meta$sample

# rownames(maint_matt)
rownames(maint_meta) <- maint_meta$sample



lact_NMDS <- NMDS_ellipse(metadata = lact_meta, OTU_table = lact_matt, grouping_set = 'condition')
lact_adon <- adonis(data = lact_meta, formula = lact_matt ~ strain + condition)

maint_NMDS <- NMDS_ellipse(metadata = maint_meta, OTU_table = maint_matt, grouping_set = 'condition')
maint_adon <- adonis(data = maint_meta, formula = maint_matt ~ strain + condition)


lact_adon$aov.tab %>% broom::tidy() %>% write_csv('./results/iTRAQ_lact_PERMANOVA.csv')
maint_adon$aov.tab %>% broom::tidy() %>% write_csv('./results/iTRAQ_maint_PERMANOVA.csv')

###

# vegdist(lact_matt)

lact_meta <- cbind(lact_meta, cmdscale(vegdist(lact_matt)))
maint_meta <- cbind(maint_meta, cmdscale(vegdist(maint_matt)))

colnames(lact_meta)[4:5] <- c('MDS1', 'MDS2')
colnames(maint_meta)[4:5] <- c('MDS1', 'MDS2')

## fig 1 ##
lact_meta %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=condition)) +
  geom_point(shape=21, size=3) + theme_bw() + geom_text_repel(aes(label=strain)) + 
  ggtitle('Lactation Diet: similarity of global protein expression')

## fig 2 ##
maint_meta %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=condition)) +
  geom_point(shape=21, size=3) + theme_bw() + geom_text_repel(aes(label=strain)) + 
  ggtitle('Maintenance Diet: similarity of global protein expression')


lact_disper <- vegan::betadisper(vegdist(lact_matt), lact_meta$condition)
maint_disper <- vegan::betadisper(vegdist(maint_matt), maint_meta$condition)

plot(lact_disper)
plot(maint_disper)


permutest(lact_disper)
permutest(maint_disper)

anova(maint_disper)
anova(lact_disper)


TukeyHSD(maint_disper)
TukeyHSD(lact_disper)



lact_meta$disper_dist <- lact_disper$distances
maint_meta$disper_dist <- maint_disper$distances



## fig 3 ##
lact_meta %>% ggplot(aes(x=condition, y=disper_dist)) +
  stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, colour = "red")+
  geom_point(aes(fill=condition), shape=21, size=3) +geom_text_repel(aes(label=strain)) + 
  theme_bw() + ylab('distance to group centroid') + 
  ggtitle('Lactation diet: Global protein expression dispersion', 'how much does each strain differ from the group centroid (typical expression for the group)')


## fig 4 ##
maint_meta %>%
  ggplot(aes(x=condition, y=disper_dist)) +
  stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, colour = "red")+
  geom_point(aes(fill=condition), shape=21, size=3) +geom_text_repel(aes(label=strain)) + 
  theme_bw() + ylab('distance to group centroid') + 
  ggtitle('Maintenance diet: Global protein expression dispersion', 'how much does each strain differ from the group centroid (typical expression for the group)')


### Normalized iTRAQ reporter intensities were used to construct a Bray-Curtis dissimilarity matrix describing the similarities of the global
### protein expression between all samples within each diet (LC-MSMS run).  WIthin the lactation diet, a significant difference in 
### global protein expression patterns was detected between the two culture conditions, (PERMANOVA p=0.05, F=5.6, R2=0.572) but no significant
### difference was detected between the different strains (PERMANOVA p=0.65) Figure 1A. In the Maintenance diet, the difference in global protein expression between the culture
### conditions was less evident (PERMANOVA p=0.08, F=3.7, R2=0.51) and no significant strain effect was detectable (PERMANOVA p=0.67) Figure 1B.
### In both diets, protein expression patterns in the in-vivo culture conditions appeared to have more variable protein expression patterns, 
### However, statistical testing failed to uncover sufficient evidence for this (PERMDISP2 p=0.30 for both maintenance and lactation diets) Figure 2.




t.test(data=lact_meta, disper_dist~condition)
t.test(data=maint_meta, disper_dist~condition)

wilcox.test(data=lact_meta, disper_dist~condition)
wilcox.test(data=maint_meta, disper_dist~condition)

# maint_disper$distances

###########

# strain diffs


lact %>%
  mutate(across(.cols = c(5:10), .fns = as.numeric)) %>% 
  pivot_longer(cols = c(5:10), names_to = 'strain' ) %>% 
  group_by(accno) %>% 
  summarise(VAR=var(value)) %>% 
  ggplot(aes(x=VAR)) +
  geom_histogram()

LACT <- 
  lact %>% 
  mutate(across(.cols = c(5:10), .fns = as.numeric), 
         DIET='lact') %>% 
  pivot_longer(cols = c(5:10), names_to = 'strain' ) %>% 
  mutate(COND=sub('.*_([vivovitro]+)','\\1',strain))
  

MAINT <- 
  maint %>% 
  mutate(across(.cols = c(5:10), .fns = as.numeric), 
         DIET='maint') %>% 
  pivot_longer(cols = c(5:10), names_to = 'strain' ) %>% 
  mutate(COND=sub('.*_([vivovitro]+)','\\1',strain))



ALLL <- bind_rows(LACT, MAINT) %>% 
  mutate(comb=paste(accno, DIET, COND, sep = '_')) %>% 
  group_by(comb) %>% 
  mutate(val2=value - min(value))

ALLL_var <- ALLL %>% 
  group_by(accno, DIET, COND, comb) %>% 
  summarise(VAR=var(value))

ALLL_var %>% 
  ggplot(aes(x=VAR)) + 
  geom_histogram() + 
  facet_wrap(~DIET + COND) + 
  geom_vline(xintercept = .5) +
  ggtitle('histogram of variance of genes within conditions') + 
  scale_y_log10()

these <- 
  ALLL_var %>% 
  filter(!grepl('CON', accno)) %>% 
  filter(VAR > .5) %>%
  pull(comb) %>%
  unique()

library(cowplot)
ALLL %>% 
  filter(DIET == 'lact' & COND == 'vivo') %>% 
  filter(comb %in% these) %>% 
  ggplot(aes(x=accno, y=val2, fill=strain)) + 
  geom_col(position =position_dodge(), width = .75, color = 'white') +
  coord_flip()+
  facet_wrap(~COND + DIET, scales = 'free', ncol = 1) + 
  theme_cowplot()



ALLL %>% 
  filter(DIET == 'lact' & COND == 'vitro') %>% 
  filter(comb %in% these) %>% 
  ggplot(aes(x=accno, y=value, fill=strain)) + 
  geom_col(position =position_dodge(), width = .75, color = 'white') +
  coord_flip()+
  facet_wrap(~COND + DIET, scales = 'free', ncol = 1) + 
  theme_cowplot()


ALLL %>% 
  filter(DIET == 'maint' & COND == 'vivo') %>% 
  filter(comb %in% these) %>% 
  ggplot(aes(x=accno, y=value, fill=strain)) + 
  geom_col(position =position_dodge(), width = .75, color = 'white') +
  coord_flip()+
  facet_wrap(~COND + DIET, scales = 'free', ncol = 1) + 
  theme_cowplot()

ALLL %>% 
  filter(DIET == 'maint' & COND == 'vitro') %>% 
  filter(comb %in% these) %>% 
  ggplot(aes(x=accno, y=value, fill=strain)) + 
  geom_col(position =position_dodge(), width = .75, color = 'white') +
  coord_flip()+
  facet_wrap(~COND + DIET, scales = 'free', ncol = 1) + 
  theme_cowplot()




#




####


lact$pval <- signif(as.numeric(sub('< ','',lact$pval)),3)
maint$pval <- signif(as.numeric(sub('< ','',maint$pval)),3)


lact_sig <- lact %>% filter(pval <= 0.05 & abs(lfc) >=.5) %>% filter(!grepl('CON_', accno))
maint_sig <- maint %>% filter(pval <= 0.05 & abs(lfc) >=.5) %>% filter(!grepl('CON_', accno))
nrow(lact_sig)

sum(lact_sig$lfc >= 0) #111
sum(lact_sig$lfc <=0)  # 43
nrow(lact)

# 154 total sigs

lact_up_vivo <- lact_sig %>% 
  dplyr::select(accno, pval, lfc) %>%
  filter(lfc > 0)
lact_up_vitro <- lact_sig %>%
  dplyr::select(accno, pval, lfc) %>% 
  filter(lfc < 0)


maint_up_vivo <- maint_sig %>%
  dplyr::select(accno, pval, lfc) %>%
  filter(lfc > 0)
maint_up_vitro <- maint_sig %>% 
  dplyr::select(accno, pval, lfc) %>%
  filter(lfc < 0)


### Write sig diffs here ###
# join with psort # 




# Supplemental tables?


####### ADD PSORT IN HERE ########
# test1 <- lact_sig %>% dplyr::select(accno, pval, lfc) %>% left_join(psort) %>% left_join(all_descriptions)
# test2 <- maint_sig %>% dplyr::select(accno, pval, lfc) %>% left_join(psort)



# vitro is reference category #
lact_sig %>% dplyr::select(accno, pval, lfc) %>%
  left_join(all_descriptions) %>% 
  left_join(psort) %>%
  write_csv('./results/iTRAQ_lact_sigs.csv')

maint_sig %>% dplyr::select(accno, pval, lfc) %>% 
  left_join(all_descriptions) %>% 
  left_join(psort) %>% 
  write_csv('./results/iTRAQ_maint_sigs.csv')

#### THESE MAY BE MISSING SOME GO ANNOTATIONS ####

lact_sig_GO <- lact_sig %>%
  left_join(GO_all) %>%
  filter(!is.na(GO)) %>%
  separate_rows(GO, sep = ',') %>%
  mutate(condition=ifelse(lfc>0, 'vivo', 'vitro'), 
         GOID=GO) %>% dplyr::select(-GO) %>%  
  left_join(GO_terms) %>% 
  dplyr::select(c('accno', 'pval', 'lfc', 'condition','GOID', 'TERM')) %>% 
  left_join(all_descriptions) %>% 
  write_csv('./results/iTRAQ_lact_sigs_GO_annotation.csv')

maint_sig_GO <- maint_sig %>%
  left_join(GO_all) %>%
  filter(!is.na(GO)) %>%
  separate_rows(GO, sep = ',') %>%
  mutate(condition=ifelse(lfc>0, 'vivo', 'vitro'), 
         GOID=GO) %>% dplyr::select(-GO) %>%  
  left_join(GO_terms) %>% 
  dplyr::select(c('accno', 'pval', 'lfc', 'condition','GOID', 'TERM')) %>% 
  left_join(all_descriptions) %>% 
  write_csv('./results/iTRAQ_maint_sigs_GO_annotation.csv')


##### Summary ##### 

# Comparing iTRAQ reporter intensities we detected many proteins that were differentially expressed between culture conditions in each of the diets.  
# More proteins were different between the in-vitro and in-vivo conditions in the lactation diet (156) compared to the maintenence diet (81).
# In the lactation diet 43 proteins were more highly expressed in the in-vitro condition and 113 were more highly expressed in the in-vivo condition.
# In the maintenance diet 25 proteins were more highly expressed in the in-vitro condition and 56 were more highly expressed in the 
# in-vivo condition.



######### GO Term enrichment for iTRAQ diffs ##########

GO_lact <- GO_all[GO_all$accno %in% lact$accno,]  # These mappings need to contain all proteins not just sigs
GO_maint <- GO_all[GO_all$accno %in% maint$accno,]

write_delim(GO_lact, delim = '\t', 'lact_gene2GO.txt')
write_delim(GO_maint, delim = '\t', 'maint_gene2GO.txt')


topGO_wrapper <- function(myInterestingGenes, # VECTOR OF GENE NAMES/IDS, MUST MATCH THOSE IN MAPPING FILE
                          mapping_file, # two column file, first column geneIDs, second column ',' delimited GOTerms
                          ont='BP',
                          algor = 'elim',
                          statistic='Fisher', 
                          nodeSize=10){
  
  require(topGO)
  
  # coreGenes <- Int_genes
  
  geneID2GO <- readMappings(mapping_file)
  geneNames <- names(geneID2GO)
  
  # Get the list of genes of interest
  # myInterestingGenes <- coreGenes$accno
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  # head(geneList)
  
  GOdata <- new("topGOdata", ontology = ont, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO, 
                nodeSize=nodeSize)
  # Run topGO with elimination test
  resultTopGO.elim <- runTest(GOdata, algorithm = algor, statistic = statistic )
  allRes <- GenTable(GOdata, pval = resultTopGO.elim,
                     orderBy = "pval", 
                     topNodes = length(GOdata@graph@nodes), #include all nodes
                     numChar=1000)
  allRes <- allRes %>% mutate(ont=ifelse(ont=='BP', 'Biological Process', 
                                         ifelse(ont=='MF', 'Molecular Function', "Cellular Component"))) %>% 
    mutate(GO_aspect = ont, 
           algorithm = algor, 
           statistic = statistic) %>% dplyr::select(-ont)
  return(allRes)
  #write.table(allRes, file = "Lact_vivo_topGO_BP_results.txt", sep = "\t", quote = F, col.names = T, row.names = F)
  
}


###

# these are limited to pval < 0.1

lact_vivo_all <- rbind(topGO_NonModel(Int_genes = lact_up_vivo, mapping_file ='lact_gene2GO.txt', ont = 'BP') ,
                       topGO_NonModel(Int_genes = lact_up_vivo, mapping_file ='lact_gene2GO.txt', ont = 'MF'),
                       topGO_NonModel(Int_genes = lact_up_vivo, mapping_file ='lact_gene2GO.txt', ont = 'CC')) %>% 
  filter(pval < 0.1) %>% dplyr::select(-algorithm, -statistic) %>% write_csv('./results/Lact_vivo_GO.csv')


lact_vitro_all <- rbind(topGO_NonModel(Int_genes = lact_up_vitro, mapping_file ='lact_gene2GO.txt', ont = 'BP') ,
                       topGO_NonModel(Int_genes = lact_up_vitro, mapping_file ='lact_gene2GO.txt', ont = 'MF'),
                       topGO_NonModel(Int_genes = lact_up_vitro, mapping_file ='lact_gene2GO.txt', ont = 'CC')) %>% 
  filter(pval < 0.1) %>% dplyr::select(-algorithm, -statistic) %>% write_csv('./results/Lact_vitro_GO.csv')

#

maint_vivo_all <- rbind(topGO_NonModel(Int_genes = maint_up_vivo, mapping_file ='maint_gene2GO.txt', ont = 'BP') ,
                       topGO_NonModel(Int_genes = maint_up_vivo, mapping_file ='maint_gene2GO.txt', ont = 'MF'),
                       topGO_NonModel(Int_genes = maint_up_vivo, mapping_file ='maint_gene2GO.txt', ont = 'CC')) %>% 
  filter(pval < 0.1) %>% dplyr::select(-algorithm, -statistic) %>% write_csv('./results/maint_vivo_GO.csv')


maint_vitro_all <- rbind(topGO_NonModel(Int_genes = maint_up_vitro, mapping_file ='maint_gene2GO.txt', ont = 'BP') ,
                        topGO_NonModel(Int_genes = maint_up_vitro, mapping_file ='maint_gene2GO.txt', ont = 'MF'),
                        topGO_NonModel(Int_genes = maint_up_vitro, mapping_file ='maint_gene2GO.txt', ont = 'CC')) %>% 
  filter(pval < 0.1) %>% dplyr::select(-algorithm, -statistic) %>% write_csv('./results/maint_vitro_GO.csv')




#############  iBAQ from MaxQuant ##########

list.files(pattern = 'proteinGroups.txt', recursive = T, full.names = T)
MQ_prot_groups <- read_delim('../O157_first/combined/txt/proteinGroups.txt', delim = '\t') %>% filter(`Q-value` < 0.05, Peptides > 1)
# MQ_prot_groups <- read_delim('../combined/txt/proteinGroups.txt', delim = '\t') %>% filter(`Q-value` < 0.05, Peptides > 1)
# MQ_iBAQ <- MQ_iBAQ %>% filter(Peptides >1) %>%
#   select(-starts_with('Fraction'), -starts_with('Reporter intensity count'))

maj_accnos <- strsplit(MQ_prot_groups$`Majority protein IDs`, split = ';') %>% map(1) %>% unlist()

intensities <- c('Intensity', 'Intensity 1', 'Intensity 2', 'iBAQ', 'iBAQ 1', 'iBAQ 2')
info <- c('accno', 'Protein names', 'Gene names','Peptides', 'Peptides 1', 'Peptides 2', 'Q-value')


MQ_prot_groups <- MQ_prot_groups %>% mutate(accno = maj_accnos) %>% dplyr::select(accno, everything())


MQ_iBAQ <- MQ_prot_groups %>% dplyr::select(info, intensities)

MQ_iBAQ <- MQ_iBAQ %>% filter(Peptides >1) %>% 
  filter(!grepl('CON', accno)) %>% filter(iBAQ !=0)  # remove contaminant proteins and those with no iBAQ info


# riBAQ
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3946283/

MQ_iBAQ$rel_iBAQ1 <- MQ_iBAQ$`iBAQ 1`/sum(MQ_iBAQ$`iBAQ 1`)*100
MQ_iBAQ$rel_iBAQ2 <- MQ_iBAQ$`iBAQ 2`/sum(MQ_iBAQ$`iBAQ 2`)*100


# MQ_iTRAQ <-MQ_prot_groups %>% dplyr::select(info, intensities, starts_with('Reporter intensity corrected'))
# MQ_iTRAQ <- MQ_iTRAQ[,colnames(MQ_iTRAQ)[c(1:4, 12:13, 22:37)]]
# 
# MQ_iTRAQ <- MQ_iTRAQ[,colnames(MQ_iTRAQ)[-c(13,14,21,22)]]
# 
# 
# 
# 
# colnames(MQ_iTRAQ)[7:18] <- c('86-24_vitro_lact', 'EDL933_vitro_lact', 'SS17_vitro_lact',
#                               '86-24_vivo_lact', 'EDL933_vivo_lact', 'SS17_vivo_lact',
#                               '86-24_vitro_maint', 'EDL933_vitro_maint', 'SS17_vitro_maint',
#                               '86-24_vivo_maint', 'EDL933_vivo_maint', 'SS17_vivo_maint')
# 
# 
# MQ_iTRAQ_gath <- MQ_iTRAQ %>% gather(key='type', value='value', -c(1:3))
# 
# 
# MQ_iTRAQ_gath %>% filter(!(type %in% c('Peptides', 'iBAQ 1', 'iBAQ 2'))) %>%
#   mutate(value=log10(value)) %>% 
#   ggplot(aes(x=type, y=value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90))
# 
# 
# MQ_iTRAQ_gath %>% filter(!(type %in% c('Peptides', 'iBAQ 1', 'iBAQ 2'))) %>% group_by(type) %>% 
#   mutate(log_intensity=log(value+1), 
#          tot_intensity=sum(log_intensity), 
#          rel_intensity=log_intensity/tot_intensity)  %>% 
#   ggplot(aes(x=type, y=rel_intensity)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90))
# 
# 
# 


#######################################
# I dont think the same amount of protein was loaded for these two runs #
# it seems to me like a more concentrated solution was loaded for the maintenance diet...


### Figure S1A
MQ_iBAQ %>% mutate(lact_iBAQ = `iBAQ 1`,
                   maint_iBAQ = `iBAQ 2`) %>%
  dplyr::select(accno, ends_with('_iBAQ')) %>% gather(key = 'diet', value = 'iBAQ', -accno) %>%
  mutate(diet=sub('_iBAQ', '', diet)) %>% 
  ggplot(aes(x=diet, y=iBAQ)) + geom_col(aes(fill=diet), color='black') + ylab('iBAQ (total intensity)') + theme_bw()


### Figure S1B
MQ_iBAQ %>% mutate(lact_iBAQ = rel_iBAQ1,
                   maint_iBAQ = rel_iBAQ2) %>%
  dplyr::select(accno, ends_with('_iBAQ')) %>% gather(key = 'diet', value = 'iBAQ', -accno) %>%
  mutate(diet=sub('_iBAQ', '', diet)) %>% 
  ggplot(aes(x=diet, y=iBAQ)) + geom_col(aes(fill=diet), color='black') + ylab('relative iBAQ') + theme_bw()


MQ_iBAQ <- MQ_iBAQ %>% mutate(FC=rel_iBAQ1/rel_iBAQ2, 
                   L2FC=log2(FC)) %>% dplyr::select(accno, `Protein names`, `Gene names`, contains('iBAQ'), L2FC)


# GO analysis on these 4 

iBAQ_up_lact <- MQ_iBAQ %>% filter(L2FC > 1) %>% dplyr::select(-iBAQ) %>% left_join(psort) %>% left_join(all_descriptions)
iBAQ_up_maint <- MQ_iBAQ %>% filter(L2FC < -1) %>% dplyr::select(-iBAQ) %>% left_join(psort) %>% left_join(all_descriptions)

iBAQ_only_maint <- MQ_iBAQ %>% filter(rel_iBAQ1 ==0) %>% dplyr::select(-iBAQ)
iBAQ_only_lact <- MQ_iBAQ %>% filter(rel_iBAQ2 ==0) %>% filter(iBAQ !=0) %>% dplyr::select(-iBAQ)


#########

iBAQ_diff_all <- rbind(iBAQ_up_lact, iBAQ_up_maint) 



iBAQ_diff_all_GO_ANNO <- iBAQ_diff_all %>%
  left_join(GO_all) %>%
  filter(!is.na(GO)) %>%
  separate_rows(GO, sep = ',') %>%
  mutate(diet=ifelse(L2FC > 0, 'lactation', 'maintenance'), 
         GOID=GO) %>% dplyr::select(-GO) %>%  
  left_join(GO_terms) %>% 
  dplyr::select(c('accno', 'L2FC', 'diet','GOID', 'TERM', 'psort_loc', 'ONTOLOGY', 'description')) %>% 
  # left_join(all_descriptions) %>% 
  write_csv('./results/iBAQ_difs_GO_annotation.csv')



#############







sum(MQ_iBAQ$`iBAQ 1` == 0) # 179 proteins unique to maintenance diet
sum(MQ_iBAQ$`iBAQ 2` == 0) # 11 proteins unique to lactation diet
sum(MQ_iBAQ$`iBAQ 2` != 0 & MQ_iBAQ$`iBAQ 1` != 0) # 566 proteins found in both diets

nrow(MQ_iBAQ) # 756 total proteins detected

nrow(iBAQ_up_lact)
nrow(iBAQ_up_maint)

write_csv(iBAQ_up_lact, './results/iBAQ_up_lact.csv')
write_csv(iBAQ_up_maint, './results/iBAQ_up_maint.csv')


### VENNS ###

library(ggforce)

# LACT VS MAINT IBAQ unique proteins

df.venn <- data.frame(x = c(0.766, -0.766),
                      y = c(-0, -0),
                      Diet = c('Lactation', 'Maintenance'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = Diet)) +
  geom_circle(alpha = .3, size = .5, colour = 'black') +
  annotate('text', x=-1.5, y=0, label='179', size=10)+
  annotate('text', x=1.5, y=0, label='11', size=10)+
  annotate('text', x=0, y=0, label='566', size=10)+
  annotate('text', x=0, y=-2, label='756 total proteins', size=10)+
    coord_fixed() +
  theme_void() + ggtitle('Proteins unique to each diet')


# LACT VS MAINT IBAQ diff proteins

nrow(iBAQ_up_lact) # 101 up in lact
nrow(iBAQ_up_maint) # 348 up in maint
756 - (101+348) # 307 not different between diets


df.venn <- data.frame(x = c(0.766, -0.766),
                      y = c(-0, -0),
                      Diet = c('Lactation', 'Maintenance'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = Diet)) +
  geom_circle(alpha = .3, size = .5, colour = 'black') +
  annotate('text', x=-1.5, y=0, label='348', size=10)+
  annotate('text', x=1.5, y=0, label='101', size=10)+
  annotate('text', x=0, y=0, label='307', size=10)+
  annotate('text', x=0, y=-2, label='756 total proteins', size=10)+
  coord_fixed() +
  theme_void() + ggtitle('Proteins with L2FC expression > 1 between diets')



##################

# iTRAQ Vitro VS vivo MAINT


nrow(maint) # 667 lact proteins
nrow(maint_up_vitro) # 25 up in vitro
nrow(maint_up_vivo) # 56 up vivo
667-(25+56) # 586 not diff


df.venn <- data.frame(x = c(0.766, -0.766),
                      y = c(-0, -0),
                      Condition = c('in-vitro', 'in-vivo'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = Condition)) +
  geom_circle(alpha = .3, size = .5, colour = 'black') +
  annotate('text', x=-1.5, y=0, label='56', size=10)+
  annotate('text', x=1.5, y=0, label='25', size=10)+
  annotate('text', x=0, y=0, label='586', size=10)+
  annotate('text', x=0, y=-2, label='667 total proteins', size=10)+
  coord_fixed() +
  theme_void() + 
  scale_fill_manual(values = c('blue', 'red'))+
  # scale_fill_brewer(palette = 'Set2')+
  ggtitle('Differentially expressed proteins, Maintenance diet')


# iTRAQ Vitro VS vivo LACT


nrow(lact) # 473 lact proteins
nrow(lact_up_vitro) # 43 up in vitro
nrow(lact_up_vivo) # 111 up vivo
473-(43+111) # 319 not diff


df.venn <- data.frame(x = c(0.766, -0.766),
                      y = c(-0, -0),
                      Condition = c('in-vitro', 'in-vivo'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = Condition)) +
  geom_circle(alpha = .3, size = .5, colour = 'black') +
  annotate('text', x=-1.5, y=0, label='111', size=10)+
  annotate('text', x=1.5, y=0, label='43', size=10)+
  annotate('text', x=0, y=0, label='319', size=10)+
  annotate('text', x=0, y=-2, label='473 total proteins', size=10)+
  coord_fixed() +
  theme_void() + 
  scale_fill_manual(values = c('blue', 'red'))+
  # scale_fill_brewer(palette = 'Set2') +
  ggtitle('Differentially expressed proteins, Lactation diet')

############################



# words
#### Because the our iTRAQ experiments were designed to give the highest resolution of differences between culture conditions within 
#### each diet (all iTRAQ labeled Maintenence diet samples pooled in one MS run and all iTRAQ labeled Lactation diet samples in annother).
#### we were unable to effectively use iTRAQ reporter intensity to compare protein expression between the diets.  However, because all samples pooled in 
#### each MS run belonged to the same diet, we are able to use iBAQ to get a rough estimate of the difference in protein expression between the diets.
#### These reported iBAQ values represent pooled protein expression for all three strains in both the in-vitro and in-vivo culture conditions for each diet.
#### As we suspected, the strains did not thrive under the conditions for the lactation diet and as such a much lower total iBAQ value was detected (FIGURE S1).
#### We converted these iBAQ values to riBAQ values to look at relative protein abundance, that is, within each diet what proportion of the total protein pool 
#### does each protein comprise? (Citation # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3946283/ ) FIGURE S2. However, caution should be used when interpreting
#### these results as compositional effects may lead to some erroneous conclusions.(more proteins under the limit of detection in the lactation diet
#### may lead to an artificial inflation of the proportions of the detected proteins)
#### We detected a total of 756 unique proteins across both diets, 566 proteins were detected in both diets, 179 were unique to the Maintenance diet and 11
#### were unique to the lactation diet.  We calculated log2FoldChange values between the diets using riBAQ values and considered any protein
#### with an absolute L2FC value of 1 or greater to be enriched in the respective diet.  Using these criteria 347 proteins were enriched in
#### the maintenance diet (Table S1A) and 101 were enriched in the lactation diet(Table S1B).





#######


# iBAQ GO mappings (from interproscan)

GO_iBAQ <- GO_all[GO_all$accno %in% MQ_iBAQ$accno,]

write_delim(GO_iBAQ, delim = '\t', 'iBAQ_gene2GO.txt')  ### THESE NEED TO CONTAIN THE ALL GENES DETECTED NOT JUST SIGS



iBAQ_GO_enrich_lact <- rbind(topGO_NonModel(Int_genes = iBAQ_up_lact, mapping_file ='iBAQ_gene2GO.txt', ont = 'BP', algor = 'elim', statistic = 'fisher'),
                             topGO_NonModel(Int_genes = iBAQ_up_lact, mapping_file ='iBAQ_gene2GO.txt', ont = 'CC', algor = 'elim', statistic = 'fisher'),
                             topGO_NonModel(Int_genes = iBAQ_up_lact, mapping_file ='iBAQ_gene2GO.txt', ont = 'MF', algor = 'elim', statistic = 'fisher'))


iBAQ_GO_enrich_maint <- rbind(topGO_NonModel(Int_genes = iBAQ_up_maint, mapping_file ='iBAQ_gene2GO.txt', ont = 'BP', algor = 'elim', statistic = 'fisher'),
                              topGO_NonModel(Int_genes = iBAQ_up_maint, mapping_file ='iBAQ_gene2GO.txt', ont = 'CC', algor = 'elim', statistic = 'fisher'),
                              topGO_NonModel(Int_genes = iBAQ_up_maint, mapping_file ='iBAQ_gene2GO.txt', ont = 'MF', algor = 'elim', statistic = 'fisher'))



iBAQ_GO_enrich_lact %>%
  filter(pval < 0.10) %>% 
  dplyr::select(-algorithm, -statistic) %>% 
  write_csv("./results/IBAQ_GO_lact.csv")



iBAQ_GO_enrich_maint %>%
  filter(pval < 0.10) %>% 
  dplyr::select(-algorithm, -statistic) %>% 
  write_csv("./results/IBAQ_GO_maint.csv")


##### 
write_lines(list.files('./results/'), path = './results/file_descriptions.csv')

