################################################################################
#  __  __ _____ _____  ___  _____  ___  _   _ __   _  ___
# |  \/  | ____|  ___|/ _ \|  ___|/ _ \| | | |  \ | |/ _ \
# | |\/| |  _| | | __  |_| | |_  | |_| | | | |   \| | |_| |
# | |  | | |___| |_| | | | |  _| | | | | |_| | |\   | | | |
# |_|  |_|_____|_____|_| |_|_|   |_| |_|_____|_| \__|_| |_|

################################################################################
##                              ARCHAEOME                     ##
################################################################################
## Preparation of your working space with the library packaging ####
rm(list=ls())
cran_packages   <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny", "microbiome",
                     "tidyverse", "miniUI", "caret", "pls", "e1071", "ggplot2",
                     "randomForest","entropart", "vegan", "plyr", "dplyr","here", "ggrepel", "nlme", "R.utils", "gridExtra", "googledrive",
                     "googlesheets", "phangorn", "devtools", "rmarkdown", "sys",
                     "reshape2", "devtools", "PMA","structSSI","ade4", "ape",
                     "Biostrings", "igraph", "ggnetwork", "intergraph", "ips",
                     "scales", "kableExtra", "pgirmess", "treemap", "microbiome")
github_packages <- c("jfukuyama/phyloseqGraphTest")
bioc_packages   <- c("phyloseq", "genefilter", "impute", "dada2", "DECIPHER")
inst <- cran_packages %in% installed.packages()
if (any(! inst)) {
  install.packages(cran_packages[!inst], repos = "http://cran.rstudio.com/") }
inst <- github_packages %in% installed.packages()
if (any(! inst)) {
  devtools::install_github(github_packages[!inst]) }
inst <- bioc_packages %in% installed.packages()
if (any(! inst)) {
  source("http://bioconductor.org/biocLite.R")
  BiocManager::install(bioc_packages[!inst]) }
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)
library(dada2)
library(stringr)


### Files preparation ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/Presentations/Archaeome")
dir.create("./Physeq_objects", recursive = T)
dir.create("./figure", recursive = T)
#load("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/03_assign_taxonomy/ps_all.RData")

ps_all <- subset_samples(ps_all, region != "Mayotte")
save(ps_all, file = "./Physeq_objects/ps_all.RData")
ps_gut <- subset_samples(ps_all, type == "gut")
ps_gut <- prune_taxa(names(which(colSums(ps_gut@otu_table)>0)), ps_gut)
save(ps_gut, file = "./Physeq_objects/ps_gut.RData")
load("./Physeq_objects/ps_all.RData") ; load("./Physeq_objects/ps_gut.RData")

## ** Archae ----
arch <- subset_taxa(ps_all, Kingdom == "Archaea")
arch <- prune_samples(sample_sums(arch) >0 , arch)
sum(sample_sums(arch))/sum(sample_sums(ps_all)) *100
save(arch, file = "./Physeq_objects/arch.RData")
arch_gut <- subset_samples(arch, type == "gut")
arch_gut <- prune_taxa(names(which(colSums(arch_gut@otu_table)>0)), arch_gut)
save(arch_gut, file ="./Physeq_objects/arch_gut.RData")
sum(sample_sums(arch_gut))/sum(sample_sums(ps_gut)) *100
arch_env <- subset_samples(arch, type != "gut")
arch_env <- prune_taxa(names(which(colSums(arch_env@otu_table)>0)), arch_env)
save(arch_env, file ="./Physeq_objects/arch_env.RData")
load("./Physeq_objects/arch.RData") ; load("./Physeq_objects/arch_gut.RData") ;  load("./Physeq_objects/arch_env.RData")

DAT_arch <- sample_data(arch) %>% as.data.frame()
DAT_all <- sample_data(ps_all) %>% as.data.frame()
table(DAT_arch$diet3)/table(DAT_all$diet3)*100
## _____ Figure Pies gut----
## __ A) Pie kingdom ----
pdf(file = "./figure/piechart.pdf", he = 7, wi = 7)
pie(c(99.915, 0.039,0.046), 
    labels = c("Bacteria","Archaea", "Eukarya") , 
    border="red", 
    col= c("darkblue", "darkred", "yellow"))
dev.off()

## __ B) Treemap gut ----
physeq_class_arch <- arch_gut %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Class)

physeq_class_arch$Class = as.character(physeq_class_arch$Class) # Avoid error message with factor for next step
group <- physeq_class_arch$Phylum
subgroup <- physeq_class_arch$Class
value <- physeq_class_arch$Abundance
data <- data.frame(group,subgroup, value)  

# treemap
pdf(file = "arch_treemap.pdf" , he = 7, wi = 7)
treemap(data,
        index=c("group","subgroup"), vSize = "value", type = "index",
        fontcolor.labels=c("white","black"),
        fontsize.labels=c(12),bg.labels=c("transparent"),
        fontface.labels=c(2,3),
        border.col=c("black","white"), border.lwds=c(4,2), 
        align.labels=list(c("center", "center"),c("left", "bottom")),
        title="Enteric Archaeome Treemap",fontsize.title=12,
        fontfamily.title ="serif")

dev.off()
## _____ PERMANOVA ----
## __A) Nature sample ----
arch_rel <- transform_sample_counts(arch, function(x) x / sum(x) )
save(arch_rel , file ="./Physeq_objects/core_algae_rel.RData")

otu <- vegdist(arch_rel@otu_table, method = "bray")
pcoa.sub <- pcoa(otu)
save(pcoa.sub,otu, file = "./Physeq_objects/beta_objects.RData")

pcoa_coord <- pcoa.sub$vectors[,1:3]
library(stringr)
DAT_arch <- data.frame(sample_data(arch))
write.table(DAT_arch, file = "./Physeq_objects/DAT_arch.txt", row.names = F, sep ="\t")
hull <- cbind(pcoa_coord, DAT_arch)
hull$type <- str_replace_all(hull$type, c("gut"= "Gut" , "sand"="Sediment", "seawater" = "Seawater"))
save(hull, file ="./Physeq_objects.RData")

pcoa <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  geom_point(data = hull, aes(x=Axis.1, y=Axis.2, color = type), alpha = 0.7, size = 5, shape = 16) +
  scale_color_manual(values = c("Seawater" = "darkblue","Gut"="darkred","Sediment"="gray"))+
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(family = "serif", size=14),
        axis.title.y = element_text(family = "serif",size=14),
        axis.text.x = element_text(family = "serif", size = 13), 
        axis.text.y = element_text(family = "serif", size = 13), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 14),
        legend.text = element_text(size = 11, family = "serif"),
        legend.title = element_text(size = 11,family = "serif")) +
  labs(size = 14, colour = "Type of sample") 

pdf(file ="./figure/pcoa.pdf", he = 7, wi= 7)
pcoa
dev.off()

DAT_arch <- read.table(file = "./Physeq_objects/DAT_arch.txt", sep ="\t", header= T)
adonis2(otu_tst ~ type, data = DAT_tst)

## __B) Gut traits ----
gut_rel <- transform_sample_counts(arch_gut, function(x) x / sum(x) )
save(gut_rel , file ="./Physeq_objects/gut_rel.RData")

otu <- vegdist(gut_rel@otu_table, method = "bray")
pcoa.sub <- pcoa(otu)
save(pcoa.sub,otu, file = "./Physeq_objects/beta_objects_gut.RData")
pcoa_coord <- pcoa.sub$vectors[,1:3]
DAT_arch_gut <- data.frame(sample_data(arch_gut))
write.table(DAT_arch, file = "./Physeq_objects/DAT_arch_gut.txt", row.names = F, sep ="\t")

hull <- cbind(pcoa_coord, DAT_arch_gut)
hull$type <- str_replace_all(hull$type, c("gut"= "Gut" , "sand"="Sediment", "seawater" = "Seawater"))
save(hull, file ="./Physeq_objects.RData")


pcoa <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  geom_point(data = hull, aes(x=Axis.1, y=Axis.2, color =diet3), alpha = 0.7, size = 5, shape = 16) +
  scale_color_manual(values = c("Carnivorous" = "darkred","Herbivorous"="darkgreen","Omnivorous"="gray"))+
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(family = "serif", size=14),
        axis.title.y = element_text(family = "serif",size=14),
        axis.text.x = element_text(family = "serif", size = 13), 
        axis.text.y = element_text(family = "serif", size = 13), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 14),
        legend.text = element_text(size = 11, family = "serif"),
        legend.title = element_text(size = 11,family = "serif")) +
  labs(size = 14, colour = "Trophic guild")

pdf(file = "./figure/diet.pdf", he = 7, wi = 7)
pcoa
dev.off()

adonis2(otu ~ diet3, data = DAT_arch_gut)
adonis2(otu ~ family, data = DAT_arch_gut)
adonis2(otu ~ genus, data = DAT_arch_gut)
adonis2(otu ~ tax1, data = DAT_arch_gut)
adonis2(otu ~ region, data = DAT_arch_gut)

## **  DSR ----
tax  <- ps_all %>% tax_table() %>% as.data.frame()
dsr <- tax %>% filter(str_detect(Genus, "Desulf") | 
                        str_detect(Family,"Desulf") |
                        str_detect(Genus, "desulf") |
                        str_detect(Family, "desulf"))

exclude <- c("Lawsonia", "Bilophila")
dsr2 <- dsr %>% filter(!Genus %in% exclude)
ps_dsr <- subset_taxa(ps_all , taxa_names(ps_all) %in% rownames(dsr2))
ps_dsr <- prune_samples(sample_sums(ps_dsr) >0 , ps_dsr)
ps_dsr <- prune_taxa(names(which(colSums(ps_dsr@otu_table)>0)), ps_dsr)
save(ps_dsr, file ="./Physeq_objects/ps_dsr.RData")

ps_dsr_gut <- subset_samples(ps_dsr, type == "gut")
ps_dsr_gut <- prune_taxa(names(which(colSums(ps_dsr_gut@otu_table)>0)), ps_dsr_gut)
save(ps_dsr_gut, file = "./Physeq_objects/ps_dsr_gut.RData")
load(file = "./Physeq_objects/ps_gut.RData")
sum(sample_sums(ps_dsr_gut))/sum(sample_sums(ps_gut))*100

## ** Venn ----
##__ A) Arch----
merged = merge_samples(arch, "type") #merge samples for herbivores, carnivores and the macroalgae
rownames(merged@otu_table@.Data)
min(rowSums(merged@otu_table@.Data)) #how many reads per sample
set.seed(10000)
venn = rarefy_even_depth(merged, sample.size = 84)

gut <- colnames(merged@otu_table@.Data)[merged@otu_table@.Data["gut",] > 0]
sw <- colnames(merged@otu_table@.Data)[merged@otu_table@.Data["seawater",] > 0]
sed <- colnames(merged@otu_table@.Data)[merged@otu_table@.Data["sand",] > 0]
source("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
pdf(file = "./figure/Venn_arch.pdf", he = 7, wi = 7)
venn_comp <- venn_diagram3(gut,sw,sed, "Gut","Seawater","Sediment", colors= c("darkred","darkblue","gray"), euler=FALSE)
dev.off()

## __ B) SRB----
merged = merge_samples(ps_dsr, "type") #merge samples for herbivores, carnivores and the macroalgae
rownames(merged@otu_table@.Data)
min(rowSums(merged@otu_table@.Data)) #how many reads per sample
set.seed(10000)
venn = rarefy_even_depth(merged, sample.size = 243)

gut <- colnames(venn@otu_table@.Data)[venn@otu_table@.Data["gut",] > 0]
algae <- colnames(venn@otu_table@.Data)[venn@otu_table@.Data["algae",] > 0]
sw <- colnames(venn@otu_table@.Data)[venn@otu_table@.Data["seawater",] > 0]
sed <- colnames(venn@otu_table@.Data)[venn@otu_table@.Data["sand",] > 0]
source("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")
pdf(file = "./figure/Venn_dsr.pdf", he = 7, wi = 7)
venn_comp <- venn_diagram4(gut,algae,sw,sed, "Gut", "Turf", "Seawater","Sediment", colors= c("darkred","darkgreen","darkblue","gray"), euler=FALSE)
dev.off()

## ** Abondance relative ----
#__ A) SRB type ----
ps_all2 <- subset_samples(ps_all, sample_names(ps_all) %in% sample_names(ps_dsr))
dsr_ab <- cbind(sample_data(ps_dsr), "Abundance" = sample_sums(ps_dsr), "Relative"= sample_sums(ps_dsr)/sample_sums(ps_all2)*100)
dsr_ab$type <- str_replace_all(dsr_ab$type ,c("gut"="Gut", "sand"="Sediment","seawater"="Seawater", "algae"="Turf"))
dsr_ab$diet3 <- str_replace_all(dsr_ab$diet3, c("turf" = "algae", "macroalgae"= "algae"))
write.table(dsr_ab, file = "./Physeq_objects/dsr_ab.txt", sep = "\t", row.names = F)
gut_dsr_ab <- dsr_ab %>% filter(type == "Gut") %>% summarise(Mean = mean(Relative) , SD = sd(Relative), Min= min(Relative), Max=max(Relative))
sed_dsr_ab <- dsr_ab %>% filter(type == "Sediment") %>% summarise(Mean = mean(Relative) , SD = sd(Relative), Min= min(Relative), Max=max(Relative))
sw_dsr_ab <- dsr_ab %>% filter(type == "Seawater") %>% summarise(Mean = mean(Relative) , SD = sd(Relative), Min= min(Relative), Max=max(Relative))
turf_dsr_ab <- dsr_ab %>% filter(type == "Turf") %>% summarise(Mean = mean(Relative) , SD = sd(Relative), Min= min(Relative), Max=max(Relative))
dsr_ab_resume <- cbind(rbind(round(gut_dsr_ab,2),round(sed_dsr_ab,2),round(sw_dsr_ab,2),round(turf_dsr_ab,2)), "Type"= c("Gut","Sediment","Seawater","Turf"))
write.table(dsr_ab_resume, file = "dsr_ab_resume.txt", sep = "\t", row.names=F)

kruskal.test(dsr_ab$Relative,dsr_ab$type)
dunn <- dunn.test(dsr_ab$Relative,dsr_ab$type, method = "bonferroni")

#__ B) SRB family ----
dsr_ab_gut <- dsr_ab %>% filter(type == "Gut")
write.table(dsr_ab_gut, file = "./Physeq_objects/dsr_ab_gut.txt", sep = "\t", row.names = F)
kruskal.test(dsr_ab_gut$Relative,dsr_ab_gut$diet3) # No difference
kruskal.test(dsr_ab_gut$Relative,dsr_ab_gut$family)
dunn <- dunn.test(dsr_ab_gut$Relative,dsr_ab_gut$family, method = "bonferroni")
sign <- dunn$comparisons[which(dunn$P.adjusted <= 0.05)]

plot2 <- ggplot(dsr_ab_gut, aes(x = family , y = Relative)) + 
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) +
  ylab("SRB Relative abundance (%)")+
  theme_bw()+
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")
pdf(file = "./figure/rel_srb_fam.pdf", he = 5, wi = 7)
plot2
dev.off()

#** T4F----
#__ A) dsrB ----
load("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Carnivory/Functions/ps_fun.RData")
tax_fun <- ps_fun %>% tax_table() %>% as.data.frame()
dsrb_ko <- tax_fun %>% filter(str_detect(ta6,"dissimilatory") | 
                        str_detect(ta6,"adenylylsulfate"))

ps_dsrB <- subset_taxa(ps_fun, taxa_names(ps_fun) %in% rownames(dsrb_ko))
ps_dsrB <- prune_samples(sample_sums(ps_dsrB) >0 , ps_dsrB)
ps_dsrB <- prune_taxa(names(which(colSums(ps_dsrB@otu_table)>0)), ps_dsrB)
save(ps_dsrB, file ="./Physeq_objects/ps_dsrB.RData")

ps_fun2 <- subset_samples(ps_fun, sample_names(ps_fun) %in% sample_names(ps_dsrB))  
sum(sample_sums(ps_dsrB))/sum(sample_sums(ps_fun2)) *100
dsrB_ab <- cbind(sample_data(ps_dsrB), "Abundance" = sample_sums(ps_dsrB), "Relative"= sample_sums(ps_dsrB)/sample_sums(ps_fun2)*100)
write.table(dsrB_ab, file = "./Physeq_objects/dsrB_ab.txt", sep = "\t", row.names = F)
kruskal.test(dsrB_ab$Relative,dsrB_ab$diet3) # No difference
kruskal.test(dsrB_ab$Relative,dsrB_ab$family)
dunn_dsrB <- dunn.test(dsrB_ab$Relative,dsrB_ab$family, method = "bonferroni")
sign_dsrB <- dunn_dsrB$comparisons[which(dunn_dsrB$P.adjusted <= 0.05)]
gut_dsrB_ab <- dsrB_ab %>%summarise(Mean = mean(Relative) , SD = sd(Relative), Min= min(Relative), Max=max(Relative))

#__ B) mcrA ----
tax_fun <- ps_fun %>% tax_table() %>% as.data.frame()
mcrA_ko <- tax_fun %>% filter(str_detect(ta6,"methyl-coenzyme M reductase "))

ps_mcrA <- subset_taxa(ps_fun, taxa_names(ps_fun) %in% rownames(mcrA_ko))
ps_mcrA <- prune_samples(sample_sums(ps_mcrA) >0 , ps_mcrA)
ps_mcrA <- prune_taxa(names(which(colSums(ps_mcrA@otu_table)>0)), ps_mcrA)
save(ps_mcrA, file ="./Physeq_objects/ps_mcrA.RData")
# --> Only 3 samples on 378 ! 
ps_fun2 <- subset_samples(ps_fun, sample_names(ps_fun) %in% sample_names(ps_mcrA))  
sum(sample_sums(ps_mcrA))/sum(sample_sums(ps_fun2)) *100
mcrA_ab <- cbind(sample_data(ps_mcrA), "Abundance" = sample_sums(ps_mcrA), "Relative"= sample_sums(ps_mcrA)/sample_sums(ps_fun2)*100)

##__ C) Both----
dsrB_ab <- read.table(file = "./Physeq_objects/dsrB_ab.txt", sep="\t",  header=T)
fun_ab <- left_join(dsrB_ab, mcrA_ab %>% select(tax3, Relative), by = "tax3")

fun_ab <- cbind("Gene" = rep("Sulfato-reduction",nrow(dsrB_ab)), dsrB_ab)
fun2_db <- cbind("Gene" = rep("Methanogenesis",nrow(mcrA_ab)), mcrA_ab)
gene_interest <- rbind(fun_ab,fun2_db)

gene = ggplot(gene_interest, aes(x = family , y = Relative, color = Gene)) + 
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + 
  theme_bw() +
  scale_color_manual(values = c("Sulfato-reduction"="black", "Methanogenesis"="red")) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title = element_text(family = "serif"),
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.text.y = element_text(family = "serif",size = 12),
        legend.position = c(.05, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 13, family = "serif", face= "bold"))+
  ylab("Predictional gene contribution (%)")+
  labs(colour = "Metabolism" )

pdf(file = "./figure/dsrb_mcrA_fam.pdf", he = 5, wi = 7)
gene
dev.off()

## qPCR dsrB results----
metadata_dsrb <- read.table(file = "metadata_dsrb.txt", sep = "\t", header = T)
data <- read.table(file = "gut_data.txt",header= T, sep ="\t")
names(data)[c(3,5)] <- c("ID" , "tax2")
gut_pcr <- metadata_dsrb %>% filter(sample_type == "gut") %>% select(Nom, Good)
names(gut_pcr)[1] <- "ID"
data <- data %>% select(ID, tax2, family, diet1, diet2, diet3)
gut_data_dsrb <- left_join(gut_pcr , data, by = "ID")
write.table(gut_data_dsrb , file ="gut_data_dsrb.txt", sep ="\t", row.names = F)
gut_data_dsrb <- read.table(file ="gut_data_dsrb.txt", sep ="\t", header=T)
dsr_p <- gut_data_dsrb %>% filter(Good == "Y") %>% unique()
table(dsr_p$family)
dsr_a <- gut_data_dsrb %>% filter(!ID %in% dsr_p$ID) %>% unique()
dsr_data_pcr <- rbind(cbind(dsr_p , "Amplified" = rep("Amplified")),
                      cbind(dsr_a , "Amplified" = rep("No amplified")))

fam_success <- as.data.frame(table(dsr_data_pcr$family, dsr_data_pcr$Amplified))
table(dsr_data_pcr$family, dsr_data_pcr$Amplified)[,1]/(table(dsr_data_pcr$family, dsr_data_pcr$Amplified)[,1]+table(dsr_data_pcr$family, dsr_data_pcr$Amplified)[,2])

plot1 <- ggplot(fam_success, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Relative Success of dsrB gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

fam_success2 <- fam_success
fam_success2[which(fam_success2$Var2 == "No amplified"),]$Freq <- -fam_success2[which(fam_success2$Var2 == "No amplified"),]$Freq 
plot1_B <- ggplot(fam_success2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="identity")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Success of dsrB gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

pdf(file ="./figure/plot_dsrB_success_fam_bis.pdf", he = 7, wi = 7)
plot1_B
dev.off()

diet_success <- as.data.frame(table(dsr_data_pcr$diet3, dsr_data_pcr$Amplified))
plot2 <- ggplot(diet_success, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Relative Success of dsrB gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 0,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

pdf(file ="./figure/plot_dsrB_success_diet3.pdf", he = 7, wi = 7)
plot2
dev.off()

diet_success2 <- diet_success
diet_success2[which(diet_success2$Var2 == "No amplified"),]$Freq <- -diet_success2[which(diet_success2$Var2 == "No amplified"),]$Freq 
plot2_B <- ggplot(diet_success2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="identity")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Success of dsrB gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

pdf(file ="./figure/plot_dsrB_success_diet3_bis.pdf", he = 7, wi = 7)
plot2_B
dev.off()

## qPCR mcrA results----
metadata_mcra <- read.table(file = "metadata_mcra.txt", sep = "\t", header = T)
data <- read.table(file= "gut_data.txt", sep ="\t", header= T)
gut_pcr_mcra <- metadata_mcra %>% filter(sample_type == "gut") %>% select(Nom, Good)
names(gut_pcr_mcra)[1] <- "ID"
gut_data_mcra <- left_join(gut_pcr_mcra , data, by = "ID")
na <- gut_data_mcra[which(is.na(gut_data_mcra$tax2), T),]
data_bis <- read.table(file= "data.txt", sep ="\t", header= T)
names(data_bis)[3] <- "ID"
na_data <- left_join(na[,c(1,2)], data_bis[,c(3,5,13,14,15,16)], by ="ID")
write.table(na_data, file = "na_data.txt", sep = "\t" , row.names=F)
na_data <- read.table(file = "na_data.txt", sep = "\t", header = T)
gut_data_mcra2 <- rbind(gut_data_mcra %>% filter(!ID %in% na_data$ID), na_data)
write.table(gut_data_mcra2, file = "gut_data_mcra.txt", sep = "\t" , row.names=F)
gut_data_mcra <- read.table(file = "gut_data_mcra.txt", sep = "\t", header = T)

mcra_p <- gut_data_mcra %>% filter(Good == "Y") %>% unique()
table(mcra_p$family)
mcra_a <- gut_data_mcra %>% filter(!ID %in% mcra_p$ID) %>% unique()
mcra_data_pcr <- rbind(cbind(mcra_p , "Amplified" = rep("Amplified")),
                      cbind(mcra_a , "Amplified" = rep("No amplified")))

fam_success <- as.data.frame(table(mcra_data_pcr$family, mcra_data_pcr$Amplified))
mean(table(mcra_data_pcr$family, mcra_data_pcr$Amplified)[,1]/(table(mcra_data_pcr$family, mcra_data_pcr$Amplified)[,1]+table(mcra_data_pcr$family, mcra_data_pcr$Amplified)[,2]))

plot1 <- ggplot(fam_success, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Relative Success of mcrA gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

fam_success2 <- fam_success
fam_success2[which(fam_success2$Var2 == "No amplified"),]$Freq <- -fam_success2[which(fam_success2$Var2 == "No amplified"),]$Freq 
plot1_B <- ggplot(fam_success2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="identity")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Success of mcrA gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

pdf(file ="./figure/plot_mcrA_success_fam.pdf", he = 7, wi = 7)
plot1
plot1_B
dev.off()

diet_success <- as.data.frame(table(mcra_data_pcr$diet3, mcra_data_pcr$Amplified))
plot2 <- ggplot(diet_success, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Relative Success of mcrA gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 0,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

diet_success2 <- diet_success
diet_success2[which(diet_success2$Var2 == "No amplified"),]$Freq <- -diet_success2[which(diet_success2$Var2 == "No amplified"),]$Freq 
plot2_B <- ggplot(diet_success2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat="identity", position="identity")+
  scale_fill_manual(values = c("darkred", "gray"))+
  ylab("Success of mcrA gene amplification")+
  theme_bw()+
  labs(fill = "Success")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "serif",size = 15),
        axis.title = element_text(family = "serif",size = 15),
        axis.text = element_text(family = "serif",size = 14), 
        axis.text.x = element_text(family = "serif",size = 12, angle = 45,hjust= 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(family = "serif", size = 15),
        legend.text = element_text(family = "serif", size = 12))

pdf(file ="./figure/plot_mcra_success_diet3.pdf", he = 7, wi = 7)
plot2
plot2_B
dev.off()



