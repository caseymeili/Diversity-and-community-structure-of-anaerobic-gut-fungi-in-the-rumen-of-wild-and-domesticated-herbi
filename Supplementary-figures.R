# assorted supplementary rumen figures 

setwd("~/Desktop/R/")
library(ggplot2)
library(vegan)
library(ape)
library(reshape2)
library(readxl)
library(ggpubr)

# plots contain AGF genera above 1% abundance in at least 50% of one animal species
# inputs contain host identity (family, species etc) in column one and a column 
  # for each AGF genus and abundance data per sample

# box plots of all sequences
all<-read_excel("~/Desktop/R/abund_1000.xlsx")
alldata <- melt(all,id.vars='Group', measure.vars=c("Neocallimastix", "Caecomyces", "Orpinomyces", "Cyllamyces", "Piromyces", "Anaeromyces", "Capellomyces", "Feramyces", "Liebetanzomyces", "Khyollomyces", "AL3", "SK3", "RH4", "NY42", "NY9", "NY8", "NY5", "NY11", "NY4", "NY6"))
allseq<-ggplot(alldata, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.size = 0.2, fill="#7BC9EF") + 
  xlab(NULL) + ylab("Percent Abundance") +
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(angle=90,hjust=1))
allseq

ggsave("all_seq_box.jpg", plot=last_plot(), height=1280, units="px")


# grouped box plot by host animal family
all<-read_excel("~/Desktop/R/allfamily.xlsx")
melted <- melt(all,id.vars='Family', measure.vars=c("Neocallimastix", "Caecomyces", "Orpinomyces", "Cyllamyces", "Piromyces", "Anaeromyces", "Capellomyces", "Feramyces", "Liebetanzomyces", "Khyollomyces", "AL3", "SK3", "RH4", "NY42", "NY9", "NY8", "NY5", "NY11", "NY4", "NY6"))

allbox<-ggplot(melted, aes(x=variable, y=value, fill=Family)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_y_continuous(trans='log10') +
  xlab(NULL) + ylab("Percent Abundance") +
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0")) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) 
allbox

ggsave("all_family_box.jpg", plot=last_plot(), height=1280, units="px")


# grouped box plot by host animal species
all<-read_excel("~/Desktop/R/allspecies.xlsx")
melted <- melt(all,id.vars='Species', measure.vars=c("Neocallimastix", "Caecomyces", "Orpinomyces", "Cyllamyces", "Piromyces", "Anaeromyces", "Capellomyces", "Feramyces", "Liebetanzomyces", "Khyollomyces", "AL3", "SK3", "RH4", "NY42", "NY9", "NY8", "NY5", "NY11", "NY4", "NY6"))

allbox<-ggplot(melted, aes(x=variable, y=value, fill=Species)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_y_continuous(trans='log10') +
  xlab(NULL) + ylab("Percent Abundance") +
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2C397B", "#418DB8", "#77A737", "#488870")) +
  theme(axis.text.x=element_text(angle=90,hjust=1))
allbox

ggsave("all_species_box.jpg", plot=last_plot(), height=1280, units="px")
