# alpha diversity box plots by family and species + rank abundance plot
# only samples with >1000 sequences were used for analysis

# all alpha diversity values were calculated using the summary.single command 
  # in mothur (sobs, shannon, and simpson)

# alpha diversity .xlsx inputs are one column of animal identity 
  # (family, age, sex, species, etc) followed by columns of genera observed, 
  # shannon, and simpson


library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(readxl)
library(ggh4x)
setwd("~/Desktop/R/")

# alpha diversity box plots
# observed number of genera by family
# calculating
SOBS<-read_excel("~/Desktop/R/alpha.div.1000.family.xlsx")
SOBdata <- melt(SOBS,id.vars='Family', measure.vars='sobs')
compare_means(value ~ Family,  data = SOBdata, exact=FALSE)
my_comparisons<-list( c("Antilocapridae", "Bovidae"), c("Antilocapridae", "Cervidae"), c("Bovidae", "Cervidae") )

#plotting
SOBbox<-ggplot(SOBdata) + geom_boxplot(aes(x=Family, y=value, fill=Family)) + 
  ggtitle("Observed Number of Genera") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Family, y=value), 
              map_signif_level = TRUE, 
              y_position = c(74, 84, 79))

SOBbox
ggsave("alpha_sobs_1000_family.jpg", plot=last_plot())


# shannon by family
# calculating
SHAN<-read_excel("~/Desktop/R/alpha.div.1000.family.xlsx")
SHANdata <- melt(SHAN,id.vars='Family', measure.vars='shannon')
compare_means(value ~ Family,  data = SHANdata)
my_comparisons<-list( c("Antilocapridae", "Bovidae"), c("Antilocapridae", "Cervidae"), c("Bovidae", "Cervidae") )

#plotting
SHANbox<-ggplot(SHANdata) + geom_boxplot(aes(x=Family, y=value, fill=Family)) + 
  ggtitle("Shannon") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0"))+
  geom_signif(comparisons=my_comparisons, 
              aes(x=Family, y=value), 
              map_signif_level = TRUE,
              y_position = c(3, 3.5, 3.25)) 
SHANbox
ggsave("alpha_shannon_1000_family.jpg", plot=last_plot())


# simpson by family
# calculating
SIMP<-read_excel("~/Desktop/R/alpha.div.1000.family.xlsx")
SIMPdata <- melt(SIMP,id.vars='Family', measure.vars='simpson')
compare_means(value ~ Family,  data = SIMPdata)
my_comparisons<-list( c("Antilocapridae", "Bovidae"), c("Antilocapridae", "Cervidae"), c("Bovidae", "Cervidae") )

#plotting
SIMPbox<-ggplot(SIMPdata) + geom_boxplot(aes(x=Family, y=value, fill=Family)) + 
  ggtitle("Simpson") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0"))+
  geom_signif(comparisons=my_comparisons, 
              aes(x=Family, y=value), 
              map_signif_level = TRUE,
              y_position = c(1.1, 1.25, 1.17)) 
SIMPbox
ggsave("alpha_simpson_1000_family.jpg", plot=last_plot())


# observed genera by species
# calculating
SOBSS<-read_excel("~/Desktop/R/alpha.div.1000.species.xlsx")
SOBSdata <- melt(SOBSS,id.vars='Species', measure.vars='sobs')
compare_means(value ~ Species,  data = SOBSdata)
my_comparisons<-list( c("Cow", "Goat"), c("Cow", "Sheep"), c("Goat", "Sheep") )

# plotting
SOBSbox<-ggplot(SOBSdata) + geom_boxplot(aes(x=Species, y=value, fill=Species)) + 
  ggtitle("Observed Number of Genera") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2C397B", "#418DB8", "#77A737")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Species, y=value), 
              map_signif_level = TRUE, 
              y_position = c(62, 72, 67))
SOBSbox
ggsave("alpha_sobs_1000_species.jpg", plot=last_plot())


# shannon by species
# calculating
SHANS<-read_excel("~/Desktop/R/alpha.div.1000.species.xlsx")
SHANSdata <- melt(SHANS,id.vars='Species', measure.vars='shannon')
compare_means(value ~ Species,  data = SHANSdata)
my_comparisons<-list( c("Cow", "Goat"), c("Cow", "Sheep"), c("Goat", "Sheep") )

# plotting
SHANSbox<-ggplot(SHANSdata) + geom_boxplot(aes(x=Species, y=value, fill=Species)) + 
  ggtitle("Shannon") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2C397B", "#418DB8", "#77A737")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Species, y=value), 
              map_signif_level = TRUE,
              y_position = c(3, 3.5, 3.25)) 
SHANSbox
ggsave("alpha_shannon_1000_species.jpg", plot=last_plot())


# simpson by species
# calculating
SIMPS<-read_excel("~/Desktop/R/alpha.div.1000.species.xlsx")
SIMPSdata <- melt(SIMPS,id.vars='Species', measure.vars='simpson')
compare_means(value ~ Species,  data = SIMPSdata)
my_comparisons<-list( c("Cow", "Goat"), c("Cow", "Sheep"), c("Goat", "Sheep") )

# plotting
SIMPSbox<-ggplot(SIMPSdata) + geom_boxplot(aes(x=Species, y=value, fill=Species)) + 
  ggtitle("Simpson") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#2C397B", "#418DB8", "#77A737")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Species, y=value), 
              map_signif_level = TRUE,
              y_position = c(1, 1.15, 1.075)) 
SIMPSbox
ggsave("alpha_simpson_1000_species.jpg", plot=last_plot())


# number of genera observed by lifestyle (wild vs domesticated)
# calculating
SOBL<-read_excel("~/Desktop/R/alpha.div.1000.lifestyle.xlsx")
SOBLdata <- melt(SOBL,id.vars='Lifestyle', measure.vars='sobs')
compare_means(value ~ Lifestyle,  data = SOBLdata)
my_comparisons<-list( c("Wild", "Domesticated"))

# plotting
SOBLbox<-ggplot(SOBLdata) + geom_boxplot(aes(x=Lifestyle, y=value, fill=Lifestyle)) + 
  ggtitle("Observed Number of Genera") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#EF873D", "#DE5034"))+
  geom_signif(comparisons=my_comparisons, 
              aes(x=Lifestyle, y=value), 
              map_signif_level = TRUE)
SOBLbox
ggsave("alpha_sobs_1000_lifestyle.jpg", plot=last_plot())


# shannon by lifestyle
# calculating
SHANl<-read_excel("~/Desktop/R/alpha.div.1000.lifestyle.xlsx")
SHANldata <- melt(SHANl,id.vars='Lifestyle', measure.vars='shannon')
compare_means(value ~ Lifestyle,  data = SHANldata)
my_comparisons<-list( c("Wild", "Domesticated"))

# plotting
SHANlbox<-ggplot(SHANldata) + geom_boxplot(aes(x=Lifestyle, y=value, fill=Lifestyle)) + 
  ggtitle("Shannon") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#EF873D", "#DE5034"))+
  geom_signif(comparisons=my_comparisons, 
              aes(x=Lifestyle, y=value), 
              map_signif_level = TRUE) 
SHANlbox
ggsave("alpha_shannon_1000_lifestle.jpg", plot=last_plot())


# simpson by lifestyle
# calculaing
SIMPL<-read_excel("~/Desktop/R/alpha.div.1000.lifestyle.xlsx")
SIMPLdata <- melt(SIMPL,id.vars='Lifestyle', measure.vars='simpson')
compare_means(value ~ Lifestyle,  data = SIMPLdata)
my_comparisons<-list( c("Wild", "Domesticated"))

# plotting
SIMPLbox<-ggplot(SIMPLdata) + geom_boxplot(aes(x=Lifestyle, y=value, fill=Lifestyle)) + 
  ggtitle("Simpson") + ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#EF873D", "#DE5034")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Lifestyle, y=value), 
              map_signif_level = TRUE) 
SIMPLbox
ggsave("alpha_simpson_1000_LIFESTYLE.jpg", plot=last_plot())





# biogeography, sex and age comparisons for cattle alpha diversity
# import cattle alpha diversity calculations
bio<-read_excel("~/Desktop/R/alpha.div.1000.cattle.xlsx")
bio

# only visualize significant comparisons on box plots
sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

# number of observed genera by country
# calculating
melted <- melt(bio,id.vars="Country", measure.vars="sobs")
melted
means<-compare_means(value ~ Country, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_country_compare_sobs.xlsx")
my_comparisons<-list( c("China", "Denmark"),  c("China", "Netherlands"), c("China", "New_Zealand"), c("China", "United_Kingdom"), c("Denmark", "Brazil"), c("Netherlands", "Brazil"), c("China", "Mexico"),  c("China", "Switzerland"),   
c("China", "Brazil"),  c("Denmark", "Mexico"),  c("Denmark", "Netherlands"),  c("Denmark", "New_Zealand"),  c("Denmark", "Switzerland"),  c("Denmark", "United_Kingdom"),   
c("Mexico", "Netherlands"),  c("Mexico", "New_Zealand"),  c("Mexico", "Switzerland"),  c("Mexico", "United_Kingdom"),  c("Mexico", "Brazil"),  c("Netherlands", "New_Zealand"),  c("Netherlands", "Switzerland"), 
c("Netherlands", "United_Kingdom"),   c("New_Zealand", "Switzerland"),  c("New_Zealand", "United_Kingdom"),  c("New_Zealand", "Brazil"),  c("Switzerland", "United_Kingdom"), 
c("Switzerland", "Brazil"),  c("United_Kingdom", "Brazil"))

# plotting
allbox<-ggplot(melted, aes(x=Country, y=value, fill=Country)) + 
  ggtitle("Observed Number of Genera") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons,
              aes(x=Country, y=value),
              map_signif_level = sigFunc,
              y_position = c(61, 65, 69, 73, 77, 81))
allbox
ggsave("alpha_cow_biogeog_sobs.jpg", plot=last_plot())


# shannon by country
# calculating
melted <- melt(bio,id.vars="Country", measure.vars="shannon")
melted
means<-compare_means(value ~ Country, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_country_compare_shannon.xlsx")
my_comparisons<-list( c("China", "Denmark"),  c("China", "Netherlands"), c("China", "New_Zealand"), c("China", "United_Kingdom"), c("China", "Brazil"), c("Denmark", "New_Zealand"),  c("Netherlands", "New_Zealand"),
                      c("Netherlands", "Brazil"), c("China", "Mexico"),  c("China", "Switzerland"),   c("Denmark", "Mexico"),  c("Denmark", "Netherlands"),  c("Denmark", "New_Zealand"),  c("Denmark", "Switzerland"),
                      c("Denmark", "United_Kingdom"),  c("Denmark", "Brazil"),  c("Mexico", "Netherlands"),  c("Mexico", "New_Zealand"),  c("Mexico", "Switzerland"),  c("Mexico", "United_Kingdom"),  c("Mexico", "Brazil"),  
                      c("Netherlands", "Switzerland"), c("Netherlands", "United_Kingdom"),   c("New_Zealand", "Switzerland"),  c("New_Zealand", "United_Kingdom"),  c("New_Zealand", "Brazil"),  c("Switzerland", "United_Kingdom"), 
                      c("Switzerland", "Brazil"),  c("United_Kingdom", "Brazil"))

# plotting
allbox<-ggplot(melted, aes(x=Country, y=value, fill=Country)) + 
  ggtitle("Shannon") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons,
              aes(x=Country, y=value),
              map_signif_level = sigFunc,
              y_position = c(2.8, 3, 3.2, 3.4, 3.6, 3.6, 2.8))
allbox
ggsave("alpha_cow_biogeog_shannon.jpg", plot=last_plot())


# simpson by country
# calculating
melted <- melt(bio,id.vars="Country", measure.vars="simpson")
melted
means<-compare_means(value ~ Country, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_country_compare_simpson.xlsx")
my_comparisons<-list( c("China", "Denmark"),  c("China", "Netherlands"),  c("China", "United_Kingdom"), c("China", "Brazil"), c("Denmark", "New_Zealand"),c("Denmark", "Switzerland"),  c("Netherlands", "New_Zealand"),  
                      c("China", "New_Zealand"), c("Denmark", "Mexico"),  c("Denmark", "Netherlands"),    c("Denmark", "Brazil"),  c("Denmark", "United_Kingdom"), c("Netherlands", "Brazil"), c("China", "Mexico"),  c("China", "Switzerland"),    
                      c("Mexico", "Netherlands"),  c("Mexico", "New_Zealand"),  c("Mexico", "Switzerland"),  c("Mexico", "United_Kingdom"),  c("Mexico", "Brazil"),   c("Netherlands", "Switzerland"), 
                      c("Netherlands", "United_Kingdom"),   c("New_Zealand", "Switzerland"),  c("New_Zealand", "United_Kingdom"),  c("New_Zealand", "Brazil"),  c("Switzerland", "United_Kingdom"), 
                      c("Switzerland", "Brazil"),  c("United_Kingdom", "Brazil"))

# plotting
allbox<-ggplot(melted, aes(x=Country, y=value, fill=Country)) + 
  ggtitle("Simpson") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons,
              aes(x=Country, y=value),
              map_signif_level = sigFunc,
              y_position = c(0.95, 1.05, 1.15, 0.85, 1.25, 1.35, 0.95))
allbox
ggsave("alpha_cow_biogeog_simpson.jpg", plot=last_plot())


# number of observed genera by sex
# calculating
melted <- melt(bio,id.vars="Sex", measure.vars="sobs")
melted
means<-compare_means(value ~ Sex, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_sex_compare_sobs.xlsx")
my_comparisons<-list( c("Male", "Female"))

# plotting
allbox<-ggplot(melted, aes(x=Sex, y=value, fill=Sex)) + 
  ggtitle("Number of Genera Observed") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Sex, y=value), 
              map_signif_level = TRUE) 
allbox
ggsave("alpha_cow_sex_sobs.jpg", plot=last_plot())


# shannon by sex
# calculating
melted <- melt(bio,id.vars="Sex", measure.vars="shannon")
melted
means<-compare_means(value ~ Sex, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_sex_compare_shannon.xlsx")
my_comparisons<-list( c("Male", "Female"))

# plotting
allbox<-ggplot(melted, aes(x=Sex, y=value, fill=Sex)) + 
  ggtitle("Shannon") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Sex, y=value), 
              map_signif_level = TRUE) 
allbox
ggsave("alpha_cow_sex_shannon.jpg", plot=last_plot())


# simpson by sex
# calculating
melted <- melt(bio,id.vars="Sex", measure.vars="simpson")
melted
means<-compare_means(value ~ Sex, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_sex_compare_simpson.xlsx")
my_comparisons<-list( c("Male", "Female"))

# plotting
allbox<-ggplot(melted, aes(x=Sex, y=value, fill=Sex)) + 
  ggtitle("Simpson") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Sex, y=value), 
              map_signif_level = TRUE) 
allbox
ggsave("alpha_cow_sex_simpson.jpg", plot=last_plot())

# number of observed genera by age
# calculating
melted <- melt(bio,id.vars="Age", measure.vars="sobs")
melted
means<-compare_means(value ~ Age, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_age_compare_sobs.xlsx")
my_comparisons<-list( c("Mature", "Old"))

# plotting
allbox<-ggplot(melted, aes(x=Age, y=value, fill=Age)) + 
  ggtitle("Number of Genera Observed") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c( "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Age, y=value), 
              map_signif_level = TRUE) 
allbox
ggsave("alpha_cow_age_sobs.jpg", plot=last_plot())


# shannon by age
# calculating
melted <- melt(bio,id.vars="Age", measure.vars="shannon")
melted
means<-compare_means(value ~ Age, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_age_compare_shannon.xlsx")
my_comparisons<-list( c("Mature", "Old"))

# plotting
allbox<-ggplot(melted, aes(x=Age, y=value, fill=Age)) + 
  ggtitle("Shannon") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c( "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Age, y=value), 
              map_signif_level = TRUE) 

allbox
ggsave("alpha_cow_age_shannon.jpg", plot=last_plot())


# simpson by age
# calculating
melted <- melt(bio,id.vars="Age", measure.vars="simpson")
melted
means<-compare_means(value ~ Age, data = melted)
means
write_xlsx(means, path="~/Desktop/R/cow_age_compare_simpson.xlsx")
my_comparisons<-list( c("Mature", "Old"))

# plotting
allbox<-ggplot(melted, aes(x=Age, y=value, fill=Age)) + 
  ggtitle("Simpson") + 
  ylab(NULL) + xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_discrete(name="") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_blank()) +
  scale_fill_manual(values=c("#B5E86A", "#7BCBB0")) +
  geom_signif(comparisons=my_comparisons, 
              aes(x=Age, y=value), 
              map_signif_level = TRUE) 
allbox
ggsave("alpha_cow_age_simpson.jpg", plot=last_plot())



# rank abundance plot
DM <-read.table("~/Desktop/R/GeneraTable.txt", header=TRUE)
library("BiodiversityR")
Rank_abundance <-rankabundance(DM)
rankabunplot(Rank_abundance, scale='abundance')
