# ordination for rumen samples including plotting NMDS and PCoA, calculating 
  # signigicance with adonis and anosim, multiple regression on matrices (MRM),
  # and Mantel

# input is phyloseq object containing samples with >1000 sequences, tree of 
  # AGF genera, and tree of samples

library(readxl)
library(phyloseq)
library(ape)
library(plyr)
library(vegan)
library(rbiom)
library(cluster)
library(ecodist)
library(ggplot2)
library(ggrepel)
library(phylosmith)

# functions
rescale_dist_mtx = function(m){
  m = m %>% as.matrix
  labs = m %>% colnames
  n_row = m %>% nrow
  n_col = m %>% ncol
  x = m %>% as.vector 
  x = scales::rescale(x) 
  m = matrix(x, nrow=n_row, ncol=n_col)
  colnames(m) = labs
  rownames(m) = labs
  m = m %>% as.dist
  return(m)
}

dist_mtx_order = function(d, x){
  m = d %>% as.matrix
  d = as.dist(m[x,x])
  return(d)
}

## create phyloseq object with the genera tree
otu_mat <-read_excel("~/Desktop/R/rumen_phyloseq2_1000.xlsx", sheet="OTU")
tax_mat<- read_excel("~/Desktop/R/rumen_phyloseq2_1000.xlsx", sheet="taxon")
Meta <-read_excel("~/Desktop/R/rumen_phyloseq2_1000.xlsx", sheet="Samples")
Meta <- Meta %>%
  tibble::column_to_rownames("Sample")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("#OTU ID")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(Meta)
samples
Tree <-ape::read.tree(file="~/Desktop/R/Genera_rooted_ny.nwk")
samples
Physeq <-phyloseq(OTU, TAX, samples, Tree)
Physeq


# create ordination plots
# bray-curtis
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="animal", shape="family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
BrayOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=animal)) +
  geom_point(size=2) + 
  facet_wrap(~method, scales="free") +
  scale_x_continuous(limits = c(-0.2, 0.25)) +
  scale_y_continuous(limits = c(-0.25, 0.2))
BrayOrd
ggsave("all_rumen_dpcoa.jpg", plot=last_plot())

# jaccard
dist="jaccard"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Animal", shape="Family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
JaccardOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=Animal, shape=Family))
JaccardOrd=JaccardOrd+geom_point(size=2)
JaccardOrd=JaccardOrd+facet_wrap(~method, scales="free")
JaccardOrd
ggsave("jaccard_rumen.jpg", plot=last_plot())

# unifrac (unweighted)
dist="unifrac"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Animal", shape="Family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
UniOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=Animal, shape=Family))
UniOrd=UniOrd+geom_point(size=2)
UniOrd=UniOrd+facet_wrap(~method, scales="free")
UniOrd
ggsave("unifrac_rumen.jpg", plot=last_plot())

#unifrac (weighted)
dist="Wunifrac"
ord_meths = c("PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="animal")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
WUniOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=animal)) +
  geom_point(size=2) +
  facet_wrap(~method, scales="free") +
  scale_color_manual(values=c("#4E67C8", "#7BC9EF", "#B5E86A", "#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8", #77A737,
                              "#488870", "#A1501D", "#8F2B19", "#7584CD", "#93D4F2", "#C4EC84", "#93D6C0"),
                     name=c("family"="Animal"))
                     
WUniOrd
ggsave("wunifrac_rumen_PCoA_animal.jpg", plot=last_plot())
samples

#ordination by family
#NMDS and PCoA using weighted unifrac
dist="Wunifrac"
ord_meths = c("NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
WUniOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=family)) +
  geom_point(size=2) +
  facet_wrap(~method, scales="free") +
  scale_color_manual(values=c("#C4EC84", "#7BC9EF", "#4E67C8"),
                     name=c("family"="Animal"))


WUniOrd
ggsave("wunifrac_rumen_NMDS_family.jpg", plot=last_plot())


# Biogeography - Cattle
# cow phyloseq object
otu_mat <-read_excel("~/Desktop/R/cow_1000_phyloseq.xlsx", sheet="OTU")
tax_mat<- read_excel("~/Desktop/R/cow_1000_phyloseq.xlsx", sheet="taxon")
Meta <-read_excel("~/Desktop/R/cow_1000_phyloseq.xlsx", sheet="Samples")
Meta <- Meta %>%
  tibble::column_to_rownames("Sample")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("#OTU ID")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(Meta)
samples
Tree <-ape::read.tree(file="~/Desktop/R/Genera_rooted_ny.nwk")
samples
Physeq <-phyloseq(OTU, TAX, samples, Tree)
Physeq

# calculating PCoA with weighted unifrac
dist="Wunifrac"
ord_meths = c("PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Age")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

# plotting
ggplot(pdataframe, aes(Axis_1, Axis_2, color=Age)) +
  geom_point(size=2) + facet_wrap(~method, scales="free") +
  scale_color_manual(values=c("#7BCBB0", "#EF873D", "#DE5034", "#2C397B", "#418DB8"),
                    labels=c("New_Zealand"="New Zealand", 
                    "United_Kingdom" = "United Kingdom"))
ggsave("wunifrac_cow_age_pcoa.jpg", plot=last_plot())





# reading in host tree, converting it to dist
Host_Sp <-ape::read.tree(file="~/Desktop/R/rumentree.nwk")

host_tree_d = Host_Sp %>% cophenetic %>% as.dist %>% rescale_dist_mtx

X = labels(host_tree_d)

## creating beta diversity matrices then ordering them
## DM is the abundance file
## factors file contains the variables
## Biom is the abundance file but flipped with genera as rows and samples as columns
DM <-read.table("~/Desktop/R/abund_1000.txt", header=TRUE, row.names=1)
Factors <-read_excel("~/Desktop/R/factors_1000.xlsx")
Biom <-read.table("~/Desktop/R/biom_1000.txt", header=TRUE, row.names=1)
Biom <-as.matrix(Biom)

Uni_w <-unifrac(Biom, weighted=TRUE, tree=Tree)

Bray <-vegdist(DM, method="bray")
Bray_o <-dist_mtx_order(Bray, X)
Euc <-vegdist(DM, method="euclidean")
Euc_o <-dist_mtx_order(Euc, X)
Jac <-vegdist(DM, method="jaccard")
Jac_o <-dist_mtx_order(Jac, X)
Unifrac_uw <-rbiom::unifrac(Biom, weighted=FALSE, tree=Tree)
Uni_uw_o <-dist_mtx_order(Unifrac_uw, X)
Unifrac_w <-unifrac(Biom, weighted=TRUE, tree=Tree)
Uni_w_o <-dist_mtx_order(Unifrac_w, X)

## Testing for significance of factors
## Anosim

Anosim1 <-anosim(Bray_o, Factors$Animal)
Anosim2 <-anosim(Euc_o, Factors$Animal)
Anosim3 <-anosim(Jac_o, Factors$Animal)
Anosim4 <-anosim(Uni_uw_o, Factors$Animal)
Anosim5 <-anosim(Uni_w_o, Factors$Animal)

Anosim6 <-anosim(Bray_o, Factors$Family)
Anosim7 <-anosim(Euc_o, Factors$Family)
Anosim8 <-anosim(Jac_o, Factors$Family)
Anosim9 <-anosim(Uni_uw_o, Factors$Family)
Anosim10 <-anosim(Uni_w_o, Factors$Family)

Anosim11 <-anosim(Bray_o, Factors$Lifestyle)
Anosim12 <-anosim(Euc_o, Factors$Lifestyle)
Anosim13 <-anosim(Jac_o, Factors$Lifestyle)
Anosim14 <-anosim(Uni_uw_o, Factors$Lifestyle)
Anosim15 <-anosim(Uni_w_o, Factors$Lifestyle)

Anosim16 <-anosim(Bray_o, Factors$Country)
Anosim17 <-anosim(Euc_o, Factors$Country)
Anosim18 <-anosim(Jac_o, Factors$Country)
Anosim19 <-anosim(Uni_uw_o, Factors$Country)
Anosim20 <-anosim(Uni_w_o, Factors$Country)

Anosim1
Anosim2
Anosim3
Anosim4
Anosim5
Anosim6
Anosim7
Anosim8
Anosim9
Anosim10
Anosim11
Anosim12
Anosim13
Anosim14
Anosim15
Anosim16
Anosim17
Anosim18
Anosim19
Anosim20


##adonis
Adonis1 <-adonis(Bray_o ~Animal, Factors)
Adonis2 <-adonis(Euc_o ~Animal, Factors)
Adonis3 <-adonis(Jac_o ~Animal, Factors)
Adonis4 <-adonis(Uni_uw_o ~Animal, Factors)
Adonis5 <-adonis(Unifrac_w ~animal, Factors)

Adonis6 <-adonis(Bray_o ~Family, Factors)
Adonis7 <-adonis(Euc_o ~Family, Factors)
Adonis8 <-adonis(Jac_o ~Family, Factors)
Adonis9 <-adonis(Uni_uw_o ~Family, Factors)
Adonis10 <-adonis(Uni_w_o ~Family, Factors)

Adonis11 <-adonis(Bray_o ~Lifestyle, Factors)
Adonis12 <-adonis(Euc_o ~Lifestyle, Factors)
Adonis13 <-adonis(Jac_o ~Lifestyle, Factors)
Adonis14 <-adonis(Uni_uw_o ~Lifestyle, Factors)
Adonis15 <-adonis(Uni_w_o ~Lifestyle, Factors)

Adonis16 <-adonis(Bray_o ~Country, Factors)
Adonis17 <-adonis(Euc_o ~Country, Factors)
Adonis18 <-adonis(Jac_o ~Country, Factors)
Adonis19 <-adonis(Uni_uw_o ~Country, Factors)
Adonis20 <-adonis(Uni_w_o ~Country, Factors)

Adonis1$aov
Adonis2$aov
Adonis3$aov
Adonis4$aov
Adonis5$aov
Adonis6$aov
Adonis7$aov
Adonis8$aov
Adonis9$aov
Adonis10$aov
Adonis11$aov
Adonis12$aov
Adonis13$aov
Adonis14$aov
Adonis15$aov
Adonis16$aov
Adonis17$aov
Adonis18$aov
Adonis19$aov
Adonis20$aov

## MRM and mantel
## creating factors distance matrices
family<-as.factor(samples$Family)
daisyfamily<-data.frame(family)
daisyfamily %>% summary
rownames(daisyfamily) = samples$Sample_ID
daisyfamily$Sample_ID= NULL
gower.family<-daisy(daisyfamily, metric="gower")
gower.family %>% summary

country<-as.factor(samples$Country)
daisycountry<-data.frame(country)
daisycountry %>% summary
rownames(daisycountry) = samples$Sample_ID
daisycountry$Sample_ID= NULL
gower.country<-daisy(daisycountry, metric="gower")
gower.country %>% summary

lifestyle<-as.factor(samples$Lifestyle)
daisylife<-data.frame(lifestyle)
daisylife %>% summary
rownames(daisylife) = samples$Sample_ID
daisylife$Sample_ID= NULL
gower.life<-daisy(daisylife, metric="gower")
gower.life %>% summary

animal<-as.factor(samples$Animal)
daisyanimal<-data.frame(animal)
daisyanimal %>% summary
rownames(daisyanimal) = samples$Sample_ID
daisyanimal$Sample_ID= NULL
gower.animal<-daisy(daisyanimal, metric="gower")
gower.animal %>% summary

Country_d <-as.dist(gower.country, diag=TRUE, upper=TRUE)
Lifestyle_d <-as.dist(gower.life, diag=TRUE, upper=TRUE)
Family_d <-as.dist(gower.family, diag=TRUE, upper=TRUE)
Animal_d <-as.dist(gower.animal, diag=TRUE, upper=TRUE)

Lifestyle_d_o <-dist_mtx_order(Lifestyle_d, X)
Country_d_o <-dist_mtx_order(Country_d, X)
Family_d_o <-dist_mtx_order(Family_d, X)
Animal_d_o <-dist_mtx_order(Animal_d, X)


Lifestyle_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Country_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Family_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Animal_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)

Uni_uw_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Uni_w_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Bray_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Euc_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Jac_o %>% lapply(function(x) x %>% as.matrix %>% dim)

MRM_Euc <-ecodist::MRM(Euc_o ~ host_tree_d + Family_d_o + Country_d_o + Lifestyle_d_o)
MRM_Euc
MRM_Jac <-ecodist::MRM(Jac_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Jac
MRM_Bray <-ecodist::MRM(Bray_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Bray
MRM_Uni_uw <-ecodist::MRM(Uni_uw_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Uni_uw
MRM_Uni_w <-ecodist::MRM(Uni_w_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Uni_w
Mantel_Euc <-ecodist::mantel(Euc_o ~ host_tree_d)
Mantel_Euc
Mantel_Euc <-ecodist::mantel(Euc_o ~ Family_d_o)
Mantel_Euc
Mantel_Euc <-ecodist::mantel(Euc_o ~ Country_d_o)
Mantel_Euc
Mantel_Euc <-ecodist::mantel(Euc_o ~ Lifestyle_d_o)
Mantel_Euc
Mantel_Jac <-ecodist::mantel(Jac_o ~ host_tree_d)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Family_d_o)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Country_d_o)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Lifestyle_d_o)
Mantel_Jac
Mantel_Bray <-ecodist::mantel(Bray_o ~ host_tree_d)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Family_d_o)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Country_d_o)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Lifestyle_d_o)
Mantel_Bray
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ host_tree_d)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Family_d_o)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Country_d_o)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Lifestyle_d_o)
Mantel_Uni_uw
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ host_tree_d)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Family_d_o)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Country_d_o)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Lifestyle_d_o)
Mantel_Uni_w
