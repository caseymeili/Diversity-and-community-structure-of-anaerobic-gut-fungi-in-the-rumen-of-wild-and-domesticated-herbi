# ordination for rumen versus feces comparison

# PCoA and NMDS using weighted unifrac and box plots of centroid distances
# calves with rumen and feces sampled on the same day and overlapping 
  # animals from global fecal mycobiome survey were used

# inputs are phyloseq objects: one for calves one for other animals

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
library(ggConvexHull)

## create phyloseq object for 8 temporal calves
otu_mat <-read_excel("~/Desktop/R/calf_phyloseq.xlsx", sheet="OTU")
tax_mat<- read_excel("~/Desktop/R/calf_phyloseq.xlsx", sheet="taxon")
Meta <-read_excel("~/Desktop/R/calf_phyloseq.xlsx", sheet="Samples")
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
Physeq <-phyloseq(OTU, TAX, samples, Tree)
Physeq


# calculate PCoA
dist = "Wunifrac"
ord_meths = c("PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

# calculate centroids
centroids <- aggregate(cbind(Axis_1, Axis_2)~sample_type,pdataframe,mean)
centroids

#plotting PCoA
ggplot(pdataframe, aes(Axis_1, Axis_2, color=sample_type)) +
  geom_point(size=2) +
  stat_ellipse(aes(x=Axis_1, y=Axis_2,group=sample_type),type = "norm") +
  facet_wrap(~method, scales="free") +
  labs(shape = "Sample Type") +
  scale_color_manual((name="Sample Type"), values=c("#4E67C8", "#7BC9EF")) + 
  geom_point(data=centroids, size=3, aes(shape=sample_type), color="black") 

ggsave("calves-with-centroids.jpg", plot=last_plot())

# use beta disper to calculate centroid distances
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}
dist_obj <- vegdist(veganotu(Physeq))
bdisper<-betadisper(dist_obj, samples$sample_type, type="centroid")

# centroid distance box plot
box<-boxplot(bdisper,  xlab="Sample Type", ylab="Distance to Centroid")


## create phyloseq object with overlapping hosts rumen/feces
otu_mat <-read_excel("~/Desktop/R/combined-rumen-and-feces-phyloseq.xlsx", sheet="OTU")
tax_mat<- read_excel("~/Desktop/R/combined-rumen-and-feces-phyloseq.xlsx", sheet="taxon")
Meta <-read_excel("~/Desktop/R/combined-rumen-and-feces-phyloseq.xlsx", sheet="Samples")
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

# calculate PCoA
dist="Wunifrac"
ord_meths = c("PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="common_name", shape="sample_type")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

# calculate centroids by sample type
centroids <- aggregate(cbind(Axis_1, Axis_2)~sample_type,pdataframe,mean)
centroids

# plotting PCoA by sample type
ggplot(pdataframe, aes(Axis_1, Axis_2)) +
  geom_point(aes(color=common_name, shape=sample_type), size=2) +
  stat_ellipse(aes(x=Axis_1, y=Axis_2,group=sample_type),type = "norm") +
  facet_wrap(~method, scales="free") +
  labs(shape = "Sample Type") +
  scale_color_manual((name="Common Name"), values=c("#4E67C8", "#7BC9EF", "#7BCBB0", 
                              "#EF873D", "#DE5034",  "#B5E86A")) + 
  geom_point(data=centroids, size=3, aes(shape=sample_type)) 

ggsave("combined-rumen-and-feces-with-centroids2.jpg", plot=last_plot())
  
# calculate centroids by animal
centroids <- aggregate(cbind(Axis_1, Axis_2)~common_name,pdataframe,mean)
centroids

# plotting PCoA by animal
ggplot(pdataframe, aes(Axis_1, Axis_2)) +
  geom_point(aes(color=common_name, shape=sample_type), size=2) +
  stat_ellipse(aes(x=Axis_1, y=Axis_2,group=common_name, color=common_name),type = "norm") +
  facet_wrap(~method, scales="free") +
  labs(shape = "Sample Type") +
  scale_color_manual((name="Common Name"), values=c("#4E67C8", "#7BC9EF", "#7BCBB0", 
                                                    "#EF873D", "#DE5034",  "#B5E86A")) + 
  geom_point(data=centroids, size=3, aes(color=common_name)) 

ggsave("combined-rumen-and-feces-with-centroids-by animal.jpg", plot=last_plot())

# beta disper centroid distances and boxplot
dist_obj <- vegdist(veganotu(Physeq))
# by sample type
bdisper<-betadisper(dist_obj, samples$sample_type, type="centroid")
boxplot(bdisper,  xlab="Sample Type", ylab="Distance to Centroid")
# by host animal
bdisper<-betadisper(dist_obj, samples$common_name, type="centroid")
boxplot(bdisper,  xlab="Animal", ylab="Distance to Centroid")
samples

