# global phylogenetic signal statistics and local indicator of phylogenetic 
  # association (LIPA) calculations

library(phylosignal)
library(ape)
library(phylobase)
library(phytools)
library(writexl)

## read the percentage abundance file. This will have samples as rows, genera as columns 
## and percent abundance as cell values
PercAbund <-read.table("~/Desktop/R/rumenabund.txt", header=TRUE)
PercAbund$Group = as.numeric(as.factor(PercAbund$Group))

## read host tree
Host_Sp <-ape::read.tree(file="~/Desktop/R/rumentree.nwk")


## read the genera tree
Taxa <-ape::read.tree(file="~/Desktop/R/Genera_rooted_ny.nwk")

## create a phylo4d object
p4dAll <-phylo4d(Host_Sp, data.frame(PercAbund))
head(PercAbund)
## calculate the Lipa index for each sample-genus pair
LipaAll <-lipaMoran(p4dAll)
LipaAll

## read the output from Lipa and copy it to excel
options(max.print = .Machine$integer.max)
LIPA<-LipaAll$lipa
LIPA<-data.frame(LIPA)
write_xlsx(LIPA,"~/Desktop/R/LIPA_rumen.xlsx")
LipaP<-LipaAll$p.value 
LIPA<-data.frame(LipaP)

## calculate other indices of phylosignal
PS <-phylosignal::phyloSignal(p4dAll, methods="all")
statdf<-data.frame(PS$stat)
write_xlsx(statdf,"~/Desktop/R/phylosignal_indices_stat.xlsx")
pvalue<-data.frame(PS$pvalue)
write_xlsx(pvalue,"~/Desktop/R/phylosignal_indices_pvalue.xlsx")


## create two trees facing one another with association in between
## the association file is two columns with animal in first column and genus in second
## Decisions on which associations to include depend on the lipa value. I used >1 as strong
## read association file
Assoc <-read.table("~/Desktop/R/strongassoc.txt", header=TRUE)
cophylo=cophylo(Host_Sp, Taxa, assoc=Assoc)
Assoc
write_xlsx(LIPA,"~/Desktop/R/LIPA_strong_assoc.xlsx")




