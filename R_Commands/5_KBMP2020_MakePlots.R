# KBMP2020_MakePlots
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Fall 2018

########################################################################################################################################################
######################## FIGURE: Abundance Correlations in Root Networks

ppi <- 300
png("~/Documents/KB_F18_FilesForUpload/Figures/Roots_AbundanceCorrelations.png", width=10*ppi, height=10*ppi, res=ppi)
par(mfrow=c(2,3),mar=c(2,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RootsVegetativeME,layout=SE_LAYOUT_RootsVegetativeME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsVegetativeME)$weight,vertex.size=2.0)
title("ME Vegetative Roots \n(60 ASVs, 55 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringME,layout=SE_LAYOUT_RootsFloweringME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsFloweringME)$weight,vertex.size=2.0)
title("ME Flowering Roots \n(180 ASVs, 261 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentME,layout=SE_LAYOUT_RootsSenescentME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsSenescentME)$weight,vertex.size=2.0)
title("ME Senescent Roots \n(257 ASVs, 696 interactions)",cex.main=2)

plot(SE_GRAPH_RootsVegetativeWW,layout=SE_LAYOUT_RootsVegetativeWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsVegetativeWW)$weight,vertex.size=2.0)
title("WW Vegetative Roots \n(29 ASVs, 20 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringWW,layout=SE_LAYOUT_RootsFloweringWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsFloweringWW)$weight,vertex.size=2.0)
title("WW Flowering Roots \n(148 ASVs, 157 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentWW,layout=SE_LAYOUT_RootsSenescentWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsSenescentWW)$weight,vertex.size=2.0)
title("WW Senescent Roots \n(286 ASVs, 490 interactions)",cex.main=2)

dev.off()

#########################################################################################################################################################

#########################################################################################################################################################
######################## FIGURE: Hub Overlap in Root Networks

# Make a list of ASVs that meet the hub criteria in more than one root network
Root_MasterHubList <- list()
Root_NetworkList <- ls(pat="SE_GRAPH_Roots")
for(n in Root_NetworkList){
  print(n)
  Network = get(n)
  HubList <- V(Network)[V(Network)$isHub=="hub"]$name
  Root_MasterHubList[[n]] <- HubList
}
Root_AllHubs<-unique(unlist(Root_MasterHubList))
Root_HubFrequencies <- table(unlist(Root_MasterHubList))
Root_RecurrentHubs <- names(Root_HubFrequencies[Root_HubFrequencies>1])

# Only one ASV is classified as a hub (90th percentile of degree and betweenness centrality) in multiple root networks:
tax_table(physeq)[rownames(tax_table(physeq)) %in% Root_RecurrentHubs,]
# 40c4afc92ffcaed17f6008cb09daa1df "Bacteria" "Actinobacteria" "Actinobacteria" "Frankiales" "Geodermatophilaceae" "Geodermatophilus" "Ambiguous_taxa"

ppi <- 300
png("~/Documents/KB_F18_FilesForUpload/Figures/Roots_RecurrentHubs.png", width=10*ppi, height=10*ppi, res=ppi)
par(mfrow=c(2,3),mar=c(2,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RootsVegetativeME,layout=SE_LAYOUT_RootsVegetativeME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsVegetativeME)$isHub=="hub" & V(SE_GRAPH_RootsVegetativeME)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsVegetativeME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsVegetativeME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsVegetativeME)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsVegetativeME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsVegetativeME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsVegetativeME)$isHub=="hub","black","grey"),vertex.label.dist=2)
title("ME Vegetative Roots \n(60 ASVs, 55 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringME,layout=SE_LAYOUT_RootsFloweringME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsFloweringME)$isHub=="hub" & V(SE_GRAPH_RootsFloweringME)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsFloweringME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsFloweringME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsFloweringME)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsFloweringME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsFloweringME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsFloweringME)$isHub=="hub","black","grey"),vertex.label.dist=2)
title("ME Flowering Roots \n(180 ASVs, 261 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentME,layout=SE_LAYOUT_RootsSenescentME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsSenescentME)$isHub=="hub" & V(SE_GRAPH_RootsSenescentME)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsSenescentME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsSenescentME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsSenescentME)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsSenescentME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsSenescentME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsSenescentME)$isHub=="hub","black","grey"),vertex.label.dist=2)
title("ME Senescent Roots \n(257 ASVs, 696 interactions)",cex.main=2)

plot(SE_GRAPH_RootsVegetativeWW,layout=SE_LAYOUT_RootsVegetativeWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub" & V(SE_GRAPH_RootsVegetativeWW)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsVegetativeWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsVegetativeWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsVegetativeWW)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub","black","grey"),vertex.label.dist=2)
title("WW Vegetative Roots \n(29 ASVs, 20 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringWW,layout=SE_LAYOUT_RootsFloweringWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsFloweringWW)$isHub=="hub" & V(SE_GRAPH_RootsFloweringWW)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsFloweringWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsFloweringWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsFloweringWW)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsFloweringWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsFloweringWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsFloweringWW)$isHub=="hub","black","grey"),vertex.label.dist=2)
title("WW Flowering Roots \n(148 ASVs, 157 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentWW,layout=SE_LAYOUT_RootsSenescentWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsSenescentWW)$isHub=="hub" & V(SE_GRAPH_RootsSenescentWW)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsSenescentWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsSenescentWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsSenescentWW)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsSenescentWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsSenescentWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsSenescentWW)$isHub=="hub","black","grey"),vertex.label.dist=2)
title("WW Senescent Roots \n(286 ASVs, 490 interactions)",cex.main=2)

dev.off()

#########################################################################################################################################################

#########################################################################################################################################################
######################## FIGURE: ASV Taxonomy in Root Networks

ppi <- 300
png("~/Documents/KB_F18_FilesForUpload/Figures/Roots_HubTaxonomy.png", width=10*ppi, height=10*ppi, res=ppi)
par(mfrow=c(3,3),mar=c(2,0,4,0),family="serif",font=1)
layout(matrix(1:9,ncol=3,nrow=3,byrow=TRUE),heights=c(3,3,1),widths=c(3,3,3))

plot(SE_GRAPH_RootsVegetativeME,layout=SE_LAYOUT_RootsVegetativeME,weighted=FALSE,vertex.label = NA,
     vertex.label.font=2,vertex.label.color="black", vertex.label.cex=0.9,
     edge.color="black",vertex.frame.color = "black")
title("ME Vegetative Roots \n(60 ASVs, 55 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringME,layout=SE_LAYOUT_RootsFloweringME,weighted=FALSE,vertex.label = NA,
     vertex.label.font=2,vertex.label.color="black", vertex.label.cex=0.9,
     edge.color="black",vertex.frame.color = "black")
title("ME Flowering Roots \n(180 ASVs, 261 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentME,layout=SE_LAYOUT_RootsSenescentME,weighted=FALSE,vertex.label = NA,
     vertex.label.font=2,vertex.label.color="black", vertex.label.cex=0.9,
     edge.color="black",vertex.frame.color = "black")
title("ME Senescent Roots \n(257 ASVs, 696 interactions)",cex.main=2)

plot(SE_GRAPH_RootsVegetativeWW,layout=SE_LAYOUT_RootsVegetativeWW,weighted=FALSE,vertex.label = NA,
     vertex.label.font=2,vertex.label.color="black", vertex.label.cex=0.9,
     edge.color="black",vertex.frame.color = "black")
title("WW Vegetative Roots \n(29 ASVs, 20 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringWW,layout=SE_LAYOUT_RootsFloweringWW,weighted=FALSE,vertex.label = NA,
     vertex.label.font=2,vertex.label.color="black", vertex.label.cex=0.9,
     edge.color="black",vertex.frame.color = "black")
title("WW Flowering Roots \n(148 ASVs, 157 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentWW,layout=SE_LAYOUT_RootsSenescentWW,weighted=FALSE,vertex.label = NA,
     vertex.label.font=2,vertex.label.color="black", vertex.label.cex=0.9,
     edge.color="black",vertex.frame.color = "black")
#legend(-0.9,-1.1,ncol=3,MyLegend[,"Categories"], fill=MyLegend[,"Palette"], col="black",cex=0.75,text.font=2, box.lwd=1.5)
title("WW Senescent Roots \n(286 ASVs, 490 interactions)",cex.main=2)

plot(0,type='n',axes=FALSE)
plot(0,type='n',axes=FALSE)
legend(0.57,2.85,ncol=2,Legend[,"Categories"], fill=Legend[,"Palette"], col="black",text.font=2, box.lwd=1.5,cex=1.2)
plot(0,type='n',axes=FALSE)

dev.off()

#########################################################################################################################################################

#########################################################################################################################################################
######################## FIGURE 1: Combines correlation data and hub overlap data for root networks

ppi <- 300
png("~/Documents/KB_F18_FilesForUpload/Figures/Figure1.png", width=10*ppi, height=10*ppi, res=ppi)
par(mfrow=c(4,3),mar=c(0,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RootsVegetativeME,layout=SE_LAYOUT_RootsVegetativeME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsVegetativeME)$weight,vertex.size=2.0)
title("ME Vegetative Roots \n(60 ASVs, 55 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringME,layout=SE_LAYOUT_RootsFloweringME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsFloweringME)$weight,vertex.size=2.0)
title("ME Flowering Roots \n(180 ASVs, 261 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentME,layout=SE_LAYOUT_RootsSenescentME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsSenescentME)$weight,vertex.size=2.0)
title("ME Senescent Roots \n(257 ASVs, 696 interactions)",cex.main=2)

par(mar=c(2,0,1,0),family="serif",font=1)

plot(SE_GRAPH_RootsVegetativeME,layout=SE_LAYOUT_RootsVegetativeME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsVegetativeME)$isHub=="hub" & V(SE_GRAPH_RootsVegetativeME)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_RootsVegetativeME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsVegetativeME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsVegetativeME)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_RootsVegetativeME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsVegetativeME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsVegetativeME)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_RootsFloweringME,layout=SE_LAYOUT_RootsFloweringME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsFloweringME)$isHub=="hub" & V(SE_GRAPH_RootsFloweringME)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_RootsFloweringME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsFloweringME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsFloweringME)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_RootsFloweringME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsFloweringME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsFloweringME)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_RootsSenescentME,layout=SE_LAYOUT_RootsSenescentME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsSenescentME)$isHub=="hub" & V(SE_GRAPH_RootsSenescentME)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_RootsSenescentME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsSenescentME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsSenescentME)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsSenescentME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsSenescentME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsSenescentME)$isHub=="hub","black","grey"),vertex.label.dist=2)

par(mar=c(0,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RootsVegetativeWW,layout=SE_LAYOUT_RootsVegetativeWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsVegetativeWW)$weight,vertex.size=2.0)
title("WW Vegetative Roots \n(29 ASVs, 20 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringWW,layout=SE_LAYOUT_RootsFloweringWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsFloweringWW)$weight,vertex.size=2.0)
title("WW Flowering Roots \n(148 ASVs, 157 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentWW,layout=SE_LAYOUT_RootsSenescentWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RootsSenescentWW)$weight,vertex.size=2.0)
title("WW Senescent Roots \n(286 ASVs, 490 interactions)",cex.main=2)

par(mar=c(2,0,1,0),family="serif",font=1)

plot(SE_GRAPH_RootsVegetativeWW,layout=SE_LAYOUT_RootsVegetativeWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub" & V(SE_GRAPH_RootsVegetativeWW)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsVegetativeWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsVegetativeWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsVegetativeWW)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsVegetativeWW)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_RootsFloweringWW,layout=SE_LAYOUT_RootsFloweringWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsFloweringWW)$isHub=="hub" & V(SE_GRAPH_RootsFloweringWW)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsFloweringWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsFloweringWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsFloweringWW)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsFloweringWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsFloweringWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsFloweringWW)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_RootsSenescentWW,layout=SE_LAYOUT_RootsSenescentWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RootsSenescentWW)$isHub=="hub" & V(SE_GRAPH_RootsSenescentWW)$name %in% Root_RecurrentHubs,V(SE_GRAPH_RootsSenescentWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RootsSenescentWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RootsSenescentWW)$name %in% Root_RecurrentHubs & V(SE_GRAPH_RootsSenescentWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RootsSenescentWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RootsSenescentWW)$isHub=="hub","black","grey"),vertex.label.dist=2)

dev.off()

#########################################################################################################################################################

#########################################################################################################################################################
######################## FIGURE S1: Abundance Correlation and Hub Overlap in Phyllosphere Networks

# Make a list of ASVs that meet the hub criteria in more than one phyllosphere network
Phyllo_MasterHubList <- list()
Phyllo_NetworkList <- c(ls(pat="SE_GRAPH_RosLeaves"),ls(pat="SE_GRAPH_Stems"))
for(n in Phyllo_NetworkList){
  print(n)
  Network = get(n)
  HubList <- V(Network)[V(Network)$isHub=="hub"]$name
  Phyllo_MasterHubList[[n]] <- HubList
}
Phyllo_AllHubs<-unique(unlist(Phyllo_MasterHubList))
Phyllo_HubFrequencies <- table(unlist(Phyllo_MasterHubList))
Phyllo_RecurrentHubs <- names(Phyllo_HubFrequencies[Phyllo_HubFrequencies>1])

ppi <- 300
png("~/Documents/KB_F18_FilesForUpload/Figures/Phyllo_RecurrentHubs.png", width=10*ppi, height=10*ppi, res=ppi)

par(mfrow=c(4,3),mar=c(0,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RosLeavesVegetativeME,layout=SE_LAYOUT_RosLeavesVegetativeME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RosLeavesVegetativeME)$weight,vertex.size=2.0)
title("ME Vegetative Leaves \n(63 ASVs, 46 interactions)",cex.main=2)

plot(SE_GRAPH_RosLeavesFloweringME,layout=SE_LAYOUT_RosLeavesFloweringME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RosLeavesFloweringME)$weight,vertex.size=2.0)
title("ME Flowering Leaves \n(90 ASVs, 75 interactions)",cex.main=2)

plot(SE_GRAPH_StemsSenescentME,layout=SE_LAYOUT_StemsSenescentME,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_StemsSenescentME)$weight,vertex.size=2.0)
title("ME Senescent Stems \n(109 ASVs, 153 interactions)",cex.main=2)

par(mar=c(2,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RosLeavesVegetativeME,layout=SE_LAYOUT_RosLeavesVegetativeME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RosLeavesVegetativeME)$isHub=="hub" & V(SE_GRAPH_RosLeavesVegetativeME)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_RosLeavesVegetativeME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RosLeavesVegetativeME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RosLeavesVegetativeME)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_RosLeavesVegetativeME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RosLeavesVegetativeME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RosLeavesVegetativeME)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_RosLeavesFloweringME,layout=SE_LAYOUT_RosLeavesFloweringME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RosLeavesFloweringME)$isHub=="hub" & V(SE_GRAPH_RosLeavesFloweringME)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_RosLeavesFloweringME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RosLeavesFloweringME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RosLeavesFloweringME)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_RosLeavesFloweringME)$isHub=="hub","red",ifelse(V(SE_GRAPH_RosLeavesFloweringME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RosLeavesFloweringME)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_StemsSenescentME,layout=SE_LAYOUT_StemsSenescentME,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_StemsSenescentME)$isHub=="hub" & V(SE_GRAPH_StemsSenescentME)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_StemsSenescentME)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_StemsSenescentME)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_StemsSenescentME)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_StemsSenescentME)$isHub=="hub","red",ifelse(V(SE_GRAPH_StemsSenescentME)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_StemsSenescentME)$isHub=="hub","black","grey"),vertex.label.dist=2)

par(mar=c(0,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RosLeavesVegetativeWW,layout=SE_LAYOUT_RosLeavesVegetativeWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RosLeavesVegetativeWW)$weight,vertex.size=2.0)
title("WW Vegetative Leaves \n(72 ASVs, 82 interactions)",cex.main=2)

plot(SE_GRAPH_RosLeavesFloweringWW,layout=SE_LAYOUT_RosLeavesFloweringWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_RosLeavesFloweringWW)$weight,vertex.size=2.0)
title("WW Flowering Leaves \n(135 ASVs, 199 interactions)",cex.main=2)

plot(SE_GRAPH_StemsSenescentWW,layout=SE_LAYOUT_StemsSenescentWW,weighted=TRUE,vertex.color="white",vertex.label=NA,edge.width=E(SE_GRAPH_StemsSenescentWW)$weight,vertex.size=2.0)
title("WW Senescent Stems \n(182 ASVs, 396 interactions)",cex.main=2)

par(mar=c(2,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RosLeavesVegetativeWW,layout=SE_LAYOUT_RosLeavesVegetativeWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RosLeavesVegetativeWW)$isHub=="hub" & V(SE_GRAPH_RosLeavesVegetativeWW)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_RosLeavesVegetativeWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RosLeavesVegetativeWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RosLeavesVegetativeWW)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_RosLeavesVegetativeWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RosLeavesVegetativeWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RosLeavesVegetativeWW)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_RosLeavesFloweringWW,layout=SE_LAYOUT_RosLeavesFloweringWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_RosLeavesFloweringWW)$isHub=="hub" & V(SE_GRAPH_RosLeavesFloweringWW)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_RosLeavesFloweringWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_RosLeavesFloweringWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_RosLeavesFloweringWW)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_RosLeavesFloweringWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_RosLeavesFloweringWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_RosLeavesFloweringWW)$isHub=="hub","black","grey"),vertex.label.dist=2)

plot(SE_GRAPH_StemsSenescentWW,layout=SE_LAYOUT_StemsSenescentWW,weighted=TRUE,vertex.label = ifelse(V(SE_GRAPH_StemsSenescentWW)$isHub=="hub" & V(SE_GRAPH_StemsSenescentWW)$name %in% Phyllo_RecurrentHubs,V(SE_GRAPH_StemsSenescentWW)$label,NA),
     vertex.label.font=2,vertex.label.color="black", edge.width=E(SE_GRAPH_StemsSenescentWW)$weight, vertex.label.cex=1.5,
     vertex.color=ifelse(V(SE_GRAPH_StemsSenescentWW)$name %in% Phyllo_RecurrentHubs & V(SE_GRAPH_StemsSenescentWW)$isHub=="hub","red",ifelse(V(SE_GRAPH_StemsSenescentWW)$isHub=="hub","black","white")),edge.color="grey",
     vertex.frame.color = ifelse(V(SE_GRAPH_StemsSenescentWW)$isHub=="hub","black","grey"),vertex.label.dist=2)

dev.off()

#########################################################################################################################################################

########################################################################################################################################################
######################## FIGURE: Abundance Correlations vs. Phylogenetic Distance in Root Networks

# Multiply phylogenetic distances by a constant so they are between 0 and 1, or rescale them if necessary.
# Make a color gradient between black (low) and red (high).
colorscale <- colorRamp(c('black','red'))
E(SE_GRAPH_RootsVegetativeME)$colorPhy = apply(colorscale(E(SE_GRAPH_RootsVegetativeME)$BLen*0.4), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
E(SE_GRAPH_RootsFloweringME)$colorPhy = apply(colorscale(E(SE_GRAPH_RootsFloweringME)$BLen*0.4), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
E(SE_GRAPH_RootsSenescentME)$colorPhy = apply(colorscale(E(SE_GRAPH_RootsSenescentME)$BLen*0.4), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
E(SE_GRAPH_RootsVegetativeWW)$colorPhy = apply(colorscale(E(SE_GRAPH_RootsVegetativeWW)$BLen*0.4), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
E(SE_GRAPH_RootsFloweringWW)$colorPhy = apply(colorscale(E(SE_GRAPH_RootsFloweringWW)$BLen*0.4), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
E(SE_GRAPH_RootsSenescentWW)$colorPhy = apply(colorscale(E(SE_GRAPH_RootsSenescentWW)$BLen*0.4), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))

ppi <- 300
png("~/Documents/KB_F18_FilesForUpload/Figures/Roots_PhylogeneticDistance.png", width=10*ppi, height=10*ppi, res=ppi)
par(mfrow=c(2,3),mar=c(2,0,4,0),family="serif",font=1)

plot(SE_GRAPH_RootsVegetativeME,layout=SE_LAYOUT_RootsVegetativeME,weighted=TRUE,vertex.color="white",vertex.label=NA,vertex.size=2.0,
     edge.width=E(SE_GRAPH_RootsVegetativeME)$weight,edge.color=E(SE_GRAPH_RootsVegetativeME)$colorPhy)
title("ME Vegetative Roots \n(60 ASVs, 55 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringME,layout=SE_LAYOUT_RootsFloweringME,weighted=TRUE,vertex.color="white",vertex.label=NA,vertex.size=2.0,
     edge.width=E(SE_GRAPH_RootsFloweringME)$weight,edge.color=E(SE_GRAPH_RootsFloweringME)$colorPhy)
title("ME Flowering Roots \n(180 ASVs, 261 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentME,layout=SE_LAYOUT_RootsSenescentME,weighted=TRUE,vertex.color="white",vertex.label=NA,vertex.size=2.0,
     edge.width=E(SE_GRAPH_RootsSenescentME)$weight,edge.color=E(SE_GRAPH_RootsSenescentME)$colorPhy)
title("ME Senescent Roots \n(257 ASVs, 696 interactions)",cex.main=2)

plot(SE_GRAPH_RootsVegetativeWW,layout=SE_LAYOUT_RootsVegetativeWW,weighted=TRUE,vertex.color="white",vertex.label=NA,vertex.size=2.0,
     edge.width=E(SE_GRAPH_RootsVegetativeWW)$weight,edge.color=E(SE_GRAPH_RootsVegetativeWW)$colorPhy)
title("WW Vegetative Roots \n(29 ASVs, 20 interactions)",cex.main=2)

plot(SE_GRAPH_RootsFloweringWW,layout=SE_LAYOUT_RootsFloweringWW,weighted=TRUE,vertex.color="white",vertex.label=NA,vertex.size=2.0,
     edge.width=E(SE_GRAPH_RootsFloweringWW)$weight,edge.color=E(SE_GRAPH_RootsFloweringWW)$colorPhy)
title("WW Flowering Roots \n(148 ASVs, 157 interactions)",cex.main=2)

plot(SE_GRAPH_RootsSenescentWW,layout=SE_LAYOUT_RootsSenescentWW,weighted=TRUE,vertex.color="white",vertex.label=NA,vertex.size=2.0,
     edge.width=E(SE_GRAPH_RootsSenescentWW)$weight,edge.color=E(SE_GRAPH_RootsSenescentWW)$colorPhy)
title("WW Senescent Roots \n(286 ASVs, 490 interactions)",cex.main=2)

dev.off()

########################################################################################################################################################
