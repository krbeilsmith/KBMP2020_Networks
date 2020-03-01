# KBMP2020_NetworkConservation
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Winter 2020

#######################################################################################################################################################
# Conservation of edges and hubs in networks from different tissues, developmental stages, and sites.
#######################################################################################################################################################

######################## CALCULATIONS

# How many edges are present in multiple networks? Are more edges conserved between networks of ASVs from the same tissue, stage, or site?
# How many hub ASVs are present in multiple networks? Are more hubs conserved between networks of ASVs from the same tissue, stage, or site?

# List network objects (use either SE_GRAPH for inverse covariance networks or CC_GRAPH for correlation networks)
NetworkList <- ls(pat="CC_GRAPH")

# Set up a dataframe for the results
shared_features <- NULL
shared_features <- data.frame(matrix(ncol = 6, nrow = 0),stringsAsFactors = FALSE)
shared_features <- setNames(data.frame(shared_features), c("Network","Vertices","Edges","% Shared Edges","Hubs","% Shared Hubs"))

# Lists for the names of edges or hubs conserved between networks
ConservedEdgeList <- list()
ConservedHubList <- list()

# For each network:
for(n in NetworkList){
  print(n)
  Network = get(n) # Get the network object
  site <- ifelse(grepl("ME", n),"ME","WW") # Get the field site
  stage <- ifelse(grepl("Flo", n),"Flowering",ifelse(grepl("Veg", n),"Vegetative","Senescent")) # Get the stage
  tissue <- ifelse(grepl("Roots", n),"Roots",ifelse(grepl("RosLeaves", n),"RosLeaves",ifelse(grepl("Stems", n),"Stems","Siliques"))) # Get the tissue
  
  vertcount <- length(V(Network)) # Count of vertices connected in the network
  edgecount <- length(E(Network)) # Count of edges in the network
  hubs <- V(Network)[V(Network)$isHub=="hub"]$name # Hubs in the network
  EdgeList <- paste(ends(Network, E(Network),names = TRUE)[,1],ends(Network, E(Network),names = TRUE)[,2],sep=" ") # Edges in the network
  
  # For every other network (conditioning on tissue, stage, or site is set three lines below):
  ConservedEdges <- c()
  HubList <- c()
  for(n2 in NetworkList[NetworkList != n & grepl(site,NetworkList)]){ # Specify other networks (of the same tissue, stage, or site)
    Network2 = get(n2) # Get the network object
    MatchingList<-c()
    for(e in EdgeList){ # For each edge in the focal network edge list, find matching edges in the other network:
      # Check for matches in the other network's edge list (there are two possible orders for vertex names in the edge)
      Matches <- ifelse(e %in% paste(ends(Network2, E(Network2),names = TRUE)[,1],ends(Network2, E(Network2),names = TRUE)[,2],sep=" ") |
                          e %in% paste(ends(Network2, E(Network2),names = TRUE)[,2],ends(Network2, E(Network2),names = TRUE)[,1],sep=" "),1,0)
      if(Matches!=0){
        MatchingList <- c(MatchingList,e) # If there is a match, add the edge to a vector
      }
    }
    ConservedEdges <- c(ConservedEdges,MatchingList) # A vector of matching edges between the focal network and other networks
    HubList <- c(HubList,V(Network2)[V(Network2)$isHub=="hub"]$name) # A vector of hubs in the other networks
  }
  
  ConservedHubList[[n]] <- hubs[hubs %in% unique(HubList)] # Add hubs from the other networks that match focal network to a list
  hub_conservation <- length(ConservedHubList[[n]])/length(hubs) # Find proportion of hubs in the focal network conserved in other networks
  
  ConservedEdgeList[[n]] <- unique(ConservedEdges) # Add edges from the focal network that match edges in other networks to a list
  edge_conservation <- length(ConservedEdgeList[[n]])/length(EdgeList) # Find the proportion of edges in the focal network conserved in other networks
  
  shared_features[nrow(shared_features)+1,] <- c(n,vertcount,edgecount,edge_conservation,length(hubs),hub_conservation) # Make dataframe entry
}

# Save dataframe
write.table(shared_features,file="~/Documents/KBMP2020_Networks/Tables/SharedFeatures_CC_Site", sep="\t", row.names=TRUE,quote=FALSE)

shared_features[,"% Shared Edges"] <- as.numeric(shared_features[,"% Shared Edges"])
mean(shared_features[,"% Shared Edges"])
range(shared_features[,"% Shared Edges"])

shared_features[,"% Shared Hubs"] <- as.numeric(shared_features[,"% Shared Hubs"])
mean(shared_features[,"% Shared Hubs"])
range(shared_features[,"% Shared Hubs"])

######################## RESULTS

# For inverse covariance networks:

# Networks from the same tissue share between 0.013 and 0.166 of their edges (mean 0.057).
# Networks from the same tissue share between 0 and 0.2 of their hubs (mean 0.044).

# Networks from the same stage share between 0.023 and 0.244 of their edges (mean 0.080).
# Networks from the same stage share between 0 and 0.571 of their hubs (mean 0.141).

# Networks from the same site share between 0.027 and 0.277 of their edges (mean 0.110).
# Networks from the same site share between 0 and 0.429 of their hubs (mean 0.187).

# For correlation networks:

# Networks from the same tissue share between 0.063 and 0.471 of their edges (mean 0.145).
# Networks from the same tissue share between 0 and 0.667 of their hubs (mean 0.204).

# Networks from the same stage share between 0.044 and 0.420 of their edges (mean 0.165).
# Networks from the same stage share between 0 and 0.571 of their hubs (mean 0.194).

# Networks from the same site share between 0.068 and 0.420 of their edges (mean 0.241).
# Networks from the same site share between 0 and 0.714 of their hubs (mean 0.346).

######################## FIGURE: Sets of edges shared across networks.

MasterEdgeList <- list()
NetworkList <- ls(pat="SE_GRAPH")
for(n in NetworkList){
  print(n)
  Network = get(n)
  EdgeList <- paste(ends(Network, E(Network),names = TRUE)[,1],ends(Network, E(Network),names = TRUE)[,2],sep=" ") # Edges in the network
  MasterEdgeList[[n]] <- EdgeList
}
AllEdges <- unique(unlist(MasterEdgeList))
EdgeFrequencies <- table(unlist(MasterEdgeList))
RecurrentEdges <- names(EdgeFrequencies[EdgeFrequencies>1])

# Create a dataframe with a row for each edge and a column for each network
EdgeTable <- NULL
EdgeTable <- data.frame(matrix(ncol=length(MasterEdgeList),nrow=length(AllEdges)),stringsAsFactors = FALSE,row.names=AllEdges)
EdgeTable <- setNames(data.frame(EdgeTable),names(MasterEdgeList))
for(h in AllEdges){ # for each ASV edge recurring at least once
  print(h)
  for(n in names(MasterEdgeList)){ # for each network's list of edges
    print(n)
    EdgeMatch <- ifelse(h %in% as.character(unlist(MasterEdgeList[n])),1,0) # if edge is in this network, add a 1 to the dataframe; otherwise, add 0
    EdgeTable[h,n] <- EdgeMatch
  }
}

# Rename columns for plot
NewNames <- c()
for(i in names(EdgeTable)){
  NetName <- as.character(strsplit(i,"GRAPH_")[[1]][2])
  NewNames <- c(NewNames, gsub("(.)([[:upper:]])", "\\1 \\2", NetName))
}
names(EdgeTable) <- NewNames

NewNames <- c()
for(i in names(EdgeTable)){
  NewNames <- c(NewNames, gsub(" ", "_", i))
}
names(EdgeTable) <- NewNames

ListForPlot <- c("Roots_Vegetative_WW","Roots_Flowering_WW","Roots_Senescen_WW", 
                 "Ros_Leaves_Vegetative_WW", "Ros_Leaves_Flowering_WW", 
                 "Stems_Senescent_WW","Siliques_Senescent_WW","Roots Vegetative ME","Roots Flowering ME","Roots Senescent ME", 
                 "Ros Leaves Vegetative ME", "Ros Leaves Flowering ME", 
                 "Stems Senescent ME","Siliques Senescent ME")

int_list <- list(list("Roots_Vegetative_WW", "Roots_Flowering_WW"), list("Roots_Flowering_WW", "Roots_Senescent_WW"),
                 list("Roots_Vegetative_WW","Ros_Leaves_Vegetative_WW"), list("Roots_Vegetative_WW","Ros_Leaves_Flowering_WW"),
                 list("Ros_Leaves_Vegetative_WW", "Ros_Leaves_Flowering_WW"),
                 list("Roots_Vegetative_WW", "Stems_Senescent_WW"), list("Roots_Vegetative_WW", "Siliques_Senescent_WW"),
                 list("Roots_Flowering_WW", "Stems_Senescent_WW"), list("Roots_Flowering_WW", "Siliques_Senescent_WW"),
                 list("Roots_Senescent_WW", "Stems_Senescent_WW"), list("Roots_Senescent_WW", "Siliques_Senescent_WW"),
                 list("Roots_Vegetative_ME", "Roots_Flowering_ME"), list("Roots_Flowering_ME", "Roots_Senescent_ME"),
                 list("Roots_Vegetative_ME","Ros_Leaves_Vegetative_ME"), list("Roots_Vegetative_ME","Ros_Leaves_Flowering_ME"),
                 list("Ros_Leaves_Vegetative_ME", "Ros_Leaves_Flowering_ME"),
                 list("Roots_Vegetative_ME", "Stems_Senescent_ME"), list("Roots_Vegetative_ME", "Siliques_Senescent_ME"),
                 list("Roots_Flowering_ME", "Stems_Senescent_ME"), list("Roots_Flowering_ME", "Siliques_Senescent_ME"),
                 list("Roots_Senescent_ME", "Stems_Senescent_ME"), list("Roots_Senescent_ME", "Siliques_Senescent_ME"))

ppi <- 300
png("~/Documents/KBMP2020_Networks/Figures/SharedEdges.png", width=8*ppi, height=8*ppi, res=ppi)

upset(EdgeTable, sets = c(ListForPlot), show.numbers=FALSE, set_size.show=FALSE,
      intersections = int_list,
      keep.order=TRUE, text.scale=1.6,
      order.by = c("degree"), matrix.dot.alpha=0,
      empty.intersections = "off", mainbar.y.label = "Edges Shared Between Networks", sets.x.label = "Edges in the Network")

dev.off()

#######################################################################################################################################################
# Where do hubs appear across tissues, developmental stages, and sites?
#######################################################################################################################################################

######################## CALCULATIONS

# Do the same ASVs meet hub criteria in different networks?
MasterHubList <- list()
NetworkList <- ls(pat="SE_GRAPH")
for(n in NetworkList){
  print(n)
  Network = get(n)
  HubList <- V(Network)[V(Network)$isHub=="hub"]$name
  #assign(paste("HubList",as.character(strsplit(n,"GRAPH")[[1]]),sep="_"),HubList)
  MasterHubList[[n]] <- HubList
}
AllHubs<-unique(unlist(MasterHubList))
HubFrequencies <- table(unlist(MasterHubList))
RecurrentHubs <- names(HubFrequencies[HubFrequencies>1])

# How many ASVs meet the hub criteria in multiple networks?
100*(length(RecurrentHubs)/length(AllHubs))

# for inverse covariance graphs
# 10.71429

# for correlation graph
# 29.62963

######################## FIGURE: Hub ASVs shared across sites, developmental stages, and tissues

# First, create a dataframe with a row for each ASV in the list of hubs and a column for each network
HubTable <- NULL
HubTable <- data.frame(matrix(ncol=length(MasterHubList),nrow=length(AllHubs)),stringsAsFactors = FALSE,row.names=AllHubs)
HubTable <- setNames(data.frame(HubTable),names(MasterHubList))
for(h in AllHubs){ # for each ASV that met the hub criteria in at least one network
  print(h)
  for(n in names(MasterHubList)){ # for each network's list of hubs
    print(n)
    HubMatch <- ifelse(h %in% as.character(unlist(MasterHubList[n])),1,0) # if the ASV is a hub in this network, add a 1 to the dataframe; otherwise, add 0
    HubTable[h,n] <- HubMatch
  }
}

# Adjust the dataframe for plotting
ForPlot <- HubTable[c(RecurrentHubs),] # select hubs that appeared in more than one network
ForPlot$Name <- rownames(ForPlot) # add a separate column for hub name
ForPlot <- gather(ForPlot,key="Network",value="Present",1:14) # transform dataframe from wide to long
ForPlot$Present<-as.factor(ForPlot$Present) # make the 1/0 code a factor rather than numeric
ForPlot$Family <- as.character(tax_table(physeq)[taxa_names(physeq) %in% ForPlot$Name,"Rank5"])

# Use the network name to assign site, developmental stage, and tissue metadata to the hubs
ForPlot$Site <- ifelse(grepl("ME", ForPlot$Network),"ME","WW")
ForPlot$Stage <- ifelse(grepl("Flo", ForPlot$Network),"Flowering",ifelse(grepl("Veg", ForPlot$Network),"Vegetative","Senescent"))
ForPlot$Tissue <- ifelse(grepl("Roots", ForPlot$Network),"Roots",ifelse(grepl("Ros", ForPlot$Network),"Rosettes",ifelse(grepl("Stems", ForPlot$Network),"Stems","Siliques")))

# Order metadata columns
ForPlot$Site <- factor(ForPlot$Site, levels=c("ME","WW"))
ForPlot$Stage <- factor(ForPlot$Stage, levels=c("Vegetative","Flowering","Senescent"))
ForPlot$Tissue <- factor(ForPlot$Tissue, levels=c("Roots","Rosettes","Stems","Siliques"))

# Heatmap plot of ASV hub status in different networks
se <- ggplot(ForPlot, aes(Name,Tissue,fill=as.factor(Present))) + 
  geom_tile(color="black")+
  #coord_equal()+
  theme_bw()+
  scale_fill_manual(breaks=c(0,1),values=c("white","black"))+
  rremove("legend")+
  facet_grid(rows=vars(Site,Stage),cols=vars(Family),switch="y",scales="free_x",space="free_x")+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 12, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 12, color = "black",face = "bold")+
  font("xlab", size = 12, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="",x="\nHub Name (16S Amplicon Sequence Variant)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=5))+
  theme(strip.text.y = element_text(size=12, face="bold"),strip.placement="outside",strip.background.y= element_rect(colour="black",fill=NA,size=1))+
  theme(strip.text.x = element_text(size = 12,angle=90),strip.background.x= element_rect(colour="black",fill=NA,size=1))+
  theme(plot.margin = unit(c(1,3,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))


# Heatmap plot of ASV hub status in different networks
cc <- ggplot(ForPlot, aes(Name,Tissue,fill=as.factor(Present))) + 
  geom_tile(color="black")+
  #coord_equal()+
  theme_bw()+
  scale_fill_manual(breaks=c(0,1),values=c("white","black"))+
  rremove("legend")+
  facet_grid(rows=vars(Site,Stage),cols=vars(Family),switch="y",scales="free_x",space="free_x")+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 12, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 12, color = "black",face = "bold")+
  font("xlab", size = 12, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="",x="\nHub Name (16S Amplicon Sequence Variant)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  scale_x_discrete(expand=c(0.025,1),label=function(x) abbreviate(x, minlength=5))+
  theme(strip.text.y = element_text(size=12, face="bold"),strip.placement="outside",strip.background.y= element_rect(colour="black",fill=NA,size=1))+
  theme(strip.text.x = element_text(size = 12,angle=90),strip.background.x= element_rect(colour="black",fill=NA,size=1))+
  theme(plot.margin = unit(c(1,3,1,1), "cm"))+
  theme(strip.placement = "outside")+ theme(panel.spacing = unit(0, "lines"))+ theme(panel.spacing.y = unit(1, "lines"))

ppi <- 300
tiff("~/Documents/KBMP2020_Networks/Figures/SE_RecurrentHubs", width=8*ppi, height=6*ppi, res=ppi)
par(margin(0,0,0,0))
se
dev.off()

ppi <- 300
tiff("~/Documents/KBMP2020_Networks/Figures/CC_RecurrentHubs", width=14*ppi, height=14*ppi, res=ppi)
par(margin(0,0,0,0))
cc
dev.off()

#######################################################################################################################################################
# How are relationships in the network changing over time?
#######################################################################################################################################################

# Is the low level of edge overlap because the ASVs present are changing or because the relationships between them are changing?

# Inverse covariance networks:
# Between 60-80% of ASVs from the previous stage are retained in the new network. 
# Between 51-87% of ASVs in a network are new.
# < 10% of edges between retained ASVs remain in the network.

# Correlation networks:
# Between 60-80% of ASVs from the previous stage are retained in the new network. 
# Between 54-60% of ASVs in a network are new.
# < 20% of edges between retained ASVs remain in the network.

intersection(SE_GRAPH_RootsVegetativeME, SE_GRAPH_RootsFloweringME, keep.all.vertices = FALSE)
# 132 new ASVs were added in the flowering network (132/180 = 73% new ASVs)
# 48 of the 60 ASVs in the vegetative network are present in the flowering network (80% of ASVs from previous stage retained)
# Among the shared ASVs, 1 of the 55 interactions from the vegetative network is present in the flowering network

intersection(SE_GRAPH_RootsFloweringME, SE_GRAPH_RootsSenescentME, keep.all.vertices = FALSE)
# 131 new ASVs were added in the senescent network (131/257 = 51% new ASVs)
# 126 of the 180 ASVs in the flowering network are present in the senescent network (70% of ASVs from previous stage retained)
# Among the shared ASVs, 9 of the 261 interactions from the flowering network are present in the senescent network

intersection(SE_GRAPH_RootsVegetativeWW, SE_GRAPH_RootsFloweringWW, keep.all.vertices = FALSE)
# 129 new ASVs were added in the flowering network (129/148 = 87% new ASVs)
# 19 of the 29 ASVs in the vegetative network are present in the flowering network (66% of ASVs from previous stage retained)
# Among the shared ASVs, 2 of the 20 interactions in the vegetative network are present in the flowering network

intersection(SE_GRAPH_RootsFloweringWW, SE_GRAPH_RootsSenescentWW, keep.all.vertices = FALSE)
# 195 new ASVs were added in the senescent network (195/286 = 68% new ASVs)
# 91 of the 148 ASVs in the flowering network are present in the senescent network (61% of ASVs from previous stage retained)
# Among the shared ASVs, 8 of the 157 interactions in the flowering network are present in the senescent network

intersection(CC_GRAPH_RootsVegetativeME, CC_GRAPH_RootsFloweringME, keep.all.vertices = FALSE)
# 83 new ASVs were added in the flowering network (83/143 = 58% new ASVs)
# 60 of the 83 ASVs in the vegetative network are present in the flowering network (72% of ASVs from previous stage retained)
# Among the shared ASVs, 12 of the 159 interactions from the vegetative network is present in the flowering network

intersection(CC_GRAPH_RootsFloweringME, CC_GRAPH_RootsSenescentME, keep.all.vertices = FALSE)
# 131 new ASVs were added in the senescent network (131/228 = 57% new ASVs)
# 97 of the 143 ASVs in the flowering network are present in the senescent network (68% of ASVs from previous stage retained)
# Among the shared ASVs, 34 of the 375 interactions from the flowering network are present in the senescent network

intersection(CC_GRAPH_RootsVegetativeWW, CC_GRAPH_RootsFloweringWW, keep.all.vertices = FALSE)
# 91 new ASVs were added in the flowering network (91/170 = 54% new ASVs)
# 79 of the 97 ASVs in the vegetative network are present in the flowering network (81% of ASVs from previous stage retained)
# Among the shared ASVs, 95 of the 448 interactions in the vegetative network are present in the flowering network

intersection(CC_GRAPH_RootsFloweringWW, CC_GRAPH_RootsSenescentWW, keep.all.vertices = FALSE)
# 156 new ASVs were added in the senescent network (156/259 = 60% new ASVs)
# 103 of the 170 ASVs in the flowering network are present in the senescent network (61% of ASVs from previous stage retained)
# Among the shared ASVs, 107 of the 1294 interactions in the flowering network are present in the senescent network

intersection(SE_GRAPH_RosLeavesVegetativeME, SE_GRAPH_RosLeavesFloweringME, keep.all.vertices = FALSE)
# 68 new ASVs were added in the flowering network (68/90 = 76% new ASVs)
# 22 of the 63 ASVs in the vegetative network are present in the flowering network (35% of ASVs from previous stage retained)
# Among the shared ASVs, 1 of the 46 interactions from the vegetative network is present in the flowering network

intersection(SE_GRAPH_RosLeavesVegetativeWW, SE_GRAPH_RosLeavesFloweringWW, keep.all.vertices = FALSE)
# 98 new ASVs were added in the flowering network (98/135 = 73% new ASVs)
# 37 of the 72 ASVs in the vegetative network are present in the flowering network (51% of ASVs from previous stage retained)
# Among the shared ASVs, 3 of the 82 interactions from the vegetative network are present in the flowering network

intersection(SE_GRAPH_RootsVegetativeME, SE_GRAPH_RosLeavesVegetativeME, keep.all.vertices = FALSE)
# 2/ 55 root interactions 

intersection(SE_GRAPH_RootsFloweringME, SE_GRAPH_RosLeavesFloweringME, keep.all.vertices = FALSE)
# 3 / 261 root interactions

intersection(SE_GRAPH_RootsSenescentME, SE_GRAPH_StemsSenescentME, keep.all.vertices = FALSE)
# 8 / 696 root interactions

intersection(SE_GRAPH_RootsVegetativeWW, SE_GRAPH_RosLeavesVegetativeWW, keep.all.vertices = FALSE)
# 0 / 20 root interactions 

intersection(SE_GRAPH_RootsFloweringWW, SE_GRAPH_RosLeavesFloweringWW, keep.all.vertices = FALSE)
# 8 / 157 root interactions

intersection(SE_GRAPH_RootsSenescentWW, SE_GRAPH_StemsSenescentWW, keep.all.vertices = FALSE)
# 13 / 490 root interactions

intersection(SE_GRAPH_RootsVegetativeME, SE_GRAPH_RootsVegetativeWW, keep.all.vertices = FALSE)
# 1 / 55 ME interactions

intersection(SE_GRAPH_RootsFloweringME, SE_GRAPH_RootsFloweringWW, keep.all.vertices = FALSE)
# 3 / 261 ME interactions

intersection(SE_GRAPH_RootsSenescentME, SE_GRAPH_RootsSenescentWW, keep.all.vertices = FALSE)
# 7 / 257 ME interactions

intersection(SE_GRAPH_RosLeavesVegetativeME, SE_GRAPH_RosLeavesVegetativeWW, keep.all.vertices = FALSE)
# 2 / 46 ME interactions

intersection(SE_GRAPH_RosLeavesFloweringME, SE_GRAPH_RosLeavesFloweringWW, keep.all.vertices = FALSE)
# 1 / 75 ME interactions

intersection(SE_GRAPH_StemsSenescentME, SE_GRAPH_StemsSenescentWW, keep.all.vertices = FALSE)
# 11 / 153 ME interactions