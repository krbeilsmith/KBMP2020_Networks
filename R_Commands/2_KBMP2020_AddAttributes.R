# 2_KBMP2020_AddAttributes
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Fall 2018

# Outline
# (pre-req) Run "ImportNetworks.R" commands to load phyloseq data, networks with edge information, and get a list of ASVs in networks
# Add abundance and prevalence of nodes to graph
# Add phylogenetic distance of edges (interactions) to graph
# Make a color legend for node taxonomy that can apply across all the networks
# Find microbial "hub" taxa
# Add network statistics (degree, centrality) to nodes
# Prune nodes with no edges from graph
# Adjust node size and color, edge size and color, and labels for plotting

# MICROBIAL HUBS
# Agler, M.T., Ruhe, J., Kroll, S., Morhenn, C., Kim, S.T., Weigel, D. and Kemen, E.M., 2016. 
# Microbial hub taxa link host and abiotic factors to plant microbiome variation. PLoS Biology, 14(1), p.e1002352.
# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002352

# Add abundance and prevalence information to ASV nodes in the networks
ntwk_graphs <- ls(pattern="_NETWORK_")
for(ng in ntwk_graphs){
  site <- ifelse(grepl("ME", ng),"ME","WW") # get the field site
  stage <- ifelse(grepl("Flo", ng),"Flowering",ifelse(grepl("Veg", ng),"Vegetative","Senescent")) # get the stage
  tissue <- ifelse(grepl("Roots",ng),"Roots",ifelse(grepl("RosLeaves", ng),"RosLeaves",ifelse(grepl("Stems", ng),"Stems","Siliques"))) # get the tissue
  print(paste(tissue,site,stage,sep=" ** "))
  MatchPhyseq <- subset_samples(physeq,PlantPart==tissue & NetworkStage== stage & Site== site) # prune to samples of the site, stage, and tissue
  MatchPhyseq <- prune_taxa(taxa_names(MatchPhyseq) %in% V(get(ng))$name, MatchPhyseq) # prune to ASVs in the network
  MatchPhyseq <- prune_samples(sample_sums(MatchPhyseq)>0, MatchPhyseq) # remove any samples that drop to zero ASVs
  Network <- get(ng)
  # Mean raw abundance of ASV node (average counts per sample) for the site, stage, and tissue
  V(Network)$RawA=c(as.numeric(as.character(rowMeans(otu_table(MatchPhyseq))[match(V(Network)$name,taxa_names(MatchPhyseq))])))
  # Mean relative abundance of ASV node (average fractional abundance) for the site, stage, and tissue
  RA_table  <- transform_sample_counts(MatchPhyseq,  function(x) x / sum(x))
  V(Network)$RelA=c(as.numeric(as.character(rowMeans(otu_table(RA_table))[match(V(Network)$name,taxa_names(RA_table))])))
  # Prevalence of ASV node (number of samples in which present) for the site, stage, and tissue
  PA_table <- transform_sample_counts(MatchPhyseq, function(x) ifelse(x>0,1,0))
  V(Network)$Prev=c(as.numeric(as.character(rowSums(otu_table(PA_table))[match(V(Network)$name,taxa_names(PA_table))])))
  assign(ng, Network)
}

# Add phylogenetic distance of edges (interactions) to graph
ntwk_graphs <- ls(pattern="_NETWORK_")
for(ng in ntwk_graphs){
  site <- ifelse(grepl("ME", ng),"ME","WW") # get the field site
  stage <- ifelse(grepl("Flo", ng),"Flowering",ifelse(grepl("Veg", ng),"Vegetative","Senescent")) # get the stage
  tissue <- ifelse(grepl("Roots",ng),"Roots",ifelse(grepl("RosLeaves", ng),"RosLeaves",ifelse(grepl("Stems", ng),"Stems","Siliques"))) # get the tissue
  print(paste(tissue,site,stage,sep=" ** "))
  MatchPhyseq <- subset_samples(physeq,PlantPart==tissue & NetworkStage== stage & Site== site) # prune to samples of the site, stage, and tissue
  MatchPhyseq <- prune_taxa(taxa_names(MatchPhyseq) %in% V(get(ng))$name, MatchPhyseq) # prune to ASVs in the network
  MatchPhyseq <- prune_samples(sample_sums(MatchPhyseq)>0, MatchPhyseq) # remove any samples that drop to zero ASVs
  Network <- get(ng)
  
  # Get the 16S tree and get all combinations of ASV pairs in the tree
  PhyTree <- phy_tree(MatchPhyseq)
  Tips <- PhyTree$tip.label
  Combos <- combn(sort(Tips),2)
  
  # Get the distance between each ASV combination from the tree
  # Store in edge list dataframe (vertex 1, vertex 2, phylogenetic distance)
  PairBranches<-NULL
  for(p in 1:ncol(Combos)){
    BLen <- cophenetic(PhyTree)[Combos[1,p], Combos[2,p]]
    NewRow <- cbind(Combos[1,p],Combos[2,p],BLen)
    PairBranches <- rbind(NewRow,PairBranches,rownames(PairBranches,"True"))
  }
  PhyloDist<-as.data.frame(PairBranches)
  PhyloDist$BLen<-as.numeric(as.character(PhyloDist$BLen))
  
  # Edge list to graph
  # Add the phylogenetic distances to each edge in the original network
  p = graph_from_data_frame(PhyloDist, directed = FALSE)
  j = intersection(Network,p,byname=TRUE,keep.all.vertices = FALSE)
  
  assign(ng, j)
}

# Filter dataset for ASVs in networks
# This uses the vector with all the ASVs present across all networks created with the ImportNetwork commands
MyTaxa <- prune_taxa(row.names(tax_table(physeq)) %in% ASVsInNetworks, physeq)

# Make color palette of right length
Categories <- as.character(unique(tax_table(MyTaxa)[,"Rank3"])[,"Rank3"])
Palette <- c(brewer.pal(n = 8, name = 'Dark2'),brewer.pal(n = 6, name = 'Set1'))
Legend <- cbind(Palette,Categories)

# Turn unclassified nodes to white
for(i in 1:nrow(Legend)){
  if(as.character(Legend[i,"Categories"])=="unclassified class"){
    Legend[i,"Palette"]<-as.character("#FFFFFFFF")
    Legend[i,"Categories"]<-as.character("unclassified")
  }
}

# For each network in the workspace...
ntwk_graphs <- ls(pattern="_NETWORK_")
for(ng in ntwk_graphs){
  g <- get(ng)
  print(paste("**********",ng,"**********",sep=" "))
  ###################################
  ################################### DESIGNATE PUTATIVE "HUB" NODES
  print("Finding hubs...")
  DegreeNum <- igraph::degree(g) # Find degree (number of edges attached to a vertex)
  Between <- igraph::betweenness(g) # Find betweenness (the number of geodesics (shortest paths) going through a vertex)
  # Get Bacterial Hubs by 90th percentile of both Degree and Betweenness
  NetStats <- merge(as.data.frame(Between),as.data.frame(DegreeNum),by=0)
  #(PlotMat2$Between-mean(PlotMat2$Between))/sd(PlotMat2$Between)
  NetStats$zScoresBetween<-scale(NetStats$Between, center = TRUE, scale = TRUE)
  NetStats$zScoresDegreeNum<-scale(NetStats$DegreeNum, center = TRUE, scale = TRUE)
  HubsDegBet<-NetStats$Row.names[which(NetStats$zScoresBetween>1.282 & NetStats$zScoresDegreeNum>1.282)]
  # Add node attribute for hub classification
  for(i in 1:length(V(g)$name)){
    if(as.character(V(g)$name)[i] %in% as.character(HubsDegBet)){
      V(g)$isHub[i]=as.character("hub")
    }else{
      V(g)$isHub[i]=as.character("no")
    }
  }
  ###################################
  ################################### ADD NETWORK STATISTICS TO NODES
  print("Adding attributes to network...")
  # Assign degree, betweenness to nodes
  V(g)$degree <- igraph::degree(g)
  V(g)$between <- igraph::betweenness(g)
  ###################################
  ################################### RESIZE NODES
  # Get degree (number of edge attachments) for each vertex, assign as vertex feature called "size" 
  # Multiply by a constant as necessary to get reasonable sizes in the plot
  V(g)$size <- (0.7*igraph::degree(g))
  ###################################
  ################################### REMOVE UNATTACHED NODES
  # Remove vertices with 0 edges
  toRemove <- V(g)[degree(g) == 0] 
  clean <- delete.vertices(g, toRemove)
  ###################################
  ################################### TAXONOMY INFORMATION
  # Add taxonomy information to ASV nodes
  V(clean)$Rank2=as.character(tax_table(MyTaxa)[,2][match(V(clean)$name,taxa_names(MyTaxa))])
  V(clean)$Rank3=as.character(tax_table(MyTaxa)[,3][match(V(clean)$name,taxa_names(MyTaxa))])
  V(clean)$Rank4=as.character(tax_table(MyTaxa)[,4][match(V(clean)$name,taxa_names(MyTaxa))])
  V(clean)$Rank5=as.character(tax_table(MyTaxa)[,5][match(V(clean)$name,taxa_names(MyTaxa))])
  ###################################
  ################################### EDGE INFORMATION
  # Color edges by correlation direction and size by magnitude (with multiplier for visibility)
  # E(clean)$color<- as.factor(E(clean)$direction)
  E(clean)$color <- ifelse(E(clean)$direction==1,"#56B489","#E69F00")
  E(clean)$weight <- abs(E(clean)$corr)*5
  ###################################
  ################################### NODE LABEL AND COLOR
  # Make short ASV node labels
  #clean<-set.vertex.attribute(clean,"label",value=paste(substr(V(clean)$Rank5,1,5),substr(V(clean)$name,1,5),sep=": "))
  clean <- set.vertex.attribute(clean,"label",value=substr(V(clean)$name,1,5))
  # Assign ASV node taxonomy color based on global legend for all networks
  V(clean)$color=as.character(Legend[,"Palette"][match(V(clean)$Rank3,Legend[,"Categories"])])
  ###################################
  ################################### SAVE GRAPH WITH ATTRIBUTES
  # Define a layout for the graph
  l <- layout_with_fr(clean)
  # Name the graph and layout
  assign(paste(strsplit(ng,"_")[[1]][1],"GRAPH",strsplit(ng,"_")[[1]][3],sep="_"),clean)
  assign(paste(strsplit(ng,"_")[[1]][1],"LAYOUT",strsplit(ng,"_")[[1]][3],sep="_"),l)
  ###################################
}
