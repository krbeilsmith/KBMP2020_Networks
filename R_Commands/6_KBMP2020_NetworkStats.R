# KBMP2020_NetworkStats
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Winter 2020

# Outline
# For each dataset, compare the inverse covariance (Cov) and correlation (Cor) networks for the following properties:
# Number of connected vertices
# Number of edges and proportion of edges with positive vs. negative correlations between the ASV abundances
# Number of vertices (ASVs) designated as hubs
# Fit of a power law to the degree distribution (to assess whether network is scale-free)
# Fit of a linear model to log10 clustering coefficient against log10 degree for each vertex (to assess whether network is hierarchical)
# Fit of a linear model of correlation absolute value vs. phylogenetic distance between the bacterial ASVs
# In addition, find the number of edges and the number of hub vertices shared between networks inferred with the covariance (Cov) and
# correlation (Cor) approaches from the same data.

# Find networks in workspace.
# Use networks pruned for vertices with 0 edges.
NetworkList <- ls(pat="_GRAPH")

# Set up data frame for results.
network_stats <- NULL
network_stats <- data.frame(matrix(ncol = 31, nrow = 0),stringsAsFactors = FALSE)
network_stats <- setNames(data.frame(network_stats), c("Site","Stage","Tissue",
                                                       "Vertices Cov","Vertices Cor",
                                                       "Edges Cov","Edges Cor","Shared Edges",
                                                       "Hubs Cov","Hubs Cor","Shared Hubs",
                                                       "Scale-free alpha Cov","Scale-free p Cov",
                                                       "Scale-free alpha Cor","Scale-free p Cor",
                                                       "Hierarchy r Cov","Hierarchy m Cov","Hierarchy p Cov",
                                                       "Hierarchy r Cor","Hierarchy m Cor","Hierarchy p Cor",
                                                       "Phy r Cov","Phy m Cov","Phy p Cov",
                                                       "Phy r Cor","Phy m Cor","Phy p Cor",
                                                       "% Pos Cov","% Pos Cor", "Modularity Z Cov", "Modularity Z Cor"))


# For each inverse covariance network:
for(l in NetworkList[grepl("SE",NetworkList)]){
  n <- get(l) # Get the network object
  
  site <- ifelse(grepl("ME", l),"ME","WW") # Get the field site
  stage <- ifelse(grepl("Flo", l),"Flowering",ifelse(grepl("Veg", l),"Vegetative","Senescent")) # Get the stage
  tissue <- ifelse(grepl("Roots",l),"Roots",ifelse(grepl("RosLeaves", l),"RosLeaves",ifelse(grepl("Stems", l),"Stems","Siliques"))) # Get the tissue
  
  vertcount_se <- length(V(n)) # Count of vertices connected in the network
  edgecount_se <- length(E(n)) # Count of edges in the network
  # Fraction of edges that represent positive correlations between ASV abundances
  prop_pos_se <- formatC(length(E(n)$direction[E(n)$direction==1])/length(E(n)$direction), format = "g", digits = 2)
  hubs_se <- V(n)[V(n)$isHub=="hub"]$name # Names of hubs in the network
  edges_se <- paste(ends(n, E(n),names = TRUE)[,1],ends(n, E(n),names = TRUE)[,2],sep=" ") # Edge names in the inverse covariance network
  
  # Fit degree distribution to power law and perform goodness of fit test; take scaling coefficient alpha and p-value from test:
  scalefree <- igraph::fit_power_law(V(n)$degree,xmin=2,implementation="plfit") # Fit power law to degree data
  scalefree_p_se <- formatC(scalefree$KS.p, format = "g", digits = 2) # Get the p value for the fit
  scalefree_a_se <- formatC(scalefree$alpha, format = "g", digits = 2) # Get the alpha for the fit
  
  # Get the degree and local clustering coefficient for each vertex, fit linear model to the data, and take r^2, slope, and p-value from model:
  V(n)$clustering <- igraph::transitivity(n,"local") # Get the clustering coefficient
  clustering_data <- as.data.frame(cbind(log10(V(n)$degree),log10(V(n)$clustering))) # Dataframe with log clustering and log degree data for vertices
  colnames(clustering_data) <- c("Degree","Clustering")
  clustering_data <- clustering_data[complete.cases(clustering_data), ] # Remove NAs
  clustering_data <- clustering_data[is.finite(clustering_data$Clustering), ] # Remove infinite values
  if(length(rownames(clustering_data))>0){
    test <- lm(Clustering ~ Degree, data = clustering_data) # Fit linear model of log10 clustering vs log 10 degree
    hierarchy_r_se <- formatC(summary(test)$adj.r.squared, format = "g", digits = 2) # Get the R^2 for the fit
    hierarchy_m_se <- formatC(coef(summary(test))["Degree","Estimate"], format = "g", digits = 2) # Get the slope for the fit
    hierarchy_p_se <- formatC(coef(summary(test))["Degree", "Pr(>|t|)"], format = "g", digits = 2) # Get the p-value for the fit
  }else{ # If no ASVs have degree and clustering != NAs or infinite log values (happens with the smallest networks), assign NA for R^2, slope, and p:
    hierarchy_r_se <- NA
    hierarchy_p_se <- NA
    hierarchy_m_se <- NA
  }
  
  # Get the 16S phylogenetic tree branch length between the ASVs connected by each edge and the edge correlation's absolute value:
  phy_data <- as.data.frame(cbind(E(n)$BLen,abs(E(n)$corr)))
  colnames(phy_data) <- c("PhyDist","Corr")
  test <- lm(Corr ~ PhyDist, data = phy_data) # Fit linear model of correlation strength vs phylogenetic distance
  phy_r_se <- formatC(summary(test)$adj.r.squared, format = "g", digits = 2) # Get the R^2 for the fit
  phy_m_se <- formatC(coef(summary(test))["PhyDist","Estimate"], format = "g", digits = 2) # Get the slope for the fit
  phy_p_se <- formatC(coef(summary(test))["PhyDist", "Pr(>|t|)"], format = "g", digits = 2) # Get the p-value for the fit
  
  # Cluster the nodes and get modularity z score compared to distribution of modularity for random networks of the same size:
  ModDist <- c()
  for(t in 1:100){
    gr <- sample_gnm(length(V(n)), length(E(n)))
    cl <- cluster_edge_betweenness(gr, weights = NULL,directed = FALSE, edge.betweenness = TRUE, merges = TRUE,
                                   bridges = TRUE, modularity = TRUE, membership = TRUE)
    ModDist <- c(ModDist,modularity(gr, weights=NULL, cl$membership))
  }
  cl <- cluster_edge_betweenness(n, weights = NULL,directed = FALSE, edge.betweenness = TRUE, merges = TRUE,
                                 bridges = TRUE, modularity = TRUE, membership = TRUE)
  z_mod_se <- (modularity(n, weights=NULL, cl$membership)-(mean(ModDist)))/sd(ModDist)
  z_mod_se <- formatC(z_mod_se, format = "g", digits = 2)
  
  # For the correlation network inferred from the same dataset:
  for(l2 in NetworkList[grepl("CC",NetworkList) & grepl(tissue,NetworkList) & grepl(stage,NetworkList) & grepl(site,NetworkList)]){ 
    
    n2 = get(l2) # Get the network object
    vertcount_cc <- length(V(n2)) # Get the field site
    edgecount_cc <- length(E(n2)) # Get the tissue
    # Fraction of edges that represent positive correlations between ASV abundances
    prop_pos_cc <- formatC(length(E(n2)$direction[E(n2)$direction==1])/length(E(n2)$direction), format = "g", digits = 2)
    hubs_cc <- V(n2)[V(n2)$isHub=="hub"]$name # Names of hubs in the network
    
    # Check for matching edges in the correlation network
    MatchingList<-c() # Make vector to store matching edges from the correlation network
    for(e in edges_se){ # For each edge in the inverse covariance network:
      # Check for matches in the correlation network edge list (there are two possible orders for vertex names in the edge)
      Matches <- ifelse(e %in% paste(ends(n2, E(n2),names = TRUE)[,1],ends(n2, E(n2),names = TRUE)[,2],sep=" ") |
                          e %in% paste(ends(n2, E(n2),names = TRUE)[,2],ends(n2, E(n2),names = TRUE)[,1],sep=" "),1,0)
      # If there is a match, add the edge to the vector
      if(Matches!=0){
        MatchingList <- c(MatchingList,e)
      }
    }
    
    # Fit degree distribution to power law and perform goodness of fit test; take scaling coefficient alpha and p-value from test
    scalefree <- igraph::fit_power_law(V(n2)$degree,xmin=2,implementation="plfit") # Fit power law to degree data
    scalefree_p_cc <- formatC(scalefree$KS.p, format = "g", digits = 2) # Get the p value for the fit
    scalefree_a_cc <- formatC(scalefree$alpha, format = "g", digits = 2) # Get the alpha for the fit
    
    # Get the degree and local clustering coefficient for each vertex, fit linear model to the data, and take r^2, slope, and p-value from model:
    V(n2)$clustering <- igraph::transitivity(n2,"local") # Get the clustering coefficient
    clustering_data <- as.data.frame(cbind(log10(V(n2)$degree),log10(V(n2)$clustering))) # Dataframe with log clustering and log degree data
    colnames(clustering_data) <- c("Degree","Clustering")
    clustering_data <- clustering_data[complete.cases(clustering_data), ] # Remove NAs
    clustering_data <- clustering_data[is.finite(clustering_data$Clustering), ] # Remove infinite values
    if(length(rownames(clustering_data))>0){
      test <- lm(Clustering ~ Degree, data = clustering_data) # Fit linear model of log10 clustering vs log 10 degree
      hierarchy_r_cc <- formatC(summary(test)$adj.r.squared, format = "g", digits = 2) # Get the R^2 for the fit
      hierarchy_m_cc <- formatC(coef(summary(test))["Degree","Estimate"], format = "g", digits = 2) # Get the slope for the fit
      hierarchy_p_cc <- formatC(coef(summary(test))["Degree", "Pr(>|t|)"], format = "g", digits = 2) # Get the p-value for the fit
    }else{ # If no ASVs have degree and clustering != NAs or infinite log values (happens with the smallest networks), assign NA for R^2, slope, and p:
      hierarchy_r_cc <- NA
      hierarchy_p_cc <- NA
      hierarchy_m_cc <- NA
    }
  }
  
  # Get the 16S phylogenetic tree branch length between the ASVs connected by each edge and the edge correlation's absolute value:
  phy_data <- as.data.frame(cbind(E(n2)$BLen,abs(E(n2)$corr)))
  colnames(phy_data) <- c("PhyDist","Corr")
  test <- lm(Corr ~ PhyDist, data = phy_data) # Fit linear model of correlation strength vs phylogenetic distance
  phy_r_cc <- formatC(summary(test)$adj.r.squared, format = "g", digits = 2) # Get the R^2 for the fit
  phy_m_cc <- formatC(coef(summary(test))["PhyDist","Estimate"], format = "g", digits = 2) # Get the slope for the fit
  phy_p_cc <- formatC(coef(summary(test))["PhyDist", "Pr(>|t|)"], format = "g", digits = 2) # Get the p-value for the fit
  
  # Cluster the nodes and get modularity z score compared to distribution of modularity for random networks of the same size:
  ModDist <- c()
  for(t in 1:100){
    gr <- sample_gnm(length(V(n2)), length(E(n2)))
    cl <- cluster_edge_betweenness(gr, weights = NULL,directed = FALSE, edge.betweenness = TRUE, merges = TRUE,
                                   bridges = TRUE, modularity = TRUE, membership = TRUE)
    ModDist <- c(ModDist,modularity(gr, weights=NULL, cl$membership))
  }
  cl <- cluster_edge_betweenness(n2, weights = NULL,directed = FALSE, edge.betweenness = TRUE, merges = TRUE,
                                 bridges = TRUE, modularity = TRUE, membership = TRUE)
  z_mod_cc <- (modularity(n2, weights=NULL, cl$membership)-(mean(ModDist)))/sd(ModDist)
  z_mod_cc <- formatC(z_mod_cc, format = "g", digits = 2)
  
  # How many edges are in the networks inferred with both approaches?
  MatchingList <- unique(MatchingList)
  edge_conservation <- length(MatchingList)
  
  # How many nodes are hubs in the networks inferred with both approaches?
  hub_conservation <- length(hubs_se[hubs_se %in% hubs_cc])
  
  # Make dataframe entry
  network_stats[nrow(network_stats)+1,] <- c(site, stage, tissue, vertcount_se, vertcount_cc, edgecount_se,edgecount_cc,edge_conservation,
                                             length(hubs_se),length(hubs_cc),hub_conservation,
                                             scalefree_a_se, scalefree_p_se, scalefree_a_cc, scalefree_p_cc,
                                             hierarchy_r_se, hierarchy_m_se, hierarchy_p_se, hierarchy_r_cc, hierarchy_m_cc, hierarchy_p_cc,
                                             phy_r_se, phy_m_se, phy_p_se, phy_r_cc, phy_m_cc, phy_p_cc,
                                             prop_pos_se, prop_pos_cc, z_mod_se, z_mod_cc)
}

# Convert categories to factors and order the dataframe
network_stats$Stage <- factor(network_stats$Stage, levels = c("Vegetative", "Flowering", "Senescent"))
network_stats$Tissue <- factor(network_stats$Tissue, levels = c("Roots", "RosLeaves", "Stems","Siliques"))
network_stats <- network_stats[with(network_stats, order(Stage,Tissue)), ]

# Return categories to characters
network_stats$Stage <- as.character(network_stats$Stage)
network_stats$Tissue <- as.character(network_stats$Tissue)

# Save dataframe
write.table(network_stats,file="~/Documents/KBMP2020_Networks/Tables/NetworkStats", sep="\t", row.names=TRUE,quote=FALSE)
# NetworkStats <- read.table("~/Documents/KBMP2020_Networks/Tables/NetworkStats", header=TRUE, sep="\t",row.names=1)

# Make table for viewing, omitting siliques because this tissue was only observed at one timepoint
z <- ztable(network_stats[network_stats$Tissue!="Siliques",])

# Group the table by stage
rgroup=c("Vegetative","Flowering","Senescent")
n.rgroup=c(4,4,4)
z <- z %>% addrgroup(rgroup=rgroup,n.rgroup=n.rgroup)

# View table
z
