# KB_F18_NetworkStructure
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Winter 2020

################################################################################################################################
# Look for C(k) consistent with hierarchical structure in the network.
# Ravasz, E. and Barabási, A.L., 2003. Hierarchical organization in complex networks. Physical review E, 67(2), p.026112.
################################################################################################################################

network_stat_plots <- NULL
network_stat_plots <- data.frame(matrix(ncol = 4, nrow = 0),stringsAsFactors = FALSE)
network_stat_plots <- setNames(data.frame(network_stat_plots), c("Site","Stage","Degree","Clustering"))

NetworkList <- ls(pat="_GRAPH")

for(l in NetworkList[grepl("SE",NetworkList)]){
  print(l)
  site <- ifelse(grepl("ME", l),"ME","WW") # Get the field site
  stage <- ifelse(grepl("Flo", l),"Flowering",ifelse(grepl("Veg", l),"Vegetative","Senescent")) # Get the stage
  n <- get(l) # Get the network object
  V(n)$clustering <- igraph::transitivity(n,"local") # Get the clustering coefficient
  #clustering_data <- as.data.frame(cbind(log10(V(n)$degree),log10(V(n)$clustering))) # Dataframe with log clustering and log degree
  clustering_data <- as.data.frame(cbind(V(n)$degree,V(n)$clustering)) # Dataframe with clustering and degree
  colnames(clustering_data) <- c("Degree","Clustering")
  clustering_data$Stage <- stage
  clustering_data$Site <- site
  network_stat_plots <- rbind(network_stat_plots, clustering_data) # Add rows for the vertices of the network to dataframe
}
# Format and order dataframe
network_stat_plots$Degree <- as.numeric(network_stat_plots$Degree)
network_stat_plots$Clustering <- as.numeric(network_stat_plots$Clustering)
network_stat_plots$Stage <- factor(network_stat_plots$Stage, levels=c("Vegetative","Flowering","Senescent"))
network_stat_plots <- network_stat_plots[complete.cases(network_stat_plots), ] # Remove NAs
network_stat_plots <- network_stat_plots[is.finite(network_stat_plots$Clustering), ] # Remove infinite values

mp <- NULL
mp <- ggscatter(data=network_stat_plots,x="Degree",y="Clustering", fill="Stage", color="Stage", shape= "Site", 
                size=2, palette= c("#009E73","#CC79A7","#E69F00"),
                position= position_jitter(width=0.25),
                add = "reg.line",conf.int = FALSE, fullrange = FALSE, rug= FALSE)+ 
  yscale("log10", .format = TRUE) + xscale("log10", .format = TRUE)+
  stat_cor(method = "pearson",aes(color=Stage, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = -.3,label.y=c(-.6,-.7,-.8),size=5)+
  geom_abline(slope = -1, color="black", linetype="dashed", size=1.5)+ 
  theme_bw()+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 16, color = "black", face = "bold")+
  labs(y="log10 Clustering coefficient\n",x="\nlog10 Degree")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),
        legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme(legend.background = element_rect(color="black", size=.5))+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1,size=12))+
  theme(strip.text.x = element_text(size = 16,angle=0),strip.background.x= element_rect(colour="black",fill=NA,size=1))

ggsave("~/Documents/KBMP2020_Networks/Figures/Hierarchy.tiff", plot = mp, device = NULL, path = NULL,
       scale = 2, width = 4, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

################################################################################################################################
# Look for P(k) consistent with a scale free network.
# Barabasi, A.-L. and Albert R. 1999. Emergence of scaling in random networks Science, 286 509–512. 
################################################################################################################################

# Get the distribution of degree and the fraction of nodes with each degree, following demonstration here:
# https://chengjunwang.com/web_data_analysis/demo2_simulate_networks/

network_stat_plots <- NULL
network_stat_plots <- data.frame(matrix(ncol = 4, nrow = 0),stringsAsFactors = FALSE)
network_stat_plots <- setNames(data.frame(network_stat_plots), c("Site","Stage","Degree","Probability"))

NetworkList <- ls(pat="_GRAPH")

for(l in NetworkList[grepl("SE",NetworkList)]){
  print(l)
  site <- ifelse(grepl("ME", l),"ME","WW") # Get the field site
  stage <- ifelse(grepl("Flo", l),"Flowering",ifelse(grepl("Veg", l),"Vegetative","Senescent")) # Get the stage
  n <- get(l) # Get the network object
  # Get degree distribution and probability
  d = degree(n, mode = "all")
  dd = degree.distribution(n, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # Delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # Make dataframe entries for this network
  degree_data = as.data.frame(cbind(degree,probability))
  degree_data$Stage <- stage
  degree_data$Site <- site
  network_stat_plots <- rbind(network_stat_plots, degree_data)
  
  # Simulate scale-free network of the same size using Barabasi-Albert model
  barabasi_graph <- sample_pa(length(V(n)), power = 1, m = 1, out.dist = NULL, out.seq = NULL,
                              out.pref = FALSE, zero.appeal = 1, directed = FALSE, algorithm = c("psumtree"))
  # Get degree distribution and probability
  d = degree(barabasi_graph, mode = "all")
  dd = degree.distribution(barabasi_graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # Delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # Make dataframe entries for this network
  degree_data = as.data.frame(cbind(degree,probability))
  degree_data$Stage <- stage
  degree_data$Site <- "Barabasi-Albert model"
  network_stat_plots <- rbind(network_stat_plots, degree_data)
}
# Format and order dataframe
network_stat_plots$degree <- as.numeric(network_stat_plots$degree)
network_stat_plots$probability <- as.numeric(network_stat_plots$probability)
network_stat_plots$Stage <- factor(network_stat_plots$Stage, levels=c("Vegetative","Flowering","Senescent"))
network_stat_plots <- network_stat_plots[complete.cases(network_stat_plots), ] # Remove NAs

mp <- NULL
mp <- ggscatter(data=network_stat_plots[network_stat_plots$Site!="Barabasi-Albert model",],x="degree",y="probability", 
                fill="Stage", color="Stage", shape= "Site", size=2, palette= c("#009E73","#CC79A7","#E69F00"), facet.by = "Stage")+ 
  yscale("log10", .format = TRUE) + xscale("log10", .format = TRUE)+
  geom_point(data=network_stat_plots[network_stat_plots$Site=="Barabasi-Albert model",],aes(x=degree,y=probability))+
  geom_smooth(data=network_stat_plots[network_stat_plots$Site!="Barabasi-Albert model",],method="lm",aes(color=Stage),se=FALSE)+
  geom_smooth(data=network_stat_plots[network_stat_plots$Site=="Barabasi-Albert model",],method="lm",color="black",linetype=2)+
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="log10 Fraction of Nodes\n",x="\nlog10 Degree")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),
        legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme(legend.background = element_rect(color="black", size=.5))+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1,size=12))+
  theme(strip.text.x = element_text(size = 16,angle=0),strip.background.x= element_rect(colour="black",fill=NA,size=1))

ggsave("~/Documents/KBMP2020_Networks/Figures/ScaleFree.tiff", plot = mp, device = NULL, path = NULL,
       scale = 2, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

################################################################################################################################
# Compare network features with survey data features. ASV Prevalence or Abundance vs Degree or Betweenness in the network.
################################################################################################################################

# Correlation of prevalence and abundance across connected vertices
assortativity(n, V(n)$Prev)
assortativity(n, V(n)$RawA)
assortativity(n, V(n)$RelA)

network_stat_plots <- NULL
network_stat_plots <- data.frame(matrix(ncol = 4, nrow = 0),stringsAsFactors = FALSE)
network_stat_plots <- setNames(data.frame(network_stat_plots), c("Site","Stage","Degree","Probability"))

NetworkList <- ls(pat="_GRAPH")

for(l in NetworkList[grepl("SE",NetworkList)]){
  print(l)
  site <- ifelse(grepl("ME", l),"ME","WW") # Get the field site
  stage <- ifelse(grepl("Flo", l),"Flowering",ifelse(grepl("Veg", l),"Vegetative","Senescent")) # Get the stage
  n <- get(l) # Get the network object
  # Get node properties
  vertex_data <- get.data.frame(n, what= c("vertices") )
  vertex_data$Stage <- stage
  vertex_data$Site <- site
  network_stat_plots <- rbind(network_stat_plots, vertex_data)
}
# Format and order dataframe
network_stat_plots$degree <- as.numeric(network_stat_plots$degree)
network_stat_plots$between <- as.numeric(network_stat_plots$between)
network_stat_plots$Prev <- as.numeric(network_stat_plots$Prev)
network_stat_plots$RawA <- as.numeric(network_stat_plots$RawA)
network_stat_plots$RelA <- as.numeric(network_stat_plots$RelA)
network_stat_plots$Stage <- factor(network_stat_plots$Stage, levels=c("Vegetative","Flowering","Senescent"))
network_stat_plots <- network_stat_plots[complete.cases(network_stat_plots), ] # Remove NAs

PrevVDegree = lm(network_stat_plots$Prev ~ network_stat_plots$degree)
summary(PrevVDegree)

RelAVDegree = lm(network_stat_plots$RelA ~ network_stat_plots$degree)
summary(RelAVDegree)

PrevVBetween = lm(network_stat_plots$Prev ~ network_stat_plots$between)
summary(PrevVBetween)

RelAVBetween = lm(network_stat_plots$RelA ~ network_stat_plots$between)
summary(RelAVBetween)

mp <- NULL
mp <- ggscatter(data=network_stat_plots,x="degree",y="Prev", 
                fill="Stage", color="Stage", shape= "Site", size=2, palette= c("#009E73","#CC79A7","#E69F00"), facet.by = "Stage")+ 
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="ASV Prevalence\n",x="\nDegree")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),
        legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme(legend.background = element_rect(color="black", size=.5))+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1,size=12))+
  theme(strip.text.x = element_text(size = 16,angle=0),strip.background.x= element_rect(colour="black",fill=NA,size=1))

ggsave("~/Documents/KBMP2020_Networks/Figures/PrevVDegree.tiff", plot = mp, device = NULL, path = NULL,
       scale = 2, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

mp <- NULL
mp <- ggscatter(data=network_stat_plots,x="between",y="Prev", 
                fill="Stage", color="Stage", shape= "Site", size=2, palette= c("#009E73","#CC79A7","#E69F00"), facet.by = "Stage")+ 
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="ASV Prevalence\n",x="\nBetweenness Centrality")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),
        legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme(legend.background = element_rect(color="black", size=.5))+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1,size=12))+
  theme(strip.text.x = element_text(size = 16,angle=0),strip.background.x= element_rect(colour="black",fill=NA,size=1))

ggsave("~/Documents/KBMP2020_Networks/Figures/PrevVBetween.tiff", plot = mp, device = NULL, path = NULL,
       scale = 2, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)


mp <- NULL
mp <- ggscatter(data=network_stat_plots,x="degree",y="RelA", 
                fill="Stage", color="Stage", shape= "Site", size=2, palette= c("#009E73","#CC79A7","#E69F00"), facet.by = "Stage")+
  yscale("log10", .format = TRUE) + xscale("log10", .format = TRUE)+
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="ASV Relative Abundance\n",x="\nDegree")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),
        legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme(legend.background = element_rect(color="black", size=.5))+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1,size=12))+
  theme(strip.text.x = element_text(size = 16,angle=0),strip.background.x= element_rect(colour="black",fill=NA,size=1))

ggsave("~/Documents/KBMP2020_Networks/Figures/RelAVDegree.tiff", plot = mp, device = NULL, path = NULL,
       scale = 2, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)

mp <- NULL
mp <- ggscatter(data=network_stat_plots,x="between",y="RelA", 
                fill="Stage", color="Stage", shape= "Site", size=2, palette= c("#009E73","#CC79A7","#E69F00"), facet.by = "Stage")+ 
  yscale("log10", .format = TRUE) + xscale("log10", .format = TRUE)+
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.background.y=element_blank(),strip.text.x = element_text(size = 20, color = "black",face = "bold"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xlab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+
  labs(y="ASV Relative Abundance\n",x="\nBetweenness Centrality")+
  theme(legend.text=element_text(size=10,color = "black",face = "bold"),
        legend.title=element_text(size=12,color = "black",face = "bold"),
        legend.position="right")+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme(legend.background = element_rect(color="black", size=.5))+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1,size=12))+
  theme(strip.text.x = element_text(size = 16,angle=0),strip.background.x= element_rect(colour="black",fill=NA,size=1))

ggsave("~/Documents/KBMP2020_Networks/Figures/RelAVBetween.tiff", plot = mp, device = NULL, path = NULL,
       scale = 2, width = 6, height = 4, units = c("in"),
       dpi = 600, limitsize = TRUE)
