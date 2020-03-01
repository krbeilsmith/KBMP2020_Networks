# KBMP2020_NetworkDynamics
# Kat Beilsmith
# University of Chicago Department of Ecology and Evolution, Bergelson Lab
# Winter 2020

#######################################################################################################################################################
# How does the relative influence of ASVs change over time? 
#######################################################################################################################################################

# Get the intersections representing nodes conserved across developmental transitions in a tissue.
# Roots vegetative --> flowering
# Roots flowering --> senescent
# Rosettes vegetative --> flowering

Comparisons <- list(c("SE_GRAPH_RootsVegetativeME", "SE_GRAPH_RootsFloweringME"), c("SE_GRAPH_RootsFloweringME", "SE_GRAPH_RootsSenescentME"),
                    c("SE_GRAPH_RootsVegetativeWW", "SE_GRAPH_RootsFloweringWW"), c("SE_GRAPH_RootsFloweringWW", "SE_GRAPH_RootsSenescentWW"),
                    c("SE_GRAPH_RosLeavesVegetativeME", "SE_GRAPH_RosLeavesFloweringME"), c("SE_GRAPH_RosLeavesVegetativeWW", "SE_GRAPH_RosLeavesFloweringWW"))

# In the intersection graphs, the attributes "degree_1" and "degree_2" represent the number of edges connected to a node before and after the transition.
# We are going to find the distribution of delta k, or the changes in node degree for conserved nodes, across different developmental transitions.

# Set up a dataframe for the results
deltaK <- NULL
deltaK <- data.frame(matrix(ncol = 4, nrow = 0),stringsAsFactors = FALSE)
deltaK <- setNames(data.frame(deltaK), c("Tissue","Transition","Site", "deltaK"))

kShift <- NULL
kShift <- data.frame(matrix(ncol = 5, nrow = 0),stringsAsFactors = FALSE)
kShift <- setNames(data.frame(kShift), c("Tissue","Transition","Site", "kShift_stat","kShift_p"))

for(cmp in Comparisons){
  print(cmp)
  tissue <- ifelse(grepl("Roots",cmp[1]), "Roots", "Rosettes") # get the tissue type
  transition <- ifelse(grepl("Senescent",cmp[1]) | grepl("Senescent",cmp[2]), "Senescence", "Flowering") # get the developmental transition
  site <- ifelse(grepl("ME",cmp[1]), "ME", "WW") # get site
  a <- get(cmp[1]) # retrieve network before transition
  b <- get(cmp[2]) # retrieve network after transition
  Graph <- intersection(a, b, keep.all.vertices = FALSE) # get conserved nodes
  GraphDF <- as_data_frame(Graph, what = c("vertices")) # put the degree of conserved nodes before and after in a dataframe
  GraphDF$deltaK <- (GraphDF$degree_2 - GraphDF$degree_1) # add a column for the difference between degree before and after (change in k)
  for(d in GraphDF$deltaK){ # for each deltaK value, add to the dataframe with metadata for plotting
    deltaK[nrow(deltaK)+1,] <- c(tissue,transition,site,d) # dataframe entry
  }
  # pairwise Wilcoxon signed rank test for shift in degree across developmental transition
  kBeforeAfter <- NULL
  kBeforeAfter <- wilcox.test(GraphDF$degree_1, GraphDF$degree_2,
                              alternative = c("two.sided"), paired = TRUE, exact = FALSE)
  kShift[nrow(kShift)+1,] <- c(tissue,transition,site,kBeforeAfter$statistic,kBeforeAfter$p.value)
}

write.table(deltaK,file="~/Documents/KBMP2020_Networks/Tables/DegreeChanges_Dist", sep="\t", row.names=TRUE,quote=FALSE)
write.table(kShift,file="~/Documents/KBMP2020_Networks/Tables/DegreeChanges_Test", sep="\t", row.names=TRUE,quote=FALSE)

ztable(kShift)

######################## FIGURE: Influence shifts across developmental transitions

deltaK$deltaK <- as.numeric(deltaK$deltaK)

deltaK$Tissue <- factor(deltaK$Tissue, levels=c("Rosettes","Roots"))

# Frequency plots
freq <- ggplot(deltaK, aes(deltaK,stat(count)),fill="black",color="black") + geom_histogram(size=1,binwidth=1)+
  facet_grid(rows=vars(Tissue,Transition),switch="y")+
  #coord_cartesian(xlim=c(0,0.1))+
  theme_bw()+
  #scale_fill_manual(values = c("#096F13","#AE6E0F","#C67D10"))+
  #scale_color_manual(values = c("#096F13","#AE6E0F","#C67D10"))+
  theme(strip.background=element_blank(),strip.text.x = element_text(size = 16, color = "black",face = "black"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  font("ylab", size = 16, color = "black",face = "bold")+
  font("xy.text", size = 12, color = "black", face = "bold")+ 
  font("xlab", size = 16, color = "black",face = "bold")+ 
  labs(x="\nChange in Degree across Developmental Transition",y="Frequency among Conserved Nodes\n")+
  theme(strip.text.y = element_text(size=12, face="bold",margin = margin(0,0,0,1, "cm"),vjust=5))+
  theme(strip.placement = "outside")+
  theme(panel.spacing = unit(1, "lines"))+
  rremove("legend")

ggsave("~/Documents/KBMP2020_Networks/Figures/DeltaK.tiff", plot = freq, device = NULL, path = NULL,
       scale = 1.5, width = 6, height = 6, units = c("in"),
       dpi = 600, limitsize = TRUE)

