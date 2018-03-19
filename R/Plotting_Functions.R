plotDependencies <- function (combined_dependencies, recurr_threshold=max(combined_dependencies$recurr)/2, suppress.plot=FALSE, seed.col="dodgerblue", target.col="goldenrod1", edge.col=c("blue", "black", "red")) {
	# TFs -> Targets Graph

	#require("igraph")
	merged <- combined_dependencies;

	# Create graph
	to_plot <- merged[abs(merged$recurr) >= recurr_threshold,]
	AdjList = data.frame(Gene=to_plot$Gene, Target=to_plot$Target, Weight=to_plot$recurr, Col = edge.col[to_plot$direction+2])
	Graph_to_plot = graph_from_data_frame(AdjList, directed=F)

	# Vertex colours
	V(Graph_to_plot)$Col <- target.col;
	V(Graph_to_plot)$Col[V(Graph_to_plot)$name %in% combined_dependencies[,1]] <- seed.col;

	# Output
	if (!suppress.plot) {plot(Graph_to_plot, edge.width=E(Graph_to_plot)$Weight, edge.color=E(Graph_to_plot)$Col, vertex.color=V(Graph_to_plot)$Col)}
	return(list(adjlist=AdjList, graph=Graph_to_plot));
}


plotInteractions<-function(combined_interactions, recurr_threshold=max(combined_interactions$pair2consistent[,3])/2, target_threshold=1, suppress.plot=FALSE, seed.col="dodgerblue", interaction.col="black") {
	#TF -> TF weighted by num interaction targets.
	#require("igraph")
	Interactions <- combined_interactions$pair2consistent
	Targets <- combined_interactions$pair2targets

	# Filters 
	Interactions <- Interactions[Interactions[,3] >= recurr_threshold,]
	Targets <- dplyr::inner_join(Targets, Interactions, by=c("TF1", "TF2"))
	Targets <- Targets[,1:4]
	Targets <- Targets[Targets[,4] >= recurr_threshold,]
	NTargets <- plyr::ddply(Targets, .(TF1, TF2), nrow);
	NTargets[,3] <- NTargets[,3]-min(NTargets[,3]);
	NTargets[,3] <- NTargets[,3]/(max(NTargets[,3])+0.01)*4+1;

	# Make and plot graph
	AdjList = data.frame(TF1=NTargets$TF1, TF2=NTargets$TF2, Weight=NTargets[,3])
	Graph_to_plot = graph_from_data_frame(AdjList, directed=F)
	if (!suppress.plot) {plot(Graph_to_plot, edge.width=E(Graph_to_plot)$Weight, edge.color=interaction.col, vertex.color=seed.col)}

	return(list(adjlist=AdjList, graph=Graph_to_plot));
}

# Include pathways?
plotInteractionsWithTargets <- function (combined_dependencies, combined_interactions, recurr_threshold=max(combined_interactions$pair2consistent[,3])/2, suppress.plot=FALSE, seed.col="dodgerblue", target.col="goldenrod1", interaction.col="black", edge.col=c("blue", "grey45", "red")) {

	# Filter Interactions
	interactions <- plotInteractions(combined_interactions, recurr_threshold, suppress.plot=TRUE, seed.col=seed.col, interaction.col=interaction.col);
	interactions <- interactions$adjlist;
	Int_Targets <- combined_interactions$pair2targets
	Int_Targets <- dplyr::inner_join(Int_Targets, interactions, by=c("TF1", "TF2"));
	# Get Directions by filtering dependencies
	dependencies <- plotDependencies(combined_dependencies, recurr_threshold, suppress.plot=TRUE, seed.col=seed.col, target.col=target.col, edge.col=edge.col);
	dependencies <- dependencies$adjlist;
	dependencies1 <- suppressWarnings(dplyr::inner_join(dependencies, data.frame(Gene=Int_Targets$TF1, Target=Int_Targets$Target), by=c("Gene", "Target")))
	dependencies2 <- suppressWarnings(dplyr::inner_join(dependencies, data.frame(Gene=Int_Targets$TF2, Target=Int_Targets$Target), by=c("Gene", "Target")))

	# adjlist
	dependencies1 <- dependencies1[,c("Gene", "Target", "Col")]
	dependencies2 <- dependencies2[,c("Gene", "Target", "Col")]
	interactions <- data.frame(Gene=interactions[,1], Target=interactions[,2], Col=rep(interaction.col, times=nrow(interactions)));
	AdjList <- rbind(dependencies1, dependencies2, interactions)
	AdjList$Width <- rep(1, times=nrow(AdjList));
	AdjList$Width[AdjList$Col == interaction.col] <- 3;
	node_cols <- c(target.col, seed.col);

	Graph = graph_from_data_frame(AdjList, directed=F)
	Node_colour <- node_cols[as.numeric(V(Graph)$name %in% as.character(AdjList[,1]))+1]
	V(Graph)$Node_colour <- Node_colour;

	if (!suppress.plot) {
	        plot(Graph, edge.color=E(Graph)$Col, edge.arrow.mode=0, vertex.color=Node_colour, edge.width=E(Graph)$Width)
	}
	return(list(adjlist=AdjList, graph=Graph))

}
