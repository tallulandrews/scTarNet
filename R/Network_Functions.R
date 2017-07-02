
process_Interactions <- function(list_of_interactionlists) {

	pair2consistent <- NA;
	pair2targets <- NA;
	for (i in 1:length(list_of_interactionlists)) {
		interact = list_of_interactionlists[[i]];
		for (j in 1:length(interact)) {
			if (!is.null(interact[[j]]$pair)) {
	                        pair = paste(interact[[j]]$pair, collapse="-")
	                        pair_mat = c(pair,i)
				if (length(pair2consistent) ==1 && is.na(pair2consistent)) {
					pair2consistent <- pair_mat
				} else {
					pair2consistent = rbind(pair2consistent, pair_mat)
				}
	                        targets = interact[[j]]$targets
	                        tar_mat = cbind(rep(pair, times=length(targets)),targets,rep(i, times=length(targets)));
				if (length(pair2targets) == 1 && is.na(pair2targets)) {
					pair2targets <- tar_mat
				} else {
					pair2targets = rbind(pair2targets, tar_mat)
				}
			}
		}
	}
	Pair2Consistent.agg = aggregate(pair2consistent, by=list(pair2consistent[,1]), function(x){length(unique(x))})
	colnames(Pair2Consistent.agg) = c("Interaction", "one", "Recurr")
	Pair2Targets.agg = aggregate(pair2targets, by=list(pair2targets[,1], pair2targets[,2]), function(x){length(unique(x))})
	colnames(Pair2Targets.agg) = c("Interaction", "Target", "one", "oneagain", "Recurr")

	return(list(pair2consistent=Pair2Consistent.agg[,c(1,3)], pair2targets=Pair2Targets.agg[,c(1,2,5)]));
}

process_Dependencies <- function(list_of_matrices) {
	merged <- NA;
	for (i in 1:length(list_of_matrices)) {
		data = list_of_matrices[[i]];
		colnames(data) = c("Gene","Target",paste("pval",i,sep="_"), paste("coeff",i,sep="_"))
		if (length(data[1,]) == 5) {
			data[,4]=data[,4]*data[,5]
		}
		cols = c(1,2,4);
		if (length(merged) == 1 && is.na(merged)) {
			merged <- data[,cols]
		} else {
			merged = merge(merged, data[,cols], by=c("Gene","Target"), all=TRUE)
		}
	}
	merged$Gene = as.character(merged$Gene)
	merged$Target = as.character(merged$Target)
	merged = merged[merged$Gene != merged$Target,];
	consistency <- function (x) { 
	        dat = as.numeric(x[3:length(x)])
	        Ndetect = sum(!is.na(dat));
		return(Ndetect);
	}
	mydirection <- function (x) { 
	        dat = as.numeric(x[3:(length(x)-1)])
	        Ndetect = sum(!is.na(dat));
	        if (sum(dat[!is.na(dat)] >= 0) == Ndetect & sum(dat[!is.na(dat)] > 0) >= 1) {
	                return(1);
	        } else if (sum(dat[!is.na(dat)] <= 0) == Ndetect & sum(dat[!is.na(dat)] < 0) >= 1) {
	                return(-1);
	        } else {
			return(0);
		}
	}
	merged$recurr = apply(merged,1,consistency)
	merged$direction  = apply(merged,1,mydirection)
	return(merged[,c("Gene", "Target", "recurr","direction")]);
}

plot_Dependencies <- function (processed_dependencies, recurr_threshold = length(list_of_matrices), adjlist_file=NA, suppress.plot=FALSE) {
	require("igraph")
	merged = processed_dependencies;
	to_plot = merged[abs(merged$recurr) >= recurr_threshold,]
	mycolours=c("blue","black","red")
	AdjList = data.frame(A=to_plot$Gene, B=to_plot$Target, Weight=to_plot$recurr, Col = mycolours[to_plot$direction+2])
	Graph_to_plot = graph_from_data_frame(AdjList, directed=F)

	if (!is.na(adjlist_file)) {write.table(AdjList, file=adjlist_file, row.names=FALSE, col.names=FALSE, quote=FALSE)}
	if (!suppress.plot) {plot(Graph_to_plot, edge.width=E(Graph_to_plot)$Weight, edge.color=E(Graph_to_plot)$Col, vertex.color="grey80")}
}

plot_Interactions<-function(threshold_interact, threshold_targets=threshold_interact, processed_interactions, output_file, adjlist_file1=NA, adjlist_file2=NA, suppress.plot=FALSE) {
	require("igraph")
	Interactions = processed_interactions$pair2consistent
	Targets = processed_interactions$pair2targets

	Interactions = Interactions[Interactions[,2] >= threshold_interact,]
	Targets = Targets[Targets[,3] >= threshold_targets & Targets[,1] %in% Interactions[,1],]

	Names_Int = as.factor(as.character(Interactions[,1]))
	palette = rainbow(n=max(as.numeric(Names_Int)))
	Colour_Int = palette[as.numeric(Names_Int)]

	tmp = strsplit(as.character(Interactions[,1]),"-");
	AdjList_Int = matrix(unlist(tmp), ncol=2, byrow=T)
	Weights_Int = Interactions[,2];

	AdjList_Tar = vector(length=2);
	Weights_Tar = -1;

	for (pair1 in unique(Targets[,1])) {
		targets1 = Targets[Targets[,1] == pair1,2];
		for (pair2 in unique(Targets[,1])) {
			targets2 = Targets[Targets[,1] == pair2,2];
			weight = sum(targets1 %in% targets2);
			if (weight > 0 & pair1 < pair2) {
				AdjList_Tar = rbind(AdjList_Tar, c(pair1, pair2));
				Weights_Tar = c(Weights_Tar, weight);
			}		
		}
	}

	stats = data.frame(N_Interact = length(AdjList_Int[,1]), N_Shared_Targets = sum(Weights_Tar)+1);
	write.table(stats, file=output_file, row.names=FALSE, col.names=TRUE, quote=FALSE)
	# Write Output #
	if (!is.na(adjlist_file1)) {write.table(cbind(AdjList_Int, Weights_Int), file=adjlist_file1, row.names=FALSE, col.names=FALSE, quote=FALSE)}
	if (length(dim(AdjList_Tar)) == 2) {
		AdjList_Tar = AdjList_Tar[2:length(AdjList_Tar[,1]),];
		Weights_Tar =  Weights_Tar[2:length(Weights_Tar)];
		if (!is.na(adjlist_file2)) {write.table(cbind(AdjList_Tar, Weights_Tar), file=adjlist_file2, row.names=FALSE, col.names=FALSE, quote=FALSE)}
	}

	if (!suppress.plot) {
	if (length(dim(AdjList_Tar)) < 2) {
		Graph_Int = graph_from_data_frame(cbind(AdjList_Int, Weights_Int,Colour_Int), directed=F)
		par(mfrow=c(1,1));
		plot(Graph_Int, edge.width=E(Graph_Int)$Weights_Int, edge.color=E(Graph_Int)$Colour_Int, vertex.color="grey80")
	} else {
		AdjList_Tar = AdjList_Tar[2:length(AdjList_Tar[,1]),];
		Weights_Tar = Weights_Tar[2:length(Weights_Tar)];
		Graph_Int = graph_from_data_frame(cbind(AdjList_Int, Weights_Int,Colour_Int), directed=F)
		Graph_Tar = graph_from_data_frame(cbind(AdjList_Tar, Weights_Tar), directed=F)
		Graph_Tar_VColour = match(V(Graph_Tar)$name,as.character(sort(Names_Int)))
		Graph_Tar_VColour = palette[Graph_Tar_VColour]
		par(mfrow=c(1,2));
		plot(Graph_Int, edge.width=E(Graph_Int)$Weights_Int, edge.color=E(Graph_Int)$Colour_Int, vertex.color="grey80")
		plot(Graph_Tar, edge.width=E(Graph_Tar)$Weights_Tar, vertex.color=Graph_Tar_VColour)
	}
	}
	obj = list(AdjListInt = AdjList_Int, AdjListTar=AdjList_Tar, Tar=Targets, Int=Interactions);
	return(obj);
}

plot_Everything <- function (processed_dependencies, processed_interactions, output_file, adjlist_file=NA, suppress.plot=FALSE, recurr_threshold = 2) {
	dependencies = processed_dependencies;
	interactions = plot_Interactions(recurr_threshold, recurr_threshold, processed_interactions, output_file, suppress.plot=TRUE);
	Targets = interactions$Tar;
	keep = vector();
	for (pair1 in unique(Targets[,1])) {
                targets1 = Targets[Targets[,1] == pair1,2];
		pair1 = gsub("-","~",pair1)
                for (pair2 in unique(Targets[,1])) {
	                targets2 = Targets[Targets[,1] == pair2,2];
			pair2 = gsub("-","~",pair2)
			if (pair1 < pair2) {
	                        shared_targets = targets1[(targets1 %in% targets2)];
				if (length(shared_targets) >=1) {
					for (t in shared_targets) {
						keep = c(keep, paste(pair1,t,sep="~"), paste(pair2,t,sep="~"));
					}
				}
			}
                }
        }
	# Colour scheme
	pos_cor = "red"
	neg_cor = "blue"
	other_cor = "grey45"
	inter_cor = "black"
	targets = "darkgoldenrod1"
	TFs = "darkmagenta"


	kept = matrix(unlist(strsplit(keep,"~")),ncol=3, byrow=T)
	int_nodes = unique(c(kept[,1],kept[,2]))
	tar_nodes = unique(kept[,3])
	Adjlist_tar = unique(rbind(kept[,c(1,3)],kept[,c(2,3)]))
	Adjlist_dir = dependencies$direction[match(paste(Adjlist_tar[,1],Adjlist_tar[,2]), paste(dependencies$Gene, dependencies$Target))]
	Adjlist_inter = unique(kept[,c(1,2)]);
	# Now need to ID and specify type of each edge ie direction of gene-target and true/false for is gene-gene.

	Adjlist_full = rbind(Adjlist_tar, Adjlist_inter)
	
	col_factor = c(neg_cor, other_cor, pos_cor)
	Edge_Colour = c(col_factor[Adjlist_dir+2], rep(inter_cor, times=length(Adjlist_inter[,1])))
	Is_Directed = c(rep(1, times = length(Adjlist_tar[,1])), rep(2, times=length(Adjlist_inter[,1])))
	Node_colour = c(rep(TFs,times=length(int_nodes)),rep(targets,times=length(unique(c(int_nodes,tar_nodes)))-length(int_nodes)))
	names(Node_colour) = unique(c(int_nodes,tar_nodes))

	if (!is.na(adjlist_file)) {write.table(cbind(Adjlist_full, Edge_Colour), file=adjlist_file, row.names=FALSE, col.names=FALSE, quote=FALSE)}
	Graph = graph_from_data_frame(cbind(Adjlist_full, Edge_Colour, Is_Directed), directed=T)
	Graph_VColour = Node_colour[match(V(Graph)$name,names(Node_colour))]
	# edge.arrow.mode = 0 (no arrows), 2=forward arrows

	if (!suppress.plot) {
	        plot(Graph, edge.color=E(Graph)$Edge_Colour, edge.arrow.mode=0, vertex.color=Graph_VColour,edge.width=E(Graph)$Is_Directed)
	}
	return(list(graph=Graph, node_col=Graph_VColour))

}

get_Targets <- function (obj, genes) {
	rows = which(obj$AdjListInt[,1] %in% genes & obj$AdjListInt[,2] %in% genes);
	if (length(rows) >1) {
		links <- apply(obj$AdjListInt[rows,],1,function(x) {paste(x[1],x[2],sep="-")});
	} else if (length(rows) == 1) {
		links <- paste(obj$AdjListInt[rows,1],obj$AdjListInt[rows,2],sep="-");
	} else {next;}
	return(obj$Tar[obj$Tar[,1] %in% links,])
}
