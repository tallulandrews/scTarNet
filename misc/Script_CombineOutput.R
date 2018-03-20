source("Network_Functions.R")
args=commandArgs(trailingOnly=T)
if (length(args) > 4) {
	prefix = args[1];
	writeadjlist = args[2];
	graphtype = args[3];
	recurr_threshold = args[4];
	IntDepend_Outputfiles = args[5:length(args)];
} else {
	stop('Please provide a prefix for outputfiles, whether to write the adjacencylist(s) for the respective graph, which graph type ("two-panel","one-panel","none"), the number of times a relationship but be seen for it to be reported, and at least one expression matrix file as arguments.');
}

inter_list <-list()
depend_list <-list()
for (i in 1:length(IntDepend_Outputfiles)) {
	load(IntDepend_Outputfiles[i])
	inter_list[[length(inter_list)+1]] <- output$interactions
	depend_list[[length(depend_list)+1]] <- output$Adj
}

interaction_data <- process_Interactions(inter_list);
dependency_data <- process_Dependencies(depend_list);

stats_file = paste(prefix,"Stats.out",sep="_")

### Calculate output & write it ###

if (grepl("one-panel", graphtype)) {
	png(paste(prefix,"_OnePanel.png",sep=""),width=10, height=10, units="in", res=300)
	if (writeadjlist) {
		adj_file = paste(prefix,"Everything_Adjlist.out",sep="_");
		plot_Everything(dependency_data, interaction_data,output_file=stats_file,adjlist_file=adj_file, suppress.plot=FALSE,recurr_threshold=recurr_threshold);
	} else {
		plot_Everything(dependency_data, interaction_data, output_file=stats_file,suppress.plot=FALSE,recurr_threshold=recurr_threshold);
	}
	dev.off()
}
if (grepl("two-panel",graphtype)) {
	png(paste(prefix,"_TwoPanel.png",sep=""),width=10*2, height=10, units="in", res=300)
	par(mfrow=c(1,2))
	if (writeadjlist) {
		adj_file1 = paste(prefix,"Interactions_Adjlist.out",sep="_");
		adj_file2 = paste(prefix,"SharedTargets_Adjlist.out",sep="_");
		plot_Interactions(threshold_interact=recurr_threshold, processed_interactions=interaction_data, output_file=stats_file, adjlist_file1=adj_file1, adjlist_file2=adj_file2);
	} else {
		plot_Interactions(threshold_interact=recurr_threshold, processed_interactions=interaction_data, output_file=stats_file);
	}
	dev.off();
}

if (grepl("none",graphtype)) {
	plot_Interactions(threshold_interact=recurr_threshold, interaction_data, output_file=stats_file, suppress.plot=TRUE);
}
