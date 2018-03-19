combineInteractions <- function(list_of_interaction_tables) {

	pair2consistent <- vector();
	pair2targets <- vector();

	for (i in 1:length(list_of_interaction_tables)) {
		if (class(list_of_interaction_tables[[i]]) == "list") {
			interact <- list_of_interaction_tables[[i]]$Int;
		} else {
			 interact <- list_of_interaction_tables[[i]]
		}
		#pairs <- plyr::ddply(interact, .(TF1,TF2), nrow);
		pairs2tar <- interact[,c("TF1", "TF2", "Target")];
		pair2consistent <- rbind(pair2consistent,unique(interact[,c("TF1", "TF2")]));
		pair2targets <- rbind(pair2targets,pairs2tar);
	}
	pair2consistent <- plyr::ddply(pair2consistent, .(TF1,TF2), nrow) # 3rd column = number of datasets seen in
	pair2targets <- plyr::ddply(pair2targets, .(TF1,TF2,Target), nrow) # 3rd column = number of datsets seen in
	colnames(pair2consistent) <- c("TF1", "TF2", "recurr")
	colnames(pair2targets) <- c("TF1", "TF2", "Target", "recurr")

	return(list(pair2consistent=pair2consistent, pair2targets=pair2targets));
}

combineDependencies <- function(list_of_correlation_tables) {
	merged <- vector();
	for (i in 1:length(list_of_correlation_tables)) {
		data <- list_of_correlation_tables[[i]];
		merged <- rbind(merged, data);
	}
	merge_directions <- function(x) { 
		val <- nrow(x);  
		dir <- x[, 5];
		if (prod(sign(dir) == 0)==1) {
			return(c(val, 0))
		} else if (prod(sign(dir) <= 0)==1) {
			return(c(val, -1))
		} else if (prod(sign(dir) >= 0)==1) {
			return(c(val, 1))
		} else {
			return(c(val, 0))
		}
	}
	merged <- plyr::ddply(merged, .(Gene, Target), merge_directions);
	colnames(merged) <- c("Gene", "Target", "recurr", "direction");

	return(merged);
}
