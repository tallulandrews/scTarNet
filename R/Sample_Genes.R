add_means <- function(old_means, new_means) {
	if (length(old_means) > 1){ 
		consistent <- names(old_means) %in% names(new_means);
		old_means <- old_means[consistent]; 
		old_means <- old_means[sort(names(old_means))];
		consistent <- names(new_means) %in% names(old_means);
		new_means <- new_means[consistent]; 
		new_means <- new_means[sort(names(new_means))];
		return(old_means+new_means);
	} else {
		return(new_means);
	}
}

overall_means <- function(list_of_expr_mats) {
	# Support matrices or rds
	gene_means<-vector()
	for (i in 1:length(list_of_expr_mats)) {
		input <- list_of_expr_mats[[i]];
		if (class(input) == "character" | is.null(dim(input)[1])) {
			if (grepl("\\.rds$", input)) {
				input <- readRDS(input) 
			} else {
				input <- read.table(input, header=T)
			}
		} 
		gene_means <- add_means(gene_means, rowMeans(input))
	}
	gene_means <- gene_means/length(list_of_expr_mats);
	return(gene_means);
}

sample_genes <- function(input_genes, background_genes="all", list_of_expr_mats) {
	orig_genes <- as.character(input_genes);
	gene_means <- overall_means(list_of_expr_mats);
	if (length(background_genes) > length(orig_genes)){
		possible_genes <- background_genes;
	} else if (background_genes[1] == "all") {
		possible_genes <- names(gene_means);
	} else {
		stop("Error: insufficient background genes");
	} 
	possible_genes <- gene_means[names(gene_means) %in% possible_genes];
	possible_genes <- names(possible_genes[possible_genes > min(gene_means) & possible_genes < max(gene_means)]);
	bin_breaks <- c(quantile(possible_genes, probs=c(0,0.2, 0.4, 0.6,0.8)), max(possible_genes))
	bins <- cut(possible_genes, bin_breaks, include.lowest=TRUE)
	new_genes <- sapply(levels(bins), function(lvl) {
			this_set <- possible_genes[bins==lvl];
			new_genes <- sample(this_set, size=sum(names(this_set) %in% orig_genes))
			return(new_genes)
			})

	return(new_genes);
}

