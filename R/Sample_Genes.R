# Read in Arguments 
args=commandArgs(trailingOnly=T)
if (length(args) > 3) {
	inputfile = args[1];
	outputfile = args[2];
	backgroundfile = args[3];
	expr_files = args[4:length(args)];
} else {
	stop('Please provide an input file of gene names,  output file name, a file of genes to sample from (or "all"), and at least one expression matrix file as arguments.');
}
#print(args)
#print(expr_files)

add_means<-function(old_means, new_means) {
	if (length(old_means) > 1){ 
		consistent = names(old_means) %in% names(new_means);
		old_means = old_means[consistent]; 
		old_means = old_means[sort(names(old_means))];
		consistent = names(new_means) %in% names(old_means);
		new_means = new_means[consistent]; 
		new_means = new_means[sort(names(new_means))];
		return(old_means+new_means);
	} else {
		return(new_means);
	}
}

# Read in files
count = 0;
expr_means = 0;
for (i in as.vector(expr_files)) {
	if (grepl("\\.rda$",i) | grepl("\\.RData$",i)) {
		# Load RData object
		objs = load(i);
		for (thing in objs) {
			mat<-try(as.matrix(thing),silent=TRUE);
			if (!(class(mat) == "try-error")){
				expr_means <- add_means(expr_means, rowMeans(mat));
				count = count+1;
			} else {
				warning(paste(i, "could not be turned into a matrix. Skipped."));
			}
		}
	} else {
		# Load text file
		thing <- read.table(i, header=T)
		mat<-try(as.matrix(thing),silent=TRUE);
		if (!(class(mat) == "try-error")){
			expr_means <- add_means(expr_means, rowMeans(mat));
			count = count+1;
		} else {
			warning(paste(i, "could not be turned into a matrix. Skipped."));
		}
	}
}
if (count == 0) {stop("No Expression Matrices could be read correctly.");}
expr_means <- expr_means/count;

# get sample
orig_genes = read.table(inputfile, header=F)
orig_genes = as.character(orig_genes[,1]);

orig_expr = expr_means[names(expr_means) %in% orig_genes];
if (!(backgroundfile == "all")) {
	possible_genes = read.table(backgroundfile, header=F);
	possible_genes = as.character(possible_genes[,1]);
} else {
	possible_genes = names(expr_means);
}
possible_genes = expr_means[names(expr_means) %in% possible_genes];

possible_genes = names(possible_genes[possible_genes > min(orig_expr) & possible_genes < max(orig_expr)]);

new_genes = sample(possible_genes, size=sum(names(expr_means) %in% orig_genes));

write.table(new_genes, file=outputfile, col.names=FALSE, row.names=FALSE, quote=FALSE)
