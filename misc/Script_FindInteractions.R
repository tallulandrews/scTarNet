source("Functions.R")

# Arguments:
# f="matrixfile"
# g="genesfile"
# p="prefix"
# Correlation Method: corr="dcor", corr="spearman", corr="pearson"
# Multiple Testing Method: mt="bon", mt="fdr=VAL", mt="str=VAL"


# Read in Arguments #
args=commandArgs(trailingOnly=T)
f  = args[1];
g  = args[2];
p  = args[3];
mt = args[4];

print(args);

#if (length(args) > 0) {
#	for(i in 1:length(args)){
#	  eval(parse(text=args[[i]]))
#	}
#}
errmessage="";
to.die = 0;
if (!exists("f")) {
	errmessage = paste(errmessage,'Please provide file containing expression matrix using: f="file.txt"\n');
	to.die=1;
}
if (!exists("g")) {
	errmessage = paste(errmessage,'Please provide file containing candidate genes using: g="file.txt"\n');
	to.die=1;
}
if (!exists("p")) {
	errmessage = paste(errmessage,'No prefix (p="prefix") provided using default: "Analyze"\n');
	p <- "Analyze";
}
#if (!exists("corr")) {
#	errmessage = paste(errmessage,'No correlation method (corr=["dcor","spearman","pearson"]) provided using default: "dcor"\n');
	corr <- "dcor";
#}
if (!exists("mt")) {
	errmessage = paste(errmessage,'No multiple testing method (mt=["bon","fdr=VAL","str=VAL"]) provided using default: "bon"\n');
	mt <- "bon";
}
if (!exists("err")) {
	errmessage = paste(errmessage,'No margin of error (err=VAL]) provided using default: 0\n');
	err <- 0;
}
if (nchar(errmessage) > 5) {
	cat(errmessage);
	if (to.die) { stop();}
}
# Check arguments are valid

# Process Arguments
if (corr == "spearman" | corr == "pearson") {
	corr = paste("cor=",corr,sep="");
}

Mat = as.matrix(read.table(f,header=T)) #expression matrix rows = genes, col = cells
Mat = Mat[rowMeans(Mat) >0,]
can_genes = read.table(g, header=F)
name = p;
candidates = as.character(can_genes[,1]) # vector of candate gene names
# First get all significant dependency relationships for all candidates.
candidates.rows = which(rownames(Mat) %in% candidates);#vector of row indicies -> eg. k-clique community from a PPI network

output = do_interactions(Mat, candidates, corr, mt, err)

#sink(paste(name,"_Interactions_Output.txt", sep=""))
#output$interactions
#sink()

#write.table(output$Adj, file=paste(name,"_Dependencies_Output.txt", sep=""))

save(output, file=paste(name,"_IntDepend_Output.Rd",sep=""))
