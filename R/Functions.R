#prevent t-test errors from crashing
my.t.test <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) {
		obj$p.value=NA
		return(obj) 
	}else{ 
		return(obj)
	}
}

cor.test.somevsall <- function(mat, rows, method) {
        N = length(mat[,1]);
        M = length(rows);
        pvals = matrix(rep(NA,times=N*M),nrow=M,ncol=N);
        strength = matrix(rep(NA,times=N*M),nrow=M,ncol=N);
        for (j in 1:M) {
                row = rows[j];
        for (i in 1:N) {
                test = cor.test(unlist(mat[row,]),unlist(mat[i,]), method=method)
                pvals[j, i] = test$p.val;
                strength[j, i] = test$estimate;
                pvals[j, i] = test$p.val;
                strength[j,i] = test$estimate;
		if (i == j) {pvals[j,i]=1} # Gene no target itself.
        }
        }
        output = list();
        output$pvals = pvals;
        output$strength = strength;
        colnames(output$strength) = rownames(mat)
        rownames(output$strength) = rownames(mat)[rows]
        colnames(output$pvals) = rownames(mat)
        rownames(output$pvals) = rownames(mat)[rows]
        output
}

# Even with or instead of and always get nothing, is something wrong with the pdcors? -> needed to unlist each row of data.
cor.test_pair_interaction <- function (x,method,Adj) {
	require("ppcor")
	gene1 = x[1]
	gene2 = x[2]
	gene1_rows = which(as.character(Adj[,1]) == gene1)
	gene2_rows = which(as.character(Adj[,1]) == gene2)

	gene1_targets = Adj[gene1_rows,2];
	gene2_targets = Adj[gene2_rows,2];

	sharedtargets = unique(gene1_targets[gene1_targets %in% gene2_targets])
	interact_targets = vector();
	for (t in sharedtargets) {
		dcor_1t = Adj[(Adj[,1] == gene1 & Adj[,2] == t),]$strength
		pdcor1t_2 = pcor.test(unlist(Mat[rownames(Mat)==gene1,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene2,]), method=method)$estimate
		dcor_2t = Adj[(Adj[,1] == gene2 & Adj[,2] == t),]$strength
		pdcor2t_1 = pcor.test(unlist(Mat[rownames(Mat)==gene2,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene1,]), method=method)$estimate
		if (abs(dcor_1t) < abs(pdcor1t_2) & abs(dcor_2t) < abs(pdcor2t_1)) {
			# INTERACTION
			if (pdcor1t_2 > 0 & pdcor2t_1 > 0) {
				interact_targets = c(interact_targets,paste("+",t,sep=""));
			} else if (pdcor1t_2 < 0 & pdcor2t_1 < 0) {
				interact_targets = c(interact_targets,paste("-",t,sep=""));
			} else {
				interact_targets = c(interact_targets,t);
			}
		}
	}
	if (length(interact_targets) >= 1) {
		return(list(pair = c(gene1,gene2), targets = interact_targets));
	} else {
		return(NA)
	}
}


dcor.test.somevsall <- function(mat, rows) {
        require("energy")
        require("pdcor2")
        N = length(mat[,1]);
        M = length(rows);
        pvals = matrix(rep(NA,times=N*M),nrow=M,ncol=N);
        strength = matrix(rep(NA,times=N*M),nrow=M,ncol=N);
        direction = matrix(rep(NA,times=N*M),nrow=M,ncol=N);
        for (j in 1:M) {
                row = rows[j];
		bins = quantile(mat[row,],c(0.25,0.75))
		low = unlist(quantile(mat[row,],c(0.25)))
		high = unlist(quantile(mat[row,],c(0.75)))
        for (i in 1:N) {
                test = dcor.ttest(unlist(mat[row,]),unlist(mat[i,]))
                pvals[j, i] = test$p.val;
                strength[j, i] = test$estimate;
		if (i == row) {pvals[j,i]=1} # Gene no target itself.
			
		test_dir = my.t.test(mat[i,(mat[row,]<=low)],mat[i,(mat[row,]>=high)])
		if (is.na(test_dir$p.value)) {
			direction[j,i] = 0;
		} else {
			if (low != high & test_dir$p.value < 0.05) {
				if (test_dir$estimate[1] < test_dir$estimate[2]) {
					direction[j,i] = 1;
				} else {
					direction[j,i] = -1;
				}
			} else {
				direction[j,i] = 0;
			}
		}
        }
        }
        output = list();
        output$pvals = pvals;
        output$strength = strength;
	output$dir = direction
        colnames(output$strength) = rownames(mat)
        rownames(output$strength) = rownames(mat)[rows]
        colnames(output$pvals) = rownames(mat)
        rownames(output$pvals) = rownames(mat)[rows]
        colnames(output$dir) = rownames(mat)
        rownames(output$dir) = rownames(mat)[rows]
        output
}

# Even with or instead of and always get nothing, is something wrong with the pdcors? -> needed to unlist each row of data.
dcor.test_pair_interaction <- function (x, Adj, margin = 0) {
        require("energy")
        require("pdcor2")
	gene1 = x[1]
	gene2 = x[2]
	gene1_rows = which(as.character(Adj[,1]) == gene1)
	gene2_rows = which(as.character(Adj[,1]) == gene2)

	gene1_targets = Adj[gene1_rows,2];
	gene2_targets = Adj[gene2_rows,2];

	sharedtargets = unique(gene1_targets[gene1_targets %in% gene2_targets])
	interact_targets = vector();
	pathway_targets = vector(); seen = 0;
	for (t in sharedtargets) {
		if (sum(as.character(Adj[,1]) == gene1 & as.character(Adj[,2]) == t) > 0 & sum(as.character(Adj[,1]) == gene2 & as.character(Adj[,2]) == t) > 0) {
		if ( gene1 != gene2 & gene1 != t & gene2 != t) {

		dcor_1t = Adj[(Adj[,1] == gene1 & Adj[,2] == t),]$strength
		pdcor1t_2 = pdcor(unlist(Mat[rownames(Mat)==gene1,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene2,]))
		dcor_2t = Adj[(Adj[,1] == gene2 & Adj[,2] == t),]$strength
		pdcor2t_1 = pdcor(unlist(Mat[rownames(Mat)==gene2,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene1,]))
		if (dcor_1t+margin < pdcor1t_2 & dcor_2t+margin < pdcor2t_1) {
			# INTERACTION
			interact_targets = c(interact_targets,t);
#		} else  if (pdcor2t_1 >= dcor_2t-margin & pdcor1t_2 <= 0+margin & pdcor1t_2 < dcor_1t-margin) {
		} else  if (pdcor2t_1 >= dcor_2t-margin & pdcor1t_2 < dcor_1t-margin) {
                        # gene1 -> gene2 -> target
                       	pathway_targets = c(pathway_targets,t);
                        if (seen) {
                                if( first != gene1 ) {seen=2} # disagreement over direction
                        } else {
                                seen = 1;
                        }
                        first <- gene1; second <- gene2;
#		} else	if (pdcor1t_2 >= dcor_1t-margin & pdcor2t_1 <= 0+margin & pdcor2t_1 < dcor_2t-margin) {
		} else	if (pdcor1t_2 >= dcor_1t-margin & pdcor2t_1 < dcor_2t-margin) {
			# gene2 -> gene1 -> target
                        interact_targets = c(interact_targets,t);
                        if (seen) {
                                if( first != gene2 ) {seen=2} # disagreement over direction
                        } else {
                                seen = 1;
                        }
                        first <- gene2; second <- gene1;
                }

		}
		}

	}
	if ((length(pathway_targets) >= 1 & seen == 1 ) & length(interact_targets) >= 1) {
		return(list(complicated = c(first,second),targets=c(pathway_targets,interact_targets)));	
	} else if (length(pathway_targets) >= 1 & seen == 1) {
		return(list(pathway = c(first, second), targets = pathway_targets));
	} else if (length(interact_targets) >= 1) {
		return(list(pair = c(gene1,gene2), targets = interact_targets));
	} else {
		return(NA)
	}
}

dcor.test_pair_pathway <- function (x, Adj, margin=0.001){
        require("energy")
        require("pdcor2")
	gene1 = x[1]
	gene2 = x[2]
	gene1_rows = which(as.character(Adj[,1]) == gene1)
	gene2_rows = which(as.character(Adj[,1]) == gene2)

	gene1_targets = Adj[gene1_rows,2];
	gene2_targets = Adj[gene2_rows,2];

	sharedtargets = unique(gene1_targets[gene1_targets %in% gene2_targets])
	interact_targets = vector();
	seen = 0;
	for (t in sharedtargets) {
		dcor_1t = Adj[(Adj[,1] == gene1 & Adj[,2] == t),]$strength
		pdcor1t_2 = pdcor(unlist(Mat[rownames(Mat)==gene1,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene2,]))
		dcor_2t = Adj[(Adj[,1] == gene2 & Adj[,2] == t),]$strength
		pdcor2t_1 = pdcor(unlist(Mat[rownames(Mat)==gene2,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==gene1,]))
		if (pdcor1t_2 >= dcor_1t-margin & pdcor2t_1 <= 0+margin & pdcor2t_1 < dcor_2t-margin) {
			# gene2 -> gene1 -> target
			interact_targets = c(interact_targets,t);
			if (seen) {
				if( first != gene2 ) {seen=2} # disagreement over direction
			} else {
				seen = 1;
			}
			first <- gene2; second <- gene1;
		}
		if (pdcor2t_1 >= dcor_2t-margin & pdcor1t_2 <= 0+margin & pdcor1t_2 < dcor_1t-margin) {
			# gene1 -> gene2 -> target
			interact_targets = c(interact_targets,t);
			if (seen) {
				if( first != gene1 ) {seen=2} # disagreement over direction
			} else {
				seen = 1;
			}
			first <- gene1; second <- gene2;
		}
	}
	if (length(interact_targets) >= 1 & seen == 1) { #only report those with consistent direction
		return(list(pair = c(first,second), targets = interact_targets));
	} else {
		return(NA)
	}
}

do_interactions <- function(Mat,candidates, cor_type, multi_test, margin=0) {
	candidates = as.character(candidates);
	candidates.rows = which(rownames(Mat) %in% candidates);
	candidates.rows = candidates.rows[!is.na(candidates.rows)];
	if (length(candidates.rows) < length(candidates)) {
		warning(paste("Only ",length(candidates.rows)," candidates have expression data.\n", sep=""));
	}
	if (length(candidates.rows) == 0) {stop("No expression data for any candidate genes. Are you sure gene IDs are consistent?");}

	# Get Pair-wise dependencies
	if (cor_type == "dcor") {
		out <- dcor.test.somevsall(Mat,candidates.rows)
	} else if (grep("cor=", cor_type))  {
		method = unlist(strsplit(cor_type,"="));
		out <- cor.test.somevsall(Mat,candidates.rows,method[2]);
	} else {
		warning("Did not recognize specified correlation method, using default (dcor)")
		out <- dcor.test.somevsall(Mat,candidates.rows)
	}
	# Tidy up result (find significant interactions
	pvals <- out$pvals
	strength <- out$strength 
	direction <- out$dir

	# Apply Multiple Testing Correction
	if (multi_test == "bon") {
		Sig <- which(pvals < 0.05/length(pvals[1,]),arr.ind=T)
	} else if (grep("fdr=", multi_test))  {
		method = unlist(strsplit(multi_test,"="));
		tmp <- pvals[p.adjust(pvals,method="fdr") < method[2]];
		Sig <- which(pvals <= max(tmp), arr.ind=T);
	} else if (grep("str=", multi_test)) {
		method = unlist(strsplit(multi_test,"="));
		Sig <- which(strength > method[2],arr.ind=T);
	} else {
		warning("Did not recognize specified multiple-testing correction, using default (bon)")
		Sig <- which(pvals < 0.05/length(pvals[1,]),arr.ind=T)
	}
	Adj <- data.frame(Gene = rownames(pvals)[Sig[,1]], Target = colnames(pvals)[Sig[,2]], pval = apply(Sig,1,function(x){pvals[x[1],x[2]]}), strength = apply(Sig,1,function(x){strength[x[1],x[2]]}), direction = apply(Sig,1,function(x){direction[x[1],x[2]]}))
	
	# Get Interactions by conditioning pairwise interactions 
	#	on each other candidate gene
	pairs <- t(combn(candidates,2))

	if (cor_type == "dcor") {
		interactions <- apply(pairs,1,dcor.test_pair_interaction,Adj=Adj,margin=margin)
	} else if (grep("cor=", cor_type))  {
		method = unlist(strsplit(cor_type,"="));
		interactions <- apply(pairs,1,method=method[2],cor.test_pair_interaction,Adj=Adj)
	} else {
		warning("Did not recognize specified partial correlation method, using default (pdcor)")
		interactions <- apply(pairs,1,dcor.test_pair_interaction,Adj=Adj)
	}

	# Clean & Return results
	interactions <- interactions[!is.na(interactions)]
	return(list(Adj = Adj,interactions = interactions));
}
