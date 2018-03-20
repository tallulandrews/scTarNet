#prevent t-test errors from crashing
getDirection <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) {
		obj$p.value=NA
		return(0) 
	}else if (obj$p.value < 0.05) {
		return(sign(obj$statistic))
	} else {
		return(0);
	}
}

dcor.test.somevsall <- function(mat, rows, n.cores) {
	# Parallelize
	cl <- parallel::makeCluster(n.cores);
	doParallel::registerDoParallel(cl);
	
        N <- length(mat[,1]);
        M <- length(rows);
        pvals <- matrix(rep(NA,times=N*M),nrow=M,ncol=N);
        strength <- matrix(rep(NA,times=N*M),nrow=M,ncol=N);
        direction <- matrix(rep(NA,times=N*M),nrow=M,ncol=N);


        x <- foreach (i = 1:M, .combine='rbind', .packages='foreach') %dopar% {
                row <- rows[i];
		bins <- quantile(mat[row,], c(0.25,0.75))
		low <- bins[1]
		high <- bins[2]
		if (low == high) {high <- high+10^-5}

        foreach (j=1:N, .combine='c') %dopar% {
		if (j == row) { # gene to self 
			return(data.frame(pv=1, st=1, d=1))
		}		
		getDirection <- function(...) {
		    obj<-try(t.test(...), silent=TRUE)
		    if (is(obj, "try-error")) {
		                obj$p.value=NA
		                return(0)
		        }else if (obj$p.value < 0.05) {
		                return(sign(obj$statistic))
		        } else {
		                return(0);
		        }
		}



                dcor_test = energy::dcor.ttest(unlist(mat[row,]),unlist(mat[j,]))

                pvals = dcor_test$p.val;
                strength = dcor_test$estimate;
		
		dir = getDirection(mat[j,(mat[row,]<=low)],mat[j,(mat[row,]>=high)])
		return(data.frame(pv=pvals, st=strength, d=dir));
        } # inner foreach (each potential target)
        } # outer foreach (each TF)
	stopCluster(cl);

        output = list();
	output$pvals<-x[,grep("pv",colnames(x))]
	output$strength<-x[,grep("st",colnames(x))]
	output$direction<-x[,grep("d",colnames(x))]

        colnames(output$strength) = rownames(mat)
        rownames(output$strength) = rownames(mat)[rows]
        colnames(output$pvals) = rownames(mat)
        rownames(output$pvals) = rownames(mat)[rows]
        colnames(output$dir) = rownames(mat)
        rownames(output$dir) = rownames(mat)[rows]
        return(output)
}

dcor_classify_interaction <- function(x, Mat, Dep, threshold.indirect, threshold.interaction) {

	tf1 <- x[1];
	tf2 <- x[2];
	threshold.interaction <- 1 + threshold.interaction

	tf1.targets <- Dep[ which(as.character(Dep[,1]) == tf1) , 2]; tf1.targets[tf1.targets != tf2]; # triplets only
	tf2.targets <- Dep[ which(as.character(Dep[,1]) == tf1) , 2]; tf1.targets[tf2.targets != tf1]; # triplets only
	sharedtargets = intersect(as.character(tf1.targets), as.character(tf2.targets))
	out <- vector()
	for (t in sharedtargets) {
		dcor_1_t <- Dep[Dep[,1] == tf1, Dep[,2] == t,]$strength;
		dcor_2_t <- Dep[Dep[,1] == tf2, Dep[,2] == t,]$strength;

		pdcor_1_t_g2 <- energy::pdcor(unlist(Mat[rownames(Mat)==tf1,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==tf2,]))
		pdcor_2_t_g1 <- energy::pdcor(unlist(Mat[rownames(Mat)==tf1,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==tf2,]))
		if (pdcor_1_t_g2 > threshold.interaction*dcor_1_t & pdcor_2_t_g1 > threshold.interaction*dcor_2_t) {
			# Interaction
			out <- rbind(out, c(sort(c(tf1, tf2)), t, "interaction"));
		} else if (pdcor_1_t_g2 < threshold.indirect*dcor_1_t) {
			# 1->2->t
			out <- rbind(out, c(tf1, tf2, t, "pathway"));
		} else if (pdcor_2_t_g1 < threshold.indirect*dcor_2_t) {
			# 2->1->t
			out <- rbind(out, c(tf2, tf1, t, "pathway"));
		}
	}
	out <- data.frame(out)
	colnames(out) <- c("TF1", "TF2", "Target", "type");
	return(out);
}

# Step 1
calculateTFstoTargets <- function(Mat, TFs, n.cores=1, mt_correction="bon=0.05"){
	# Process arguments
	undetected <- rowSums(Mat > 0) == 0;
	if (sum(undetected) > 0) {
		print(paste("Removing", sum(undetected), "undetected genes."))
		Mat <- Mat[!undetected,];
	}

	MTmethod = unlist(strsplit(mt_correction,"="));
	
	TFs <- as.character(TFs)
	TFs.rows <- which(rownames(Mat) %in% TFs);
	TFs.rows <- TFs.rows[!is.na(TFs.rows)];
	

	out <- dcor.test.somevsall(Mat, TFs.rows, n.cores)
	pvals <- as.matrix(out$pvals)
 	strength <- as.matrix(out$strength)
  	direction <- out$dir

	if (MTmethod[1] == "bon") {
		Sig <- which(pvals < as.numeric(MTmethod[2])/length(pvals[1,]),arr.ind=T)
	} else if (MTmethod[1] == "str") {
		Sig <- which(strength > method[2], arr.ind=T);
	} else {
		tmp <- max(pvals[p.adjust(pvals,method=method[1]) < method[2]]);
		Sig <- which(pvals <= tmp, arr.ind=T);
	}
	Dep <- data.frame(Gene = rownames(pvals)[Sig[,1]], Target = colnames(pvals)[Sig[,2]], pval = unlist(pvals[Sig]), strength = unlist(strength[Sig]), direction = unlist(direction[Sig]))
	return(Dep);
}

# Step 2
calculateConditionalCors <- function(Mat, TFs, Dep, n.cores=1, threshold.interaction=0.01, bidirectional=TRUE, threshold.indirect=0.5, exclude.indirect=TRUE) {
	# threshold indirect = conditional dcors < this*original decor are equivalent to zero (set to a negative to turn off)
	# threshold interaction = % increase in dcor after conditioning for a correlation to be considered an interaction
	threshold.interaction <- 1 + threshold.interaction

	undetected <- rowSums(Mat > 0) == 0;
	if (sum(undetected) > 0) {
		print(paste("Removing", sum(undetected), "undetected genes."))
		Mat <- Mat[!undetected,];
	}

	TFs <- as.character(TFs)
	TFs.rows <- which(rownames(Mat) %in% TFs);
	TFs.rows <- TFs.rows[!is.na(TFs.rows)];
	
	pairs <- t(combn(TFs,2))
	# Parallelize
	cl <- parallel::makeCluster(n.cores);
	doParallel::registerDoParallel(cl);
        inter <- foreach (i = 1:nrow(pairs), .combine='rbind', .packages='foreach') %dopar% {
		dcor_classify_interaction <- function(x, Mat, Dep, threshold.indirect, threshold.interaction) {

		        tf1 <- x[1];
		        tf2 <- x[2];

		        tf1.targets <- Dep[ which(as.character(Dep[,1]) == tf1) , 2]; tf1.targets <- tf1.targets[tf1.targets != tf2]; # triplets only
		        tf2.targets <- Dep[ which(as.character(Dep[,1]) == tf2) , 2]; tf2.targets <- tf2.targets[tf2.targets != tf1]; # triplets only
		        sharedtargets = intersect(as.character(tf1.targets), as.character(tf2.targets))
		        out <- vector()
		        for (t in sharedtargets) {
		                dcor_1_t <- Dep[Dep[,1] == tf1 & Dep[,2] == t,]$strength;
		                dcor_2_t <- Dep[Dep[,1] == tf2 & Dep[,2] == t,]$strength;

		                pdcor_1_t_g2 <- energy::pdcor(unlist(Mat[rownames(Mat)==tf1,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==tf2,]))
		                pdcor_2_t_g1 <- energy::pdcor(unlist(Mat[rownames(Mat)==tf2,]), unlist(Mat[rownames(Mat)==t,]),unlist(Mat[rownames(Mat)==tf1,]))
				print(c(dcor_1_t, pdcor_1_t_g2))
				print(c(dcor_2_t, pdcor_2_t_g1))
		                if ((pdcor_1_t_g2 > threshold.interaction*dcor_1_t & pdcor_2_t_g1 > threshold.interaction*dcor_2_t) |
					!bidirectional & (pdcor_1_t_g2 > threshold.interaction*dcor_1_t | pdcor_2_t_g1 > threshold.interaction*dcor_2_t )) {
		                        # Interaction
		                        out <- rbind(out, c(sort(c(tf1, tf2)), t, "interaction"));
		                } else if (pdcor_1_t_g2 < threshold.indirect*dcor_1_t) {
		                        # 1->2->t
		                        out <- rbind(out, c(tf1, tf2, t, "pathway"));
		                } else if (pdcor_2_t_g1 < threshold.indirect*dcor_2_t) {
		                        # 2->1->t
		                        out <- rbind(out, c(tf2, tf1, t, "pathway"));
		                }
		        }
		        out <- data.frame(out)
			if (nrow(out) > 0) {
			        colnames(out) <- c("TF1", "TF2", "Target", "type");
			        return(out);
			}
		}

		dcor_classify_interaction(pairs[i,], Mat, Dep, threshold.indirect, threshold.interaction)
	}
	stopCluster(cl);

	if (exclude.indirect) {
		indirect <- unique(inter[inter[,"type"]=="pathway", c("TF2", "Target")])
		direct_int <- dplyr::anti_join(inter, indirect, by=c("TF2", "Target"))
		colnames(Dep) <- c("TF2", "Target", "pval", "strength", "direction");
		Dep <- suppressWarnings(dplyr::anti_join(Dep, indirect, by=c("TF2", "Target"))) # get warnings if not all TFs involved in interactions/pathways
		colnames(Dep) <- c("Gene", "Target", "pval", "strength", "direction");
	} else {
		direct_int <- inter;
	}

	return(list(Dep=Dep, Int=direct_int));
}
