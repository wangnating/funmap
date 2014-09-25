#################################################################
#
# New Systems Mapping Application(SysMap1)#
# Common library
#
#    1) interpermu.execute
#    2) interpermu.summary
#    3) interpermu.plot
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################


#--------------------------------------------------------------
# interpermu.execute
# 
# do a permutation to get a p-value vs cutoff table. 
# Same as T10 test.
#
# Input 
#    dat    : data object
# Output
#    Cut_off table
#--------------------------------------------------------------
interpermu.execute<-function( dat_obj )
{
	x2 = array(0, dim=c(0,3));
	#nLoop <- as.numeric( FM2.get_value("permu_loop", def=1000) );
	nLoop=10;
	FM_sys$task_start("Execute the permutation, nCount=", nLoop, "...\n");
	
	n_cluster <- as.numeric( FM2.get_value("cluster_count", def=1) );
	if( n_cluster > 1 && suppressWarnings(require(snow,quietly=TRUE)))
	{
		cat("#Running permutations via a cluster of", n_cluster, "nodes.\n");
		cl <- makeCluster( n_cluster );
		clusterSetupRNG( cl );
		clusterEvalQ(cl, require(Funmap, quietly=TRUE) ) ;
		n_perm <- ceiling( nLoop/n_cluster);
		operm <- clusterCall(cl, fin.InterPermutation_cluster, TRUE, dat_obj, n_perm );
		stopCluster(cl);
		for(j in 1:length(operm))
			x2 <- rbind(x2, operm[[j]]);
	}
	else
		x2 <- fin.InterPermutation_cluster(FALSE, dat_obj, nLoop );
		
	FM_sys$task_stop("Permutation is done.\n");
			
		p_cut <- c( 0.9,  0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
				0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01);
	
		if (nLoop>=1000)
			p_cut <- c(p_cut, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001 );
	
		if (nLoop>=10000)
			p_cut <- c(p_cut, 0.0009, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002, 0.0001);
	
		if (nLoop>=100000)
			p_cut <- c(p_cut, 0.00001);
	
	pv2<-array(0, dim=c(length(p_cut), 2) );
	for(i in 1:length(p_cut))
	{
		pv2.order <- round (length(x2[,3]) * p_cut[i]);
		if ( pv2.order>0)
		{
			pv2[i,1] <- p_cut[i];
			pv2[i,2] <- sort(x2[,3], decreasing=TRUE )[ pv2.order ];
		}
	}
	ret<-list(
		name      = paste(dat_obj$curve_name, ".permu",sep=""),
		pheno_file= dat_obj$pheno_file,
		cross_type= dat_obj$cross_type,
		nLoop	  = nLoop,
		pv_table  = pv2,
		fullres   = x2);
		
	return (ret);
}

#--------------------------------------------------------------
# interpermu.summary
# 
# Summarize the result object of the permutation tests. 
# Used by summary( INT.ret.perm object );
#--------------------------------------------------------------
interpermu.summary<-function( object )
{
	res <- object;
	str <- "";
		str <- paste( str, "Permutation result: \n", sep="" );		
		str <- paste( str, "------------------------------------\n", sep="" );		
		str0 <- sprintf("%15s: %s\n", 	 "Cross", 	FM2.cross$name  );
		str1 <- sprintf("%15s: %d\n", 	 "Loop", 	res$nLoop);
		str <- paste( str, str0, str1, sep="" );
		str <- paste( str, "------------------------------------\n", sep="" );		
	str <- paste( str, "\n", sep="" );
	
	str<- paste( str, "p-value\tCutoff\n", "\n", sep=" " );	
	for (i in 1:length(res$pv_table[,1]) )
	{
		st0 <- sprintf("%7.5f\t%10.5f\n", res$pv_table[i,1], res$pv_table[i,2] )
		str <- paste( str, st0, sep="" );		
	}
	
	str<- paste( str, "\n", sep=" " );
	
	strFile <- fpt.plot_to_doc( res$pheno_file, FM2.get_value("plot_doctype") ) 
	fpt.plot_permutation ( res$pv_table );
	dev.off();
	
	str <- paste(str,"*1:The figure 1 for the permutation result is saved to ", strFile, ".\n", sep="");
	str<- paste( str, "\n", sep="" );
	
	return(str);
}

#--------------------------------------------------------------
# interpermu.plot
# 
# Plot the result object of the permutation tests. 
# Used by plot( INT.ret.perm object );
#--------------------------------------------------------------
interpermu.plot<-function( object )
{
	x11();
	fpt.plot_permutation ( object$pv_table );
}

#--------------------------------------------------------------
# public: fin.InterPermutation_cluster
#
# Permutation test assigned at single cluster.
#
# input 
#     dat : data object
#   ht_obj: hypothesis test object, such as XX.RIL.T10,
#           
#--------------------------------------------------------------
fin.InterPermutation_cluster<-function(bClusterUsed, dat, n_perm )
{
	x2 = array(0, dim=c(0,3));
	
	i<-1
	while (i<=n_perm)
	{
		dat0 <- dat;
		n <- length( dat$phenos_table ); 
		new_index <- sample( n );
		new_p <- dat0$phenos_table[new_index];
		dat0$phenos_table <- new_p;
		
		try( r <- INTER.model$qtlscan(dat0), FALSE);
		if (class(r)=="try-error")
			next;
			
		x2 <- rbind( x2, c( i, mean(r$full_res[,3], na.rm=TRUE ), max(r$full_res[,3], na.rm=TRUE) ) );
		cat(x2)
		
		cat("Permutation Test(", i,"): LR2=", max(r$full_res[,3], na.rm=TRUE), "\n");
		if (!bClusterUsed)
			FM_sys$task_elapsed( i/n_perm, "Loop:", i, "/", n_perm, ",", "$SYS_PROMPT$", "\n" )
		i<-i+1;
	}
	return(x2)
}
