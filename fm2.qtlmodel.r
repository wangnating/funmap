Qtlmodel.mlefunc<-function( parin, y, qtl.prob, time.std, par_fixed=NULL)
{
	len.cov <- FM2.covar$get_par_num(y);
	len.gen <- FM2.curve$par_num;
	par.covar <- c( par_fixed, parin[1:len.cov] );
	if (!FM2.covar$is_valid( par.covar ))
		return(NaN);
	
	cov <- FM2.covar$get_mat( par.covar, time.std );
	if (is.na(cov) )
	{
		cat("CORVAR ERR:", par.covar, "\n");
		browser();
		return(NaN);
	}
	
	m  <- length(y[1,]);
	n  <- length(y[,1]);

	yy0<- NULL;
	yy1<- NULL;
	yy2<- NULL;

	sd.x <- apply(y, 2, function(x){sd(x, na.rm=T);});
	mu.x <- colMeans(y, na.rm=T);

	len <- len.cov;
	if (FM2.cross$gen_QQ)
	{
		QQ2_par<- parin[(len+1):(len+len.gen)];
		mu2 <- FM2.curve$get_mu( QQ2_par, time.std  );
		yy2 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu2, nrow=1 ); 
		len <- len+len.gen;

		p0 <- pnorm( mu2,  mu.x, sd.x);
		#if (length(which(p0<0.05))>=1) return(NaN);
	}

	if (FM2.cross$gen_Qq)
	{
		Qq1_par<- parin[(len+1):(len+len.gen)];
		mu1 <- FM2.curve$get_mu( Qq1_par, time.std  );
		yy1 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu1, nrow=1 ); 
		len <- len+len.gen;
		p0 <- pnorm( mu1,  mu.x, sd.x);
		#if (length(which(p0<0.05))>=1) return(NaN);
	}
	
	if (FM2.cross$gen_qq)
	{
		qq0_par<- parin[(len+1):(len+len.gen)];
		mu0 <- FM2.curve$get_mu( qq0_par, time.std );
		yy0 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu0, nrow=1 ); 
		len <- len+len.gen;
		p0 <- pnorm( mu0,  mu.x, sd.x);
		#if (length(which(p0<0.05))>=1) return(NaN);
	}

	for ( i in 1:length(y[1,]) )
	{
		y.miss <- which( is.na(y[,i]) );
		if(!is.null(yy0)) yy0[y.miss, i]<-0;
		if(!is.null(yy1)) yy1[y.miss, i]<-0;
		if(!is.null(yy2)) yy2[y.miss, i]<-0;
	}
	
	fy0 <- NULL;
	fy1 <- NULL;
	fy2 <- NULL;
	if(!is.null(yy0))
	{
		fy0 <- try( dmvnorm(yy0, rep(0, m), cov), T);
		if (class(fy0)=="try-error") return (NaN);
	}
	if(!is.null(yy1))
	{
		fy1 <- try( dmvnorm(yy1, rep(0, m), cov), T);
		if (class(fy1)=="try-error") return (NaN);
	}
	if(!is.null(yy2))
	{
		fy2 <- try( dmvnorm(yy2, rep(0, m), cov), T);
		if (class(fy2)=="try-error") return (NaN);
	}
	
	pf <- cbind(fy0, fy1, fy2)*qtl.prob;
	A <- -sum(log(rowSums(pf)));

	FM2.set_value( "mle.pdf", length(which( pf>1 | pf<0 ) ) );
cat(A, parin, "\n");
	return (A);
}

Qtlmodel.get_est_param<-function( dat )
{
	len.cov<- FM2.covar$get_par_num(dat);
	len.mu <- FM2.curve$par_num;
	parin <- c(dat$covar$par, dat$curve$par);

	for (i in 2:FM2.cross$gen_num)
		parin <- c(parin, dat$curve$par*runif(1, min=0.99, max=1.01) );

	return(parin);
}

Qtlmodel.old_geno_pos <<- -1;

Qtlmodel.qtlscan<-function( dat )
{
	FM2.set_value( "likelihood_null", NULL );
   	phenos <- as.matrix( dat$phenos_table );
	parin  <- Qtlmodel.get_est_param( dat );

cat("INIT PARAMS:", parin, "\n");

	marker_obj <- dat$marker_obj;
	FM_sys$task_start( "Execute the hypothesis test ", FM2.cross$hp_desc, "...\n" );	
	Qtlmodel.old_geno_pos <<- -1;

  	res<-c();
	for ( grp_idx in 1:marker_obj$count )
	{	
		grp_idx=1;
		marker_grp = marker_obj$grps[[grp_idx]];
		qtls <- fin.get_qtl_scanpoints_grp(marker_grp); 
		filterqtl <- fin.filter.qtl( dat, qtls, marker_grp);
		# search for all marker for current group(chromosome)
		for ( nQtl in  filterqtl )
		{
			y <- phenos;
			mrk_idx <- qtls[nQtl,3]
			par.cross <-list(lmarker = marker_grp$start_idx + mrk_idx-1,
					 rmarker = marker_grp$start_idx + mrk_idx,
					 qtl.pos = qtls[nQtl,1], 
					 dist    = qtls[nQtl,2])
			ret <- Qtlmodel.get_est_LR2(dat, par.cross, parin)
			cat( "QTLSCAN:", grp_idx, mrk_idx, qtls[nQtl,4], ret, "\n");
			res <- rbind(res, c( chr_grp=grp_idx, pos = qtls[nQtl,4] , ret) );
		}
		FM_sys$task_elapsed(finished = grp_idx/marker_obj$count, "Group:", grp_idx, "/", marker_obj$count, ",", "$SYS_PROMPT$","\n" );
	}

	FM_sys$task_stop( "The hypothesis test is done\n");	

	#if (!is.null(ht_obj$colnames))
	#	colnames(res)<-c(  "chr_grp", "pos", ht_obj$colnames )
	ret <- fin.qtl_locate( dat, res);
	return( ret );
}

Qtlmodel.optim<-function(parin, mlefunc, y, qtl.prob =NULL, time.std=NULL )
{
	g_optim  <- FM2.get_value("optim", "BFGS");
	loop  <- 0;
	h.ret <- NA;

	while( loop < 1 )
	{
		h <- NA;
		if (is.null(qtl.prob))
		{
			if (g_optim=="BFGS" || g_optim=="Nelder-Mead" || g_optim=="SANN")
				try( h<- optim( parin, mlefunc, y = y, time.std=time.std, method =g_optim ), F)
			else if (g_optim=="CE")
				h<- optim_ce( parin, mlefunc, y = y, time.std=time.std);
				
		}
		else
		{
			if (g_optim=="BFGS" || g_optim=="Nelder-Mead" || g_optim=="SANN")
				try( h<- optim( parin, mlefunc, y = y, time.std=time.std, qtl.prob=qtl.prob, method =g_optim ), F)
			else if (g_optim=="CE")
				h<- optim_ce( parin, mlefunc, y = y, time.std=time.std, qtl.prob=qtl.prob );

		}
		
		if ( class(h) != "try-error" &&  !is.na(h) && h$value>0 )
		{
			len.cov <- FM2.covar$get_par_num();
			gen_par <- h$par[(len.cov+1):(len.cov+ FM2.curve$par_num )];
			sd0 <- apply(y, 2, function(x){sd(x, na.rm=T);});
			mu0 <- FM2.curve$get_mu( gen_par, time.std  );
			p0 <- pnorm( mu0, colMeans(y, na.rm=T), sd0);
			p1 <-c();
			if(length(h$par)>len.cov+ FM2.curve$par_num*1)
			{
				gen_par <- h$par[(len.cov+1+FM2.curve$par_num):(len.cov+ FM2.curve$par_num*2 )];
				mu1 <- FM2.curve$get_mu( gen_par, time.std  );
				p1 <- pnorm( mu1, colMeans(y, na.rm=T), sd0);
			}
			p2 <- c();
			if(length(h$par)>len.cov+ FM2.curve$par_num*2)
			{
				gen_par <- h$par[(len.cov+1+FM2.curve$par_num*2):(len.cov+ FM2.curve$par_num*3 )];
				mu2 <- FM2.curve$get_mu( gen_par, time.std  );
				p2 <- pnorm( mu2, colMeans(y, na.rm=T), sd0);
			}
			
			p <- c(p0,p1,p2);
			#if (length(which(p<0.05 | p> 0.95) )==0)
			#{
				loop <- loop + 1;
				if (is.na(h.ret)) h.ret <- h;
				if (h.ret$value > h$value) h.ret <- h;
			#}
			
cat("LOOP=", loop, h$value, h$par, "\n");
		}			 

		parin <- parin * runif( length(parin), 0.95, 1.05 );
		if (parin[1]>1 ) parin[1] <- 0.5;
	}	
		
	return(h.ret);
}

Qtlmodel.get_est_LR2<-function(dat, par.cross, parin)
{
	y <- as.matrix( dat$phenos_table );

	flank_snps <- dat$genos_table[, c(par.cross$lmarker, par.cross$rmarker)]; 
	missing <- which(flank_snps[,1] == -1 | flank_snps[,2] == -1 );

	if ( length(missing) > 0)
	{
		flank_snps <- flank_snps[ -(missing), ];
		y <- y[ -(missing), ];
	}		
	if ( par.cross$lmarker != Qtlmodel.old_geno_pos )
	{
		FM2.set_value( "likelihood_null", NULL );
		Qtlmodel.old_geno_pos <- par.cross$lmarker;
	}
			
	# probability table for QQ,Qq,qq at markers MiMjNiNj;
	# 0: homozygote(qq)
	# 1: heterozygote(Qq)
	# 2: homozygote+additive(QQ)
	allprob <- FM2.cross$get_qtl_prob( flank_snps[,1], flank_snps[,2], par.cross$dist , par.cross$qtl.pos )			
	allprob[which(allprob==0)]<-0.005;
	allprob[which(allprob==1)]<-0.995;

	bSuccess <- TRUE;
	h0 <- FM2.get_value( "likelihood_null" );
	if ( is.null( h0) )
	{
		par <- parin[1:(FM2.covar$get_par_num(dat) + FM2.curve$par_num)];
		h0 <- Qtlmodel.optim( par, FM2.curve$get_mlefunc, y = y, time.std=dat$sample_times );
		if ( any(is.na(h0)) ) 
			bSuccess <- FALSE
		else
			FM2.set_value( "likelihood_null", h0 );
	}
	else 
		if (h0$value==0 ) bSuccess <- FALSE;

	if (bSuccess)
		parin[ 1:length(h0$par) ]<- h0$par;

	h1 <- Qtlmodel.optim( parin, Qtlmodel.mlefunc, y = y, qtl.prob = allprob, time.std=dat$sample_times );
	if ( any(is.na(h1) ) ) 
		bSuccess <- FALSE
	else
	{
		bSuccess <- TRUE
		#cat(h1$par, "\n");
	}

	if (bSuccess)
	{
		#ret.x <- ret[1] * length(phenos[,1])/( length(phenos[,1])-length(missing) );
		return ( c( 2*( h0$value - h1$value ), h1$value, h1$par ) )
	}		
	else
	{
		return ( c( NaN, NaN, rep(NaN,length(parin))) );
	}
}

#--------------------------------------------------------------
# private: fin.get_qtl_scanpoints_grp
#
# Get QTL intervals between two markers. 
#
# input 
#       mrk_dist: distance between two markers.
#           step: step
# output: the vector of the intervals. 
#--------------------------------------------------------------
fin.get_qtl_scanpoints_grp <- function(marker_grp, step=NA)
{
	if (is.na(step) )
		step <- as.numeric( FM2.get_value("scan_step", def=2) );
	qtls <- c();
	
	nLastDist <- 0;
	for ( mrk_idx in 1: (marker_grp$count-1)  )
	{
		mrk_dist  <- marker_grp$dists[[ mrk_idx+1 ]] - marker_grp$dists[[ mrk_idx ]];

		if ( mrk_dist > step )
			qtls_0 <- seq(0, mrk_dist, step)
		else
			qtls_0 <- c(0,mrk_dist);

		if (!(mrk_dist %in% qtls_0))
			qtls_0 <- c(qtls_0, mrk_dist);

		for (i in 1:length(qtls_0))
		{
			nRelDist <- nLastDist + qtls_0[i];
			qtls<-rbind(qtls, c( qtls_0[i], mrk_dist, mrk_idx, nRelDist) );
		}
		nLastDist <- nLastDist + mrk_dist;
	}

	if (length(qtls[,1])>3)
		for(i in seq( length(qtls[,1])-1,2,-1))
		{
			if (qtls[i,2] == qtls[i,1])
				qtls <- qtls[-i,];
		}

	return ( matrix(qtls, ncol=4) );
}

fin.filter.qtl <- function( dat, qtls, marker_grp)
{
	stay <- c();
	for ( nQtl in  1:length(qtls[,1]) )
	{
		y <- as.matrix( dat$phenos_table );
		mrk_idx <- qtls[nQtl,3];
		par.cross <-list(lmarker = marker_grp$start_idx + mrk_idx-1,
					rmarker = marker_grp$start_idx + mrk_idx,
					qtl.pos = qtls[nQtl,1], 
					dist    = qtls[nQtl,2]);
		flank_snps <- dat$genos_table[, c(par.cross$lmarker, par.cross$rmarker)]; 
		missing <- which(flank_snps[,1] == -1 | flank_snps[,2] == -1 );
		
		if ( length(missing) > 0)
		{
			flank_snps <- flank_snps[ -(missing), ];
		}
		allprob <- FM2.cross$get_qtl_prob( flank_snps[,1], flank_snps[,2], par.cross$dist , par.cross$qtl.pos )	;
		allprob[which(allprob<0.6)]<-0;
		allprob[which(allprob>=0.6)]<-1;
		#len.gen <- FM2.cross$gen_num;
		tp1 <- which(allprob[,1]==1);
		tp2 <- which(allprob[,2]==1);
		Type1 <- apply(dat$phenos_table[tp1,], 2, function(x){mean(x, na.rm=T);});
		Type2 <- apply(dat$phenos_table[tp2,], 2, function(x){mean(x, na.rm=T);});
		sample <- apply(dat$phenos_table, 2, function(x){mean(x, na.rm=T);});
		pho1 <- cor(Type1, sample);
		pho2 <- cor(Type2, sample);
		if(abs( pho1 - pho2) >0.0003)
		{
			stay<-c(stay,nQtl);
		}
		
	}
		filterqtl <- qtls[stay,];
		return(stay);
			
}

#--------------------------------------------------------------
# private: get_qtl_scanpoints
#
# Get QTL intervals between two markers. 
#
# input 
#       mrk_dist: distance between two markers.
#           step: step
# output: the vector of the intervals. 
#--------------------------------------------------------------
fin.get_qtl_scanpoints <- function( mrk_dist, step=NA)
{
	if (is.na(step) )
		step <- as.numeric( FM2.get_value("scan_step", def=2) );
	
	if ( mrk_dist > step )
		qtls <- seq(0, mrk_dist, step)
	else
		qtls <- c(0,mrk_dist);
	
	if (!(mrk_dist %in% qtls))
		qtls<-c(qtls, mrk_dist);
		
	return (qtls);
}

fin.qtl_locate <- function( dat_obj, res)
{
	nPeakCount <- as.numeric( FM2.get_value("peak_count") );
	
	res <- fin.remove_invalid_row(res);
	if ( length( res[,1]) >0 )
	{
		nPeaks <- fin.select_peaks_by_simple_way(res[,3]);
		if ( length(nPeaks) > 0 )
		{
			nPeaks5 <-c();
			chromo5 <-c();
			for (i in 1:length(nPeaks))
			{
				if (! (res[nPeaks[i],1] %in% chromo5) )
				{
					chromo5 <- c( chromo5, res[nPeaks[i],1] );
					nPeaks5 <- c( nPeaks5, nPeaks[i] );
					if (length(chromo5)>=nPeakCount)
						break;
				}	
			}

			FM_sys$set_LR_peaks( matrix( res[nPeaks5,c(1,2,3)], ncol=3) );
			r <- which.max(res[,3]);

			wlen <- length(res[1,]);
			ret<-list(
				name      = "Qtlmodel.ret",
				ht_no     = 10,
				success   = TRUE,
				cross_type= dat_obj$cross_type,
				covar_type= dat_obj$covar_type,
				curve_name= dat_obj$curve_name,
				qtl_pos   = res[r,2],
				qtl_grp   = res[r,1],
				qtl_LR    = res[r,3],
				qtl_pvalue= 1- pchisq( res[r,3], 1 ),
				LR_peaks  = matrix( res[nPeaks5,c(1:wlen)], ncol=wlen),
				curve_par = c( res[r,c(5:wlen)] ) ,
				full_res  = res);
				class(ret) <- "FM2.Qtlmodel.ret";

				return(ret);
		}				
	}
	ret<-list(
		name       = "Qtlmodel.ret",
		ht_no      = 10,
		success	   = FALSE ,
		cross_type = dat_obj$cross_type,
		curve_name = dat_obj$curve_name,
		full_res   = res );
	class(ret) <- "FM2.Qtlmodel.ret";

	return (ret);
}

#--------------------------------------------------------------
# Qtlmodel.summary_ret
# 
# Summarize the result object of the hypothesis test. 
# Used by summary( Qtlmodel.ret object );
#--------------------------------------------------------------
Qtlmodel.summary_ret<-function( res10, dat_obj )
{
	if (!is.null(res10) && res10$name != "Qtlmodel.ret")
		stop("Not a result for hypothesis test 10");

	if (!res10$success )
	{
		str<- paste( "Failed to do hypothesis test(10).\n\n", sep="" )		
		return(str);
	}
	
	str <- "";
	str<- paste( str, "Hypothesis test 10: \n", sep="" );	
	str<- paste( str, FM2.cross$hp_desc, "\n", sep="" )
	
	str<- paste( str, "------------------------------------\n", sep="" );

	str0 <- sprintf("%15s: %s\n", 	 "Curve",	dat_obj$curve_name );
	str1 <- sprintf("%15s: %s\n", 	 "Cross", 	FM2.cross$name );
	str2 <- sprintf("%15s: %s\n", 	 "Covariance", 	FM2.covar$name );
	str <- paste( str, str0, str1, str2, sep="" );

	st0 <- sprintf( "%15s: %-5.1f (Group:%d)\n", "QTL pos.", res10$qtl_pos, res10$qtl_grp );
	str <- paste(str, st0, sep="");
	st0 <- sprintf( "%15s: %-8.3f \n", "QTL LR", res10$qtl_LR);
	str <- paste(str, st0, sep="");
	st0 <- sprintf( "%15s: %-8.3f \n", "QTL p-value", res10$qtl_pvalue);
	str <- paste(str, st0, sep="");

	st0<- sprintf("%15s: %-8.3f\n", "rho", 	res10$curve_par[1] );
	str <- paste(str, st0, sep="");
	st0<- sprintf("%15s: %-8.3f\n", "sigma2", res10$curve_par[2] );
	str <- paste(str, st0, sep="");

	wlen <- (length( res10$curve_par )-2 )/FM2.cross$gen_num;
	nstart <- 2;
	if(FM2.cross$gen_QQ)
	{
		mu_QQ_str <- fin.format_param( res10$curve_par[(nstart+1):(wlen+nstart)] );
		st0<- sprintf("%15s: mu=%s\n", "QTL para(QQ)", mu_QQ_str);
		str <- paste(str, st0, sep="");
		nstart <- nstart+wlen;
	}
	if(FM2.cross$gen_Qq)
	{
		mu_Qq_str <- fin.format_param( res10$curve_par[(nstart+1):(wlen+nstart)] );
		st0<- sprintf("%15s: mu=%s\n", "QTL para(Qq)", mu_Qq_str );
		str <- paste(str, st0, sep="");
		nstart <- nstart+wlen;
	}	
	if(FM2.cross$gen_qq)
	{
		mu_qq_str <- fin.format_param( res10$curve_par[(nstart+1):(wlen+nstart)] );
		st0<- sprintf("%15s: mu=%s\n", "QTL para(qq)", mu_qq_str );
		str <- paste(str, st0, sep="");
		nstart <- nstart+wlen;
	}

	str<- paste( str, "------------------------------------\n", sep="" );		
	str<- paste( str, "\n", sep="" );	

	st0 <- sprintf(" No.  Grp      Pos.       LR        \n" )
	str <- paste(str, st0, sep="");

	for(i in 1:length(res10$LR_peaks[,1]))
	{
		st0 <- fin.format_param( res10$LR_peaks[i,c(1:4)]);
		str <- paste(str, st0, sep="\n");
	}
	str <- paste(str, "", sep="\n");
	browser()
	figlist <- Qtlmodel.plot_ret(res10, dat_obj, bPrintOut=FALSE);
	browser()
	for (i in 1:length(figlist[,1]))
	{
		str <- paste(str,"*", i, " The figure(",figlist[i,1],") for the hypothesis test(10) is saved to ",  figlist[i,2], ".\n", sep="");
	}
	str <- paste( str, "\n", sep="" );		

	return(str);
}

#--------------------------------------------------------------
# Qtlmodel.report_ret
# 
# Summarize the result object of the hypothesis test. 
# Used by summary( NP.ret object );
#--------------------------------------------------------------
Qtlmodel.report_ret<-function( res10, dat )
{
	if (!is.null(res10) && res10$name != "Qtlmodel.ret")
		stop("Not a result for hypothesis test 10");
	
	if (!res10$success )
		stop("Failed to do hypothesis test(10).");

	len.par <- length( res10$LR_peaks[1,]);
	ret<-c();
	if (dat_obj$cross_type==CROSS_BC || dat_obj$cross_type==CROSS_RIL)
	{
		len.mu  <- (len.par-6) / 2;
		cols.mu <- c();
		for( i in 1:len.mu)
			cols.mu <- c(cols.mu, paste("mu", i, sep=""));
		ret<-res10$LR_peaks[,-4, drop=FALSE];
		colnames(ret) <- c("Grp.", "Pos.",  "LR", "rho", "s2", cols.mu, cols.mu);
	}
	else
	{
		len.mu  <- (len.par-6) / 3;
		cols.mu <- c();
		for( i in 1:len.mu)
			cols.mu <- c(cols.mu, paste("mu", i, sep=""));
		ret<-res10$LR_peaks[,-4, drop=FALSE];
		colnames(ret) <- c("Grp.", "Pos.",  "LR", "rho", "s2", cols.mu, cols.mu, cols.mu);
	}

	plotdoc.old <- FM2.get_value("plot_doctype");
	FM2.set_value("plot_doctype", "png");
	figlist <- Qtlmodel.plot_ret(res10, dat, bPrintOut=FALSE);
	FM2.set_value("plot_doctype", plotdoc.old);

	return(list(sum=ret, fig=figlist));
}

#--------------------------------------------------------------
# Qtlmodel.plot_ret
# 
# Plot the figure for the hypothesis test. 
# Used by plot( NP.ret object );
#--------------------------------------------------------------
Qtlmodel.plot_ret<-function( res10, dat_obj, bPrintOut=TRUE  )
{
	browser()
	if (!is.null(res10) && res10$name != "Qtlmodel.ret")
		stop("Not a result for hypothesis test 10");
	
	if (!res10$success )
		stop( "Failed to do hypothesis test(10).");		
	
	figlist <- c();

	#### Graph 1/3 ###########
	if (!bPrintOut) strFile <- fpt.plot_to_doc( dat_obj$pheno_file, FM2.get_value("plot_doctype") ) else X11();
  	
  	err.fig<-try ( fpt.plot_qtl_map( dat_obj, res10$full_res ) );
	if (class(err.fig)!="try-error")
		title("The LR profile for all chromosomes");

	if (!bPrintOut)
	{
		dev.off();
		figlist<- rbind(figlist, c("QTL profile",strFile));
	}
	cat("++++++++++++++++++++++++++++++++++");
	browser()
	for (i in 1:length(res10$LR_peaks[,1]))
	{
		grp <- res10$LR_peaks[i,1];
		qtls <- res10$LR_peaks[which( res10$LR_peaks[,1]==grp ), 2]; 

		#### Graph 2/3 ###########
		if (!bPrintOut) strFile <- fpt.plot_to_doc( dat_obj$pheno_file, FM2.get_value("plot_doctype") ) else X11();
		err.fig <- try ( fpt.plot_qtl_pos( res10$LR_peaks[i,1], dat_obj, res10$full_res, qtl_ps = c(qtls)  ) );
		if (class(err.fig)!="try-error")
			title(paste("The LR profile for QTL position(Group:", grp, ")", sep="") );
		
		if (!bPrintOut)
		{
			dev.off();
			figlist<- rbind(figlist, c("QTL position",strFile));
		}

		#### Graph 3/3 ###########
		if (!bPrintOut) strFile <- fpt.plot_to_doc( dat_obj$pheno_file, FM2.get_value("plot_doctype") ) else X11();
		nLen <- length( res10$LR_peaks[i,] )
		
		par0<-FM2.cross$get_gen_par(res10$LR_peaks[i,7:nLen])
		err.fig <- try ( fpt.plot_com_curve ( max(dat_obj$sample_times), 
								max(dat_obj$sample_times) + 4,
								f_curve_mu = FM2.curve$get_mu,
								dat    = dat_obj,
								QQ_par = par0$QQ, 
								Qq_par = par0$Qq, 
								qq_par = par0$qq ) );
	
		if (class(err.fig)!="try-error")
			title(paste("Curves(Group:", res10$LR_peaks[i,1], ",Pos.:",res10$LR_peaks[i,2], ")", sep="") );
		
		if (!bPrintOut)
		{
			dev.off();
			figlist<- rbind(figlist, c("QTL curve",strFile));
		}
	}
	cat("_________________________");
	browser()
	return (figlist);
}

#--------------------------------------------------------------
# public: fin.select_peaks_by_simple_way
#
#--------------------------------------------------------------
fin.select_peaks_by_simple_way<-function( curve )
{
	curv2 <- curve[2:length(curve)];
	diff  <- curv2 - curve[1:(length(curve)-1)];
	
	peaks1<-c();
	for (i in 2:length(diff) )
	{
		if ( diff[i-1]>0 && diff[i]<0 ) 
			peaks1 <- c(peaks1, i);
	}

	index <- sort(curve[peaks1], index.return=TRUE )$ix;
	peak_index <- c(index[1]);

	x <- curve[peaks1];
	o <- order(x, decreasing = TRUE);
	return ( peaks1[o] );
}	

#--------------------------------------------------------------
# public: fin.remove_invalid_row
#
#--------------------------------------------------------------
fin.remove_invalid_row<-function( res )
{
	for(i in 1:length(res[1,]) )
	{
		Nans <- which( is.nan(res[,i] ) );
		if (length(Nans)>0)
			res <- res[(-Nans),,drop=FALSE];
	}
	
	return (res);
}

#--------------------------------------------------------------
# public: fin.format_param
#
#--------------------------------------------------------------
fin.format_param<-function(par)
{
	st0<-"";
	for (i in 1:length(par))
		st0 <- paste( st0, sprintf("%8.3f ", par[i]) ) ;
	
	return(st0);
}

#--------------------------------------------------------------
# Qtlmodel.set_lr2_cutoff
#
# 
#--------------------------------------------------------------
Qtlmodel.set_lr2_cutoff<-function( res, p05, p01)
{
	newres <- res;
	newres$qtl_p05<- p05;
	newres$qtl_p01<- p01;
	if (newres$qtl_pvalue<p05)
		newres$qtl_pvalue<- ">0.05"
	else
	if (newres$qtl_pvalue<p01)
		newres$qtl_pvalue<- "<0.05"
	else
		newres$qtl_pvalue<- "<0.01"
	
	return (newres);
}

#--------------------------------------------------------------
# summary.FM2.Qtlmodel.ret
#
# used by summary( res_obj, dat_obj )
#--------------------------------------------------------------
summary.FM2.Qtlmodel.ret<-function( res_obj, dat_obj, file=NA , append = TRUE  )
{
	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	str <- Qtlmodel.summary_ret( res_obj, dat_obj )
	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) );
}

#--------------------------------------------------------------
# report.FM2.Qtlmodel.ret
#
# used by summary( FM2.ret.hp.obj )
#--------------------------------------------------------------
report.FM2.Qtlmodel.ret<-function( res_obj, dat_obj )
{
	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	rptlist<- Qtlmodel.report_ret(res_obj, dat_obj );
	return(rptlist);
}

#--------------------------------------------------------------
# plot.FM2.Qtlmodel.ret
#
# used by plot( FM2.ret.hp.obj )
#--------------------------------------------------------------
plot.FM2.Qtlmodel.ret<-function( res_obj, dat_obj )
{
	if (is.null(FM2.curve))
	{
		stop("inavlid curve type");
	}

	
	Qtlmodel.plot_ret( res_obj, dat_obj );
	invisible();
}

class_qtlmodel<<-list(
	name 	       = "QTL Model",
	get_est_param  = Qtlmodel.get_est_param,
	get_est_LR2    = Qtlmodel.get_est_LR2,
	qtlscan        = Qtlmodel.qtlscan,
	set_lr2_cutoff = Qtlmodel.set_lr2_cutoff);

FM2.model <<- class_qtlmodel;