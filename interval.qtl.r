#################################################################
#
# New Systems Mapping Application(SysMap1)#
# Common library
#
#    1) interqtl.mlefunc
#    2) interqtl.get_est_param
#    3) interqtl.optim
#    4) interqtl.get_est_LR2
#    5) interqtl.qtlscan
#    6) interqtl.summary_ret
#
#
# History:
# 12/15/2011 Version 1.1
#
###################################################################

#--------------------------------------------------------------
#  interqtl.mlefunc
#
#  MLE function, a function to be maximized 
#  
#--------------------------------------------------------------
interqtl.mlefunc<-function( parin, y, qtl.prob, par_fixed=NULL)
{
	len.y <- length(y);
	
	yy0<- NULL;
	yy1<- NULL;
	yy2<- NULL;
	
	sd <- sd(y, na.rm=T);
	mu <- mean(y, na.rm=T);
	
	len <- 1;
	
	if(FM2.cross$gen_QQ)
	{
		mu2 <- parin[(len+1)];
		len <- len+1;
		yy2 <- y - mu2;	
		p0 <- pnorm( mu2, mu, sd);	
	}
	
	if(FM2.cross$gen_Qq)
	{
		mu1 <- parin[(len+1)];
		len <- len+1;
		yy1 <- y- mu1;
		p0 <- pnorm( mu1, mu, sd);
	}
	
	if(FM2.cross$gen_qq)
	{
		mu0 <- parin[(len+1)];
		len <- len+1
		yy0 <- y - mu0;
		
		p0 <- pnorm( mu0, mu, sd);
	}
	
	for( i in 1:length(y))
	{
		y.miss <- which( is.na(y));
		if(!is.null(yy0)) yy0[y.miss]<-0;
		if(!is.null(yy1)) yy1[y.miss]<-0;
		if(!is.null(yy2)) yy2[y.miss]<-0;
	}
	
	fy0 <-NULL;
	fy1 <-NULL;
	fy2 <-NULL;
	
	sdd <- parin[1];
	
	if(!is.null(yy0))
	{
		fy0 <- try( dnorm( yy0, 0, sdd),T);
		if (class(fy0)=="try-error") return (NaN);
		
	}
	
	if(!is.null(yy1))
	{
		fy1 <- try( dnorm( yy1, 0, sdd),T);
		if (class(fy1)=="try-error") return (NaN);
		
	}
	
	if(!is.null(yy2))
	{
		fy2 <- try( dnorm( yy2, 0, sdd),T);
		if (class(fy2)=="try-error") return (NaN);
	}
	pf <- cbind(fy0, fy1, fy2)*qtl.prob;
	A <- -sum(log(rowSums(pf)));
	
	#FM2.set_value( "mle.pdf", length(which( pf>1 | pf<0 ) ) );
#cat(A, parin, "\n");

	return (A);
}

#--------------------------------------------------------------
#  interqtl.get_est_param
#
#  Used by estimating parameters 
#  
#--------------------------------------------------------------
interqtl.get_est_param <- function( dat )
{
	y <- dat$phenos_table ;
	
	sd <- sd(y, na.rm=T);
	mu <- mean(y, na.rm=T);
	parin <- c( sd, mu);
	
	for(i in 2:FM2.cross$gen_num)
		parin <- c(parin, mu*runif(1, min=0.99, max=1.01));
		
	return(parin);
}

#--------------------------------------------------------------
#  interqtl.get_est_param
#
#  Used by estimating parameters 
#  
#--------------------------------------------------------------
interqtl.optim <- function(parin, mlefunc, y, qtl.prob =NULL )
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
				try( h<- optim( parin, mlefunc, y =y, method =g_optim ), F)
			else if (g_optim=="CE")
				h<- optim_ce( parin, mlefunc, dat =dat);
		}
		else
		{
			if (g_optim=="BFGS" || g_optim=="Nelder-Mead" || g_optim=="SANN")
				try( h<- optim( parin, mlefunc, y =y, qtl.prob=qtl.prob, method =g_optim ), F)
			else if (g_optim=="CE")
				h<- optim_ce( parin, mlefunc, y =y, qtl.prob=qtl.prob );
		}
		if(class(h) != "try-error" &&  !is.na(h) && h$value>0)
		{			
			sd <- sd(y, na.rm=T);
			mu <- mean(y, na.rm=T);
			
			mu2 <- h$par[2];
			p2 <- pnorm( mu2, mu, sd);
			
			p1 <- c();
			if(length(h$par)>2)
			{
				mu1 <- h$par[3];
				p1 <- pnorm( mu1, mu, sd);
			}
			
			p0 <- c();
			if(length(h$par)>3)
			{
				mu0 <- h$par(4);
				p0 <- pnorm( mu0, mu, sd);
			}
			
			p <- c(p0, p1, p2);
			#if(length(which(p<0.05 | p> 0.95))==0)
			
				loop <- loop+1;
				if(is.na(h.ret)) h.ret <- h;
				if(h.ret$value > h$value )  h.ret <- h;
			
cat("LOOP=", loop, h$value, h$par, "\n");
		}
			parin <- parin * runif( length(parin), 0.95, 1.05);
	}
	
	return(h.ret);
}
#--------------------------------------------------------------
#  interqtl.get_est_param
#
#  Used by estimating parameters 
#  
#--------------------------------------------------------------
interqtl.get_est_LR2 <- function(dat, par.cross, parin)
{
	
	y <- unlist( dat$phenos_table );
	
	# current two markers in current inteval
	flank_snps <- dat$genos_table[, c(par.cross$lmarker, par.cross$rmarker)]; 
	missing <- which(flank_snps[,1] == -1 | flank_snps[,2] == -1 );
	
	if ( length(missing) > 0)
	{
		flank_snps <- flank_snps[ -(missing), ];
		y <- y[ -(missing)];
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
		par <- parin[1:2];
		h0 <- interqtl.optim( par, interval.mlefunc, y = y );  
		if ( any(is.na(h0)) ) 
			bSuccess <- FALSE
		else
			FM2.set_value( "likelihood_null", h0 );
	}
	else
		if (h0$value==0 ) bSuccess <- FALSE;
	
	if (bSuccess)
		parin[ 1:length(h0$par) ]<- h0$par;
	h1 <- interqtl.optim( parin, interqtl.mlefunc ,y =y, qtl.prob = allprob);	
	if ( any(is.na(h1) ) ) 
		bSuccess <- FALSE
	
	if (bSuccess)
		return ( c( 2*( h0$value - h1$value ), h1$value, h1$par ) )
	else
		return ( c( NaN, NaN, rep(NaN,length(parin))));
}

#--------------------------------------------------------------
#  interqtl.get_est_param
#
#  Used by estimating parameters 
#  
#--------------------------------------------------------------
interqtl.qtlscan <- function(dat)
{
	FM2.set_value("likelihood_null",NULL);
	phenos <- unlist( dat$phenos_table );
	parin <- interqtl.get_est_param(dat);
	
cat("INIT PARAMS:", parin, "\n");
	
	marker_obj <- dat$marker_obj;
	FM_sys$task_start( "Execute the hypothesis test ", FM2.cross$hp_desc, "...\n" );
	interqtl.old_geno_pos <<- -1;
	
	res<-c();
	
	for ( grp_idx in 1:marker_obj$count )
	{
		marker_grp = marker_obj$grps[[grp_idx]];
		qtls <- fin.get_qtl_scanpoints_grp(marker_grp); 
		
		# search for all marker for current group(chromosome)
		for ( nQtl in  1:length(qtls[,1]) )
		{
			mrk_idx <- qtls[nQtl,3];
			par.cross <-list(lmarker = marker_grp$start_idx + mrk_idx-1,
					 rmarker = marker_grp$start_idx + mrk_idx,
					 qtl.pos = qtls[nQtl,1], 
					 dist    = qtls[nQtl,2])
		
			ret <- interqtl.get_est_LR2(dat, par.cross, parin);
			cat( "QTLSCAN:", grp_idx, mrk_idx, qtls[nQtl,4], ret, "\n");
			res <- rbind(res, c( chr_grp=grp_idx, pos = qtls[nQtl,4] , ret) );
		}
	        FM_sys$task_elapsed(finished = grp_idx/marker_obj$count, "Group:", grp_idx, "/", marker_obj$count, ",", "$SYS_PROMPT$","\n" );
	 
	  }
	  FM_sys$task_stop( "The hypothesis test is done\n");
	  
	  ret <- fin.interqtl_locate( dat, res);  
	  
	  return( ret );
}

#--------------------------------------------------------------
# interqtl.summary_ret
# 
# Summarize the result object of the hypothesis test. 
# Used by summary( interqtl.ret object );
#--------------------------------------------------------------
interqtl.summary_ret <- function(res10, dat_obj)
{
	if (!is.null(res10) && res10$name != "interqtl.ret")
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
	
	str0 <- sprintf("%15s: %s\n", 	 "Cross", 	FM2.cross$name );
	str <- paste( str, str0, sep="" );
	
	st0 <- sprintf( "%15s: %-5.1f (Group:%d)\n", "QTL pos.", res10$qtl_pos, res10$qtl_grp );
	str <- paste(str, st0, sep="");
	st0 <- sprintf( "%15s: %-8.3f \n", "QTL LR", res10$qtl_LR);
	str <- paste(str, st0, sep="");
	st0 <- sprintf( "%15s: %-8.3f \n", "QTL p-value", res10$qtl_pvalue);
	str <- paste(str, st0, sep="");
	
	len <- length(res10$par);
	st0<- sprintf("%15s: %-8.3f\n", "sigma", res10$par[1] );
	str <- paste(str, st0, sep="");
	
	nstart <- 1;
	
	if(FM2.cross$gen_QQ)
	{
		mu_QQ_str <- fin.format_param( res10$par[(nstart+1)] );
		st0<- sprintf("%15s: mu=%s\n", "QTL para(QQ)", mu_QQ_str);
		str <- paste(str, st0, sep="");
		nstart <- nstart+1;
	}
	
	if(FM2.cross$gen_Qq)
	{
		mu_Qq_str <- fin.format_param( res10$par[(nstart+1)] );
		st0<- sprintf("%15s: mu=%s\n", "QTL para(Qq)", mu_Qq_str );
		str <- paste(str, st0, sep="");
		nstart <- nstart+1;
	}

	if(FM2.cross$gen_qq)
	{
		mu_qq_str <- fin.format_param( res10$par[(nstart+1)] );
		st0<- sprintf("%15s: mu=%s\n", "QTL para(qq)", mu_qq_str );
		str <- paste(str, st0, sep="");
		nstart <- nstart+1;
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
	
	return(str);
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

fin.interqtl_locate <- function(dat_obj, res)
{
	nPeakCount <- as.numeric( FM2.get_value("peak_count") );
	
	res <- fin.remove_invalid_row(res);
	if ( length( res[,1]) >0 )
	{
		nPeaks <- fin.select_peaks_by_simple_way(res[,3]);
		if(length(nPeaks) > 0)
		{
			nPeaks5 <-c();
			chromo5 <-c();
			for (i in 1:length(nPeaks))
			{
				if (! (res[nPeaks[i],1] %in% chromo5) )
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
				name      = "interqtl.ret",
				ht_no     = 10,
				success   = TRUE,
				qtl_pos   = res[r,2],
				qtl_grp   = res[r,1],
				qtl_LR    = res[r,3],
				qtl_pvalue= 1- pchisq( res[r,3], 1 ),
				LR_peaks  = matrix( res[nPeaks5,c(1:wlen)], ncol=wlen),
				par = c( res[r,c(5:wlen)] ) ,       
				full_res  = res);
				class(ret) <- "interqtl.ret";
			
				return(ret);
		}
	}
	ret<-list(
		name       = "interqtl.ret",
		ht_no      = 10,
		success	   = FALSE ,
		full_res   = res );
		class(ret) <- "interqtl.ret";
	
		return (ret);
}

#--------------------------------------------------------------
# interqtl.set_lr2_cutoff
#
# 
#--------------------------------------------------------------
interqtl.set_lr2_cutoff <- function( res, p05, p01)
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
# summary.interqtl.ret
#
# used by summary( res_obj, dat_obj )
#--------------------------------------------------------------
summary.interqtl.ret<-function( res_obj, dat_obj, file=NA , append = TRUE  )
{
	str <- interqtl.summary_ret( res_obj, dat_obj )
	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) );
}

class_interqtl<<-list(
	name 	       = "INTERQTL Model",
	get_est_param  = interqtl.get_est_param,
	get_est_LR2    = interqtl.get_est_LR2,
	qtlscan        = interqtl.qtlscan,
	set_lr2_cutoff = interqtl.set_lr2_cutoff );
	
INTER.model <<- class_interqtl;			