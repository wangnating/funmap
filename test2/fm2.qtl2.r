#--------------------------------------------------------------
# COM2.hp_test
# 
# Hypothesis Test 10 according to the paper.
#
#  input: data object and test methods,
#         such as test=c(10:1x) or test=c(10,12).
# output: result object
#--------------------------------------------------------------
COM2.hp_test<-function( dat_obj, test=c(10) )
{
	cross <- dat_obj$cross_type;
	
	ret<-list();
	for (i in 1:length(test))
	{
		tn <- test[i];

		if (tn <10 || tn>10)
			next;
		
		ret_t <- COM_BC_hp_test2( dat, COM.TEST2 );

		ret[[i]]<- ret_t;
	}
	
	return (ret);
}


#--------------------------------------------------------------
# COM.likelihood
#
#  Used by permutaion.
#
#	prob: the probility of Qq
#     y : the phenotype * times data
#    par: initial parameter for likelihood
#         including a1,b1,r1,a0,b0,r0,rho,sigma2
#--------------------------------------------------------------
COM.Test2_likelihood<-function(par, y, prob, cross)
{
	m  <- length(y[1,]);
	n  <- length(y[,1]);

	rho  <- par[1];
	s2   <- par[2];
	if (rho>0.999 || rho<0 || s2<0.001)
		return(NaN);

 	time.std <- as.numeric( colnames(y) );
	if (cross==CROSS_F2)
	{
		len    <- (length(par)-2)/3;
		QQ2_par<- par[3:(2+len)];
		Qq1_par<- par[(3+len):(2+2*len)];
		qq0_par<- par[(3+len*2):(2+3*len)];
		mu2<- FM2.model$get_mu( QQ2_par, time.std );
		mu1<- FM2.model$get_mu( Qq1_par, time.std );
		mu0<- FM2.model$get_mu( qq0_par, time.std );
	}
	else
	{
		len    <- (length(par)-2)/4;
		qq00_par<- par[3:(2+len)];
		qQ01_par<- par[(3+1*len):(2+2*len)];
		Qq10_par<- par[(3+2*len):(2+3*len)];
		QQ11_par<- par[(3+3*len):(2+4*len)];

		mu11<- FM2.model$get_mu( QQ11_par, time.std );
		mu10<- FM2.model$get_mu( Qq10_par, time.std );
		mu01<- FM2.model$get_mu( qQ01_par, time.std );
		mu00<- FM2.model$get_mu( qq00_par, time.std );
 	}

	
	yy00 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu00, nrow=1 ); 
	yy01 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu01, nrow=1 ); 
	yy10 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu10, nrow=1 ); 
	yy11 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu11, nrow=1 ); 
	if (cross==CROSS_F2)
		yy2 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu2, nrow=1 ); 

	for ( i in 1:length(y[1,]) )
	{
		y.miss <- which( is.na(y[,i]) );
		if (cross==CROSS_F2)
			yy2[y.miss, i]<-0;
		yy00[y.miss, i]<-0;
		yy01[y.miss, i]<-0;
		yy10[y.miss, i]<-0;
		yy11[y.miss, i]<-0;
	}
		

	A <-NaN;
	if ( cross == CROSS_F2 )
	{
		stop("unimplemented");
		A<- sum( log( 
	      		prob[,1]*exp(-((1-rho^2)*yy2[,m]^2+rowSums((yy2[,c(1:(m-1))]-rho*yy2[,c(2:m)])^2 ))/s2/2/(1-rho^2))
				+ prob[,2]*exp(-((1-rho^2)*yy1[,m]^2+rowSums((yy1[,c(1:(m-1))]-rho*yy1[,c(2:m)])^2 ))/s2/2/(1-rho^2)) 
				+ prob[,3]*exp(-((1-rho^2)*yy0[,m]^2+rowSums((yy0[,c(1:(m-1))]-rho*yy0[,c(2:m)])^2 ))/s2/2/(1-rho^2)) )
				-(m-1)*log(1-rho^2)/2-m*log(s2)/2-m/2*log(2*pi));
	}
	else
	{		
		A<- sum( log( prob[,1]*exp(-((1-rho^2)*yy00[,m]^2+rowSums((yy00[,c(1:(m-1))]-rho*yy00[,c(2:m)])^2 ))/2/s2/(1-rho^2)) +
				      prob[,2]*exp(-((1-rho^2)*yy01[,m]^2+rowSums((yy01[,c(1:(m-1))]-rho*yy01[,c(2:m)])^2 ))/2/s2/(1-rho^2)) +
				      prob[,3]*exp(-((1-rho^2)*yy10[,m]^2+rowSums((yy10[,c(1:(m-1))]-rho*yy10[,c(2:m)])^2 ))/2/s2/(1-rho^2)) +
				      prob[,4]*exp(-((1-rho^2)*yy11[,m]^2+rowSums((yy11[,c(1:(m-1))]-rho*yy11[,c(2:m)])^2 ))/2/s2/(1-rho^2))  )
				-(m-1)*log(1-rho^2)/2-m*log(s2)/2-m/2*log(2*pi));
	}

	return (-A);
}

#--------------------------------------------------------------
# COM.test2_likelihood_null
#
#  Used by permutaion. Same as H0 of the hypothesis test 10.
#
# input
#     y : the phenotype * times data
#    par: initial parameter for likelihood
#         including a, b, r, rho, sigma2
#--------------------------------------------------------------
COM.Test2_likelihood_null<-function(par, y )
{
	m  <- length(y[1,]);
	n  <- length(y[,1]);

	rho  <- par[1];
	s2   <- par[2];
	if (rho>0.999 || rho<0 || s2<0.001)
		return(NaN);

	len  <- length(par)-2;
	par<- par[3:(2+len)];

 	time.std <- as.numeric( colnames(y) );
	mu<- FM2.model$get_mu( par, time.std );

	yy1 <- y - matrix(rep(1,n),byrow=n)%*%matrix(mu, nrow=1); 
	for ( i in 1:length(y[1,]) )
	{
		y.miss <- which( is.na(y[,i]));
		yy1[y.miss, i]<-0;
	}

	rho.std <- rho^(time.std[c(2:m)]- time.std[c(1:(m-1))] );

	A <- sum( - ( (1-rho^2)*yy1[,m]^2+rowSums( (yy1[,1:(m-1)] - rho * yy1[,2:m])^2)) /2/s2/(1-rho^2)
		-(m-1)*log(1-rho^2)/2 - m*log(s2)/2 - m/2*log(2*pi) ); ###-10.E5*(rho>0.999)-10.E5*(s2<0.001);

	return(-A)
}


#--------------------------------------------------------------
# private: COM.Test2_EM test object
#
# H0:  a1 ==  a0, b1 == b0, c1 == c0
# H1:  (above)  !=  (above)
#       at t=(1:m)  
#--------------------------------------------------------------
COM.Test2_EM <- function( parin, y, allprob, cross  )
{
	bSuccess <- TRUE;
	
	h0 <- FM2.get_value( "likelihood_null" );
	if ( is.null( h0) )
	{
		if (cross == CROSS_F2 )
			len  <- (length(parin)-2)/3
		else
			len  <- (length(parin)-2)/2;

		try( h0<- optim( parin[1:(2+len)], COM.Test2_likelihood_null, y = y, method ="BFGS" ), TRUE);
		if ( any(is.na(h0)) ) bSuccess <- FALSE;
	}
	else 
		if (h0$value==0 ) bSuccess <- FALSE;

	h1 <- NA;
	try( h1 <- optim( parin, COM.Test2_likelihood, y = y, prob = allprob, cross=cross, method = "BFGS" ), TRUE);
	if ( any(is.na(h1) ) ) bSuccess <- FALSE;

	if (bSuccess)
		return ( c( 2*( h0$value - h1$value ), h1$value, h1$par ) )
	else
	{
		return ( c( NaN, NaN, rep(NaN,length(parin))) );
	}
}

COM.Test2_conclude <- function( dat_obj, res)
{
	nPeakCount <- as.numeric( FM2.get_value("peak_count") );
	
	res <- remove_invalid_row(res);
	if ( length( res[,1]) >0 )
	{
		r <- which.max( res[,5] );
		if (r>1 && r<=length(res[,5]) )
		{
			nPeaks <- select_peaks_by_simple_way(res[,5]);
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
			FM_sys$set_LR_peaks( matrix( res[nPeaks5,c(1:5)], ncol=5) );
			
			wlen <- length(res[1,]);
			ret<-list(
				name      = "COM.ret.T10",
				ht_no     = 10,
				success   = TRUE,
				cross_type= dat_obj$cross_type,
				model_name= dat_obj$model_type,
				np.order  = as.numeric( FM2.get_value("np.order") ),
				qtl_pos1   = res[r,4],
				qtl_grp1   = res[r,3],
				qtl_pos2   = res[r,2],
				qtl_grp2   = res[r,1],
				qtl_LR    = res[r,5],
				qtl_pvalue= 1- pchisq( res[r,5], 1 ),
				LR_peaks  = matrix( res[nPeaks5,c(1:wlen)], ncol=wlen),
				curve_par = c( res[r,c(7:wlen)] ) ,
				full_res  = res);

			return(ret);
	   	}
	}

	ret<-list(
		name       = "COM.ret.T10",
		ht_no      = 10,
		success	   = FALSE ,
		cross_type = dat_obj$cross_type,
		model_name = dat_obj$model_name );

	return (ret);
}

COM.Test2_Pre<-function(parin, y, allprob)
{
	if ( is.null( FM2.get_value( "likelihood_null" ) ) )
	{
		len  <- (length(parin)-2)/2;
		r0<-NA;
		try( r0<- optim( parin[1:(2+len)], COM.Test2_likelihood_null, y = y, method ="BFGS" ), TRUE);
		if ( any(is.na(r0) ) )
		{
			cat("Failed to estimate likelihood NULL hypothesis.\n");					
			r0 <- list(value=0)
		}

		FM2.set_value( "likelihood_null", r0 );
	}
}

COM.TEST2<-list(
	name         ="COM.Test2",
	em_pre       = COM.Test2_Pre,
	em_f         = COM.Test2_EM,
	sum_f        = COM.Test2_conclude,
	chromo_range = TRUE)




#--------------------------------------------------------------
# public: COM_BC_hp_test2
#
# Hypothesis test for backcross based on the paper of Changxing Ma 
# et al.
#
# input 
#     dat : data object
#   ht_obj: hypothesis test object, such as XX.BC.T10, XX.BC.T11,
#           XX.BC.T12....
#--------------------------------------------------------------
COM_BC_hp_test2<-function( dat, ht_obj )
{
	marker_obj <- recalc_marker( dat$marker_table );
	parin <- COM_get_init_value_test2( dat);
	
	FM_sys$task_start( "Execute the hypothesis test ", ht_obj$name, "...\n" );	
	FM2.set_value( "likelihood_null", NULL );
	
	phenos<-as.matrix( dat$phenos_table );
   	res<-c();
	
	old_geno_pos = -1;
	for ( grp_idx in 1:marker_obj$count )
	{
		marker_grp = marker_obj$grps[[grp_idx]];
		qtls <- get_qtl_scanpoints_grp(marker_grp); 

		# search for all marker for current group(chromosome)
		for ( nQtl in  1:length(qtls[,1]) )
		{
			y <- phenos;
			mrk_idx <- qtls[nQtl,3]

			# current two markers in current inteval
			flank_snps <- dat$genos_table[, c( marker_grp$start_idx + mrk_idx-1, marker_grp$start_idx + mrk_idx)]; 
			missing <- which(flank_snps[,1] == -1 | flank_snps[,2] == -1 );

			if ( length(missing) > 0)
			{
				flank_snps <- flank_snps[ -(missing), ];
				y <- phenos[ -(missing), ];
			}
			
			if ( marker_grp$start_idx + mrk_idx-1 != old_geno_pos )
			{
				FM2.set_value( "likelihood_null", NULL );
				old_geno_pos = marker_grp$start_idx + mrk_idx-1;
			}

			if ( !is.null(ht_obj$em_pre) )
				ht_obj$em_pre(parin, y, allprob);

			flank_snps[ which(flank_snps[,1] == 2), 1] <- 0;
			flank_snps[ which(flank_snps[,2] == 2), 2] <- 0;

			# probability table for Qq at markers 1 1, 1 0, 0 1, 0 0;
			# 0: homozygote(qq)
			# 1: heterozygote(Qq)
			theta <- ( qtls[nQtl,1]/100 ) /( qtls[nQtl,2]/100  );
			prob <- c(1, 1-theta, theta, 0); 
			allprob<- t( prob [4 - ( flank_snps[,1] * 2 + flank_snps[,2]) ] );
			
			res_inner <- COM_BC_hp_test2_inner(dat, ht_obj, parin, allprob, grp_idx, qtls[nQtl,4]);
			res <- rbind(res, res_inner);

		}

		FM_sys$task_elapsed(finished = grp_idx/marker_obj$count, "Group:", grp_idx, "/", marker_obj$count, ",", "$SYS_PROMPT$","\n" );
	}

	if (!is.null(ht_obj$colnames) )
		colnames(res)<-c( "chr_grp", "pos", ht_obj$colnames );

	if ( !is.null(ht_obj$sum_f) )
		ret <- ht_obj$sum_f( dat, res );

	FM_sys$task_stop( "The hypothesis test is done\n");	
		
	return( ret );
}


COM_BC_hp_test2_inner<-function( dat, ht_obj, parin, allprob1, grp_idx1, pos1)
{
	marker_obj <- recalc_marker( dat$marker_table );
	phenos<-as.matrix( dat$phenos_table );

	res <- c();
	old_geno_pos = -1;
	for ( grp_idx in 1:marker_obj$count )
	{
		marker_grp = marker_obj$grps[[grp_idx]];
		qtls <- get_qtl_scanpoints_grp(marker_grp); 

		# search for all marker for current group(chromosome)
		for ( nQtl in  1:length(qtls[,1]) )
		{

			y <- phenos;
			mrk_idx <- qtls[nQtl,3]

			# current two markers in current inteval
			flank_snps <- dat$genos_table[, c( marker_grp$start_idx + mrk_idx-1, marker_grp$start_idx + mrk_idx)]; 
			missing <- which(flank_snps[,1] == -1 | flank_snps[,2] == -1 );

			if ( length(missing) > 0)
			{
				flank_snps <- flank_snps[ -(missing), ];
				y <- phenos[ -(missing), ];
			}
			
			if ( marker_grp$start_idx + mrk_idx-1 != old_geno_pos )
			{
				FM2.set_value( "likelihood_null", NULL );
				old_geno_pos = marker_grp$start_idx + mrk_idx-1;
			}

			if ( !is.null(ht_obj$em_pre) )
				ht_obj$em_pre(parin, y, allprob);

			flank_snps[ which(flank_snps[,1] == 2), 1] <- 0;
			flank_snps[ which(flank_snps[,2] == 2), 2] <- 0;

			# probability table for Qq at markers 1 1, 1 0, 0 1, 0 0;
			# 0: homozygote(qq)
			# 1: heterozygote(Qq)
			theta <- ( qtls[nQtl,1]/100 ) /( qtls[nQtl,2]/100  );
			prob <- c(1, 1-theta, theta, 0); 
			allprob2<- t( prob [4 - ( flank_snps[,1] * 2 + flank_snps[,2]) ] );

			qtls_prob <- calc_qtl2_prob(allprob1, allprob2 );

			if (ht_obj$chromo_range || FM_sys$is_LR_peak( grp_idx, qtls[nQtl,4] ) )
			{
				ret <- ht_obj$em_f( parin, y, qtls_prob, dat$cross_type );
				ret.x <- ret[1] * length(phenos[,1])/( length(phenos[,1])-length(missing) );
				ret[1] <- ret.x;
				cat(grp_idx1, pos1, grp_idx, mrk_idx, qtls[nQtl,4], ret, "\n");
			
				res <- rbind(res, c( grp_idx1=grp_idx1, pos1= pos1, chr_grp2=grp_idx, pos2 = qtls[nQtl,4] , ret) );
			}

		}

		FM_sys$task_elapsed(finished = grp_idx/marker_obj$count, "Group:", grp_idx, "/", marker_obj$count, ",", "$SYS_PROMPT$","\n" );
	}


	return(res);

}



#--------------------------------------------------------------
# public: COM_get_init_value_test2
#
#--------------------------------------------------------------
COM_get_init_value_test2<-function( dat, cross = -1 )
{
	if (!is.null(FM2.model$get_estvalue))
	{
		####keep safe call begin
		pin <- FM2.model$get_estvalue( dat )
		####keep safe call end
	}
	else
		pin <- COM_likelihood_estvalue(dat);

	if (cross==0)
		return(pin);

	if ( dat$cross == CROSS_BC || dat$cross == CROSS_RIL)
	{
		parin <- c(pin[1], pin[2]);
		p_len <- length(pin)-2;
		parin <- c(parin, pin[3:(p_len+2)]*0.95 );
		parin <- c(parin, pin[3:(p_len+2)]*0.98 );
		parin <- c(parin, pin[3:(p_len+2)]*1.01 );
		parin <- c(parin, pin[3:(p_len+2)]*1.04 );
	}
	else
	{
		parin <- c(pin[1], pin[2]);
		p_len <- length(pin)-2;
		parin <- c(parin, pin[3:(p_len+2)]*1.02 );
		parin <- c(parin, pin[3:(p_len+2)]*1.00 );
		parin <- c(parin, pin[3:(p_len+2)]*0.98 );
	}
	
	return(parin);
}


calc_qtl2_prob<-function(allprob1, allprob2 )
{
	prob <- c();
	for (i in 1:length(allprob1))
	{
		thet1 <- allprob1[i];
		thet2 <- allprob2[i];

		prob<-rbind( prob, c( (1-thet1)*(1-thet2), thet1*(1-thet2), (1-thet1)*thet2, thet1*thet2 ) );
	}


	return(prob);
}
