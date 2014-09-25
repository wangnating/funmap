g.browser<<-0;

#--------------------------------------------------------------
# curve.curve_mle
#
#  Used by permutaion.
#
#	prob: the probility of Qq
#     y : the phenotype * times data
#    par: initial parameter for likelihood
#         including a1,b1,r1,a0,b0,r0,rho,sigma2
#--------------------------------------------------------------
curve.leastsq<-function( par, y, time.std, loop)
{
	n <- dim(y)[1];
	
	curve <- FM2.curve$get_mu( par, time.std  );
	mu0 <- colMeans(y, na.rm=T)
	sd0 <- apply(y, 2, function(x){sd(x, na.rm=T);});
	p0  <- pnorm( curve, mu0, sd0 );
	
	prob0 <- FM2.get_value("fitting.prob", 0.05);
	prob1 <- 1-FM2.get_value("fitting.prob", 0.05);
	unfit.max <- floor(100/loop);

	#if ( length(which( p0 < prob0 | p0 > prob1 ) ) >= unfit.max ) 

	# Optim will try to fit many points, but for only one point, optim can not do it well.
	if ( length(which( p0 < prob0 | p0 > prob1 ) ) ==1 ) 
	{
		#cat(loop, length(which( p0 < prob0 | p0 > prob1 ) ), "\n");
		return(NaN);
	}

	yy0 <- y - matrix( rep(1,n),byrow=n )%*%matrix( curve, nrow=1 ); 
	
	for ( i in 1:length(y[1,]) )
	{
		y.miss <- which( is.na(y[,i]) );
		if( length(y.miss)>0 ) yy0[y.miss, i]<-0;
	}

	A <- sum(yy0^2);
	
	return (A);
}


curve.mlefunc<-function( par, y, time.std)
{
	len.cov <- FM2.covar$get_par_num(y);
	par.covar <- par[1:len.cov];
	if (!FM2.covar$is_valid(par.covar))
		return(NaN);

	cov <- FM2.covar$get_mat(par.covar, 1:length(time.std) );

	m  <- length(y[1,]);
	n  <- length(y[,1]);

	gen_par<- par[(len.cov+1):(len.cov+ FM2.curve$par_num )];
	mu0 <- FM2.curve$get_mu( gen_par, time.std  );
	yy0 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu0, nrow=1 ); 

	sd0 <- apply(y, 2, function(x){sd(x, na.rm=T);});
	p0 <- pnorm( mu0, colMeans(y, na.rm=T), sd0);
	#if(length(which(p0<0.01 || p0>0.99))>=1)
	#	return(NaN);
	
	for ( i in 1:length(y[1,]) )
	{
		y.miss <- which( is.na(y[,i]) );
		if( length(y.miss)>0 ) yy0[y.miss, i]<-0;
	}

	hh <- try( fy0<-  dmvnorm( yy0, rep(0, m), cov, log=T), T );
	if (class(hh)=="try-error") return(NaN);

	A <- -sum( fy0 );

	return (A);
}


curve.get_est_param<-function(dat, f_init_rand)
{
	phenos  <- as.matrix( dat$phenos_table );
	par.num <- FM2.curve$par_num;
	prob0   <- FM2.get_value("fitting.prob", 0.05);

	val <- Inf;
	loop <-0;
	loop.Failed <- 0
	rand.seed <- runif(1);
	rset.seed <- 0;

       
	if ( !is.null(f_init_rand))
		par0 <- f_init_rand(dat, dat$sample_times)
	else	
		par0 <- runif(par.num)*mean(phenos, na.rm=T);

	par <- par0	

cat("CURVE.INIT=", par, "Prob.=", prob0, "\n");		
	while( loop <= 100 )
	{
		if (loop==0)
			parin <- par
		else
		{
			if (round(rand.seed%%3)==0)
				parin <- runif( par.num, 0.98, 1.02) *par
			else	
				parin <- runif( par.num, 1/10+0.8*loop/100, 10 - 8.9*loop/100 )* par
		}
		
		r0 <- NA;
		bFailed <- FALSE;
		optim.method<-ifelse(round(rand.seed%%2)==0, "Nelder-Mead", "BFGS" )
		try( r0<- optim( parin, curve.leastsq, y = phenos, time.std=dat$sample_times, loop=loop, method=optim.method), T);
		if (class(r0)=="try-error" || is.na(r0) )
			bFailed <- TRUE
		else
		{
			curve <- FM2.curve$get_mu( r0$par, dat$sample_times  );
			mu <- colMeans(phenos, na.rm=T)
			sd <- apply(phenos, 2, function(x){sd(x, na.rm=T);});
			p0  <- pnorm( curve, mu, sd );

			prob1 <- 1 - prob0;
			unfit.max <- round( length(dat$sample_times)*(80-loop)/80 );
			if (unfit.max<1) unfit.max <- 1;
			
			if ( length(which( p0 < prob0 | p0 > prob1 ) ) >= unfit.max ) 
				bFailed <- TRUE;
		
		}
		
		if(bFailed)
		{
cat(".");
			rand.seed <- (rand.seed + runif(1, 1, 5) )%%1000;
			loop.Failed <- loop.Failed + 1
			if ( loop.Failed> 50 )
			{
cat("\n");	
				#par <- runif(par.num, 0.95, 1.05)*par0;
				loop <- 1;
				rset.seed <- rset.seed + 1;
				loop.Failed <- 0;
		
				if ( rset.seed > 20 )
				{
					cat("\nFailed to estimate curve!.\n");
					return(list(val=NA, par=NA, method="LeastSq",rset.seed=rset.seed, fitting.prob=prob0));
				}	
			}
			
			set.seed( rand.seed^2 );
			next;
		}
		
		if (r0$value < val )
		{
			par <- r0$par;
			val <- r0$value;
		}
		
		loop.Failed <- 0;
		loop <- loop + 1;
cat("o");
	}
	
	if(1)
	{
cat("\n");	
		curve <- FM2.curve$get_mu( par, dat$sample_times  );
		mu <- colMeans(phenos, na.rm=T)
		sd <- apply(phenos, 2, function(x){sd(x, na.rm=T);});
		p0  <- pnorm( curve, mu, sd );
		r.tmp <- rbind(curve, mu, sd, p0);

		show(r.tmp);
		min.prob <- min(p0, 1-p0);
	}	
	
	cat("\nCurve Parameters:", "val(L)=", val, "PAR=", par, "Prob=", min.prob, "\n");

	return(list(val=val, par=par, method="LeastSq", fitting.prob=prob0, min.prob=min.prob));
}




class_curve<<-list(
	name 	 	= "Base Curve",
	get_est_param 	= curve.get_est_param,
	get_mlefunc     = curve.mlefunc);
