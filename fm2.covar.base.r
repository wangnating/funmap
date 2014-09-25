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
covar.mlefunc<-function( par, y, time.std, curve.par)
{
	if (!FM2.covar$is_valid(par))
		return(NaN);
	
	cov <- FM2.covar$get_mat(par, 1:length(time.std) );

	m  <- length(y[1,]);
	n  <- length(y[,1]);

	mu0 <- FM2.curve$get_mu( curve.par, time.std );
	yy0 <- y - matrix( rep(1,n),byrow=n )%*%matrix( mu0, nrow=1 ); 
	
	for ( i in 1:length(y[1,]) )
	{
		y.miss <- which( is.na(y[,i]) );
		if( length(y.miss)>0 ) yy0[y.miss, i]<-0;
	}

	hh <- try( fy0 <- dmvnorm(yy0, rep(0, m), cov), T );
	if (class(hh)=="try-error") return(NaN);

	FM2.set_value( "mle.pdf", length(which( fy0>1 | fy0<0 ) ) );

	A <- -sum(log(fy0));

	return (A);
}


covar.get_est_param<-function(dat, curve.par)
{
	phenos <- as.matrix( dat$phenos_table );
	par.num <- FM2.covar$get_par_num(dat);

	val <- Inf;
	loop <-0;
	par <- FM2.covar$get_init_rand(dat);

cat("COVAR.INIT=", par, "\n");		
		
	while( loop < 20 )
	{
		if(loop<=10)
			parin <- runif(par.num) * 2 * par
		else	
			parin <- runif(par.num, 0.8, 1.2 )* par;
		
		r0 <- NA;
cat(".");
		try( r0<- optim( parin, covar.mlefunc, y = phenos, time.std=dat$sample_times, curve.par=curve.par$par,  method ="Nelder-Mead" ), T);
		if (class(r0)=="try-error" || is.na(r0) )
			next;
		
		if(any(is.na(r0)))
			next;
		
		if (r0$value < val )
		{
			par <- r0$par;
			val <- r0$value;
		}
cat("o");
		loop <- loop + 1;
	}
	
	cat("\nCovariance Parameters:", "log(L)=", val, "PAR=", par, "\n");
	return( list(par=par, val=val, method=~"Likelihood"));
}

class_covar<<-list(
	name 	 	= "Base Covar",
	get_est_param 	= covar.get_est_param,
	get_mlefunc     = covar.mlefunc);

