#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  SAD1
#  
#  Variance : SAD1
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SAD1.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	phi   <- par[1];
	sigma <- par[2]
	m <- length( times );

	cov.1 <- array(0, dim=c( m , m ));
	for (k0 in 1:m)
	for (k1 in k0:m)
	{
		w <- sigma^2 * phi^abs(times[k1]-times[k0]) * (1-phi^(2*times[k0]))/(1-phi^2)				
		cov.1[k0, k1] <- w
		cov.1[k1, k0] <- w
	}
        
	return(cov.1);
}

SAD1.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	dat.t0<-dat$phenos_table[ , 1] ;
	par.s2 <- sd(dat.t0, na.rm=T);

	dat.t1<-dat$phenos_table[ , 2 ];
	par.rho <- cor(dat.t1, dat.t0, use="complete.obs");
	par.phi <- par.rho/par.s2^2;
	
	cat("cov.1",par.phi,par.s2,"\n");
	return( c( par.phi, par.s2));
}

SAD1.is_valid<-function(par)
{
	return( par[1] != 1);
}

SAD1.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.phi",  covar$par[1] );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par[2] );
	
	return(paste(str0, str1, str2, sep=""));
}

SAD1.get_par_num<-function(dat)
{
	return(2);
}

SAD1.guess_simu_sigma<-function( sim.mu, prob, H2 )
{
}

covar_SAD1<-list(
	name 	= "SAD1",
	desc 	= "SAD model",
	get_par_num  = SAD1.get_par_num, 
	is_valid     = SAD1.is_valid,
	get_mat      = SAD1.get_mat,
	get_init_rand= SAD1.get_init_rand,
	get_summary  = SAD1.get_summary,
	get_est_param= covar.get_est_param,
	guess_simu_sigma  = SAD1.guess_simu_sigma);

ZZZ.regcovar_sad1<-function()
{
	COVAR_SAD1  <<- FM2.reg_covar( covar_SAD1 );
}