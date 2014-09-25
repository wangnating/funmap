#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Autoregressive model(1)
#  
#  Variance : AR1
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

AR1.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);
        phi<--log(par[1]);
	m <- length( times );
	Ar.1 <- array(0, dim=c(m,m));
	for (k0 in 1:m)
		for (k1 in 1:m)
			Ar.1[k0,k1] <- par[2]^2 * exp(-phi*abs( times[k0] - times[k1] ));
			
	return(Ar.1);
}

AR1.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2]; 
	dat.t0<-dat$phenos_table[ , m ] ;
	par.sigma <- sd(dat.t0, na.rm=T);

	dat.t1<-dat$phenos_table[ , m-1 ];
	par.rho <- cor(dat.t1, dat.t0, use="complete.obs" );
	
	
	return( c( par.rho, par.sigma));
}

AR1.is_valid<-function(par)
{
	return( par[1] >= 0 && par[1] <= +1&&par[2]>=0);
}

AR1.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.rho",  covar$par[1] );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par[2] );
	
	return(paste(str0, str1, str2, sep=""));
}

AR1.get_par_num<-function(dat)
{
	return(2);
}

AR1.guess_simu_sigma<-function( sim.mu, prob, H2 )
{
	m <- dim(dat$phenos_table)[2];

	s1<- sim.mu[1, ] - sim.mu[2, ]; 
	s2<- sim.mu[2, ] - sim.mu[3, ];
	phe.max <- sim.mu[, which.max(s1*s1+s2*s2) ];

	mu <- (sim.mu[1,] + sim.mu[3,])/2
	a  <- phe.max[1,] - mu;
	d  <- phe.max[2,] - mu;
	q <- 1-prob;

	alpha <- a+(q-prob)*d;
	sig.a <- 2*prob*q*alpha;
	var_a <- 2*prob*q*alpha^2;
	var_d <- (2*prob*q*d)^2;
	var_g <- var_a + var_d;
	var_e <- (1/H2-1)*var_g;

	return(sqrt(var_e));
}

covar_AR1<-list(
	name 	= "AR1",
	desc 	= "parametric stationary autoregressive(AR) model",
	get_est_param= covar.get_est_param,
	get_par_num  = AR1.get_par_num, 
	get_summary  = AR1.get_summary,
	is_valid     = AR1.is_valid,
	get_mat      = AR1.get_mat,
	get_init_rand= AR1.get_init_rand,
	guess_simu_sigma  = AR1.guess_simu_sigma);

ZZZ.regcovar_ar1<-function()
{
	COVAR_AR1  <<- FM2.reg_covar( covar_AR1 );
}
