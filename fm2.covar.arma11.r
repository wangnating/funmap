#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  ARMA(1,1)
#  
#  Variance : ARMA11
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

ARMA11.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	t_len <- length( times );
	Arma11 <- array(0, dim=c(t_len,t_len));
	phi<--log(par[1]);
	for (k0 in 1:t_len)
		for (k1 in 1:t_len)
		{       if(k0==k1)
			Arma11[k0,k1] <- par[2]^2
			else
			Arma11[k0,k1] <- par[2]^2 *par[3]* exp(-phi*abs( times[k0] - times[k1] ));}
		
	return(Arma11);
	
}

ARMA11.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	dat.t0<-dat$phenos_table[ , m ] ;
	dat.t1<-dat$phenos_table[ , 1 ] ;
	dat.t2<-dat$phenos_table[ , 2 ] ;
	dat.t3<-dat$phenos_table[ , 3 ] ;
	par.s2 <- sd(dat.t0, na.rm=T);

	
	par.rho12 <- cor(dat.t1, dat.t2 );
	par.rho13 <- cor(dat.t1, dat.t3 );
	par.rho <- par.rho13/par.rho12;
	par.phi <- par.rho12/par.rho;
		
	return( c( par.rho, par.s2 ,par.phi));
}

ARMA11.is_valid<-function(par)
{
	return( par[1] >= -1 && par[1] <= +1);
}

ARMA11.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.rho",  covar$par[1] );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par[2] );
	str3 <- sprintf(format, "-Covar.phi",  covar$par[3] );
	
	return(paste(str0, str1, str2, str3,  sep=""));
}

ARMA11.get_par_num<-function(dat)
{
	return(3);
}

ARMA11.guess_simu_sigma<-function( sim.mu, prob, H2 )
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

covar_ARMA11<-list(
	name 	= "ARMA11",
	desc 	= "parametric stationary ARMA(1,1) model",
	get_est_param= covar.get_est_param,
	get_par_num  = ARMA11.get_par_num, 
	is_valid     = ARMA11.is_valid,
	get_mat      = ARMA11.get_mat,
	get_init_rand= ARMA11.get_init_rand,
	get_summary  = ARMA11.get_summary,
	guess_simu_sigma  = ARMA11.guess_simu_sigma);

ZZZ.regcovar_arma11<-function()
{
	COVAR_ARMA11 <<- FM2.reg_covar( covar_ARMA11 );
}