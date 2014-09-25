#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  
#  
#  Variance : Scaled Identity
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

SI.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	t_len <- length( times );
	
	SI <-par[1]^2*diag(t_len)
	
	return(SI);
}
SI.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	dat.t0<-dat$phenos_table[ , m ] ;
	par.s2 <- sd(dat.t0, na.rm=T);
        
	return( par.s2);
}
SI.is_valid<-function(par)
{
	return(TRUE);
}


SI.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.rho",  0 );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par );
	
	return(paste(str0, str1, str2, sep=""));
}


SI.get_par_num<-function(dat)
{
	return(1);
}

SI.guess_simu_sigma<-function( sim.mu, prob, H2 )
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

covar_SI<-list(
	name 	= "SI",
	desc 	= "Scaled Identity model",
	get_est_param= covar.get_est_param,
	get_par_num  = SI.get_par_num, 
	is_valid     = SI.is_valid,
	get_mat      = SI.get_mat,
	get_init_rand= SI.get_init_rand,
	get_summary  = SI.get_summary,
	guess_simu_sigma  = SI.guess_simu_sigma);

ZZZ.regcovar_si<-function()
{
	COVAR_SI  <<- FM2.reg_covar( covar_SI );
}


