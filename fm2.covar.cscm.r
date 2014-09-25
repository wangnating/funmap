#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Compound Symmetry:Correlation Metric
#  
#  Variance : CSCM  
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

CSCM.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	t_len <- length( times );
	CSCM <- array(0, dim=c(t_len,t_len));
	for (k0 in 1:t_len)
		for (k1 in 1:t_len)
		{
		        if(k0==k1)
		        	CSCM[k0,k1] <- par[2]^2
		        else
				CSCM[k0,k1] <- par[2]^2*par[1];
		}
        
	return(CSCM);
}

CSCM.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	dat.t0<-dat$phenos_table[ , m ] ;
	par.s2 <- sd(dat.t0, na.rm=T);

	dat.t1<-dat$phenos_table[ , m-1 ];
	par.rho <- cor(dat.t1, dat.t0 );
	cat("CSCM",par.rho,par.s2,"\n");
	return( c( par.rho, par.s2));
}

CSCM.is_valid<-function(par)
{
	return( par[1] >= -1 && par[1] <= +1);
}

CSCM.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.rho",  covar$par[1] );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par[2] );
	
	return(paste(str0, str1, str2, sep=""));
}

CSCM.get_par_num<-function(dat)
{
	return(2);
}

CSCM.guess_simu_sigma<-function( sim.mu, prob, H2 )
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

covar_CSCM<-list(
	name 	= "CSCM",
	desc 	= "parametric stationary Compound Symmetry:Correlation Metric(CSCM) model",
	get_est_param= covar.get_est_param,
	get_par_num  = CSCM.get_par_num, 
	is_valid     = CSCM.is_valid,
	get_mat      = CSCM.get_mat,
	get_init_rand= CSCM.get_init_rand,
	get_summary  = CSCM.get_summary,
	guess_simu_sigma  = CSCM.guess_simu_sigma);

ZZZ.regcovar_cscm<-function()
{
	COVAR_CSCM  <<- FM2.reg_covar( covar_CSCM );
}