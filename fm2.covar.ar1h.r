#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  AR(1): Heterogenous
#  
#  Variance : AR1H
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

AR1H.get_mat <- function(par0, times, options=list())
{      
	
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	t_len <- length( times );	
	AR1H <- array(0, dim=c(t_len,t_len));
	phi<--log(par[1]);
	for (k0 in 1:t_len)
		for (k1 in 1:t_len)
			AR1H[k0,k1] <- par[k0+1]*par[k1+1] * exp(-phi*abs( times[k0] - times[k1] ));
       
	return(AR1H);
}

AR1H.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	
	dat.t0<-c();
	par.s2<-c();
	for(i in 1:m){
	dat.t0<-dat$phenos_table[ , i ] ;
	par.s2 [i]<- sd(dat.t0, na.rm=T);}

	dat.t1<-dat$phenos_table[ , m-1 ];
	par.rho <- cor(dat.t1, dat.t0,use="complete.obs");
	
	return( c( par.rho, par.s2));                     
}

AR1H.is_valid<-function(par)
{
	return( par[1] >= -1 && par[1] <= +1);
}

AR1H.get_par_num<-function(dat)
{
        m<-length(dat$sample_times);
	return(m+1);
}

AR1H.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.rho",  covar$par[1] );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par[-1] );
	
	return(paste(str0, str1, str2, sep=""));
}

AR1H.guess_simu_sigma<-function( sim.mu, prob, H2 )
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

covar_AR1H<-list(
	name 	= "AR1H",
	desc 	= "parametric stationary autoregressive(AR) heterogenous model",
	get_est_param= covar.get_est_param,
	get_par_num  = AR1H.get_par_num, 
	is_valid     = AR1H.is_valid,
	get_mat      = AR1H.get_mat,
	get_init_rand= AR1H.get_init_rand,
	get_summary  = AR1H.get_summary,
	guess_simu_sigma  = AR1H.guess_simu_sigma);

ZZZ.regcovar_ar1h<-function()
{
	COVAR_AR1H  <<- FM2.reg_covar( covar_AR1H );
}

