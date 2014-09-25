#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Factor Analytic: First-Order
#  
#  Variance : FA1
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

FA1.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	t_len <- length( times );
	FA.1 <- array(0, dim=c(t_len,t_len));
	for (k0 in 1:t_len)
		for (k1 in 1:t_len)
		{
		        if(k0==k1)
		         FA.1[k0,k1]<-par[k0+1]*par[k1+1]+par[1]
		        else
		         FA.1[k0,k1]<-par[k0+1]*par[k1+1];
		  }        
		
	return(FA.1);
}

FA1.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	
	dat.t0<-c();
	par.s2<-c();
	
	for(i in 1:m){
	             dat.t0<-dat$phenos_table[ , i ] ;
	             par.s2 [i]<- sd(dat.t0, na.rm=T);}
	                         
	dat.t1<-dat$phenos_table[ , 1 ];
	dat.t2<-dat$phenos_table[ , 2 ];
	par.rho <- cor(dat.t1, dat.t2,use="complete.obs" );
	
	par.d<-(par.s2[1]+par.s2[2]-sqrt((par.s2[1]-par.s2[2])^2+4*par.rho^2))/2;
	
	par.lambda<-c();
	for(i in 1:m)
	{
	             par.lambda[i]<-sqrt(par.s2[i]-par.d);
	}
	
	return( c( par.d,par.lambda));
}

FA1.is_valid<-function(par)
{
	return( TRUE);
}

FA1.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.d",  covar$par[1] );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par[-1] );
	
	return(paste(str0, str1, str2, sep=""));
}

FA1.get_par_num<-function(dat)
{      
        m <- dim(dat$phenos_table)[2];
	return(m+1);
}

FA1.guess_simu_sigma<-function( sim.mu, prob, H2 )
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

covar_FA1<-list(
	name 	= "FA1",
	desc 	= "parametric stationary Factor Analytic: First-Order(FA1) model",
	get_est_param= covar.get_est_param,
	get_par_num  = FA1.get_par_num, 
	is_valid     = FA1.is_valid,
	get_mat      = FA1.get_mat,
	get_init_rand= FA1.get_init_rand,
	get_summary  = FA1.get_summary,
	guess_simu_sigma  = FA1.guess_simu_sigma);

ZZZ.regcovar_fa1<-function()
{
	COVAR_FA1  <<- FM2.reg_covar( covar_FA1 );
}