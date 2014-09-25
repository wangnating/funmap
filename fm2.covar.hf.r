#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Huynh-Feldt
#  
#  Variance : HF
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

HF.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);

	t_len <- length( times );	
	HF <- array(0, dim=c(t_len,t_len));
	for (k0 in 1:t_len)
		for (k1 in 1:t_len)
		{       if(k0==k1)
		        HF[k0,k1]<-par[k0+1]^2
		        else
		        HF[k0,k1]<-(par[k0+1]^2+par[k1+1]^2)/2-par[1];
		 }
		        
	return(HF);	
}

HF.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	
	dat.t0<-c();
	par.s2<-c();
	for(i in 1:m){
	dat.t0<-dat$phenos_table[ , i ] ;
	par.s2[i]<- sd(dat.t0, na.rm=T);}
 
        par.rho<-c();
        lambda<-c();
        lambda[1]<-0;
        dat.t1<-dat$phenos_table[ , 1 ];
        for(i in 2:m){
	dat.t2<-dat$phenos_table[ , i ];
	par.rho[i] <- cor(dat.t1, dat.t2,use="complete.obs");
	lambda[i]  <- (par.s2[1]^2+par.s2[i]^2)/2-par.rho[i];  }
	
	par.lambda0 <- sum(lambda)/(m-1);
	cat("HF",par.lambda0,par.s2,"\n");
	return( c( par.lambda0, par.s2));                     
}

HF.is_valid<-function(par)
{
	return( TRUE);
}

HF.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
		
	str1 <- sprintf(format, "-Covar.lambda",  covar$par[1] );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par[-1] );
	return(paste(str0, str1, str2, sep=""));
}

HF.get_par_num<-function(dat)
{
        m<-length(dat$time.std);
	return(m+1);
}

HF.guess_simu_sigma<-function( sim.mu, prob, H2 )
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

covar_HF<-list(
	name 	= "HF",
	desc 	= "parametric stationary Huynh-Feldt(HF) model",
	get_est_param= covar.get_est_param,
	get_par_num  = HF.get_par_num, 
	is_valid     = HF.is_valid,
	get_mat      = HF.get_mat,
	get_init_rand= HF.get_init_rand,
	get_summary  = HF.get_summary,
	guess_simu_sigma  = HF.guess_simu_sigma);

ZZZ.regcovar_hf<-function()
{
	COVAR_HF  <<- FM2.reg_covar( covar_HF );
}