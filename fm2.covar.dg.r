#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Diagonal
#  
#  Variance : DG
#  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

source("fm2.covar.base.r");

DG.get_mat <- function(par0, times, options=list())
{
	par<-par0;
	if (class(par0)=="list")
		par <- unlist(par0);
		
	sigma<-par*par;
	DG<-diag(sigma);
	
	return(DG);
}

DG.get_init_rand<-function(dat)
{
	m <- dim(dat$phenos_table)[2];
	par.s2<-c();
	dat.t0<-c();
	for(i in 1:m)
	{
	          dat.t0<-dat$phenos_table[ , i ] ;
	          par.s2[i] <- sd(dat.t0, na.rm=T);
	 }
        cat("DG",par.s2,"\n");
	return( par.s2);
}

DG.is_valid<-function(par)
{
	return( par[1] >= -1 && par[1] <= +1);
}

DG.get_summary <-function( covar, format)
{
	str0 <- "";
	if (!is.null(covar$val))
		str0 <- sprintf(format, "-Covar.Lval",  covar$val );
	str1 <- sprintf(format, "-Covar.rho",  0 );
	str2 <- sprintf(format, "-Covar.sigma",  covar$par);
	
	return(paste(str0, str1, str2, sep=""));
}

DG.get_par_num<-function(dat)
{         
        m <- dim(dat$phenos_table)[2];
	return(m);
}

DG.guess_simu_sigma<-function( sim.mu, prob, H2 )
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

covar_DG<-list(
	name 	= "DG",
	desc 	= "parametric stationary Diagonal(DG) model",
	get_est_param= covar.get_est_param,
	get_par_num  = DG.get_par_num, 
	is_valid     = DG.is_valid,
	get_mat      = DG.get_mat,
	get_init_rand= DG.get_init_rand,
	get_summary  = DG.get_summary,
	guess_simu_sigma  = DG.guess_simu_sigma);

ZZZ.regcovar_dg<-function()
{
	COVAR_DG  <<- FM2.reg_covar( covar_DG );
}