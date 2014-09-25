#-----------------------------------------------------------------
# CURVE_GOMP:  Gompertz Curve
#-----------------------------------------------------------------

GOMP.get_mu <- function(par, times, options=list())
{
       
	return ( par[[1]]*(par[[2]]/par[[1]])^exp(-par[[3]]*times) );
	
}

GOMP.get_init_rand<-function(dat, times)
{
	m <- dim( dat$phenos_table)[2];
	um <- mean(dat$phenos_table[,m], na.rm=T);
	um.1 <- mean(dat$phenos_table[,m-1], na.rm=T);
	u2 <- mean(dat$phenos_table[,2], na.rm=T);
	u1 <- mean(dat$phenos_table[,1], na.rm=T);
	if (u2>u1) 
		u0 <- u1*3/4 
	else
		u0 <- u1*4/3;
	
	par.m <- u0;
	par.k <- um * um / um.1;
	par.r <- try( -(log((log(um/par.k))/log(par.m/par.m)))/par.k );
	if (is.na(par.r)) par.r <- 0.6;
	if(class(par.r)=="try-error")
		par.r <- 1
	return(c(par.k, par.m, par.r));
}



GOMP.get_summary <-function( curve, format)
{
	str0 <- sprintf(format, "-Curve.Lval",  curve$val );
	str1 <- sprintf(format, "-Curve.k",  curve$par[1] );
	str2 <- sprintf(format, "-Curve.m",  curve$par[2] );
	str3 <- sprintf(format, "-Curve.r",  curve$par[3] );
	
	return(paste(str0, str1, str2, str3, sep=""));
}

class_GOMP<-list(
	type 	 	= 9,
	name 	 	= "GOMP",
	desc 	 	= "Gompertz Curve",
	par_num    	= 3, 
	par_name 	= c( "k", "m", "r"),
	get_mu      	= GOMP.get_mu,
	get_summary    	= GOMP.get_summary,
	get_init_rand   = GOMP.get_init_rand);
	


par_GOMP<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1,2,4,5,7,9,10),
	time.std = c(1,2,4,5,7,9,10),
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#qtl position
	simu_qtlpos = 50,			

	simu_cross = list(),
	
	simu_covar = list(),
	
	#AR1H model
	covar.AR1H = list(
		#rho
		rho   = 0.5,
		s2    = c(1,1.2,1.3,1.4,1.5,1.5,1.6) ),
		
	#AR1 model
	covar.AR1 = list(
		#rho
		rho   = 0.5,
		s2    = 1.22 ),
		
	#SI  model
	covar.SI = list(
			
	        s2    = 1.56 ),
	        
	#CSCM  model
	covar.CSCM = list(
		#rho
		rho   = 0.759,
		s2    = 1.5 ),	
			
	#CSH   model
	covar.CSH = list(
		#rho
		rho   = 0.7543,
		s2    = c(1,1.2,1.3,1.4,1.5,1.5,1.6) ),
		
	#DG   model
	covar.DG = list(
			
		s2    = c(1,1.2,1.3,1.4,1.5,1.5,1.6) ),
			
	#ARMA11
	covar.ARMA11 = list(
	        #rho
		rho   = 0.5,
		s2    = 1.22,
		phi   = 0.52  ),
			
			
	#FA1
	covar.FA1 = list(
	         d       = 1.7543,
		 lambda    = c(1,1.2,1.3,1.4,1.5,1.5,1.8) ),
		
	#HF
	covar.HF = list(
				
	  	lambda    = 0.01,
		s2    = c(1,1.02,1.03,1.04,1.05,1.05,1.08) ),
		
	#CS
	covar.CS = list(
			       
		sigma1   = 1.1,
		sigma      =1.5        ),
		
	#SAD1
	covar.SAD1 = list(
				       
		 phi   = 1.1,
		 s2       =1.5        ),
	
	
	
	# 	QQ2         
	QQ2 = list(
		simu_k   = 50,
		simu_m   = 5,
		simu_r   = 0.5),

	# 	Qq1         
	Qq1 = list(
		simu_k   = 45,
		simu_m   = 4,
		simu_r   = 0.3),
	
	#	qq0
	qq0 = list(
		simu_k   = 30,
		simu_m   = 3,
		simu_r   = 0.2)
);

ZZZ.regcovar_gomp<-function()
{
	CURVE_GOMP <<- FM2.reg_curve( class_GOMP );
}

