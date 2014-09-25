#-----------------------------------------------------------------
# CURVE_KORF:  Korf's Curve
#-----------------------------------------------------------------

KORF.get_mu <- function(par, times, options=list())
{
       
	return ( par[[1]]*exp(-par[[2]]*times^(-par[[3]])) ) ;
	
}


KORF.get_init_rand<-function(dat, times)
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
	u11 <- mean(c(u1,u2));
	
	par.a <- um * um / um.1;
	par.b <- -log(u11/par.a);
	par.c <- try( (log(m))/log((log(um/par.a))/(-par.b)) );
	if (is.na(par.c)) par.c <- 0.6;
	if(class(par.c)=="try-error")
		par.c <- 1
	return(c(par.a, par.b, par.c));
}


KORF.get_summary <-function( curve, format)
{
	str0 <- sprintf(format, "-Curve.Lval",  curve$val );
	str1 <- sprintf(format, "-Curve.a",  curve$par[1] );
	str2 <- sprintf(format, "-Curve.b",  curve$par[2] );
	str3 <- sprintf(format, "-Curve.r",  curve$par[3] );
	
	return(paste(str0, str1, str2, str3, sep=""));
}

class_KORF<-list(
	type 	 	= 8,
	name 	 	= "KORF",
	desc 	 	= "Korf's Curve",
	par_num    	= 3, 
	par_name 	= c( "a", "b", "r"),
	get_mu      	= KORF.get_mu,
	get_summary    	= KORF.get_summary,
	get_init_rand   = KORF.get_init_rand);


par_KORF<-list(
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
		simu_a   = 21.9824,
		simu_b   = 9.7768,
		simu_r   = 0.4699),

	# 	Qq1         
	Qq1 = list(
		simu_a   = 19.9810,
		simu_b   = 8.776,
		simu_r   = 0.47),
	
	#	qq0
	qq0 = list(
		simu_a   = 15.95,
		simu_b   = 7.58,
		simu_r   = 0.48)
);

ZZZ.regcovar_korf<-function()
{
	CURVE_KORF <<- FM2.reg_curve( class_KORF );
}