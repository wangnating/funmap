#-----------------------------------------------------------------
# CURVE_LINE:  Linear model
#-----------------------------------------------------------------

LINE.get_mu <- function(par, times, options=list())
{
       
	return ( par[[1]]+times*par[[2]] );
	
}

LINE.get_init_rand<-function(dat, times)
{
	m <- dim( dat$phenos_table)[2];
	um <- mean(dat$phenos_table[,m], na.rm=T);
	u2 <- mean(dat$phenos_table[,2], na.rm=T);
	u1 <- mean(dat$phenos_table[,1], na.rm=T);
	if (u2>u1) 
		u0 <- u1*3/4 
	else
		u0 <- u1*4/3;
	
	par.m <- u0;
	par.r <- (um-par.m)/m;
	
	return(c(par.m, par.r));
}


LINE.get_summary <-function( curve, format)
{
	str0 <- sprintf(format, "-Curve.Lval",  curve$val );
	str1 <- sprintf(format, "-Curve.m",  curve$par[1] );
	str2 <- sprintf(format, "-Curve.r",  curve$par[2] );
	
	return(paste(str0, str1, str2, sep=""));
}

class_LINE<-list(
	type 	 	= 4,
	name 	 	= "LINE",
	desc 	 	= "Linear model",
	par_num    	= 2, 
	par_name 	= c( "m", "r"),
	get_mu      	= LINE.get_mu,
	get_summary    	= LINE.get_summary,
	get_init_rand   = LINE.get_init_rand);
	

par_LINE<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1,2,4,5,7,9,10),
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
		simu_m   = 1,
		simu_r   = 10),

	# 	Qq1         
	Qq1 = list(
		simu_m   = 1,
		simu_r   = 5),
	
	#	qq0
	qq0 = list(
		simu_m   = 1,
		simu_r   = 1)
);

ZZZ.regcovar_line<-function()
{
	CURVE_LINE <<- FM2.reg_curve( class_LINE );
}