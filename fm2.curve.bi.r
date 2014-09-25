#-----------------------------------------------------------------
#
#
# CURVE_BI:  Bi-exponential Curve
#
#
#-----------------------------------------------------------------

BI.get_mu <- function(par, times, options=list())
{
	return(  par[[1]]*exp(-par[[2]]*times) + par[[3]]*exp(-par[[4]]*times) );
}

class_BI<-list(
	type 	 	= 1,
	name 	 	= "BI",
	desc 	 	= "Bi-exponential Curve",
	par_num    	= 4, 
	par_name 	= c( "p1", "d1", "p2", "d2"),
	get_mu      	= BI.get_mu );

BI.get_summary <-function( curve, format)
{
	str0 <- sprintf(format, "-Curve.Lval",  curve$val );
	str1 <- sprintf(format, "-Curve.p1",  curve$par[1] );
	str2 <- sprintf(format, "-Curve.d1",  curve$par[2] );
	str3 <- sprintf(format, "-Curve.p2",  curve$par[3] );
	str4 <- sprintf(format, "-Curve.d2",  curve$par[4] );
	
	return(paste(str0, str1, str2, str3, str4, sep=""));
}

class_BI<-list(
	type 	 	= 2,
	name 	 	= "BI",
	desc 	 	= "Bi-exponential Curve",
	par_num    	= 4, 
	par_name 	= c( "p1", "d1", "p2", "d2"),
	get_mu      	= BI.get_mu
	get_summary    	= BI.get_summary
	#get_init_rand   = BI.get_init_rand);

par_BI<-list(
	#sample size
	simu_N      = 100,
	#trait measured in time 1 to 7
	simu_times  = c(1:7),
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#qtl position
	simu_qtlpos = 50,			

	simu_cross = list(),
	
	simu_covar = list(),
	
	#AR1H model
	covar.AR1H = list(
			rho   = 0.5,
			s2    = c(1,1.2,1.3,1.4,1.5,1.5,1.6) ),
	
	#AR1 model	
	covar.AR1 = list(
			rho   = 0.5,
			s2    = 1.22 ),
	
	#SI  model
	covar.SI = list(	
		        s2    = 1.56 ),
		        
	#CSCM  model
	covar.CSCM = list(
			rho   = 0.759,
			s2    = 1.5 ),	
				
	#CSH   model
	covar.CSH = list(
			rho   = 0.7543,
			s2    = c(1,1.2,1.3,1.4,1.5,1.5,1.6) ),
			
	#DG   model
	covar.DG = list(	
			s2    = c(1,1.2,1.3,1.4,1.5,1.5,1.6) ),
				
	#ARMA11 model
	covar.ARMA11 = list(
			rho   = 0.5,
			s2    = 1.22,
			phi   = 0.52  ),
					
	#FA1  model
	covar.FA1 = list(
		        d       = 1.7543,
			lambda    = c(1,1.2,1.3,1.4,1.5,1.5,1.8) ),
			
	#HF  model
	covar.HF = list(		
		  	lambda    = 0.01,
			s2    = c(1,1.02,1.03,1.04,1.05,1.05,1.08) ),
			
	#CS  model
	covar.CS = list(	       
			sigma1   = 1.1,
			sigma      =1.5        ),
		
	#SAD1  model
	covar.SAD1 = list(		       
			phi   = 1.1,
			s2       =1.5        ),
	
	# 	QQ2         
	QQ2 = list(
		simu_p1  = 19.9824,
		simu_d1  = 0.4699,
		simu_p2  = 8.7768,
		simu_d2  = 1.4699),

	# 	Qq1         
	Qq1 = list(
		simu_p1  = 17.9824,
		simu_d1  = 0.0699,
		simu_p2  = 9.7768,
		simu_d2  = 1.0699),
	
	#	qq0
	qq0 = list(
		simu_p1  = 15.9507,
		simu_d1  = 0.1836,
		simu_p2  = 10.5737,
		simu_d2  = 1.8836)
);

ZZZ.regcovar_bi<-function()
{
	CURVE_BI <<- FM2.reg_curve( class_BI );
}