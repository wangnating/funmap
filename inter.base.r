#################################################################
#
# New Systems Mapping Application(SysMap1)#
# Common library
#
#    1) interval.mlefunc
#    2) inter_get_init
#    3) interval.get_est_param
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################


#--------------------------------------------------------------
#  interval.mlefunc
#
#  MLE function, a function to be maximized 
#  
#--------------------------------------------------------------
interval.mlefunc <- function( par, y)
{
	hh <- try( fy0 <- dnorm( y, par[2], par[1], log=T), T );
	if (class(hh)=="try-error") return(NaN);
	A <- -sum( fy0 );
	
	return (A);

}

#--------------------------------------------------------------
#  inter_get_init
#
#  Used by getting an initial value 
#  
#--------------------------------------------------------------
inter_get_init <- function( dat )
{
	y <- as.matrix(dat$phenos_table);
	sd <- sd(y, na.rm=T);
	mu <- mean(y, na.rm=T);
	parin <- c( sd, mu);
	
	return(parin);
}

#--------------------------------------------------------------
#  interval.get_est_param
#
#  Used by estimating parameters 
#  
#--------------------------------------------------------------
interval.get_est_param<-function(dat, parin )
{
	y <- unlist( dat$phenos_table );
	
	bSuccess <- TRUE;
	h0 <- FM2.get_value( "likelihood_null");
	if ( is.na( h0) )
	{
		par <- parin;
		h0 <- interqtl.optim( par, interval.mlefunc, y = y );  
		if ( any(is.na(h0)) ) 
			bSuccess <- FALSE
		else
			FM2.set_value( "likelihood_null", h0 );
	}
	else
		if (h0$value==0 ) bSuccess <- FALSE;	
cat("\nINTERVAL Parameters:", "val(L)=", h0$value, "Par", h0$par, "\n");
		
	return(list(val=h0$value, sigma=par[1], mu=par[2], method="mlefunc"));
}