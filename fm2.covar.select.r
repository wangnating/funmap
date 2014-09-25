curve.select <- function(dat)
{
  len.curve <- length(FM2.curves);
  Aic <- c();
  str2 <- "";
 
  	for(i in 1:len.curve)
  	{
  	        curve <- FM2.curves[[i]];
  	        FM2.curve<<-curve;
  	        cat("NAME",FM2.curve$name);
  	        q <- curve$par_num;
  	        
  	        curve.obj <- curve$get_est_param(dat,curve$f_init_rand);
  	       
  	        if(is.na(curve.obj)||is.na(curve.obj$val))
  	        return(NAN);
  	        
  	        lr <- curve.obj$val;
  	        aic <- 2*lr+2*q;
  
                str1 <- sprintf("%10s: %10s: %-8.3f %10s: %-2.0f %5s: %-8.3f\n", 
  		                 curve$name, "-LOG(L)", lr, "Para Num", q, "AIC", aic);
  		str2 <- paste(str2, str1, sep=""); 
  		
  		Aic <- c(Aic,aic);
  	}
  	
  	   curve <- FM2.curves[[which.min(Aic)]];
  	  
  	   FM2.curve <<- curve;
	   Curve.obj <- FM2.curve$get_est_param(dat,curve$f_init_rand);
	   
	   if(is.na(Curve.obj)||is.na(Curve.obj$par))
	  		return(NAN);
	   par <-Curve.obj$par;
	   
	   str3 <- sprintf("%10s\n","Min(AIC):");
	   str4 <- sprintf("%10s: %10s: %-10.3f %10s ", curve$name ,"AIC" ,min(Aic), "Paras");      
	
	   str <- "";
	   str <- paste( str,"----------------------------------------------------------------------\n",sep="");
	   str <- paste( str,"Curve Selection: \n",sep="");
	   str <- paste( str, str2,sep="");
	   
	   str <- paste( str,str3,str4,sep="");
	   cat(str,capture.output(show(par)),"\n");
	   strr <- "----------------------------------------------------------------------\n";
	   cat(strr);
	   
   return(curve$reg_no);	 
  	
}

covar.select <- function(dat,curve)
{
 
  len.cov <- length(FM2.covars);
  Aic <- c();
  str2 <- "";
  
  	for(i in 1:len.cov)
  	{
  	       
  		cov <- FM2.covars[[i]];
  		FM2.covar<<-cov;
  		cat("NAME",FM2.covar$name);
  		q <- cov$get_par_num(dat);
  		
  		cov.obj <- cov$get_est_param(dat, curve);
  		if(is.na(cov.obj)||is.na(cov.obj$val))
  		return(NAN);
  		
  		lr <- cov.obj$val;
  		aic <- 2*lr+2*q;
  		
  		str1 <- sprintf("%10s: %10s: %-8.3f %10s: %-2.0f %5s: %-8.3f\n", 
  		                 cov$name, "-LOG(L)", lr, "Para Num", q, "AIC", aic);
  		str2 <- paste(str2, str1, sep="");
  		
  		Aic <- c(Aic,aic);
  	 }
  	
   
   cov <- FM2.covars[[which.min(Aic)]];
   Cov.obj <- cov$get_est_param(dat, curve);
   
   if(is.na(Cov.obj)||is.na(Cov.obj$par))
  		return(NAN);
   par <- Cov.obj$par;
   
   str3 <- sprintf("%10s\n","Min(AIC):");
   str4 <- sprintf("%10s: %10s: %-10.3f %10s ", cov$name ,"AIC" ,min(Aic), "Paras");      

   str <- "";
   str <- paste( str,"----------------------------------------------------------------------\n",sep="");
   str <- paste( str,"Covariance Selection: \n",sep="");
   str <- paste( str, str2,sep="");
   
   str <- paste( str,str3,str4,sep="");
   cat(str,capture.output(show(par)),"\n");
   strr <- "----------------------------------------------------------------------\n";
   cat(strr);
   return(cov);	 
   
}




#FM2.auto_select_test<-function( par_obj, curve_type, cross_type, covar_type= COVAR_AR1 )