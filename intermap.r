###############################################################
# 
# New Systems Mapping Application(SysMap1)#
# 
#
# Routine:
#  1) inter.param
#  2) summary.inter.par
#  3) inter.simulate
#  4) inter.data_est
#  5) summary.inter.dat
#  6) summary.FM2.par 
#  7) interdat.get_simuparam				          
#  8) interdat.get_simuparam
#  9) interdat.summary_par
# 10) interpar.get_summary
# 11) interdat.simulate
# 12) interdat.summary		      
# 13) INTER.qtlscan
# 14) inter.permutation
# 15) summary.INT.ret.perm
# 16) plot.INT.ret.perm
#
# History:
# 03/17/2010 Version 0.1
#
###############################################################

#--------------------------------------------------------------
# inter.param
#
# Create a param object for the simulation(BC,F2) 
#--------------------------------------------------------------
inter.param <- function( par_obj , cross_type )
{
	if (is.null( FM_sys ) ) 
		FM2.start();
	
	cross <- FM2.get_cross(cross_type);
	if (is.null(cross))
	{
		stop("inavlid cross type");
	}
	FM2.cross <<- cross;
	
	par <- interdat.get_simuparam( par_obj ,cross_type );
	
	class( par )<- "inter.par";
	return(par);
}

#--------------------------------------------------------------
# summary.inter.par
#
# used by summary( FM2.par.obj )
#--------------------------------------------------------------
summary.inter.par <- function( par_obj , file=NA , append = TRUE )
{
	str <- interdat.summary_par(par_obj );
	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) )
}	

#--------------------------------------------------------------
# inter.simulate
#
# Create a simulation data set by parameter object
#--------------------------------------------------------------
inter.simulate <- function(par_obj)
{
	dat <- interdat.simulate( par_obj );
	class(dat) <- "inter.dat";
	
	return(dat);	
}

#--------------------------------------------------------------
# inter.data_est
#
# Load a real data, phenotype file, genotype file, marker file 
# is necessary.
#--------------------------------------------------------------
inter.data_est <- function( dat,cross_type, prob=0.05)
{
	if (is.null( FM_sys ) )
		FM2.start();
	
	if(class(dat)!="inter.dat")
	{
		stop("inavlid data");
	}
	
	cross.obj <- FM2.get_cross(cross_type);
	if (is.null(cross.obj))
	{
		stop("inavlid cross type");
	}
	FM2.cross <<- cross.obj;	

	FM2.set_value("fitting.prob", prob);
	parin <- inter_get_init( dat );
cat("INIT", parin ,"\n");
	param <- interval.get_est_param(dat, parin );
	if (is.na(param)||is.na(param$sigma)||is.na(param$mu))
		return(NA);
	
	dat$param<-list( sigma=param$sigma, mu= param$mu, val=param$val, method=param$method );
	
	return( dat);
}

#--------------------------------------------------------------
# summary.inter.dat
#
# used by summary( FM2.dat.obj )
#--------------------------------------------------------------
summary.inter.dat <- function( dat_obj, file=NA , append = TRUE  )
{	
	str<- interdat.summary( dat_obj );
	
	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) )
}

#--------------------------------------------------------------
# interdat.get_simuparam
#
# Create a parameter object for BC, F2, RIL or NP.
#--------------------------------------------------------------
interdat.get_simuparam <- function( par_obj , cross_type )
{
	par <- par_obj;
	par$cross_type <- cross_type;
	par$par_num    <- FM2.cross$gen_num+1;
	par$crossname       <- FM2.cross$name;
	
	if (cross_type==CROSS_BC) par$QQ2   <- NULL;
	if (cross_type==CROSS_RIL)par$Qq1   <- NULL;
	
	return( par );
}

#--------------------------------------------------------------
# interdat.summary_par
# 
# Summarize the parameter object for the F2 simulation
# Used by summary( XX.F2.par object)
#--------------------------------------------------------------
interdat.summary_par<-function( par_obj )
{	
	cross  <- par_obj$cross_type;
	parname<- FM2.cross$name;
	
	if ( par_obj$crossname != parname )
	{
		sErrMsg<- paste( "Error: Not a parameter object for intermap model. par$name=", par_obj$name, sep="");
		stop( sErrMsg );
	}
	marker_s <- paste( cumsum( par_obj$simu_mrkdist ), collapse=",", sep="");
	
	strt <- sprintf("The parameter for intermap model:\n");
	stru <- sprintf("------------------------------------\n");
	str0 <- sprintf("%15s: %s\n", 	   "Date", 		Sys.time() );
	str1 <- sprintf("%15s: %-10.0f\n", "Sample size", 	par_obj$simu_N );
	str2 <- sprintf("%15s: %s\n", 	   "Marker pos.", 	marker_s );
	str3 <- sprintf("%15s: %-10.0f\n", "QTL pos.", 	  	par_obj$simu_qtlpos );
	str4 <- "";
	if (!is.null(interpar.get_summary))
		str4 <- interpar.get_summary(par_obj, "%15s: %-10.5f\n")
	strd <- sprintf("------------------------------------\n\n");
	
	str <- paste(strt,  stru, str0, str1, str2, str3,  str4, strd, sep="" );
	return ( str );
}

#--------------------------------------------------------------
# inter.load_data
#
# Load a real data, phenotype file, genotype file, marker file 
# is necessary.
#--------------------------------------------------------------
inter.load_data <- function( pheno_file, geno_file, marker_file )
{
	if (is.null( FM_sys ) )
		FM2.start();
		
	dat <- interdat.load(pheno_file, geno_file, marker_file );
	
	class(dat)<-"inter.dat";
	
	return(dat);
}

#--------------------------------------------------------------
# private: interdat.summar_dat
#
# summarize the important information for the simulate data or 
# real data.
#
# input : data object.
#--------------------------------------------------------------
interdat.summary <- function( dat_obj )
{
	if ( class(dat_obj) != "inter.dat" )
	{
		sErrMsg<- "Error: Not a dat set for InterMap model.";
		stop( sErrMsg );
	}
	strt <- sprintf("The data set for InterMap model:\n");
	stru <- sprintf("-----------------------------------\n");
	str0 <- sprintf("%15s: %s\n", 	 "Date", 	   Sys.time() );
	str1 <- sprintf("%15s: %s\n", 	 "Pheno. file",	   dat_obj$pheno_file );
	str2 <- sprintf("%15s: %-10.0f\n", "Sample size",  dat_obj$sample_N );
	str3 <- sprintf("%15s: %-10.6f\n", "val",  dat_obj$param$val );
	str4 <- sprintf("%15s: %-10.6f\n", "sigma",  dat_obj$param$sigma );
	str5 <- sprintf("%15s: %-10.6f\n", "mu",  dat_obj$param$mu );
	strd <- sprintf("------------------------------------\n\n");
	
	str <- paste(strt, stru, str0, str1, str2, str3, str4, str5, strd, sep="" );
	
	return (str);	
}

#--------------------------------------------------------------
# INTER.qtlscan
#
# Hypothesis test for data object, the test methods should be 
# (10,11,12,13,14,15)
#
# 1) scan_step, default=2, an interval distance used to scan flanking 
#    marker, default is 2cm.
# 2) peak_count, default=5, a number shows how many significant QTL will 
#    be selected.
# 3) plot_doctype, default=pdf, the figure output type for summary command.
#
#--------------------------------------------------------------	
INTER.qtlscan <- function( dat_obj, options=list() ) 
{
	if(class(dat_obj)!="inter.dat")
	{
		stop("Invalid data object.");
	}
	
	set.scan_step <- FALSE;
	if (!(is.null(options$scan_step) || is.na(options$scan_step) ))
	{	
		set.scan_step 	<- TRUE;
		old.scan_step   <- FM2.set_value("scan_step", options$scan_step);
	}
	
	set.peak_count 	<- FALSE;
	
	if (!(is.null(options$peak_count) || is.na(options$peak_count) ) )
	{	
		set.peak_count 	<- TRUE;
		old.peak_count  <- FM2.set_value("peak_count", options$peak_count);
	}
	ret<- INTER.model$qtlscan( dat_obj );

	if (set.scan_step)
			FM2.set_value("scan_step", old.scan_step);
		if (set.peak_count)
			FM2.set_value("peak_count", old.peak_count);
	
		return (ret);
}

#--------------------------------------------------------------
# INTER.permutation
#
# Permutation tests.
#--------------------------------------------------------------
inter.permutation <- function( dat_obj, options=list() )
{
	old.debug <- NULL;
	set.debug <- FALSE;
	if (!is.null(options$debug) && !is.na(options$debug) )
	{
	        old.debug <- FM2.set_value("debug", as.numeric( options$debug));
		set.debug <- TRUE;
	}
	
	old.cluster_count <- NULL;
	set.cluster_count <- FALSE;
	if (!is.null(options$cluster_count) && !is.na(options$cluster_count) )
	{
	        old.cluster_count <- FM2.set_value("cluster_count", as.numeric( options$cluster_count));
		set.cluster_count <- TRUE;
	}
	
	old.permu_loop <- NULL;
	set.permu_loop <- FALSE;
	if (!is.null(options$permu_loop) && !is.na(options$permu_loop))
	{
	        old.permu_loop <- FM2.set_value("permu_loop", as.numeric( options$permu_loop) );
		set.permu_loop <- TRUE;
	}
		
	#loops <- as.numeric( FM2.get_value("permu_loop", def=1000) );
	loops <- 10;
	if (loops<100)
	{
		warning("Permutation loop is too few(<100).\n");
	}
	
	ret<- interpermu.execute( dat_obj );
	class( ret ) <- "INT.ret.perm";
	
	if (set.debug)
			FM2.set_value("debug", old.debug);
		if (set.cluster_count)
			FM2.set_value("cluster_count", old.cluster_count);
		if (set.permu_loop)
			FM2.set_value("permu_loop", old.permu_loop);
	
	return( ret );
}

#--------------------------------------------------------------
# summary.INT.ret.perm
#
# used by summary( INT.ret.perm.obj )
#--------------------------------------------------------------
summary.INT.ret.perm <- function( perm_obj, file=NA, append=TRUE )
{
	str <- interpermu.summary(perm_obj);
	
	if (is.na(file))
		return ( cat( str ) )
	else
		return ( cat( str, file= file, append = append) )
}

#--------------------------------------------------------------
# plot.INT.ret.perm
#
# used by plot( FM.ret.perm object )
#--------------------------------------------------------------
plot.INT.ret.perm <- function( perm_obj )
{
	interpermu.plot( perm_obj )
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function: inter.one_step
#
# The result is saved into [pheno_csv.rdata]
#
# Options:
#1) cluster_count, default=1, the cluster count for parallel permutation.
#2) permu_loop, default=1000, the count of permutation loop.
#3) file.summary, default=NA, a filename where the summary information 
#   which displayed in the console will be saved.
#4) file.rdata, default=NA, a filename for RDATA format file where the 
#   data object, result object and permutation object will be saved.
#5) file.report, default=NA, a filename for PDF report file where the 
#   summaries and figures will be saved.
#6) scan_step, default=2, an interval distance used to scan flanking 
#   marker, default is 2cm.
#7) peak_count, default=5, a number shows how many significant QTL will 
#   be selected.
#8) plot_doctype, default=pdf, the figure output type for summary command.
#9) np.order, default=6, the order of Legendre polynomial for nonparametric 
#   method.
#10) sg.order, default=3, the order of segmental lines.
#11) CM.model1, vector, default=c(CURVE_LC, 3), the 1st model(curve) type 
#    and count of parameters for composite model.
#12) CM.model2, vector, default=c(CURVE_NP, 4) ), the 2nd model(curve) type 
#    and count of parameters for composite model.
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AM.one_step<-function(pheno_csv, geno_csv, marker_csv, cross_type, options=list(prob=0.05))
{
	if (is.null( FM_sys ) )
		FM2.start();

	dat_os<- inter.load_data( pheno_csv, geno_csv, marker_csv );
	
	if (is.null(dat_os))
		stop("STOP:Failed to load the phenotype.\n");
		
	est <- inter.data_est(dat_os, cross_type, options$prob);	
	if (is.na(est) || is.null(est))
		cat("Failed to estimate the data.\n")
	
	file.summary<-""
	if (!is.null(options$file.summary))
		file.summary <- options$file.summary;
	cat("\nSummary report:\n", file=file.summary);
	
	try( summary( dat_os, file=file.summary ) );
	
	ret_os<- INTER.qtlscan( dat_os,options=options );
	summary(ret_os, dat_os);
	##Save 1;
	if (is.null( options$file.rdata) )
		save(dat_os, ret_os, file=paste(pheno_csv, ".rdata", sep="") )
	else
		save(dat_os, ret_os, file=options$file.rdata );
	
	nLoop <- FM2.get_value("permu_loop", def=1000);
	if ( !is.null(options$permu_loop) && !is.na(options$permu_loop) )
		nLoop <- as.numeric( options$permu_loop );

	if (nLoop>0)
	{
		ret_perm_os<- inter.permutation( dat_os, options=options);
		summary( ret_perm_os, file=file.summary );

		p05 <- 0;
		p01 <- 0;
		idx05 <- which(ret_perm_os$pv_table[,1]==0.05);
		idx01 <- which(ret_perm_os$pv_table[,1]==0.01);
		if (length(idx05)>0 )
			p05 <- ret_perm_os$pv_table[idx05[1],2];
		if (length(idx01)>0 )
			p01 <- ret_perm_os$pv_table[idx01[1],2];

		if (p05>0 && p01>0)
			ret_os <- INTER.model$set_lr2_cutoff( ret_os, p05, p01);

		try( summary(ret_os, dat_os, file=file.summary) );

		##Save 2;
		if (is.null( options$file.rdata) )
			save(dat_os, ret_os, ret_perm_os, file=paste(pheno_csv, ".rdata", sep="") )
		else
			save(dat_os, ret_os, ret_perm_os, file=options$file.rdata );
	}
	else
		try( summary(ret_os, dat_os, file=file.summary) );

	return(ret_os);
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AM.simu_test
#
# Abstract: any model, any cross by simulation data
# 
# The result is saved into [LC_simu_test_XX_XX.rdata]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AM.simu_test <- function( cross_type, par_obj=inter_par )
{
	if (is.null( FM_sys ) )
		FM2.start();
	
	par_as <- inter.param( par_obj, cross_type );
	summary(par_as);
	
	dat_as<- inter.simulate( par_as);
	
	if (is.null(dat_as))
		stop("STOP:Failed to load the phenotype.\n");
	
	
	est <- inter.data_est(dat_as, cross_type);	
	if (is.null(est))
		stop("STOP:Failed to estimate the data.\n");
	
	dat_as$param <- est$param;
	summary(dat_as);
	
	ret_as<- INTER.qtlscan( dat_as);
	summary(ret_as, dat_as);
	##Step 1
	rdata.name <- paste( "simu_test_", FM2.cross$name,".rdata", sep="");
	save(dat_as, ret_as, file = rdata.name );
	
	##Step 2(Ppermutation) 
	FM2.set_value("debug", TRUE)
	ret_perm_as<- inter.permutation( dat_as, options=list(permu_loop=100) );
	summary( ret_perm_as );
	
	##Step 3(Save All Data)
	save(dat_as, ret_as, ret_perm_as, file = rdata.name );
	cat("\nNotice: The results(dat_as, ret_as and ret_perm_as) are stored into ", rdata.name, "\n");
	       
	return( list(dat_as, ret_as, ret_perm_as) );
		
}