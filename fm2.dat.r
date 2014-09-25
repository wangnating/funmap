#################################################################
#
# New Systems Mapping Application(SysMap1)#
# Common library
#
#    1) dat.get_simuparam()
#    2) dat.summary_par()
#    3) dat.simulate()
#    4) dat.load()
#    5) dat.summary()
#    6) dat.report()
#    7) dat.plot()
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################

#--------------------------------------------------------------
# dat.get_simuparam
#
# Create a parameter object for BC, F2, RIL or NP.
#--------------------------------------------------------------
dat.get_simuparam<-function( par_obj, cross_type, covar_type )
{
	par <- par_obj;
	par$cross_type <- cross_type;
	par$covar_type <- covar_type;
	par$curve_type <- FM2.curve$curve_type;
	par$curve_name <- FM2.curve$name;
	par$par_num    <- FM2.curve$par_num;
	par$trait_num  <- FM2.curve$trait_num;
	#par$time.std <- FM2.curve$time.std;
	par$name       <- paste(par$curve_name, FM2.covar$name, FM2.cross$name, "par", sep=".");
		
	if (cross_type==CROSS_BC) par$QQ2   <- NULL;
	if (cross_type==CROSS_RIL)par$Qq1   <- NULL;
	#if (covar_type==COVAR_MAR1)  par$simu_covar <- par_obj$covar_ar1;
	
	covar.exp <- paste("par_obj$covar", FM2.covar$name, sep=".");
	par$simu_covar <- eval( parse(text=covar.exp) )

	return( par );
}
	
#--------------------------------------------------------------
# dat.summary_par
# 
# Summarize the parameter object for the F2 simulation
# Used by summary( XX.F2.par object)
#--------------------------------------------------------------
dat.summary_par<-function( par_obj )
{
	covar  <- par_obj$covar_type;
	cross  <- par_obj$cross_type;
	curve  <- par_obj$curve_type;
	
	parname<- paste( par_obj$curve_name, FM2.covar$name, FM2.cross$name, "par", sep=".");

	if ( par_obj$name != parname )
	{
		sErrMsg<- paste( "Error: Not a parameter object for Funmap model. par$name=", par_obj$name, sep="");
		stop( sErrMsg );
	}
	marker_s <- paste( cumsum( par_obj$simu_mrkdist ), collapse=",", sep="");

	strt <- sprintf("The parameter for FunMap model:\n");
	stru <- sprintf("------------------------------------\n");
	str0 <- sprintf("%15s: %s\n", 	   "Date", 		Sys.time() );
	str1 <- sprintf("%15s: %-10.0f\n", "Sample size", 	par_obj$simu_N );
	str2 <- sprintf("%15s: %-10.0f\n", "Sample times",	length(par_obj$simu_times) );
	str3 <- sprintf("%15s: %s\n", 	   "Marker pos.", 	marker_s );
	str4 <- sprintf("%15s: %-10.0f\n", "QTL pos.", 	  	par_obj$simu_qtlpos );
	str5 <- sprintf("%15s: %s\n", 	   "Curve",		par_obj$curve_name );
	str6 <- "";
	
	if (!is.null(FM2.curve$get_summary ))
	{
		str6.1 <- sprintf("%15s: %s\n", 	   "QQ", 		"" );
		curve.par <- list(par=unlist(par_obj$QQ2))	
		str6.2 <- FM2.curve$get_summary( curve.par, "%15s: %-10.5f\n" );
		str6.3 <- sprintf("%15s: %s\n", 	   "Qq", 		"" );
		curve.par <- list(par=unlist(par_obj$Qq1))	
		str6.4 <- FM2.curve$get_summary( curve.par, "%15s: %-10.5f\n" );
		str6.5 <- sprintf("%15s: %s\n", 	   "qq", 		"" );
		curve.par <- list(par=unlist(par_obj$qq0))	
		str6.6 <- FM2.curve$get_summary( curve.par, "%15s: %-10.5f\n" );
		str6 <- paste( str6.1, str6.2, str6.3, str6.4, str6.5, str6.6, sep="");
	}
	
	str7 <- sprintf("%15s: %s\n", 	   "Covariance",	FM2.covar$name );
	
	covar.par <- list(par = unlist(par_obj$simu_covar) );
	str8 <- FM2.covar$get_summary(covar.par, "%15s: %-10.5f\n");

	strd <- sprintf("------------------------------------\n\n");
	str <- paste(strt,  stru, str0, str1, str2, str3, 
				 str4,  str5, str6, str7, str8, strd, sep="" );
	
	return ( str );
}

#--------------------------------------------------------------
# dat.simulate
#
# Simulate a data set for backcros, including phenotype, genotype, 
# marker table.
#
# Input  
#       parameter object for Backcross
# output:
#       data object
#--------------------------------------------------------------
dat.simulate<-function( par_obj )
{
	dat<-list(
		curve_name   = par_obj$curve_name,
		curve_type   = par_obj$curve_type,
		cross_type   = par_obj$cross_type,
		covar_type   = par_obj$covar_type,
		sample_N     = par_obj$simu_N,
		#trait_num    = 1,
		sample_times = par_obj$simu_times,
		time.std = par_obj$simu_times,
		name         = paste( par_obj$curve_name, "dat", sep="."),
		pheno_file   = paste( "simu.pheno",  par_obj$curve_name, FM2.covar$name, FM2.cross$name, sep="."),
		geno_file    = paste( "simu.geno",   par_obj$curve_name, FM2.covar$name, FM2.cross$name, sep="."),
		marker_file  = paste( "simu.marker", par_obj$curve_name, FM2.covar$name, FM2.cross$name, sep="."), 
		phenos_table = NULL,
		genos_table  = NULL,
		marker_table = NULL,
		marker_obj   = NULL);
        
	dat$genos_table <- FM2.cross$get_simu_marker( par_obj$simu_N, par_obj$simu_mrkdist, par_obj$simu_cross );

	mk_name  <- c();
	mk_dist  <- c()
	mk_index <- c();
	mk_group <- c();
	mrkplace = cumsum(par_obj$simu_mrkdist);
	
	for (i in 1:length(par_obj$simu_mrkdist) )
	{
		mk_name  <- c( mk_name, paste('marker',i) );
		mk_dist  <- c( mk_dist, mrkplace[i])
		mk_index <- c( mk_index, 1);
		mk_group <- c( mk_group, 'G1' );
	}
	
 	dat$marker_table <- data.frame(Marker=mk_name, Dist=mk_dist,grp_idx=mk_index, Group=mk_group);
 	dat$marker_obj   <- fin.get_markerobj( dat$marker_table );
 	
	#generate traits
	sim.mu  <-  fin.get_traits_mu( par_obj );
	sim.covar<- FM2.covar$get_mat( par_obj$simu_covar, par_obj$simu_times );
	idx     <- max(which( mrkplace < par_obj$simu_qtlpos ));
	qtlmrk1 <- idx[1];
	qtlmrk2 <- idx[1]+1;

	dat$phenos_table<- array(0, dim=c(par_obj$simu_N, length(par_obj$simu_times)  ) )
	gen.qtl <- FM2.cross$get_simu_qtl( par_obj$simu_N, dat$genos_table[,qtlmrk1], dat$genos_table[,qtlmrk2],
									  par_obj$simu_qtlpos, mrkplace[qtlmrk1], mrkplace[qtlmrk2], par_obj$simu_cross );
         
	print(sim.mu);
	
	cat("\nH2=0.1\n");
	#if (FM2.cross$gen_num>2)
	#	show(FM2.covar$guess_simu_sigma(sim.mu,  0.7, 0.1));

	cat("\nH2=0.4\n");
	#if (FM2.cross$gen_num>2)
	#	show(FM2.covar$guess_simu_sigma(sim.mu,  0.7, 0.4));

	for (i in 1:par_obj$simu_N)
	{
		 y <- rmultnorm(1, sim.mu[ gen.qtl[i],], sim.covar );
		 dat$phenos_table[i,] <- y;
	}
	colnames( dat$phenos_table ) <- sort(rep(dat$sample_times )) ;

	cat("Data simulation is done!\n");
	return(dat);
}

#--------------------------------------------------------------
# public: dat.load
#
# load a real data set for backcross and F2 experiment.
# 
# input:
#  file : pheno_file, geno_file, marker_file
#  cross: cross type, BC or F2
#--------------------------------------------------------------
dat.load<-function( pheno_file, geno_file, marker_file, head=TRUE )
{
	dat<-list(
		name         = paste(FM2.curve$name, ".dat",sep=""),
		sample_N     = 0,
		sample_times = NULL,
		pheno_file   = pheno_file,
		geno_file    = geno_file,
		marker_file  = marker_file,
		phenos_table = NULL,
		genos_table  = NULL,
		marker_table = NULL,
		marker_obj   = NULL,
		phenos_info  = NULL);
	
	tb0 <- read.csv( pheno_file, sep=",", header=head);
	tb1 <- read.csv( geno_file, sep=",", header=head);

   	#if the ids in two files are not consistent, it would be bad data.
   	#if ( any(tb1[,1] != tb0[,1]) )
   	{
   		warnings("Warnings: the ids in phenotype file and genotype fils are not consistent!");
   		tb1.ord <- c();
   		tb0.ord <- c();
   		tb1.miss <- c();
   		tb1.multi <- c();
   		for( i in 1:length(tb0[,1]) )
   		{
   			tb1.mat <- which( tb0[i,1]==tb1[,1] )
   			if (length(tb1.mat)<0) tb1.miss <- c( tb1.miss, tb0[i,1] );
   			if (length(tb1.mat)>1) tb1.multi <- c( tb1.multi, tb0[i,1] );
   			if (length(tb1.mat)==1) { tb0.ord <-c(tb0.ord, i); tb1.ord <- c(tb1.ord, tb1.mat[1]); };
   		}
   		
   		if (length(tb1.miss)>0) warning("Warnings: missing genotype data for id, ", tb1.miss, "\n");
   		if (length(tb1.multi)>0) warning("Warnings: multiple genotype data for id, ", tb1.multi, "\n");

   		tb0 <- tb0[tb0.ord,];
   		tb1 <- tb1[tb1.ord,];
   	}
   	
	dat$phenos_info$ids <- tb0[,1];
   	
	tb0 <- tb0[,-1];
	tb1 <- tb1[,-1];

   	tb0.miss <- c();
   	for( i in 1:dim(tb0)[1] )
   		if (all(is.na(tb0[i,]))) tb0.miss<-c( tb0.miss, i );

	if( length(tb0.miss)>0 )
	{
		tb0 <- tb0[-tb0.miss,];	
		tb1 <- tb1[-tb0.miss,];	
	}

	dat$phenos_info$removed_missing <- dat$phenos_info$ids[ tb0.miss ];
	dat$phenos_info$ids <- dat$phenos_info$ids[ -tb0.miss ];
   	dat$phenos_table <- tb0;
	dat$sample_N <- dim(tb0)[1];
	
	time.str <- colnames(tb0);
	if ( substr(time.str[1],1,1)=="X" )
	{
		time.str<-substring(time.str, 2 );
	}
	time.std <- as.numeric( time.str );
	if ( any( is.na(time.std) ) )
		time.std <- c(1:length(time.std))
	if ( max(time.std) <=1)
		time.std <- time.std* length(time.std);
	
	dat$time.std <- time.std;
	colnames( dat$phenos_table ) <- time.std;
   	
	if (as.numeric(FM2.get_value("log_pheno", def="0")) != 0 )
	{
		dat$phenos_table <- log(dat$phenos_table);
		dat$log <- TRUE;
	}

   	dat$genos_table  <- tb1;
   	tb2              <- read.csv(file=marker_file,sep=",",head=TRUE);
   	dat$marker_table <- tb2[,-1];
   	colnames(dat$marker_table) <- c("Marker", "Dist", "grp_idx", "Group");
   	dat$marker_obj   <- fin.get_markerobj(dat$marker_table);
   	dat$sample_N     <- length( dat$phenos_table[,1] )
   	dat$sample_times <- time.std; ###c(1:length( dat$phenos_table[1,] ) )
	
cat("DATA LOADING", pheno_file, "Sample:", dat$sample_N, "Time", dat$sample_times , "\n");
	
  	return( dat );
  	
}
   	
#--------------------------------------------------------------
# private: dat.summar_dat
#
# summarize the important information for the simulate data or 
# real data.
#
# input : data object.
#--------------------------------------------------------------
dat.summary<-function( dat_obj )
{
	if ( class(dat_obj) != "FM2.dat" )
	{
		sErrMsg<- "Error: Not a dat set for FunMap model.";
		stop( sErrMsg );
	}
	
	strt <- sprintf("The data set for FunMap model:\n");
	stru <- sprintf("-----------------------------------\n");
	str0 <- sprintf("%15s: %s\n", 	 "Date", 	   Sys.time() );
	str1 <- sprintf("%15s: %s\n", 	 "Pheno. file",	   dat_obj$pheno_file );
	str3 <- sprintf("%15s: %-10.0f\n", "Sample size",  dat_obj$sample_N );
	str4 <- sprintf("%15s: %-10.0f\n", "Sample times", length(dat_obj$phenos_table[1,]) );
	str7 <- sprintf("%15s: %s\n", 	 "Curve Model",	   dat_obj$curve$name );
	str8 <- sprintf("%15s: %s\n", 	 "Covariance", 	   dat_obj$covar$name);
	
	stra <- "";
	if (!is.null( FM2.curve$get_summary ) && !is.null( dat_obj$curve ) )
		stra <- FM2.curve$get_summary(dat_obj$curve, "%15s: %-10.5f\n" );

	strb <- "";
	if (!is.null( FM2.covar$get_summary ) && !is.null( dat_obj$covar ) )
		strb <- FM2.covar$get_summary(dat_obj$covar, "%15s: %-10.5f\n" );
		
	strd <- sprintf("------------------------------------\n\n");

	str <- paste(strt, stru, str0, str1, str3, str4, 
				 str7, stra, str8, strb, strd, sep="" );

	# Figure 1
	strFile <- fpt.plot_to_doc( dat_obj$pheno_file, FM2.get_value("plot_doctype", def="pdf") )
	err.fig1 <- try( fpt.plot_tiled_curves( dat_obj, max_curves=8*8 ) );
	if ( class(err.fig1) != "try-error" )
		title(paste("The ",dat_obj$curve$name," for all individuals.", sep=""));
	dev.off();
	
	str <- paste(str,"\n*1:The figure 1 for all individuals is saved to ", strFile, ".\n", sep="");
	
	# Figure 2
	strFile <- fpt.plot_to_doc( dat_obj$pheno_file, FM2.get_value("plot_doctype", def="pdf") )
	err.fig2 <- try( fpt.plot_overlapping_curves( dat_obj ) );			 
	if ( class(err.fig2) != "try-error" )
		title(paste("The ", dat_obj$curve$name," for all individuals.", sep=""));
	dev.off();

	str <- paste(str,"*2:The figure 2 for all individuals is saved to ", strFile, ".\n\n", sep="");

	return (str);
}

#--------------------------------------------------------------
# private: dat.report_dat
#
# summarize the important information for the simulate data or 
# real data.
#
# input : data object.
#--------------------------------------------------------------
dat.report<-function( dat_obj )
{
	str1 <- sprintf("%15s: %s\n", 	 "Pheno. file",	 dat_obj$pheno_file );
	str2 <- sprintf("%15s: %s\n", 	 "Geno. file", 	 dat_obj$geno_file );
	str3 <- sprintf("%15s: %s\n", 	 "Maker file", 	 dat_obj$marker_file );
	str4 <- sprintf("%15s: %-10.0f\n", "Sample size",  dat_obj$sample_N );
	str5 <- sprintf("%15s: %-10.0f\n", "Sample times", length(dat_obj$phenos_table[1,]) );
	str6 <- sprintf("%15s: %-10.0f\n", "Marker count", length(dat_obj$marker_table[,1]) );
	str7 <- sprintf("%15s: %s\n", 	 "Cross", 	 FM2.cross$name );
	str8 <- sprintf("%15s: %s\n", 	 "Model",	 dat_obj$curve_name );
	str9 <- sprintf("%15s: %s\n", 	 "Model",	 dat_obj$curve$name );

	stra <- "";
	if (!is.null( FM2.curve$get_summary ) && !is.null( dat_obj$curve ) )
		stra <- FM2.curve$get_summary(dat_obj$curve, "%15s: %-10.5f\n" );

	str10 <- sprintf("%15s: %s\n", 	 "Covariance", 	 dat_obj$covar$name );
	strb <- "";
	if (!is.null( FM2.covar$get_summary ) && !is.null( dat_obj$covar ) )
		strb <- FM2.covar$get_summary(dat_obj$covar, "%15s: %-10.5f\n" );

	str <- paste(str1, str2, str3, str4, str5, str6, str7, str8, str9, stra, str10, strb, sep="" );

	# Figure 1
	plot.old <- FM2.get_value("plot_doctype", def="pdf");
	FM2.set_value("plot_doctype", "png");
	fig1 <- fpt.plot_to_doc( dat_obj$pheno_file, FM2.get_value("plot_doctype") )
	err.fig1 <- try( fpt.plot_tiled_curves( dat, max_curves=8*8 ) );
	if ( class(err.fig1) != "try-error" )
		title(paste("The ",dat_obj$curve_name," for all individuals.", sep=""));
	dev.off();
	
	# Figure 2
	fig2 <- fpt.plot_to_doc( dat_obj$pheno_file, FM2.get_value("plot_doctype") )
	err.fig2 <- try( fpt.plot_overlapping_curves( dat ) );			 
	if ( class(err.fig2) != "try-error" )
		title(paste("The ", dat_obj$curve_name," for all individuals.", sep=""));
	dev.off();

	FM2.set_value("plot_doctype", plot.old);

	return (list(str, fig1, fig2) );
}

#--------------------------------------------------------------
# dat.plot
# 
# Plot the Pharmacology curve for data object. 
# Used by plot( EXP.dat object );
#--------------------------------------------------------------
dat.plot<- function( dat_obj, plot_type=NA  )
{
	if( is.na(plot_type) || plot_type==1)
	{
		x11();
		err.fig1 <- try( fpt.plot_tiled_curves( dat_obj ) );
		if ( class(err.fig1) != "try-error" )
			title(paste("The ", dat_obj$curve_name, " for all individuals.", sep=""));
	}
	
	if( is.na(plot_type) || plot_type==2 )
	{
		x11();
		err.fig2 <- try( fpt.plot_overlapping_curves( dat_obj ) );			 
		if ( class(err.fig2) != "try-error" )
			title(paste("The ", dat_obj$curve_name," for all individuals.", sep=""));
	}
	
	if( is.na(plot_type) || plot_type==2 )
	{
		x11();
		err.fig3 <- try( fm.genetic.map( dat_obj ) );			 
		if ( class(err.fig3) != "try-error" )
			title(paste("The genetic mapping"));
	}
	
	if( is.na(plot_type) || plot_type==2 )
	{
		x11();
		err.fig4 <- try( fpt_plot_genotype( dat_obj ) );			 
		if ( class(err.fig3) != "try-error" )
			title(paste("The marker mapping "));
	}		
}

#--------------------------------------------------------------
# public: fin.get_traits_mu
#
#--------------------------------------------------------------
fin.get_traits_mu<-function( par )
{
	mu    <- c();
	if (!is.null(par$qq0) )	  
		mu    <- rbind( mu, FM2.curve$get_mu( par$qq0, par$simu_times ) )
	else
		mu    <- rbind( mu, rep(NA, length(par$simu_times) ) );
		
	if (!is.null(par$Qq1) )
		mu    <- rbind( mu, FM2.curve$get_mu( par$Qq1, par$simu_times ) )
	else
		mu    <- rbind( mu, rep(NA, length(par$simu_times) ) );

	if (!is.null(par$QQ2) )
		mu    <- rbind( mu, FM2.curve$get_mu( par$QQ2, par$simu_times ) )
	else
		mu    <- rbind( mu, rep(NA, length(par$simu_times) ) );
	
	return(mu);
}
		
#--------------------------------------------------------------
# private: fin.recalc_marker
# 
# make a marker object according to maker table of the data set
#
# input : marker table
# output: marker object
#--------------------------------------------------------------
fin.get_markerobj<-function( marker_table )
{
	marker_obj<-list(
			count = NULL,
			grps  = list() );
	
	marker_grp<-list(
			index     = -1,
			id        = "",
			count     = NULL,
			start_idx = NULL,
			dists     = list(),
			names     = list() );

	grp_idx <- -1;
	index   <- 1;
	index2  <- 1;
	dists   <- list();
	names   <- list();
	
	for (i in 1:length(marker_table[,3]))
	{
		if (marker_table[i,3]>grp_idx)
		{
			grp_idx <- marker_table[i,3];

			marker_obj$grps[[index]]       <- list();
			marker_obj$grps[[index]]$index <- marker_table[i,3];
			marker_obj$grps[[index]]$id    <- marker_table[i,4];
			marker_obj$grps[[index]]$start_idx <- i;

			if (index>1)
			{
				marker_obj$grps[[index-1]]$dists <- dists;
				marker_obj$grps[[index-1]]$names <- names;
				marker_obj$grps[[index-1]]$count <- index2-1;
				dists  <- list();
				names  <- list();
				index2 <- 1;
			}
			index <- index+1;
		}

		dists[[index2]] <- marker_table[i,2];
		names[[index2]] <- marker_table[i,1];
		index2 <- index2+1;
	}

	marker_obj$grps[[index-1]]$dists <- dists;
	marker_obj$grps[[index-1]]$names <- names;
	marker_obj$grps[[index-1]]$count <- index2-1;
	marker_obj$count <- index-1;
	
	return (marker_obj);
}
