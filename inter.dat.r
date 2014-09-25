#################################################################
#
# New Systems Mapping Application(SysMap1)#
# Common library
#
#    1) inter_par
#    2) interpar.get_summary
#    3) interdat.simulate
#    4) interdat.load
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################

#--------------------------------------------------------------
# inter_par
#
# parameter object for Backcross
#
# Used by simulating a data 
#
#--------------------------------------------------------------
inter_par <- list(
		#sample size
		simu_N      = 100,
		#marker distance
		simu_mrkdist= c(0, 20, 20, 20, 20, 20),
		#qtl position
		simu_qtlpos = 50,
		
		simu_cross = list(),
		
		sigma = list(
			1.2  ),
		
		# QQ2         
		QQ2 = list(
			mu = 13.8823 ),
		
		# Qq1         
		Qq1 = list(
			mu = 13.2345 ),
			
		# qq0
		qq0 = list(
			mu =  5.3245 )
);

#--------------------------------------------------------------
# interpar.get_summary
#
# Used by summary( parameter)
#--------------------------------------------------------------
interpar.get_summary <-function( par_obj, format)
{
	str0 <- "";
	if (!is.null(par_obj$val))
		str0 <- sprintf(format, "-Lval",  par_obj$val );
	
	str1 <- sprintf(format, "sigma", par_obj$sigma );
	str2 <- sprintf(format, "QQmu2",  par_obj$QQ2 );
	str3 <- sprintf(format, "Qqmu1",  par_obj$Qq1 );
	str4 <- sprintf(format, "qqmu0",  par_obj$qq0 );
	
	return(paste(str0, str1,str2, str3, str4, sep=""));
}

#--------------------------------------------------------------
# interdat.simulate
#
# Simulate a data set for backcross, including phenotype, genotype, 
# marker table.
#
# Input  
#       parameter object for Backcross
# output:
#       data object
#--------------------------------------------------------------
interdat.simulate <- function( par_obj )
{
	dat<-list(
		cross_type   = par_obj$cross_type,
		sample_N     = par_obj$simu_N,
		#trait_num   = 1,
		name         = par_obj$crossname,
		pheno_file   = paste( "simu.pheno", FM2.cross$name, sep="."),
		geno_file    = paste( "simu.geno",  FM2.cross$name, sep="."),
		marker_file  = paste( "simu.marker", FM2.cross$name, sep="."), 
		phenos_table = NULL,
		genos_table  = NULL,
		marker_table = NULL,
		marker_obj   = NULL);
	
	dat$genos_table <- FM2.cross$get_simu_marker( par_obj$simu_N, par_obj$simu_mrkdist, par_obj$simu_cross );
	
	mk_name  <- c();
	mk_dist  <- c();
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
 	
 	sim.mu <- inter.get_traits_mu( par_obj );
 	sim.sigma <- unlist( par_obj$sigma );
 	
 	idx     <- max(which( mrkplace < par_obj$simu_qtlpos ));
	qtlmrk1 <- idx[1];
	qtlmrk2 <- idx[1]+1;
	
	dat$phenos_table<- array(0, dim=c(par_obj$simu_N, 1) );
	gen.qtl <- FM2.cross$get_simu_qtl( par_obj$simu_N, dat$genos_table[,qtlmrk1], dat$genos_table[,qtlmrk2],
						par_obj$simu_qtlpos, mrkplace[qtlmrk1], mrkplace[qtlmrk2], par_obj$simu_cross );

	print(sim.mu);
	for (i in 1:par_obj$simu_N)
	{
		 y <- rnorm(1,sim.mu[ gen.qtl[i]], sim.sigma );
		 dat$phenos_table[i] <- y;
	}
	
	cat("Data simulation is done!\n");
	return(dat);
}

#--------------------------------------------------------------
# public: interdat.load
#
# load a real data set for backcross and F2 experiment.
# 
# input:
#  file : pheno_file, geno_file, marker_file
#  cross: cross type, BC or F2
#--------------------------------------------------------------
interdat.load <- function( pheno_file, geno_file, marker_file, head=TRUE )
{

	interdat<-list(
		#name         = paste(FM2.curve$name, ".dat",sep=""),
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
   	
   	interdat$phenos_info$ids <- tb0[,1];
   	tb0 <- tb0[,-1,drop=F];
	tb1 <- tb1[,-1];
	
	tb0.miss <- c();
	
	for( i in 1:length(tb0) )
   		if (all(is.na(tb0[i]))) tb0.miss<-c( tb0.miss, i );
   		
   	if( length(tb0.miss)>0 )
	{
		tb0 <- tb0[-tb0.miss];	
		tb1 <- tb1[-tb0.miss,];	
		interdat$phenos_info$removed_missing <- interdat$phenos_info$ids[ tb0.miss ];
		interdat$phenos_info$ids <- interdat$phenos_info$ids[ -tb0.miss ];
	}
	
   	interdat$phenos_table <- unlist(tb0);
   	interdat$sample_N <-length(tb0);
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
	
	interdat$time.std <- time.std;
	#colnames( interdat$phenos_table ) <- time.std;
	
	if (as.numeric(FM2.get_value("log_pheno", def="0")) != 0 )
	{
		interdat$phenos_table <- log(interdat$phenos_table);
		interdat$log <- TRUE;
	}
	
	interdat$genos_table  <- tb1;
	tb2  <- read.csv(file=marker_file,sep=",",head=TRUE);
	interdat$marker_table <- tb2[,-1];
	colnames(interdat$marker_table) <- c("Marker", "Dist", "grp_idx", "Group");
	interdat$marker_obj   <- fin.get_markerobj(interdat$marker_table);
	interdat$sample_N     <- length( interdat$phenos_table );
   	interdat$sample_times <- time.std; ###c(1:length( interdat$phenos_table[1,] ) );
   	cat("DATA LOADING", pheno_file, "Sample:", interdat$sample_N, "Time", interdat$sample_times , "\n");
		
	return( interdat );
	  	
}

#--------------------------------------------------------------
# public: fin.get_traits_mu
#
#--------------------------------------------------------------
inter.get_traits_mu <- function( par )
{
	mu <- c();
	if (!is.null(par$qq0) )
		mu	<- rbind( mu,par$qq0$mu )
	else
		mu	<- rbind( mu, NA );
	
	if (!is.null(par$Qq1) )
		mu	<- rbind( mu,par$Qq1$mu )
	else
		mu	<- rbind( mu, NA );
		
	if (!is.null(par$qq0) )
		mu 	<- rbind( mu,par$QQ2$mu )
	else
		mu	<- rbind( mu, NA );
		
	return(mu);
}