source("..//sys-mgr.r");
source("..//Plot.r");
source("..//COM.model.r");
source("..//LC.model.r");
source("..//BI.model.r");
source("..//PC.model.r");
source("..//NP.model.r");
source("..//SG.model.r");
source("..//CM.model.r");
source("..//Funmap.r");
source("..//report.r");


test_ril<-function( pheno_file )
{
	FM2.set_value("peak_count", 5);
	FM2.set_value("log_pheno", 1)
	
	dat<-FM.load_data(pheno_file, "simurima.ril.geno.csv", "simurima.ril.map.csv", MODEL_LC, cross_type=CROSS_RIL );
	#summary(dat);
	ret<-FM.hp_test(dat, test=c(10));
	#summary(ret, dat);
	
	rdata <- sub(".csv", paste("-r", ".rdata", sep=""), pheno_file, fixed = TRUE); 
	save(dat, ret, file=rdata );

	#make a report
	FM.report( dat, ret );
}


permu_ril<-function(pheno_file)
{
	FM2.set_value("peak_count", 5);
	
	dat<-FM.load_data(pheno_file, "simurima.ril.geno.csv", "simurima.ril.map.csv", MODEL_LC, cross_type=CROSS_RIL );
	#summary(dat);
	ret<-FM.permutation(dat, options=list(permu_loop=100, debug=1) );

	rdata <- sub(".csv", paste("-permu-r", np.order, ".rdata", sep=""), pheno_file, fixed = TRUE) 
	save(dat, ret, file=rdata );

	summary(ret);
	plot(ret);
	
	return(ret);
}

#ril1 <-test_ril("grainwt.ril.phe.csv");
ril2 <-test_ril("stemweight.ril.phe.csv");
#perm5 <-permu_ril( "grainwt.ril.phe.csv" );
