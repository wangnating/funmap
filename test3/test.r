rm(list=ls(all=TRUE));

source("..//com.sys.r");
source("..//com.report.r");
source("..//fm2.api.r");
source("..//fm2.dat.r");
source("..//fm2.cross.r");
source("..//fm2.covar.ar1.r");
source("..//fm2.covar.sad2.r");
source("..//fm2.plot.r");
source("..//fm2.report.r");
source("..//fm2.curve.base.r");
source("..//fm2.curve.ext.r");
source("..//fm2.qtlmodel.r");
source("..//fm2.permu.r");
source("..//fm2.eval.r");
source("..//fm2.optim.ce.r");
source("..//zzz.r");

.onAttach();

phe.csv  <- "Pheno 1997 circonf ss na.csv";
gen.csv  <- "Geno URO P1997 LG ref ss na.csv";
mark.csv <- "Carte URO P1997 LG ref.csv";

test_bc<-function()
{
	FM2.set_value("peak_count", 5);
	FM2.set_value("optim", "CE");
	#FM2.set_value("log_pheno", 0)
	
	dat<-FM2.load_data(phe.csv, gen.csv, mark.csv, CURVE_NP, cross_type=CROSS_BC );
	#summary(dat);
	try( plot(dat) );

	ret<- FM2.qtlscan( dat );
	summary(ret, dat);
	try( plot(ret, dat) );

	##Step1
	save(dat, ret, file="test3.rdata");

	FM2.report("test3.pdf", dat, ret);

	ret_perm<- FM2.permutation( dat, options=list(permu_loop=100) );
	summary(ret_perm);
	try(plot(ret_perm));

	save(dat, ret, ret_perm, file="test3.rdata");

	return(ret);

}

ret.bc <-test_bc();

