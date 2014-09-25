rm(list=ls(all=TRUE));

source("com.sys.r");
source("com.report.r");
source("fm2.api.r");
source("fm2.dat.r");
source("fm2.cross.r");
source("fm2.covar.mar1.r");
source("fm2.plot.r");
source("fm2.report.r");
source("fm2.curve.base.r");
#source("fm2.curve.ext.r");
source("fm2.curve.lc.r");
source("fm2.curve.mlog.r");
source("fm2.qtlmodel.r");
source("fm2.permu.r");
source("fm2.eval.r");
source("zzz.r");

.onAttach();

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Test code 2: LC_BC_file_test
#
# Abstract: LC, Backcross, real data(populus)
# 
# The result is saved into [LC_BC_file_test.rdata]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FM2.LC_BC_file_test<-function()
{
	dat<- FM2.load_data( "test/LC.populus.pheno.csv", 
			     "test/LC.populus.geno.csv", 
			     "test/LC.populus.marker.csv", CURVE_LC, CROSS_BC, COVAR_AR1 );
	summary(dat);
	#try( plot(dat) );

	ret<- FM2.qtlscan( dat );
	summary(ret, dat);
	try( plot(ret, dat) );

	FM2.report( "LC_BC_file_test.pdf", dat, ret);

	##Step1
	save(dat, ret, file="LC_BC_file_test.rdata");

	ret_perm<- FM2.permutation( dat, options=list(permu_loop=100) );
	summary(ret_perm);
	try(plot(ret_perm));

	save(dat, ret, ret_perm, file="LC_BC_file_test.rdata");

	return(ret);
}

#ret1 <- FM2.LC_BC_file_test();
#ret1 <- FM2.simu_test( par_LC, CURVE_LC, CROSS_BC, COVAR_AR1);
#ret1 <- FM2.simu_test( par_LC, CURVE_LC, CROSS_F2, COVAR_AR1);
#ret1 <- FM2.simu_test( par_BI, CURVE_BI, CROSS_BC, COVAR_AR1);
#ret1 <- FM2.simu_test( par_BI, CURVE_BI, CROSS_F2, COVAR_AR1);
#ret1 <- FM2.simu_test( par_PC, CURVE_PC, CROSS_BC, COVAR_AR1);
#ret1 <- FM2.simu_test( par_EXP, CURVE_EXP, CROSS_F2, COVAR_AR1);
#ret1 <- FM2.simu_test( par_NP, CURVE_NP, CROSS_F2, COVAR_AR1);
#ret1 <- FM2.simu_test( par_LC, CURVE_LC, CROSS_RIL, COVAR_AR1);
#ret1 <- FM2.qtlscan("LC.populus.pheno.csv", "LC.populus.geno.csv", "LC.populus.marker.csv", CURVE_LC, CROSS_BC, COVAR_AR1)

FM2.set_value("optim", "Nelder-Mead");
ret1 <- FM2.simu_test( par_MLC, CURVE_MLC, CROSS_BC, COVAR_MAR1)

