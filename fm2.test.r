rm(list=ls(all=TRUE));

source("com.sys.r");
source("com.report.r");
source("fm2.api.r");
source("fm2.dat.r");
source("fm2.plot.r");
source("fm2.report.r");
source("fm2.cross.r");
source("fm2.covar.ar1.r");
#source("fm2.covar.sad1.r");
#source("fm2.covar.ar1h.r");
#source("fm2.covar.si.r");
#source("fm2.covar.cscm.r");
#source("fm2.covar.csh.r");
#source("fm2.covar.dg.r");
#source("fm2.covar.arma11.r");
#source("fm2.covar.fa1.r");
#source("fm2.covar.hf.r");
#source("fm2.covar.cs.r");
#source("fm2.covar.sad1.r");
source("fm2.covar.base.r");
source("fm2.curve.base.r");
#source("fm2.curve.ext.r");
source("fm2.qtlmodel.r");
source("fm2.permu.r");
source("zzz.r");
source("fm2.covar.select.r");
source("fm2.curve.lc.r");
#source("fm2.curve.bi.r");
#source("fm2.curve.exp.r");
#source("fm2.curve.line.r");
#source("fm2.curve.pl.r");
#source("fm2.curve.hoss.r");
#source("fm2.curve.mono.r");
#source("fm2.curve.weib.r");
#source("fm2.curve.korf.r");
#source("fm2.curve.gomp.r");
source("inter.dat.r");
source("interpermu.r");
source("interval.qtl.r");
source("intermap.r");
source("inter.base.r");


.onAttach();


#ret1 <- FM2.simu_test( CROSS_BC );
#ret1 <- AM.simu_test( CROSS_BC );

ret1 <- FM2.one_step("test/LC.populus.pheno.csv", "test/LC.populus.geno.csv", "test/LC.populus.marker.csv", CURVE_LC, CROSS_BC, COVAR_AR1)
#ret1 <- AM.one_step("test/Inter.populus.pheno.csv", "test/LC.populus.geno.csv", "test/LC.populus.marker.csv", CROSS_BC);
#ret1 <- FM2.one_step("test/LC.mouse.2.pheno.csv", CURVE_LC, COVAR_SAD1)
#ret1 <- FM2.one_step("soy-leaf/2006-leafonly.csv", CURVE_LC, COVAR_SAD1)
#ret1 <- FM2.one_step("test2/grainwt.ril.phe.csv", CURVE_LC, COVAR_SAD1, options=list(prob=0.02) )
#ret1 <- FM2.one_step("test2/stemweight.ril.phe.csv", CURVE_LC, COVAR_SAD1, options=list(prob=0.02) )
#ret1 <- FM2.one_step("test3/Pheno 1997 circonf ss na.csv", CURVE_LC, COVAR_SAD1, options=list(prob=0.01) )