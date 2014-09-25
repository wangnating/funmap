#################################################################
#
# New Systems Mapping Application(SysMap1)
#
# New Curve demostration:
#
# MODEL_FOO:  FOO method 
#
# History:
# 12/15/2011 Version 1.1
#
##################################################################

source("com.sys.r");
source("com.report.r");
source("fm2.api.r");
source("fm2.com.r");
source("fm2.inner.r");
source("fm2.plot.r");
source("fm2.report.r");

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Foo Curve
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
simu_Foo<-list(
	#rho & sigma^2
	simu_rho = 0.7543,
	simu_s2  = 1.2260,

	#sample size
	simu_N      = 100,
	#marker distance
	simu_mrkdist= c(0, 20, 20, 20, 20, 20),
	#treat measured in time 1 to 7
	simu_times  = c(1:7),
	#qtl position
	simu_qtlpos = 50,			
	
	# 	QQ2         
	QQ2 = list(
		simu_a   = -0.8,
		simu_b   = 0.9,
		simu_c   = 38),

	# 	Qq1         
	Qq1 = list(
		simu_a   = -0.75,
		simu_b   = 0.92,
		simu_c   = 38.1),
	
	#	qq0
	qq0 = list(
		simu_a   = -0.84,
		simu_b   = 0.85,
		simu_c   = 37)
);

Foo_get_mu <- function(par, times, options=list())
{
	return ( par[1]*times^2+ par[2]*times+par[3] );
}

model_FOO<-list(
	type 	= 1,
	name 	= "FOO",
	desc 	= "FOO Curve",
	simu_par= simu_Foo,
	par_num = 3, 
	get_mu  = Foo_get_mu,
	get_est_par = NULL);


MODEL_FOO <<- SM.reg_model( model_FOO );


par<-FM2.param(MODEL_FOO, CROSS_F2);
summary(par);
dat<-FM2.simulate(par);
summary(dat);
ret<-FM2.hp_test(dat);
summary(ret, dat);
FM2.report("foo.test.pdf", dat, ret);
ret.eval2 <- FM2.evaluate(par, 10);
summary(ret.eval2);
