#  File src/library/grDevices/R/zzz.R

.noGenerics <- TRUE

msg <- function(...)
{
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
 
    cat("##\n## Funmap Package v.2.2-1\n")
    cat("## Build date: ", date(), "\n")
    cat("## Copyright (C) 2011-", yr, ", http://statgen.psu.edu\n", sep="")
    cat("## Written by Zhong Wang(zzw2@psu.edu)\n\n")
}

.onAttach <- function(...)
{
	msg();
	#ZZZ.regcovar_bi();
	ZZZ.regcovar_lc();
	#ZZZ.regcovar_exp();
	#ZZZ.regcovar_line();
	#ZZZ.regcovar_pl();
	#ZZZ.regcovar_hoss();
	#ZZZ.regcovar_mono();
	#ZZZ.regcovar_weib();
	#ZZZ.regcovar_korf();
	#ZZZ.regcovar_gomp();
	ZZZ.regcross();
	ZZZ.regcovar_ar1();
	#ZZZ.regcovar_sad1();
	#ZZZ.regcovar_ar1h();
	#ZZZ.regcovar_si();
	#ZZZ.regcovar_cscm();
	#ZZZ.regcovar_csh();
	#ZZZ.regcovar_dg();
	#ZZZ.regcovar_arma11();
	#ZZZ.regcovar_fa1();
	#ZZZ.regcovar_hf();
	#ZZZ.regcovar_cs();
	
}


.onLoad <- function(libname, pkgname)
{
}


