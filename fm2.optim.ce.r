optim_ce<-function(parin, mlefunc, y, time.std, qtl.prob=NULL )
{
	option.N <- 30*length(parin);
	option.beta <- 0.9;
	option.alpha <- 0.3
	option.eps <- 0.1;
	option.h <- 10;
	option.q <- 6;
	option.inject <- 3;
	option.N.elite <- option.N*0.2;
	
	par.n   <- length(parin);
	par.mu  <- array(0, dim=c(0, par.n))
	par.sg2 <- array(0, dim=c(0, par.n))
	par.mu  <- rbind( par.mu, parin);
	par.sg2 <- rbind( par.sg2, c(0.5, abs(parin[-1]*10) ) );
	par.minL <- c();
	level.t  <- 1;
	inject.no<- 0;
	inject_loop<-0;

	sample.Ps<- array(0, dim=c(0, par.n))
	sample.Ls<- c();
	
	h.ret <- list(value=Inf, par=parin);

	stop.iteration <- FALSE;
	while( !stop.iteration )
	{
		smooth.mu  <- par.mu[ level.t, ]
		smooth.sg2 <- par.sg2[ level.t, ]
		if (level.t>1)
		{
			smooth.mu  <- option.alpha*par.mu[ level.t, ] + (1 - option.alpha)*par.mu[ level.t-1, ];
			beta.t <- option.beta - option.beta*( 1 - 1/level.t )^option.q;
			smooth.sg2 <- beta.t * par.sg2[ level.t, ] + (1 - beta.t) * par.sg2[ level.t-1, ];
		}

		if (max(smooth.sg2)<=option.eps*0.5^inject.no || inject_loop>=100)
		{

cat("[", inject.no, "]", smooth.sg2, "\n");			
			if ( inject.no >= option.inject )
			{
				stop.iteration <- TRUE;
				h.ret$par <- par.mu[ level.t, ]
				h.ret$value <- 0;
				if (is.null(qtl.prob))
					h.ret$value<- mlefunc(h.ret$par, y, time.std)
				else
				{
					h.ret$value<- mlefunc(h.ret$par, y, qtl.prob, time.std);
					browser();				
				}					

				FM2.set_value( "mle.pdf", 0);	
				return(h.ret);
			}

			smooth.sg2 <- smooth.sg2 + abs(par.minL[level.t-1]-par.minL[level.t-2])*option.h;
			inject.no <- inject.no + 1;
			inject_loop <- 0;
		}

		#sample.Ps<- array(0, dim=c(0, par.n))
		#sample.Ls<- c();
		i<-1;
		sample.min <- Inf;
		while (i<option.N)
		{
			sample.p <- rnorm( par.n, smooth.mu, smooth.sg2 );
			if (sample.p[1] >= 1.0) sample.p[1] <- 0.9;
			if (sample.p[1] <= 0.0) sample.p[1] <- 0.1;
			sample.p[2] <- abs(sample.p[2]);
			
			L.value <- Inf;
			if (is.null(qtl.prob))
				L.value<- mlefunc(sample.p, y, time.std)
			else
				L.value<- mlefunc(sample.p, y, qtl.prob, time.std, );
			
			if ( is.nan(L.value) || L.value == Inf)
				next;

			if (L.value <= 0 || FM2.get_value( "mle.pdf",0) >= 1 )
				next;
			
			sample.Ps <- rbind(sample.Ps, sample.p);
			sample.Ls <- c(sample.Ls, L.value);
			if (sample.min>L.value ) sample.min <- L.value

			i <- i+1;
		}

		ord.Ls0 <- order( sample.Ls, decreasing = FALSE);
		ord.Ls  <- ord.Ls0[1:option.N.elite]
		new.mu  <- colMeans( sample.Ps[ord.Ls, ] );
		new.sg2 <- apply( sample.Ps[ord.Ls, ], 2, function(x){sd(x, na.rm=T)} );
		par.mu  <- rbind(par.mu, new.mu);
		par.sg2 <- rbind(par.sg2, new.sg2);
		par.minL<- c(par.minL, sample.min );

		sample.Ps <- sample.Ps[ord.Ls[1:5],, drop=F];
		sample.Ls <- sample.Ls[ord.Ls[1:5]];
		
cat("\nLEVEL", level.t, sample.Ls, "--", sort(new.sg2,decreasing=T),"\n");
		level.t <- level.t + 1;
		inject_loop <- inject_loop+1;
	}

	return(h.ret);
}
