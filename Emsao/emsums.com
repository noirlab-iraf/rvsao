common /emsums/ sumvel, sumerr, sumpix, sumdwl, spx, swl,
		evel, epix, ewl, nap
double  sumvel, sumerr, sumpix, sumdwl, spx, swl
pointer	evel, epix, ewl
int	nap	# Number of apertures for which shift has been computed
