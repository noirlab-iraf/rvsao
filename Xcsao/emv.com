# File rvsao/lib/emv.com
# January 26, 2007

#  Common blocks for emission line velocity pacakge

common/emv/ emparams,wlrest,pxobs,wlobs,wtobs,wspan,zsig,disperr,gvel,mincont,
	    nfound,nfit,npfit,nlcont,esmooth,override,linedrop,nmobs,
	    vplot,lfit,vfit,sfit,xfit,skyspec,dispem,dispabs,plotcorr

double	emparams[10,2,MAXREF]	# Parameters of emission line fits
			# 1,2,3 are polynomial coefficients
			# 4,5,6 are Gaussian coefficients
			# 7 is equivalent width
			# 8 is chi^2 and degrees of freedom
			# 9 is velocity in z
			# 10 is acceptance flag (0 or velocity error)
double	wlrest[MAXREF]	# Rest wavelengths of observed emission lines
double	pxobs[MAXREF]	# Pixel centers of observed emission lines
double	wlobs[MAXREF]	# Observed wavelengths of emission lines
double	wtobs[MAXREF]	# Weights for observed emission lines
double	wspan		# Delta wavelength  to search for line
double	zsig		# Minimum number of sigma for emission line
double	disperr		# Mean dispersion error in km/sec
double	gvel		# Guessed velocity
double	mincont		# Minimum continuum for equivalent width
int	nfound		# Number of emission lines found
int	nfit		# Number of emission lines fit
int	npfit		# Number of points to fit for line peak
int	nlcont		# Number of points to fit for line continuum
int	esmooth		# Number of times to smooth data for fit
int	override[MAXREF] # Override swiches for observed lines
int	linedrop[MAXREF] # Error flags for dropped observed lines
char	nmobs[SZ_ELINE,MAXREF]	# Names of observed emission lines
int	vplot		# Switch to select which velocity to plot
bool	lfit		# Switch to retry line search and fit
bool	vfit		# Switch to redo velocity fit
bool	sfit		# Switch to redo velocity search, line search, and fit
bool	xfit		# Switch to rerun cross-correlation (used in XCSAO)
bool	skyspec		# True if sky spectrum is present
bool	dispem		# True if emission lines are to be marked
bool	dispabs		# True if absorption lines are to be marked
bool	plotcorr	# True to plot correlation after labelled spectrum


# Emission lines for velocity guess

common/search/ restwave, bluelim, redlim, nsearch, nmsearch
double	restwave[MAXSEARCH]	# Rest wavelengths
double	bluelim[MAXSEARCH]	# Blue wavelength limit
double	redlim[MAXSEARCH]	# Red wavelength limit
int	nsearch			# Number of emission lines to try
char	nmsearch[SZ_ELINE,MAXSEARCH]	# Names of emission lines


# Reference emission lines
common/emlines/ wlref, bcont, rcont, wgfit, nref, nmref, plotem, lfound
double	wlref[MAXREF]		# Wavelengths of reference emission lines
double	bcont[MAXREF]		# Blue limit of region to sample continuum RMS
double	rcont[MAXREF]		# Red limit of region to sample continuum RMS
double	wgfit[MAXREF]		# Half-width of region to fit in angstroms
int	nref			# Number of reference emission lines
char	nmref[SZ_ELINE,MAXREF]		# Names of reference emission lines
bool	plotem			# True to always label emission lines
int	lfound[MAXREF]		# Index of line in list of lines found and fit


# Lines to be fit together, if possible 
common/combo/ elines, eht, edwl, numcom, ncombo, edrop
double	elines[MAXCLINES,MAXCOMB]	# line center wavelengths
double	eht[MAXCLINES,MAXCOMB]		# line heights
double	edwl[MAXCOMB]			# Wavelength to fit beyond each edge 
int	numcom[MAXCOMB]			# Number of lines in tuple 
int	ncombo				# Number of combination lines
bool	edrop[MAXCLINES,MAXCOMB]	# Flag to allow dropping unfound lines


# Reference absorption lines
common/ablines/ wlabs, nabs, nmabs
double	wlabs[MAXABS]		# Wavelengths of absorption lines
int	nabs			# Number of reference absorption lines
char	nmabs[SZ_ELINE,MAXABS]		# Names of absorption lines

# Aug  3 1994	Add dispersion error as velocity (for emission error)
# Aug  4 1994	Parameterize length of line names

# Sep 19 1995	Add correlation plot flag PLOTCORR
# Oct  2 1995	Add inital velocity guess GVEL

# Feb 14 1997	Add unfound combo lines flag
# Feb 19 1997	Add flag for why lines were dropped
# May  9 1997	Add minimum continuum for equivalent width

# May 12 2004	Add plotem to /emlines/ common to always label emission lines

# Jan 26 2007	Add lfound to emlines common to index lines which were found
