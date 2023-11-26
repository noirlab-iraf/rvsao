# File rvsao/lib/rvsao.com
# By Jessica Mink, Smithsonian Astrophysical Observatory
# May 1, 2014

# Common task parameters for XCSAO, EMVEL, LINESPEC, SUMSPEC, and EQWIDTH
 
# Parameters that usually stay fixed
 
common /fixp/ pkfrac, tshift, minvel, maxvel, cvel, c0, npts, han, zpad,
	      lo,toplo,nrun,topnrn, filter,xcrdif,xcr0,sfchop,schop,svcor,
	      pkmode, pkmode0,npkmode,pksrch,contfit,contdim,nmspec,specnums,
	      nsmooth,correlate,vinit,taskname,cortemp
 
double	pkfrac		# <1: fraction of peak, >1: number of peak points to fit
double	tshift		# Night to night zero point velocity shift
double	minvel		# Minimum velocity in km/sec
double	maxvel		# Maximum velocity in km/sec
double	cvel		# Initial value for 1+z
double	c0		# Speed of light
int	npts		# Number of points in the spectrum and transform
real	han		# Fraction of spectrum windowed by a cosine bell
bool	zpad		# Pad transform with an equal number of zeroes
int	lo		# Lowest frequency included in the fit
int	toplo		# Upper frequency of bottom ramp
int	nrun		# Initial number of bins included
int	topnrn		# Upper frequency final cutoff
int	filter		# Filter type
real	xcrdif		# Half-range for velocity in final plot
real	xcr0		# Center velocity in final plot (0=CZ)
int	sfchop		# Flag for emission line chopping in spectrum
bool	schop		# True to chop emission lines from current spectrum
int	svcor		# Velocity correction for spectrum
int	pkmode0		# Peak fit mode (1=parab., 2=quartic, 3=cosx/x, 0=all 3)
int	pkmode		# Current peak fit mode
int	npkmode		# Number of peak fit modes
int	pksrch		# Maximum distance of peak from guess
bool	contfit		# True to use IRAF interactive continuum fitting
int	contdim		# Polynomial divided out of the spectrum if not IRAF
int	vinit		# Switch to select source of initial v guess
int	nsmooth		# Number of times to smooth data for search and plot
int	correlate	# Space for correlation (velocity pixel wavelength no)
int	nmspec		# Number of apertures to read
char	specnums[SZ_LINE]	# List of multispec apertures to read
char	taskname[SZ_LINE]	# Name of IRAF task
char	cortemp[SZ_FNAME]	# Name of template used if VINIT is CORTEMP
 

# Parameters that govern user interaction
 
common /users/ pltspec, plttemp, pltcon, pltapo, pltfft, pltuc, pltcor, plttft,
	       debug, nlogfd, logfd
 
bool	pltspec		# Yes to plot object spectrum
bool	plttemp		# Yes to plot template spectrum
bool	pltcon		# Yes to plot continuum-subtracted spectra
bool	pltapo		# Yes to plot apodized spectra
bool	pltfft		# Yes to plot transforms
bool	pltuc		# Yes to plot unfiltered cross-correlation
bool	pltcor		# Yes to plot cross-correlation
bool	plttft		# Yes to plot transformed transforms
bool	debug		# Yes to print program information
int	nlogfd		# Number of simultaneous log files
int	logfd[MAXLOG]	# File descriptors for log files

 
# Parameters that describe the object spectrum
 
common /fitp/ spstrt,spfnsh,spechcv,spvel,sperr,spxvel,spxerr,spxr,
	      spevel,speerr,spwl0,spdwl,tvel,spmean,waverest,
	      spvqual,specdc,specchop,specpix,spnl,spnlf,specref,
	      z0,sig0,gam0,z,ze, mext,specvb,savevel,renorm,
	      maxlength,specname,specid,skyname,edited,savespec
 
double	spstrt		# Starting log wavelength
double	spfnsh		# Endin log wavelength
double	spechcv		# Spectrum heliocentric velocity
double	spvel		# Spectrum velocity
double	sperr		# Spectrum velocity error
double	spxvel		# Spectrum cross-correlation velocity
double	spxerr		# Spectrum cross-correlation velocity error
double	spxr		# Spectrum cross-correlation r-value
double	spevel		# Spectrum emission line velocity
double	speerr		# Spectrum emission line velocity error
double	spwl0		# Wavelength (or log wavelength) of first pixel
double	spdwl		# Wavelength (or log wavelength) per pixel
double	tvel		# Known part of velocity
double	spmean		# Mean value of input spectrum
double	waverest	# Rest wavelength of line used to set velocity
int	spvqual		# Quality 0=none 1=bad 2,3=questionable 4=good
int	specdc		# 1 = log, 0 = linear
int	specchop	# 1 = emmission lines chopped, 0 = not
int	specpix		# Number of points in spectrum
int	spnl		# Number of emission lines found
int	spnlf		# Number of emission lines in velocity fit
int	specref		# Reference line for spectrum
real	z0		# Initial guess at the redshift
real	sig0		# Initial guess at the velocity dispersion
real	gam0		# Initial guess at the line strength parameter
real	z[NPAR]		# Fitted values of the above
real	ze[NPAR]	# Fitted errors
int	mext		# Image extension number containing spectra
bool	specvb		# True if barycentric velocity from file
bool	savevel		# Save velocity, error, and R in data file header
bool	renorm		# Renormalize spectrum
int	maxlength	# Max. length of printed file name, use default if 0
char	specname[SZ_PATHNAME] # Name of current spectrum
char	specid[SZ_PATHNAME] # Current spectrum file, line, and aperture
char	skyname[SZ_PATHNAME] # Name of sky for current spectrum
bool	edited		# TRUE if spectrum has been changed
bool	savespec	# TRUE to save an edited spectrum in place

# Parameters that describe the template spectra
 
common /ptemp/ testrt,tefnsh,tempvel,temphcv,tempshift,tempwl1,tempwl2,
		twl1,twl2,ntemp,
		tempdc,temppix,tempfilt,tvcor,tfchop,tchop,tcont,npcor,
		overlap,tempname,tempid
 
double	testrt[MAXTEMPS]	# Start wavelength in log space
double	tefnsh[MAXTEMPS]	# Finish wavelength in log space
double	tempvel[MAXTEMPS]	# Observed velocities of each template
double	temphcv[MAXTEMPS]	# Heliocentric velocity of templates
double	tempshift[MAXTEMPS]	# Night to night zero point velocity shift
double	tempwl1[MAXTEMPS]	# Starting wavelengths or pixels of templates
double	tempwl2[MAXTEMPS]	# Ending wavelengths or pixels of templates
double	twl1[MAXTEMPS]		# Starting wavelengths of overlap regions
double	twl2[MAXTEMPS]		# Ending wavelengths of overlap regions
int	ntemp			# Total number of cached templates
int	tempdc[MAXTEMPS]	# 1 = log, 0 = linear
int	temppix[MAXTEMPS]	# Number of points in spectrum
int	tempfilt[MAXTEMPS]	# 1: do not filter template
				# 2: do not high-filter spectrum or template
				# 3: do not high-filter spectrum;
				#    do not filter template at all
int	tvcor			# Velocity correction for templates
int	tfchop			# Flag for emission line chopping in templates
bool	tchop			# True to chop emission lines from this template
bool	tcont			# True to remove continuum from template
int	npcor			# Number of points in correlation transform
bool	overlap			# True to use full overlap pf spectra
char	tempname[SZ_PATHNAME,MAXTEMPS]	# Names of the template objects
char	tempid[SZ_PATHNAME,MAXTEMPS]	# Template filename[aperture]

 
# Parameters that describe the wavelength scales
 
common /wvscl/ dlogw, logw0, delwav, wave0, delpix, deltav, specsh
 
double	dlogw		# Spacing of each pixel in log wavelength space
double	logw0		# Its base 10 log 
double	delwav		# Spacing of each rebinned pixel in wavelength space
double	wave0		# Wavelength of first point of temp file
double	delpix		# Spacing of each rebinned pixel in original pixel space
real	deltav		# Corresponding velocity width
pointer	specsh		# Object spectrum header

common/qpl/ qplot, newresults
bool	qplot		# True if only qplotting
bool	newresults	# True if new results have been gotten

# May 24 1993	Add MWCS pointers for spectrum and templates
# Jun  2 1993	Add multispec line number for spectrum and templates
# Jun 30 1993	Drop MWCS pointer for templates
# Jul  7 1993	Drop specsh and skysh from common

# Feb 11 1994	Add tsig1 and tsig2 for special type 6 report
# Mar 23 1994	Add template file id's
# May 23 1994	Add switch to save information in spectrum header
# Jun 22 1994	Add all previous velocity information for object spectrum
# Aug  3 1994	Replace emchop with separate parameters for spectrum & template
# Aug  4 1994	Add nmspec and specnums to print out history once per spectrum
# Aug 10 1994	Add IRAF task name
# Aug 17 1994	Add ZPAD to common
# Aug 19 1994	Make TWL1 and TWL2 double
# Nov 16 1994	Add SPVQUAL quality flag
# Dec  7 1994	Make TEMPFILT integer instead of boolean

# May 10 1995	Change filename lengths from sz_fname to sz_pathname
# Jun 26 1995	Add OVERLAP parameter to override limit parameters
# Aug 18 1995	Change MINSHIFT and MAXSHIFT to MINVEL and MAXVEL
# Sep 21 1995	Add QPLOT to indicate that a quality flag should be written

# Jan 24 1997	Drom PLTPHS; add PLTTFT
# Feb  5 1997	Add C0, NSMOOTH, VINIT, and CVEL from emv.com
# Mar 14 1997	Add spectrum header
# May 19 1997	Add correlation template specification for initial velocity

# Apr  6 1999	Add template header and spectrum structures

# Jul  5 2000	Add renorm and spmean for renormalization listing
# Sep 19 2000	Add wavelength/pixel rebinned spacing
# Sep 26 2000	Add template limiting wavelengths

# Jul 31 2001	Add delpix to provide shifts in pixels as well as wavelength

# Jan 31 2007	Add waverest and specref to /fitp/ common
# Aug  2 2007	Drop no-longer-used tmp* template caching buffers and nctemp

# Jul  8 2010	Add mext to /fitp/ for image extension

# May 24 2012	Add setlength to allow for long file names in output

# Aug 22 2013	Declare nmspec

# May  1 2014	Add edited to denote an edited spectrum
# May  1 2014	Add savespec to save an edited spectrum
