# File lib/results.com
# August 18, 1999

# Cross-correlation answers from XCSAO
 
common /answr/ zvel,cz,czerr,czr,xcrmax, tcent,thght,twdth,tvw,tarms,tsrms,
		taa, tvshift,tsig1,tsig2,tpcent,tphght,tpwdth,
		tnpfit,itmax,npmax,tschop,tachop, tscont,tconproc,toverlap
 
double	zvel[MAXTEMPS]		# True barycentric velocities per template
double	cz[MAXTEMPS]		# Cz for template correlations
double	czerr[MAXTEMPS]		# Cz error for template correlations
double	czr[MAXTEMPS]		# R-values for template correlations
double	xcrmax			# Cz for template with maximum R value
double	tcent[MAXTEMPS]		# Centers for template correlations
double	thght[MAXTEMPS]		# Fit heights for template correlations
double	twdth[MAXTEMPS]		# Peak width(pixel)  for template correlations
double	tvw[MAXTEMPS]		# Peak width(velocity) for template correlations
double	tarms[MAXTEMPS]		# Antisymmetric rms for template correlations
double	tsrms[MAXTEMPS]		# Symmetric rms for template correlations
double	taa[MAXTEMPS]		# Pixels per log wavelength for each template
double	tvshift[MAXTEMPS]	# Velocity shift from successive guesses
double	tsig1[MAXTEMPS]		# Sigma of spectrum transforms
double	tsig2[MAXTEMPS]		# Sigma of template transforms
double	tpcent[MAXTEMPS]	# Measured centers for template correlations
double	tphght[MAXTEMPS]	# Measured heights for template correlations
double	tpwdth[MAXTEMPS]	# Measured widths for template correlations
int	tnpfit[MAXTEMPS]	# Number of points fit for each template
int	itmax			# Template with maximum R value
int	npmax			# Number of points fit for max R template
bool	tschop[MAXTEMPS]	# Emission line chopping per template
bool	tachop[MAXTEMPS]	# Absorption line chopping per template
bool	tscont[MAXTEMPS]	# Subtract continuum from template
int	tconproc[MAXTEMPS]	# Method for continuum removal
bool	toverlap[MAXTEMPS]	# Use overlapping spectrum region, not one set

common /answrq/ sig,sigerr,gam,gamerr,perdeg

double	sig[MAXTEMPS]
double	sigerr[MAXTEMPS]
double	gam[MAXTEMPS]
double	gamerr[MAXTEMPS]
real	perdeg[MAXTEMPS]

# Parameters for error testing

common/errtest/ tczerr1, tczerr2, czerr1, czerr2

double	tczerr1			# Error from peak counting
double	tczerr2			# Error from peak width
double	czerr1[MAXTEMPS]	# Error from peak counting
double	czerr2[MAXTEMPS]	# Error from peak width

common/templ/ itr
int	itr[MAXTEMPS]	# Template indices sorted by R-value

# Mar 13 1995	Add TACHOP as a flag for absorption line chopping
# Jun 26 1995	Add TOVERLAP to bypass limiting wavelength parameters
# Jul  3 1995	Add TPHGHT to save measured height

# Apr  7 1997	Add ITR for template indices sorted by R-value

# Aug 18 1999	Change bool tdivcon to int tcontproc
