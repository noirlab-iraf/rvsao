# Common for equivalent width computation

common/eqw/ magzero,z1,exptime,dnrows,bindim,rmode,norm,net,mag,torest,byexp,
	    bypix,fitcont, readfilt, contname, fspec

double	magzero		# Magnitude zeropoint for magnitude output
double	z1		# 1+z redshift (1.0 = use rest wavelength)
double	exptime		# Number of seconds to be normalized out
double	dnrows		# Number of rows to be normalized out
int	bindim		# Number of rows in unbinned (bin by 1) spectrum
int	rmode		# Report mode 1=old 2=line names in heading
bool	norm		# Normalize bands by response?
bool	net		# Output net flux in main region, not equivalent width
bool	mag		# Output magnitudes instead of fluxes?
bool	torest		# If true, shift from observed velocity to rest
bool	byexp		# If true, divide by exposure time in seconds
bool	bypix		# If true, divide by number of rows added (sky)
bool	fitcont		# If true, use fit continuum from contpars.par
bool	readfilt	# If true, filter name follows wavelength limits
bool	contname	# If true, continuum name precedes wavelength limits
char	fspec[16]	# Format for spectrum name in output set in eq_header()
