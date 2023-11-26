# Common for IRAF ICFIT continuum normalization
# July 17, 2009

common/contin/ sample, confunc, function, naverage, order, niterate,
		lowrej, hirej, abrej, emrej, grow, interact, conproc,
		congt, congp

char	sample[SZ_FNAME]
char	confunc[SZ_FNAME,2]	# Continuum fitting function type
int	function[2]		# Continuum fitting function code
int	naverage		# number of points to average
int	order[2]		# order of fit
int	niterate		# number of iterations for fit
real	lowrej[2]		# fit rejection limit in sigma below continuum
real	hirej[2]		# fit rejection limit in sigma above continuum
real	abrej[2]		# line rejection limit in sigma below continuum
real	emrej[2]		# line rejection limit in sigma above continuum
real	grow			# radius around rejected points to replace
bool	interact		# Use interactive fitting (true or false)
int	conproc[2]		# Continuum removal method (none,sub,div,zero)
pointer	congt			# pointer to gtools structure
pointer	congp			# pointer to graphics structure

# May 15 1995	Add continuum division option

# Aug 18 1999	Replace divcon with continuum removal type conproc

# Jul 17 2009	Add second value to function and order
