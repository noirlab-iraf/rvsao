# File rvsao/lib/emv.h
# March 18, 2008

# Parameters for EMSAO

# Search this fractional amount beyond the wavelength limits
define	WEXTRA	0.02  

# Fit +- this amount around emission line center wavelength
define	WGFIT	15.0

#  Maximum number of reference emission and absorption lines
define	MAXREF		500
define	MAXABS		10

#  Maximum number of emission line combinations
define	MAXCOMB		20
define	MAXCLINES	5

#  Maximum number of emission lines for velocity from single line
define	MAXSEARCH	5

define	SZ_ELINE	9

# Profile types.
define  PTYPES  "|gaussian|lorentzian|voigt|"
define  GAUSS           1       # Gaussian profile
define  LORENTZ         2       # Lorentzian profile
define  VOIGT           3       # Voigt profile

# Aug  4 1994	New file
# Feb  5 1997	Move velocity flags to rvsao.h
# Feb 27 1997	Add line profile types

# Mar 14 2008	Increase maximum number of reference lines from 50 to 500
# Mar 18 2008	Increase maximum number of combination lines from 5 to 20
