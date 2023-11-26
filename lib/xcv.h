# File lib/xcv.h
# April 7, 1997

# Parameters for XCSAO

#  Filter function flags
define	XC_FTYPES	"|square|ramp|welch|hanning|cos-bell|"
define  SQUARE		1	# Filter type
define	RAMP		2	# Filter type
define	WELCH		3	# Filter type
define	HANNING		4	# Filter type
define	COSBELL		5	# Filter type
define	NOFILT		6	# No filter

#  Emission line chopping flags
define	XC_CTYPES	"|no|yes|tempfile|specfile|"
define  NOCHOP		1	# no correction
define  CHOP		2	# Always remove emission lines
define	TEMPFILE	3	# Remove emission lines if template CHOPEM=T
define	SPECFILE	4	# Remove emission lines if spectrum CHOPEM=T

#  Zero padding flags
define	XC_ZTYPES	"|no|yes|tempfile|"
define  NOZPAD		1	# No zero padding
define  ZPAD		2	# Zero-pad object and template transforms
#define	TEMPFILE	3	# Zero-pad lines if template ZEROPAD=T

#  XCSAO initial velocity flags
define	XC_VTYPES	"|zero|guess|correlation|emission|combination|"
define	IVZERO		1	# Set to zero
define	IVGUESS		2	# Use guess from parameter file
define  IVXC		3	# Use previous cross-correlation velocity
define	IVEM		4	# Use previous emission line fit velocity
define	IVCOMB		5	# Use previous combined velocity and

define	MAXCOR		8192	# Maximum number of points in correlation

# Apr  6 1995	New common
#
# Apr  7 1997	Add maximum length of correlation vector
