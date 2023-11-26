# File rvsao/lib/contin.h
# August 18, 1999

# Continuum fitting functions
define  CN_INTERP_MODE  "|spline3|legendre|chebyshev|spline1|"
define  CN_SPLINE3              1
define  CN_LEGENDRE             2
define  CN_CHEBYSHEV            3
define  CN_SPLINE1              4

# Continuum Subtraction Parameter Commands
define	CONT_KEYWORDS	"|c_interactive|c_sample|naverage|c_function\
			 |cn_order|s_low_reject|s_high_reject\
			 |t_low_reject|t_high_reject\
			 |s_abs_reject|s_em_reject\
			 |t_abs_reject|t_em_reject|niterate|grow|"

# Continuum computation parameters
define	CNT_INTERACTIVE		1	# Do it interactively?
define	CNT_SAMPLE		2	# Sample string to use
define	CNT_NAVERAGE		3	# Npts to average in sample
define	CNT_FUNCTION		4	# Fitting function
define	CNT_CN_ORDER		5	# Order of function
define	S_LOW_REJECT		6	# Low rejection in sigma of fit
define	S_HIGH_REJECT		7	# High rejection in sigma of fit
define	T_LOW_REJECT		8	# Low rejection in sigma of fit
define	T_HIGH_REJECT		9	# High rejection in sigma of fit
define	S_ABS_REJECT		10	# Low rejection in sigma of fit
define	S_EM_REJECT		11	# High rejection in sigma of fit
define	T_ABS_REJECT		12	# Low rejection in sigma of fit
define	T_EM_REJECT		13	# High rejection in sigma of fit
define	CNT_NITERATE		14	# Number of rejection iterations
define	CNT_GROW		15	# Rejection growing radius

# Continuum removal method
define	CONT_TYPE		"|none|subtract|divide|zerodiv|"
define	NOCONT			1
define	SUBCONT			2
define	DIVCONT			3
define	ZEROCONT		4

# Feb  5 1997	Add continuum removal options

# Aug 18 1999	Add ZERODIV continuum removal option
