# File rvsao/Makespec/contsum.x
# March 20, 1998
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# CONTSUM  - Support routines for the 'contsum' named external pset.  

# 	This file include routines for filling the /contin/ common as well
# as command handling.  Command handling is limited to changing the parameter
# values or resetting them to the default values.  Routines included here are
# as follows:
# 
# 		  csum_get_pars ()
# 		 csum_parupdate ()
# 
# 	The 'cmd_' prefix indicates that the routine is called from a colon 
# command to either print the current value or set the new value for that
# field.  Other routines should be self-explanatory

include "contin.h"

# Default values for the XCONPARS pset
define	DEF_INTERACTIVE		FALSE		# Fit continuum interactively?
define	DEF_TYPE		DIFF		# Type of output(fit|diff|ratio)
define	DEF_SAMPLE		"*"		# Sample of points to use in fit
define	DEF_NAVERAGE		1		# Npts in sample averaging
define	DEF_FUNCTION		CN_SPLINE3	# Fitting function
define	DEF_CONFUNC		"spline3"	# Fitting function
define	DEF_ORDER		1		# Order of fitting function
define	DEF_S_LOW_REJECT	2.		# Low rejection in sigma--fit
define	DEF_S_HIGH_REJECT	2.		# High rejection in sigma--fit
define	DEF_T_LOW_REJECT	2.		# Low rejection in sigma--fit
define	DEF_T_HIGH_REJECT	2.		# High rejection in sigma--fit
define	DEF_NITERATE		10		# Number of rejection iterations
define	DEF_GROW		1.		# Rejection growing radius


# CSUM_GET_PARS - Get the continuum fitting parameters from the pset.

procedure csum_get_pars (ctype)

char ctype[ARB]		# Type of continuum removal
			# (|subtract|divide|zerodiv|no|s2|s3|s4|s5|d2|d3|d4|d5|)

pointer	pp, clopset()
int	strdic(), clgpseti()
real	clgpsetr()
bool	clgpsetb(), streq()
include "contin.com"

begin
	# Get continuum parameters.
	iferr (pp = clopset("contsum"))
	    call error (0, "CONTSUM: Error opening `contsum' pset")

	call clgpset (pp, "c_function", confunc, SZ_LINE)
	if (streq(confunc,"") || streq(confunc," "))
	    call error (0,"Continpars.function specified as empty string.")
	function = strdic (confunc, confunc, SZ_LINE, CN_INTERP_MODE)
	if (function == 0) 
	    call error (0, "Unknown fitting function type")

	call clgpset (pp, "c_sample", sample, SZ_LINE)
	if (streq(sample,"") || streq(sample," "))
	    call strcpy ("*", sample, SZ_FNAME)

	order = clgpseti (pp, "order")
	niterate = clgpseti (pp, "niterate")
	naverage = clgpseti (pp, "naverage")
	grow = clgpsetr (pp, "grow")
	lowrej[1] = clgpsetr (pp, "s_low_reject")
	hirej[1] = clgpsetr (pp, "s_high_reject")
	lowrej[2] = clgpsetr (pp, "t_low_reject")
	hirej[2] = clgpsetr (pp, "t_high_reject")
	interact = clgpsetb(pp, "c_interactive")

	call clcpset (pp)				# Close pset

	# Find out if the spectrum is being split into pieces for fit
	if (ctype[2] == '2')
	    ncfit = 2
	else if (ctype[2] == '3')
	    ncfit = 3
	else if (ctype[2] == '4')
	    ncfit = 4
	else if (ctype[2] == '5')
	    ncfit = 5
	else
	    ncfit = 1

	# Get continuum parameters for second part of spectrum
	if (ncfit > 1) {
	    iferr (pp = clopset("contsum"))
		call error (0, "CONTSUM: Error opening `contsum2' pset")

	    call clgpset (pp, "c_function", confunc, SZ_LINE)
	    if (streq(confunc,"") || streq(confunc," "))
		call error (0,"Continpars.function specified as empty string.")
	    function[2] = strdic (confunc, confunc, SZ_LINE, CN_INTERP_MODE)
	    if (function[2] == 0) 
		call error (0, "Unknown fitting function type")
	    call clgpset (pp, "c_sample", sample, SZ_LINE)
	    if (streq(sample,"") || streq(sample," "))
		call strcpy ("*", sample, SZ_FNAME)

	    order[2] = clgpseti (pp, "order")
	    niterate[2] = clgpseti (pp, "niterate")
	    naverage[2] = clgpseti (pp, "naverage")
	    grow[2] = clgpsetr (pp, "grow")
	    lowrej[2] = clgpsetr (pp, "t_low_reject")
	    hirej[2] = clgpsetr (pp, "t_high_reject")
	    interact[2] = clgpsetb(pp, "c_interactive")

	    }
end


# CSUM_PARUPDATE - Update the pset with the current values of the struct.

procedure csum_parupdate ()

pointer	sp, b1
pointer	pp, clopset()
errchk  clopset
include "contin.com"

begin
	# Update contin params
	iferr (pp = clopset ("contsum")) {
	    call printf ("CONTSUM: Error opening `contsum' pset.")
	    return
	}

	call smark (sp)
	call salloc (b1, SZ_LINE, TY_CHAR)

	call clppseti (pp, "order", order)
	call clppseti (pp, "naverage", naverage)
	call clppseti (pp, "niterate", niterate)

	call clppsetr (pp, "s_low_reject", lowrej[1])
	call clppsetr (pp, "s_high_reject", hirej[1])
	call clppsetr (pp, "t_low_reject", lowrej[2])
	call clppsetr (pp, "t_high_reject", hirej[2])
	call clppsetr (pp, "grow", grow)

	call clppsetb (pp, "c_interactive", interact)

	call clppset (pp, "c_function", confunc)

	call clppset (pp, "c_sample", sample)

	call clcpset (pp)
	call sfree (sp)
end

# Feb  3 1997	New subroutine in rvsao/Sumtemp

# Mar 20 1998	Fix error messages
