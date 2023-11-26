# File rvsao/Util/contpars.x
# January 30, 2007
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# CONTPARS  - Support routines for the 'contpars' named external pset.  

#	This file include routines for filling the /contin/ common as well
# as command handling.  Command handling is limited to changing the parameter
# values or resetting them to the default values.  Routines included here are
# as follows:
#
#		  cont_get_pars ()
#		 cont_parupdate ()
#		   cont_unlearn ()
#		      cont_show ()
#		     cont_colon (cmdstr)
#		cmd_interactive ()
#		     cmd_sample ()
#		   cmd_naverage ()
#		   cmd_function ()
#		   cmd_cn_order ()
#	       cmd_s_low_reject ()
#	      cmd_s_high_reject ()
#	       cmd_t_low_reject ()
#	      cmd_t_high_reject ()
#		   cmd_niterate ()
#		       cmd_grow ()
#
#	The 'cmd_' prefix indicates that the routine is called from a colon 
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


# CONT_GET_PARS - Get the continuum fitting parameters from the pset.

procedure cont_get_pars ()

pointer	pp, clopset()
int	strdic(), clgpseti()
real	clgpsetr()
bool	clgpsetb(), streq()
include "contin.com"

begin
	# Get continuum parameters.
	iferr (pp = clopset("contpars"))
	    call error (0, "XCSAO: Error opening `contpars' pset")

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
end


# CONT_PARUPDATE - Update the pset with the current values of the struct.

procedure cont_parupdate ()

pointer	sp, b1
pointer	pp, clopset()
errchk  clopset
include "contin.com"

begin
	# Update contin params
	iferr (pp = clopset ("contpars")) {
	    call printf ("XCSAO: Error opening `contpars' pset.")
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


# CONT_UNLEARN - Unlearn the pset and replace with the default values.

procedure cont_unlearn ()

include "contin.com"

begin
	order = DEF_ORDER
	lowrej[1] = DEF_S_LOW_REJECT
	hirej[1] = DEF_S_HIGH_REJECT
	lowrej[2] = DEF_T_LOW_REJECT
	hirej[2] = DEF_T_HIGH_REJECT
	niterate = DEF_NITERATE
	grow = DEF_GROW
	naverage = DEF_NAVERAGE
	interact = DEF_INTERACTIVE

	call strcpy (DEF_SAMPLE, sample, SZ_FNAME)
	call strcpy (DEF_CONFUNC, confunc,SZ_FNAME)
	function = DEF_FUNCTION
end


# CONT_SHOW - Show the current contin parameters

procedure cont_show (fd)

pointer	fd			#I output file descriptor

pointer	sp, str, str1
include "contin.com"

begin
	if (fd == NULL)
	    return

	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)
	call salloc (str1,SZ_LINE, TY_CHAR)

	call fprintf (fd, "%21tProcesspars PSET Values\n")
	call fprintf (fd, "%21t-----------------------\n\n")

	# Print the contpars info
	call fprintf (fd, "CONTINUUM parameters:\n")

	call fprintf (fd, "c_interactive%15t= %b\n")
	    call pargb (interact)
	call fprintf (fd, "c_sample%15t= '%.10s'\n")
	    call pargstr (sample)
	call fprintf (fd, "naverage%15t= %d\n")
	    call pargi (naverage)
	call fprintf (fd, "c_function%15t= '%.10s'\n")
	    call pargstr (confunc)
	call fprintf (fd, "order%15t= %d\n")
	    call pargi (order)
	call fprintf (fd, "s_low_reject%15t= %g\n")
	    call pargr (lowrej[1])
	call fprintf (fd, "s_high_reject%15t= %g\n")
	    call pargr (hirej[1])
	call fprintf (fd, "t_low_reject%15t= %g\n")
	    call pargr (lowrej[2])
	call fprintf (fd, "t_high_reject%15t= %g\n")
	    call pargr (hirej[2])
	call fprintf (fd, "niterate%15t= %d \n")
	    call pargi (niterate)
	call fprintf (fd, "grow%15t= %g\n")
	    call pargr (grow)

	call fprintf (fd, "\n\n")
	call sfree (sp)
end


# CONT_COLON -- Process the contpars task colon commands.

procedure cont_colon (cmdstr)

char	cmdstr[SZ_LINE]			#I command string

pointer	sp, cmd
int	strdic()
include "contin.com"

begin
	call smark (sp)
	call salloc (cmd, SZ_LINE, TY_CHAR)

	call sscan (cmdstr)
	call gargwrd (Memc[cmd], SZ_LINE)

	# Unpack the keyword from the string and look it up in the
	# dictionary.  Switch on command and call the appropriate routines.

	switch (strdic(Memc[cmd], Memc[cmd], SZ_FNAME, CONT_KEYWORDS)) {

	    case CNT_INTERACTIVE:
		call cmd_interactive ()
	    case CNT_SAMPLE:
		call cmd_sample ()
	    case CNT_NAVERAGE:
		call cmd_naverage ()
	    case CNT_FUNCTION:
		call cmd_cnfunc ()
	    case CNT_CN_ORDER:
		call cmd_cn_order ()
	    case S_LOW_REJECT:
		call cmd_s_low_reject ()
	    case S_HIGH_REJECT:
		call cmd_s_high_reject ()
	    case T_LOW_REJECT:
		call cmd_t_low_reject ()
	    case T_HIGH_REJECT:
		call cmd_t_high_reject ()
	    case CNT_NITERATE:
		call cmd_niterate ()
	    case CNT_GROW:
		call cmd_grow ()
	    default:
	    }

	call sfree (sp)
end


# CMD_INTERACTIVE - Set/Show the interactive continuum subtraction flag.

procedure cmd_interactive ()

int	nscan()
bool	bval
include "contin.com"

begin
	call gargb (bval)
	if (nscan() == 2) {
	    interact = bval
	    }
	else {
	    call printf ("contpars.c_interactive = %b")
		call pargb (interact)
	    }
end


# CMD_SAMPLE - Set/Show the sample regions for continuum fitting.

procedure cmd_sample ()

pointer	sp, buf
bool	streq()
include "contin.com"

begin
	call smark (sp)
	call salloc (buf, SZ_LINE, TY_CHAR)

	call gargstr (Memc[buf], SZ_FNAME)
	if (Memc[buf] != EOS) {
	    if (streq(Memc[buf],"") || streq(Memc[buf]," "))
	        call error (0, "contpars.c_sample specified as empty string.")
	    call strcpy (Memc[buf+1], sample, SZ_LINE)
	    }
	else {
	    call printf ("contpars.c_sample = '%s'")
	        call pargstr (sample)
	    }

	call sfree (sp)
end


# CMD_NAVERAGE - Set/Show the number of points to average in the fit.

procedure cmd_naverage ()

int	ival, nscan()
include "contin.com"

begin
	call gargi (ival)
	if (nscan() == 2) {
	    naverage = ival
	    }
	else {
	    call printf ("contpars.naverage = %d")
		call pargi (naverage)
	    }
end


# CMD_CNFUNC - Set/Show the fitting function used.

procedure cmd_cnfunc ()

pointer	sp, buf, bp
int	strdic()
include "contin.com"

begin
	call smark (sp)
	call salloc (buf, SZ_LINE, TY_CHAR)
	call salloc (bp, SZ_LINE, TY_CHAR)

	call gargstr (confunc, SZ_FNAME)
	if (Memc[buf] != EOS) {
	    function = strdic (confunc, confunc, SZ_LINE, CN_INTERP_MODE)
	    }
	else {
	    call printf ("contpars.c_function = '%s'")
		call pargstr (confunc)
	    }

	call sfree (sp)
end


# CMD_CN_ORDER - Set/Show the order of the function fit.

procedure cmd_cn_order ()

int	ival, nscan()
include "contin.com"

begin
	call gargi (ival)
	if (nscan() == 2) {
	    order = ival
	    }
	else {
	    call printf ("contpars.order = %d")
		call pargi (order)
	    }
end


# SET_CN_ORDER - Set the order of the function fit.

procedure set_cn_order (ival)

int	ival

include "contin.com"

begin
	order = ival
end


# CMD_S_LOW_REJECT - Set/Show the lower sigma rejection limit.

procedure cmd_s_low_reject ()

real	rval
int	nscan()
include "contin.com"

begin
	call gargr (rval)
	if (nscan() == 2) {
	    lowrej[1] = rval
	    }
	else {
	    call printf ("contpars.s_low_reject = %g")
		call pargr (lowrej[1])
	    }
end


# CMD_S_HIGH_REJECT - Set/Show the upper sigma rejection limit.

procedure cmd_s_high_reject ()

real	rval
int	nscan()
include "contin.com"

begin
	call gargr (rval)
	if (nscan() == 2) {
	    hirej[1] = rval
	    }
	else {
	    call printf ("contpars.s_high_reject = %g")
		call pargr (hirej[1])
	    }
end


# CMD_T_LOW_REJECT - Set/Show the lower sigma rejection limit.

procedure cmd_t_low_reject ()

real	rval
int	nscan()
include "contin.com"

begin
	call gargr (rval)
	if (nscan() == 2) {
	    lowrej[2] = rval
	    }
	else {
	    call printf ("contpars.t_low_reject = %g")
		call pargr (lowrej[2])
	    }
end


# CMD_T_HIGH_REJECT - Set/Show the upper sigma rejection limit.

procedure cmd_t_high_reject ()

real	rval
int	nscan()
include "contin.com"

begin
	call gargr (rval)
	if (nscan() == 2) {
	    hirej[2] = rval
	    }
	else {
	    call printf ("contpars.t_high_reject = %g")
		call pargr (hirej[2])
	    }
end


# CMD_NITERATE - Set/Show the number of iterations in the fit.

procedure cmd_niterate ()

int	ival, nscan()
include "contin.com"

begin
	call gargi (ival)
	if (nscan() == 2) {
	    niterate = ival
	    }
	else {
	    call printf ("contpars.niterate = %d")
		call pargi (niterate)
	    }
end


# CMD_GROW - Set/Show the rejection growing radius.

procedure cmd_grow ()

real	rval
int	nscan()
include "contin.com"

begin
	call gargr (rval)
	if (nscan() == 2) {
	    grow = rval
	    }
	else {
	    call printf ("contpars.grow = %g")
		call pargr (grow)
	    }
end
# Apr 20 1994	Eliminate 4th argument to CLPPSET
# Apr 21 1994	Use SSCAN as function instead of subroutine

# Mar 29 1995	Remove emission and absorption rejection for output spectra

# Feb  5 1997	Change help statements at start to comments

# Jan 30 2007	Add set_cn_order() to set order in common
