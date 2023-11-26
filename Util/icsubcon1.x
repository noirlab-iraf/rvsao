# File rvsao/Util/icsubcon.x
# July 21, 2011
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After IRAF RV.CONTINUUM subroutine

include <pkg/gtools.h>
include	<error.h>
include "../lib/contin.h"
include	"../lib/rvsao.h"
include	"../lib/emv.h"


# ICSUBCON - Do the continuum normalization, with optional interaction

procedure icsubcon (n, fdata, wltemp, name, spectype, nsm, fit)

int	n			# Number of pixels to fit
real	fdata[ARB]		# Spectrum to be fit
real	wltemp[ARB]		# Wavelengths for spectrum
char	name[ARB]		# Name of spectrum
int	spectype		# 1=spectrum, 2=template
int	nsm			# Number of times to smooth spectrum
real	fit[ARB]		# Fit array (Returned)

int	i,newrej,nrej,ns,nsp
char	title_line[SZ_LINE]
char	wlab[SZ_LINE]
double	tfit
real	fitmean			# Mean continuum value for scaling
bool	chop
real	cveval()
pointer	spec0
pointer	ic			# ICFIT ptr
pointer	gt			# GTOOLS pointer
pointer	gp			# graphics pointer
pointer	cv			# pointer for curve fitting
pointer	sp
pointer	x, wts, rejpts		# buffers for curve fitting
pointer gopen(), gt_init()
errchk  gopen, gt_init

include "../lib/contin.com"
include	"../lib/rvsao.com"
include	"../lib/emv.com"
include "../lib/xwt.com"

begin
	call strcpy ("Wavelength in Angstroms", wlab, SZ_LINE)

# If the order to be fit is zero, return spectrum unchanged
	if (order < 1) {
	    do i = 1, n {
		fit[i] = fdata[i]
		}
	    return
	    }

# Check to see if spectrum is all zeroes
	nsp = 0
	do i = 1, n {
	    if (fdata[i] != 0.0)
		nsp = nsp + 1
	    }

# If spectrum is all zeroes, return it unchanged with fit all zeroes, too
	if (nsp < 1) {
	    do i = 1, n {
		fit[i] = fdata[i]
		}
	    return
	    }

# Set the ICFIT pointer structure.
	call ic_open (ic)
	call ic_pstr (ic, "sample", sample)
	call ic_pstr (ic, "function", confunc)
	call ic_puti (ic, "naverage", naverage)
	call ic_puti (ic, "order", order)
	call ic_puti (ic, "niterate", niterate)
	call ic_putr (ic, "low", lowrej[spectype])
	call ic_putr (ic, "high", hirej[spectype])
	call ic_putr (ic, "grow", grow)
	call ic_pstr (ic, "ylabel", "")

#  Allocate memory for curve fitting.
	call smark (sp)
	call salloc (x, n, TY_REAL)
	call salloc (wts, n, TY_REAL)
	call salloc (rejpts, n, TY_INT)

#  Initlialize WTS array
	call amovkr (1., Memr[wts], n)

#  Avoid fitting zeroes at ends of spectrum
	do i = 1, n {
	    if (fdata[i] == 0.) {
		Memr[wts+i-1] = 0.
#		if (debug) {
#		    call printf ("ICSUBCON: %d = 0\n")
#			call pargi (i)
#		    call flush (STDOUT)
#		    }
		}
	    else
		break
	    }
	do i = n, 1, -1 {
	    if (fdata[i] == 0.) {
		Memr[wts+i-1] = 0.
#		if (debug) {
#		    call printf ("ICSUBCON: %d = 0\n")
#			call pargi (i)
#		    call flush (STDOUT)
#		    }
		}
	    else
		break
	    }
	call amovki (NO, Memi[rejpts], n)

#  Initlialize X array
	do i = 1, n {
	    Memr[x+i-1] = real (i)
#	    if (debug) {
#		call printf ("ICSUBCON: %d=%f: %f (%f)\n")
#		    call pargi (i)
#		    call pargr (Memr[x+i-1])
#		    call pargr (fdata[i])
#		    call pargr (Memr[wts+i-1])
#		call flush (STDOUT)
#		}
	    }

	if (spectype == 1)
	    chop = schop
	else if (spectype == 2)
	    chop = tchop

#  Update icfit struct
	call ic_putr (ic, "xmin", 1.)
	call ic_putr (ic, "xmax", real(n))

#	if (debug) {
#	    call printf ("ICSUBCON: %d points being fit from spectrum %d\n")
#		call pargi (n)
#		call pargi (spectype)
#	    call flush (STDOUT)
#	    }

#  If the interactive flag is set then use icg_fit to set the fitting
#  parameters.  Only done if task is run interactively as well.

	if (interact) {
	    iferr {
		gp = gopen ("stdgraph", NEW_FILE, STDGRAPH)
		gt = gt_init ()
		call gt_sets (gt, GTTYPE, "line")
		} then
		    call error (0, "Error opening `stdgraph'.")

	    call icg_fit (ic,gp,"cursor",gt,cv,Memr[x],fdata,Memr[wts],n)
	
	# Now recover any parameters that were changed
	    call recover_icfit_pars (ic, spectype)
	    call gt_free (gt)
	    call gclose (gp)
	    }

#  Do the fit non-interactively using ic_fit.
	else {
	    call ic_fit (ic,cv,Memr[x],fdata,Memr[wts],n,YES,YES,YES,YES)
	    }

#  Save original spectrum if plotting and chopping emission lines
	spec0 = NULL
	if (pltcon && chop) {
	    call salloc (spec0,n,TY_REAL)
	    do i = 1, n {
		Memr[spec0+i-1] = fdata[i]
		}
	    }

#  Find deviant points
	if (chop && (abrej[spectype] > 0. || emrej[spectype] > 0.)) {
	    nrej = 0
	    newrej = 0
	    call ic_deviantr (cv, Memr[x],fdata,Memr[wts], Memi[rejpts],
			      n,abrej[spectype],emrej[spectype],grow,NO,
			      nrej,newrej)

#  Replace points out of bounds with continuum
	    do i = 1, n {
		if (Memi[rejpts+i-1] == YES)
		    fdata[i] = cveval (cv, Memr[x+i-1])
		}

# If plot enabled, plot the lines that were removed
	    if (pltcon && newrej > 0) {
		do i = 1, n {
		    Memr[spec0+i-1] = Memr[spec0+i-1] - fdata[i]
		    }
		call strcpy (name, title_line, SZ_LINE)
		if (emrej[spectype] > abrej[spectype])
		    call strcat (" absorption lines", title_line, SZ_LINE)
		else
		    call strcat (" emission lines", title_line, SZ_LINE)
		ns = nsm
		call plotspec (n,Memr[spec0],title_line,wltemp, wlab, ns)
		}
	    }

#  Now remove the fit continuum
	call cvvector (cv, Memr[x], fit, n)
	if (conproc == DIVCONT || conproc == ZEROCONT) {
	    tfit = 0.d0
	    do i = 1, n {
		tfit = tfit + fit[i]
		}
	    fitmean = tfit / double (n)
	    do i = 1, n {
		if (Memr[wts+i-1] > 0.) {
		    if (fit[i] <= 0.0)
			fdata[i] = 1.0
		    else
			fdata[i] = fdata[i] / fit[i]
		    }
		}
#	    call adivr (fdata, fit, fdata, n)
	    if (conproc == ZEROCONT)
		call asubkr (fdata, 1.0, fdata, n)
	    do i = 1, n {
		fdata[i] = fdata[i] * fitmean
		}
	    }
	else if (conproc == SUBCONT) {
	    do i = 1, n {
		if (Memr[wts+i-1] > 0.)
		    fdata[i] = fdata[i] - fit[i]
		}
#	    call asubr (fdata, fit, fdata, n)
	    }
	else if (conproc == CONTFIT) {
	    do i = 1, n {
		fdata[i] = fit[i]
		}
	    }

# If plot enabled, show the object spectrum with continuum removed
	if (pltcon) {
	    call strcpy (name, title_line, SZ_LINE)
	    if (conproc == NOCONT)
		call strcat (" continuum NOT removed", title_line, SZ_LINE)
	    else if (conproc == CONTFIT)
		call strcat (" continuum", title_line, SZ_LINE)
	    else
		call strcat (" continuum removed", title_line, SZ_LINE)
	    ns = nsm
	    call plotspec (n,fdata,title_line,wltemp, wlab, ns)
	    }

# If weighting spectrum is to be used, apply it now
	if (wtsh != NULL && !NOCONT && !CONTFIT) {
	    do i = 1, n {
		if (Memr[wts+i-1] > 0.) {
		    if (wtspec[i] <= 0.0)
			fdata[i] = 1.0
		    else
			fdata[i] = fdata[i] / wtspec[i]
		    }
		}
	    if (pltcon) {
		call strcat (" weighted", title_line, SZ_LINE)
		ns = nsm
		call plotspec (n,fdata,title_line,wltemp, wlab, ns)
		}
	    }

	call ic_closer (ic)
	call cvfree (cv)
	call sfree (sp)
	return
end


# ICSUBCONY - Do the continuum normalization, with optional interaction
#             Set range of counts/flux in graph

procedure icsubcony (n, fdata, wltemp, name, spectype, nsm, fit, ymin, ymax)

int	n			# Number of pixels to fit
real	fdata[ARB]		# Spectrum to be fit
real	wltemp[ARB]		# Wavelengths for spectrum
char	name[ARB]		# Name of spectrum
int	spectype		# 1=spectrum, 2=template
int	nsm			# Number of times to smooth spectrum
real	fit[ARB]		# Fit array (Returned)
real	ymin, ymax		# Minimum, maximum flux/count values to plot

int	i,newrej,nrej,ns,nsp
char	title_line[SZ_LINE]
char	wlab[SZ_LINE]
double	tfit
real	fitmean			# Mean continuum value for scaling
bool	chop
real	cveval()
pointer	spec0
pointer	ic			# ICFIT ptr
pointer	gt			# GTOOLS pointer
pointer	gp			# graphics pointer
pointer	cv			# pointer for curve fitting
pointer	sp
pointer	x, wts, rejpts		# buffers for curve fitting
pointer gopen(), gt_init()
errchk  gopen, gt_init

include "contin.com"
include	"rvsao.com"
include	"emv.com"

begin
	call strcpy ("Wavelength in Angstroms", wlab, SZ_LINE)

# If the order to be fit is zero, return spectrum unchanged
	if (order < 1) {
	    do i = 1, n {
		fit[i] = fdata[i]
		}
	    return
	    }

# Check to see if spectrum is all zeroes
	nsp = 0
	do i = 1, n {
	    if (fdata[i] != 0.0)
		nsp = nsp + 1
	    }

# If spectrum is all zeroes, return it unchanged with fit all zeroes, too
	if (nsp < 1) {
	    do i = 1, n {
		fit[i] = fdata[i]
		}
	    return
	    }

# Set the ICFIT pointer structure.
	call ic_open (ic)
	call ic_pstr (ic, "sample", sample)
	call ic_pstr (ic, "function", confunc)
	call ic_puti (ic, "naverage", naverage)
	call ic_puti (ic, "order", order)
	call ic_puti (ic, "niterate", niterate)
	call ic_putr (ic, "low", lowrej[spectype])
	call ic_putr (ic, "high", hirej[spectype])
	call ic_putr (ic, "grow", grow)
	call ic_pstr (ic, "ylabel", "")

#  Allocate memory for curve fitting.
	call smark (sp)
	call salloc (x, n, TY_REAL)
	call salloc (wts, n, TY_REAL)
	call salloc (rejpts, n, TY_INT)

#  Initlialize WTS array
	call amovkr (1., Memr[wts], n)

#  Avoid fitting zeroes at ends of spectrum
	do i = 1, n {
	    if (fdata[i] == 0.) {
		Memr[wts+i-1] = 0.
#		if (debug) {
#		    call printf ("ICSUBCON: %d = 0\n")
#			call pargi (i)
#		    call flush (STDOUT)
#		    }
		}
	    else
		break
	    }
	do i = n, 1, -1 {
	    if (fdata[i] == 0.) {
		Memr[wts+i-1] = 0.
#		if (debug) {
#		    call printf ("ICSUBCON: %d = 0\n")
#			call pargi (i)
#		    call flush (STDOUT)
#		    }
		}
	    else
		break
	    }
	call amovki (NO, Memi[rejpts], n)

#  Initlialize X array
	do i = 1, n {
	    Memr[x+i-1] = real (i)
#	    if (debug) {
#		call printf ("ICSUBCON: %d=%f: %f (%f)\n")
#		    call pargi (i)
#		    call pargr (Memr[x+i-1])
#		    call pargr (fdata[i])
#		    call pargr (Memr[wts+i-1])
#		call flush (STDOUT)
#		}
	    }

	if (spectype == 1)
	    chop = schop
	else if (spectype == 2)
	    chop = tchop

#  Update icfit struct
	call ic_putr (ic, "xmin", 1.)
	call ic_putr (ic, "xmax", real(n))

#	if (debug) {
#	    call printf ("ICSUBCON: %d points being fit from spectrum %d\n")
#		call pargi (n)
#		call pargi (spectype)
#	    call flush (STDOUT)
#	    }

#  If the interactive flag is set then use icg_fit to set the fitting
#  parameters.  Only done if task is run interactively as well.

	if (interact) {
	    iferr {
		gp = gopen ("stdgraph", NEW_FILE, STDGRAPH)
		gt = gt_init ()
		call gt_sets (gt, GTTYPE, "line")
		} then
		    call error (0, "Error opening `stdgraph'.")

	    call icg_fit (ic,gp,"cursor",gt,cv,Memr[x],fdata,Memr[wts],n)
	
	# Now recover any parameters that were changed
	    call recover_icfit_pars (ic, spectype)
	    call gt_free (gt)
	    call gclose (gp)
	    }

#  Do the fit non-interactively using ic_fit.
	else {
	    call ic_fit (ic,cv,Memr[x],fdata,Memr[wts],n,YES,YES,YES,YES)
	    }

#  Save original spectrum if plotting and chopping emission lines
	spec0 = NULL
	if (pltcon && chop) {
	    call salloc (spec0,n,TY_REAL)
	    do i = 1, n {
		Memr[spec0+i-1] = fdata[i]
		}
	    }

#  Find deviant points
	if (chop && (abrej[spectype] > 0. || emrej[spectype] > 0.)) {
	    nrej = 0
	    newrej = 0
	    call ic_deviantr (cv, Memr[x],fdata,Memr[wts], Memi[rejpts],
			      n,abrej[spectype],emrej[spectype],grow,NO,
			      nrej,newrej)

#  Replace points out of bounds with continuum
	    do i = 1, n {
		if (Memi[rejpts+i-1] == YES)
		    fdata[i] = cveval (cv, Memr[x+i-1])
		}

# If plot enabled, plot the lines that were removed
	    if (pltcon && newrej > 0) {
		do i = 1, n {
		    Memr[spec0+i-1] = Memr[spec0+i-1] - fdata[i]
		    }
		call strcpy (name, title_line, SZ_LINE)
		if (emrej[spectype] > abrej[spectype])
		    call strcat (" absorption lines", title_line, SZ_LINE)
		else
		    call strcat (" emission lines", title_line, SZ_LINE)
		ns = nsm
		call plotsum (n,Memr[spec0],title_line,wltemp, wlab, ns, ymin, ymax)
		}
	    }

#  Now remove the fit continuum
	call cvvector (cv, Memr[x], fit, n)
	if (conproc == DIVCONT || conproc == ZEROCONT) {
	    tfit = 0.d0
	    do i = 1, n {
		tfit = tfit + fit[i]
		}
	    fitmean = tfit / double (n)
	    do i = 1, n {
		if (Memr[wts+i-1] > 0.) {
		    if (fit[i] <= 0.0)
			fdata[i] = 1.0
		    else
			fdata[i] = fdata[i] / fit[i]
		    }
		}
#	    call adivr (fdata, fit, fdata, n)
	    if (conproc == ZEROCONT)
		call asubkr (fdata, 1.0, fdata, n)
	    do i = 1, n {
		fdata[i] = fdata[i] * fitmean
		}
	    }
	else if (conproc == SUBCONT) {
	    do i = 1, n {
		if (Memr[wts+i-1] > 0.)
		    fdata[i] = fdata[i] - fit[i]
		}
#	    call asubr (fdata, fit, fdata, n)
	    }
	else if (conproc == CONTFIT) {
	    do i = 1, n {
		fdata[i] = fit[i]
		}
	    }

# If plot enabled, show the object spectrum with continuum removed
	if (pltcon) {
	    call strcpy (name, title_line, SZ_LINE)
	    if (conproc == NOCONT)
		call strcat (" continuum NOT removed", title_line, SZ_LINE)
	    else if (conproc == CONTFIT)
		call strcat (" continuum", title_line, SZ_LINE)
	    else
		call strcat (" continuum removed", title_line, SZ_LINE)
	    ns = nsm
	    call plotsum (n,fdata,title_line,wltemp, wlab, ns, ymin, ymax)
	    }

# If weighting spectrum is to be used, apply it now
	if (wtsh != NULL && !NOCONT && !CONTFIT) {
	    do i = 1, n {
		if (Memr[wts+i-1] > 0.) {
		    if (wtspec[i] <= 0.0)
			fdata[i] = 1.0
		    else
			fdata[i] = fdata[i] / wtspec[i]
		    }
		}
	    if (pltcon) {
		call strcat (" weighted", title_line, SZ_LINE)
		ns = nsm
		call plotsum (n,fdata,title_line,wltemp, wlab, ns, ymin, ymax)
		}
	    }

	call ic_closer (ic)
	call cvfree (cv)
	call sfree (sp)
	return
end


# RECOVER_ICFIT_PARS - Since the ICFIT parameters may have been changed in
# an interactive operation, we need to get the new values from the ICFIT
# structure.

procedure recover_icfit_pars (ic,intype)

pointer	ic		# ICFIT pointer
int	intype		# 1=spectrum, 2=template

int	strdic(), ic_geti()
real	ic_getr()
include "contin.com"

begin
	naverage = ic_geti (ic, "naverage")
	order = ic_geti (ic, "order")
	niterate = ic_geti (ic, "niterate")
	lowrej[intype] = ic_getr (ic, "low")
	hirej[intype] = ic_getr (ic, "high")
	grow = ic_getr (ic, "grow")

	call ic_gstr (ic, "sample", sample, SZ_LINE)
	call ic_gstr (ic, "function", confunc, SZ_LINE)
	function = strdic (confunc, confunc, SZ_LINE, CN_INTERP_MODE)
	if (function == 0)
	    call error (0, "Unknown fitting function type")
	return
end
# Nov 18 1991	Use malloc instead of salloc
# Nov 20 1991	Free buffers before allocating them
# Dec 16 1991	Plot spectrum and emission lines removed
#		Fix bug by filling in rejected points
#		Chop emission lines only if chop is set
# Dec 17 1991	Don't plot emission lines if none have been removed

# Oct 22 1992	Always free spec0 before returning

# Apr 13 1994	Cleanup code after ftnchek
# Apr 26 1994	Clean up code
# May  3 1994	Add smoothed vector argument to PLOTSPEC call
# Jun 15 1994	Add smoothing argument to PLOTSPEC call
# Aug  3 1994	Use separate template and spectrum emission line chopping flags
# Aug  3 1994	Change common and header from fquot to rvsao

# Mar 13 1995	Allow line chopping in absorption as well as emission
# Mar 15 1995	Add argument for name of spectrum
# May 15 1995	Add template-driven option to divide, not subtract, continuum
# Oct 16 1995	Try salloc instead of malloc

# Dec 13 1996	Remove zeroes and negative numbers from continuum if dividing
# Dec 16 1996	Subtract one after dividing continuum

# May  6 1997	Add smoothing parameter as argument
# Dec 22 1997	Add units to wavelength graph label

# Aug 16 1999	Comment out subtraction of 1 when dividing out continuum
# Aug 18 1999	Add ZEROCONT for continuum divided - 1
# Aug 18 1999	Use int conproc instead of bool divcon

# Jan 25 2000	Add option to not remove continuum

# Apr 26 2001	Drop zeroed pixels at spectrum ends from continuum fit

# May  7 2004	Add option to return continuum fit in place of spectrum

# Apr 25 2006	Return fit and spectrum as all zeroes if spectrum is all zeroes
# Jul  7 2006	If dividing out continuum and it is zero or negative, set to 1
# Sep 18 2006	Add ICSUBCONY to plot continuum with flux/count limits

# Jul 16 2009	Renormalize by mean continuum if dividing out spectrum
# Jul 16 2009	Return spectrum unchanged if order is zero

# Jul  7 2011	Divide by optional weighting function, such as normalized blaze
