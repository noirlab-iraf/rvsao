# File rvsao/Util/icsubcon.x
# July 7, 2006
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After IRAF RV.CONTINUUM subroutine

include <pkg/gtools.h>
include	<error.h>
include "contin.h"
include	"rvsao.h"
include	"emv.h"


# ICSUBCONY - Do the continuum normalization, with optional interaction
#             Set range of counts/flux in graph

procedure icsubcon (n, fdata, wltemp, name, spectype, nsm, fit, ymin, ymax)

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

	call ic_closer (ic)
	call cvfree (cv)
	call sfree (sp)
	return
end
