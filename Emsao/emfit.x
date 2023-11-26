# File Emsao/emfit.x
# October 7, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1997-2008 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Compute radial velocity from emission line shift

include <smw.h>
include "rvsao.h"
include "emv.h"

procedure emfit (spectrum,wlspec,sky,specfile,mspec,specim,pix1,pix2,rmode)

real	spectrum[ARB]	# Unsmoothed object spectrum (wavelength-binned)
real	wlspec[ARB]	# Wavelengths for object spectrum
real	sky[ARB]	# Sky spectrum (wavelength-binned)
char	specfile[ARB]	# Name of spectrum file
int	mspec		# Spectrum number if multispec file, else 0
pointer	specim		# Header structure for spectrum
int	pix1		# Blue limit of spectrum in pixels
int	pix2		# Red limit of spectrum in pixels
int	rmode		# Report format (1=normal,2=one-line)

double	tcz		# Cz for object (returned)
double	tczrms		# 1+Z rms error for object (returned)
double	tczchi		# 1+Z  chi^2 for object (returned)
double	wblue		# Blue limit of spectrum in wavelength
double	wred		# Red limit of spectrum in wavelength
int	i
double	drms, meanwidth
double	speerr0
double	czvel		# Initial velocity guess in km/sec
bool	pixfill
int	nz
int	npix
int	itemp
double	sptot, spmin, spmax, specd
int	nsp

pointer	sp
pointer	smplot		# Smoothed object spectrum (wavelength-binned)
pointer	smspec		# Spectrum smoothed for searching (continuum removed)
pointer	smcont		# Continuum smoothed for searching and plotting
pointer	sfspec		# Spectrum smoothed for fitting
pointer	sfcont		# Continuum smoothed for fitting

pointer	smused		# Smoothed object spectrum within wavelength limits
pointer	smcusd		# Smoothed continuum spectrum within wavelength limits
pointer	sfused		# Object spectrum for fitting within wavelength limits
pointer	sfcusd		# Continuum smoothed for fitting within limits

bool	linefit		# True to fit line profiles (centers)
bool	velfit		# True to fit velocity for multiple lines
int	spectype	# =1 for object spectrum continuum substraction
int	iord
int	mspec0
double	dindef
char	xlab[SZ_LINE+1]
bool	clgetb()
double	clgetd()
real	clgetr()
int	clgeti()
int	clscan()
int	strlen(), strcmp()
double	wcs_p2w()

include	"emv.com"
include	"rvsao.com"
include	"contin.com"
include "results.com"

define	restart_  10
define	fitlines_ 20
define	fitvel_   30

begin

	newresults = FALSE
	dindef = INDEFD
	drms = 0.d0
	drms = clgetd ("disperr")
	meanwidth = 0.d0
	iord = 0
	spectype = 1
	do i = 1, MAXREF {
	    override[i] = 0
	    }
	call strcpy (LABEL(specsh),xlab,SZ_LINE)
	if (strlen (UNITS(specsh)) > 0) {
	    call strcat (" in ",xlab,SZ_LINE)
	    call strcat (UNITS(specsh),xlab,SZ_LINE)
	    }

#  Set wavelength limts for line-labelled plot
	call wcs_set (specsh)
	wred = wcs_p2w (double (pix2))
	wblue = wcs_p2w (double (pix1))
	if (wred < wblue) {
	    wblue = wcs_p2w (double (pix2))
	    wred = wcs_p2w (double (pix1))
	    }


# Line peak must have highest value in +/- this number of pixels
	npfit = 2
	npfit = clgeti ("npfit")

# Number of coefficients to fit for line continuum
	nlcont = 2
	nlcont = clgeti ("nlcont")

# Line search wavelength half-width in angstroms
	wspan = 10.d0
	wspan = clgeti ("wspan")

# Line peak must be this number of standard deviations above continuum
	zsig = 2.d0
	zsig = clgetd ("linesig")

	tczrms = 0.d0
	tczchi = 0.d0

	debug = clgetb ("debug")

#  Spectrum smoothing for line fitting (Use with care!)
	esmooth = 0
	esmooth = clgeti ("esmooth")

#  Continuum parameters
	call cont_get_pars()
	abrej[1] = 0.
	emrej[1] = 0.
	if (clscan("mincont") != EOF)
	    mincont = clgetd ("mincont")
	else
	    mincont = 0.d0

	sp = NULL
	sfspec = NULL
	sfcont = NULL
	smcont = NULL
	smspec = NULL

#  Allocate unsmoothed spectrum to be fit
	call smark (sp)
	call salloc (sfspec,specpix,TY_REAL)
	if (sfspec == NULL)
	    call printf (" EMFIT:  cannot allocate sfspec\n")

#  Allocate unsmoothed continuum for fit
	call salloc (sfcont, specpix, TY_REAL)
	if (sfcont == NULL)
	    call printf (" EMFIT:  cannot allocate sfcont\n")

#  Allocate smoothed spectrum for plotting
	call salloc (smplot,specpix,TY_REAL)
	if (smplot == NULL)
	    call printf (" EMFIT:  cannot allocate smplot\n")

#  Allocate smoothed spectrum for line search
	call salloc (smspec,specpix,TY_REAL)
	if (smspec == NULL)
	    call printf (" EMFIT:  cannot allocate smspec\n")

#  Allocate smoothed continuum for line search
	call salloc (smcont, specpix, TY_REAL)
	if (smcont == NULL)
	    call printf (" EMFIT:  cannot allocate smcont\n")

#  Set pointers for wavelength region actually used
	npix = pix2 - pix1 + 1
	smused = smspec + pix1 - 1
	smcusd = smcont + pix1 - 1
	sfused = sfspec + pix1 - 1
	sfcusd = sfcont + pix1 - 1

	nsmooth = 10
	nsmooth = clgeti ("nsmooth")

# If plot enabled, show the object and spectra.     
	if (pltspec) {
	    call plotspec (npix,spectrum[pix1],specname,wlspec[pix1],xlab,nsmooth)
	    if (skyspec) {
		call plotspec (npix,sky[pix1],skyname,wlspec[pix1],xlab,nsmooth)
		}
	    }

# Eliminate bad lines from spectrum
	if (clscan("fixbad") != EOF) {
	    if (clgetb ("fixbad")) {
		call filllist (npix, spectrum, wlspec, iord, pixfill, debug)
		if (pltspec)
		    call plotspec (npix,spectrum[pix1],specname,
				   wlspec[pix1],xlab,nsmooth)
		}
	    }

# Find spectrum maximum and minimum values
        spmax = spectrum[pix1]
        spmin = spectrum[pix2]
        do i = pix1, pix2 {
	    specd = spectrum[i]
	    if (specd > spmax)
		spmax = specd
	    if (specd < spmin)
		spmin = specd
	    }
	if (debug) {
	    call printf ("EMFIT: %8g < counts < %8g from %.3fA(%d) - %.3fA(%d)\n")
		call pargd (spmin)
		call pargd (spmax)
		call pargd (wblue)
		call pargi (pix1)
		call pargd (wred)
		call pargi (pix2)
	    call flush (STDOUT)
	    }

# Data renormalization
	renorm = clgetb ("renormalize")

# If RENORM is not set, set it if maximum counts are less than 1
        if (!renorm) {
            if (spmax < 1.0)
                renorm = TRUE
            else
                renorm = FALSE
            }

# Renormalize spectrum if requested
	if (renorm) {

#	If any pixels are negative, add a floor so fits will work
	    if (spmin < 0.0) {
		do i = 1, specpix {
		    spectrum[i] = spectrum[i] - spmin
		    }
		}
	    sptot = 0.d0
	    nsp = 0
	    do i = pix1, pix2 {
		if (spectrum[i] != 0.) {
		    sptot = sptot + double (spectrum[i])
		    nsp = nsp + 1
		    }
		}
	    if (nsp > 0) {
		spmean = 0.001d0 * sptot / double (nsp)
		if (spmean == 0.d0)
		    spmean = 1.d0
		do i = 1, specpix {
		    if (spectrum[i] != 0.) {
			spectrum[i] = spectrum[i] / spmean
			}
		    }
		}
	    else {
		call eprintf ("*** Spectrum is all zeroes\n")
		call close_image (specim, specsh)
		call sfree (sp)
		return
		}
	    if (pltspec)
		call plotspec (npix,spectrum[pix1],specname,wlspec[pix1],
			       xlab,nsmooth)
	    }

# Smooth spectrum vector for plotting
	call amovr (spectrum,Memr[smplot],specpix)
	call smooth (Memr[smplot],specpix,nsmooth)
#	call gsmooth (Memr[smplot],spectrum,specpix,2,nsmooth, 1.0d0)

	linefit = TRUE
	linefit = clgetb ("linefit")
	if (linefit) {
	    velfit = TRUE
	    spvqual = 0
	    }
	else {
	    velfit = FALSE
	    savevel = FALSE
	    }
	czvel = dindef

#  Make copy of spectrum from which to subtract continuum for line finding
restart_

#  Fit and subtract continuum for search
	call amovr (Memr[smplot],Memr[smspec],specpix)
	call icsubcon (npix,Memr[smused],wlspec[pix1],specname,spectype,
		       nsmooth,Memr[smcusd])
	if (debug) call printf ("EMFIT: Continuum subtracted from spectrum\n")

#  Set initial value of Cz before fitting lines
	if (linefit) {
	    switch (vinit) {
		case VSEARCH:
		    call emguess (Memr[smspec],Memr[smcont],specpix,
				  wblue,wred,czvel,debug)
		case VGUESS:
		    czvel = clgetd ("czguess")
		case VCORREL:
		    czvel = spxvel
		case VEMISS:
		    czvel = spevel
		case VCOMB:
		    czvel = spvel
		case ZGUESS:
		    czvel = c0 * clgetd ("czguess")
		case VCORTEMP:
		    if (ntemp > 0) {
			do itemp = 1, ntemp {
			    if (strcmp (cortemp,tempid[1,itemp]) == 0)
				czvel = zvel[itemp]
			    }
			}
		default:
		}
	    if (czvel == dindef)
		czvel = 0.d0
	    }
	cvel = 1.d0 + (czvel / c0)
	if (debug) {
	    call printf ("Initial velocity guess is %8.2f -> 1+z= %6.4f\n")
		call pargd (czvel)
		call pargd (cvel)
	    }
	call flush (STDOUT)

#  Find emission lines
fitlines_
	if (linefit) {
	    call emsrch (Memr[smspec],Memr[smcont],specpix,
			 wblue,wred,cvel,debug)

#	Smooth spectrum for line fitting
	    call amovr (spectrum,Memr[sfspec],specpix)
	    call smooth (Memr[sfused],npix,esmooth)
#	    call gsmooth (Memr[sfused],spectrum,npix,2,esmooth,1.d0)

#	Set up continuum to be subtracted before fitting lines
	    call icsubcon (npix,Memr[sfused],wlspec[pix1],specname,spectype,
			   esmooth,Memr[sfcusd])

#	Fit profiles to lines
	    call emlfit (Memr[sfspec],Memr[sfcont],sky,specpix,debug)
	    }

#  Compute velocity combining results for fit lines
fitvel_
	if (velfit) {
	    newresults = TRUE
	    if (nfound > 0)
	        call emvfit (wblue,wred,drms,meanwidth,debug,
			     tcz,tczrms,tczchi,nz)
	    else {
		spevel = dindef
		nz = 0
		tcz = 0.d0
		tczrms = 0.d0
		}
	    if (tcz != 0.d0)
		spevel = (tcz - 1.d0) * c0 + spechcv
	    else if (nfound > 0)
		spevel = spechcv
	    speerr = tczrms * c0
	    nfit = nz
	    if (nfit == 1) {
		speerr0 = 0.d0
		speerr0 = clgetr ("sigline")
		if (speerr0 > 0)
		    speerr = speerr0
		}

#  Combine emission with cross-correlation velocity
	    spnl = nfound
	    spnlf = nfit
	    call vcombine (spxvel,spxerr,spxr,spevel,speerr,spnlf,spvel,sperr,debug)
	    }

#  Print table of results
	if (debug) {
	    call printf ("EMFIT: Ready to call EMRSLTS\n")
	    call flush (STDOUT)
	    }
	call emrslts (specfile,mspec,specim,rmode)

	if (nfound == 0 && spmin < 0.0) {
	    do i = 1, specpix {
		spectrum[i] = spectrum[i] + spmin
		}
	    }

#  Plot spectrum and label emission and absorption lines
	mspec0 = mspec
	call emplot (specfile,specim,mspec,specpix,spectrum,Memr[smplot],
		     Memr[sfspec], Memr[smspec], wlspec,pix1,pix2)
	if (mspec != mspec0) {
	    return
	    }

#  Start with new initial velocity if switch was set in interactive cursor mode
	if (sfit) {
	    linefit = TRUE
	    velfit = TRUE
	    goto restart_
	    }

#  Refit lines if switch was set in interactive cursor mode
	if (lfit) {
	    linefit = TRUE
	    velfit = TRUE
	#  Fit and subtract continuum from smoothed spectrum for line search
	    call amovr (Memr[smplot],Memr[smspec],specpix)
	    call icsubcon (npix,Memr[smused],wlspec[pix1],specname,spectype,
			   nsmooth, Memr[smcusd])
	    goto fitlines_
	    }

#  Refit velocity if switch was set in interactive cursor mode
	if (vfit) {
	    velfit = TRUE
	    goto fitvel_
	    }
	call sfree (sp)

	if (debug) {
	    call printf ("EMFIT: All done!\n")
	    call flush (STDOUT)
	    }

	return

end
# Oct 24 1991	Subtract continuum using icfit
# Nov 20 1991	Free buffers before malloc'ing them
# Dec  3 1991	Pass continuum fit to line finding and fitting routines
# Dec 12 1991	Set velocity error for single line based on parameter SIGLINE
# Dec 19 1991	Pass through wavelength vector

# Mar 23 1992	Use emission line information from file
# Mar 26 1992	Get number of points to fit from parameter file
# Apr 22 1992	Always use spwl0 and spdwl, not w0 and wpc
# May 29 1992	Change true to TRUE; don't change value of vel until combining
# Aug 11 1992	Set nlfit before calling vcombine which now uses it
# Aug 12 1992	Drop npcont; it's not used
# Aug 13 1992	For single line, fix error bug
# Oct 22 1992	Free smplot before returning from this subroutine
# Nov 19 1992	Get nlcont (#coeffs for line continuum fit) from parameter file
# Nov 24 1992	Zero MAXREF instead of 12 locations in override
# Dec  1 1992	Put spectrum resmoothing into cursor control subroutine

# Feb 10 1993	Use spechcv instead fbcv
# Jun  2 1993	Implement MWCS for pixel<-> wavelength transformations
# Jun 16 1993	Pass wavelength limits in as arguments
# Dec  2 1993	Pass multispec spectrum number to output subroutines

# Mar 23 1994	Pass number of sigma limit in labelled common
# Apr  7 1994	Allow fit spectrum to be smoothed from cursor mode
# Apr 11 1994	Restart more completely after deleting points
# Apr 12 1994	Pass all smoothing parameters in labelled common
# Apr 14 1994	Pass red and blue wavelength limits to EMPLOT
# Apr 22 1994	Drop image structure from EMPLOT argument list
# Apr 22 1994	Use stack pointer salloc instead of malloc
# May  3 1994	Pass smoothed spectrum for plotting as argument
# May  5 1994	Pass blue and red pixel limits as arguments
# May  5 1994	Optionally plot object spectrum before fitting
# May  9 1994	Move smoothing after optional object plot
# May 16 1994	Add report mode switch so 1-line reports can be generated
# Jun  9 1994	Fit continuum only over specified wavelength range
# Jun 15 1994	Use SPECDC instead of DCFLAG to check for log-lambda file
# Jun 23 1994	Eliminate use of getim labelled common
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug  8 1994	Drop specfile from emhead; change filename to specfile
# Aug 15 1994	Read dispersion fit RMS in angstroms as parameter
# Aug 15 1994	Read previous emission line fit
# Aug 16 1994	Add spectrum header to EMPLOT argument list

# Mar 15 1995	Pass specname to ICSUBCON
# Jul 13 1995	Add debugging argument to VCOMBINE
# Jul 13 1995	Drop image header reading with EMRHEAD
# Jul 21 1995	Do not recompute combined velocity if VELFIT not true
# Aug  7 1995	Do not initialize NFIT and NFOUND in this subroutine
# Sep 21 1995	Drop parameter VELFIT; set NEWRESULT
# Sep 25 1995	If not LINEFIT, set SAVEVEL false
# Oct  2 1995	Change initial guess velocity GVEL to CZVEL
# Oct  3 1995	Use GVEL for initial velocity parameter

# Feb 22 1996	Set quality flag to 0 when doing new fit

# Feb  3 1997	Add code to optionally cut bad lines such as night sky lines
# Feb 26 1997	Try Gaussian smoothing
# Mar 14 1997	Set X axis label from spectrum header
# Apr 25 1997	Add zguess as option for initial velocity
# May  2 1997	Test against dindef instead of INDEFD
# May  6 1997	Add NSMOOTH argument to ICSUBCON
# May  9 1997	Add MINCONT for continuum limit for equivalent width
# May 19 1997	Add CORTEMP to get initial velocity from specific template
# May 22 1997	Set NZ to zero if NFOUND is zero; it stayed at last setting
# Sep 30 1997	upgrade comments
# Nov 13 1997	Pass continuum-subtracted spectrum to plotting subroutine

# Feb 12 1998	Fix bug which occurs when single line sigma is zero

# Dec  2 2002	Move renormalization to this subroutine and fix test
# Dec  2 2002	Fix bad pixels before renormalization
# Dec  2 2002	Compute renormalization factor over selected spectrum portion
# Dec  2 2002	Add a floor to make all pixels positive if renorm is yes
# Dec  2 2002	Plot renormalized spectrum

# Nov  3 2005	Take out floor before plotting if no lines found

# Feb  7 2006	Label 1+z correctly

# Jan 30 2007	Drop wavelength limits from input argument list
# Jan 30 2007	Send pixel limits, not wavelength limits to emplot()
# Jan 31 2007	Print both count and wavelength limits on one line
# Mar 30 2007	Keep pixel limit order separate from wavelength limit order

# Mar  5 2008	Add pixfill argument to filllist() call
# May  8 2008	Add order argument to filllist() call
# Oct  7 2008	If mspec is changed by emplot, return immediately
