# File rvsao/Emsao/emplot.x
# May 2, 2014
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1991-2014 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  EMPLOT uses IRAF plot routines to plot a spectrum and label its
#  emission and absorption lines

###########################################################################
 
include	<gset.h>
include	<gio.h>
include	<pkg/gtools.h>
include <ttyset.h>
include <fset.h>
include <imio.h>
include <imhdr.h>
include	"rvsao.h"
include	"emv.h"
 
procedure emplot (specfile,specim,mspec0,npoints,spectrum,pix,cspec,cpix,wpix,
		  ipix01,ipix02)
 
char	specfile[ARB]	# Name of spectrum file
pointer	specim		# IRAF image descriptor structure
int	mspec0		# Spectrum number if multispec, else 0
int	npoints		# Number of points in spectrum
real	spectrum[ARB]	# Points in spectrum (unsmoothed)
real	pix[ARB]	# Points to plot (smoothed)
real	cspec[ARB]	# Continuum-subtracted spectrum (unsmoothed)
real	cpix[ARB]	# Continuum-subtracted points to plot (smoothed)
real	wpix[ARB]	# Vector of corresponding wavelengths
int	ipix01		# First pixel
int	ipix02		# Last pixel
#pointer	gt		# Plotting structure
 
int	ipix1,ipix2	# First and last pixels to plot
int	npix		# Number of points to plot
int	npix0		# Number of points in full plot
int	mspec		# Spectrum number if multispec, else 0
bool	dfirst,zfirst, newgraph
bool	heading		# True to plot filename and velocity in heading
char	command[SZ_FNAME]
int	wc, key, key1, it, it0, ipx1, ipx2, ix, itemp
int	iw,i,llnam,llref,i1,ip, dispmode, mspec1, nmax
real	x, y, x1, y1
double	px, wx, wr, bcz, dtemp
double	diff, dindef
pointer	gt, gfd, gfp
char	device[SZ_FNAME] # Display on which to plot data
bool	cursor		# true if waiting for cursor after plotting
bool	hardcopy	# true to make automatic hard copy and plot
char	plotter[SZ_FNAME] # Device on which to make hard copies
char	lchar		# single character line identification code
char	lname[SZ_LINE]	# name of emission or absorption line
double	velp		# Velocity for plotted shift
double	wblue		# Blue wavelength limit in angstroms
double	wred		# Red wavelength limit in angstroms
double	splvel		# Redshift velocity from single line
int	ndim
int	clgcur()
int	strlen(),strncmp(), ctod(), ctoi()
bool	clgetb()
int	clgeti(), clscan()
int	imaccf()
bool	displot
bool	contsub
double	vcomb,vcerr
bool	reversed
double	wcs_p2w(), wcs_w2p()
int	gt_gcur()
pointer	gopen(), gt_init()
char	linebuf[SZ_LINE]
int	getline()
int	rcompe()
extern	rcompe()

common/emp/ shift
double	shift
 
include	"rvsao.com"
include	"emv.com"
include	"ansum.com"
include	"results.com"
 
begin
	dindef = INDEFD
	heading = TRUE
	plotcorr = FALSE
	c0 = 299792.5d0
	vcomb = spvel
	vcerr = sperr
	mspec = mspec0
	splvel = dindef
	lname[0] = EOS
	ndim = IM_NDIM(specim)
	if (ndim < 2)
	    nmax = 1
	else if (ndim < 3)
	    nmax = IM_LEN(specim,2)
	else if (ndim < 4) {
	    if (IM_LEN(specim,2) > 1)
		nmax = IM_LEN(specim,2)
	    else
		nmax = IM_LEN(specim,3)
	    }
	if (plotem)
            dispem = TRUE
	specref = 0
	waverest = 0.d0
	if (imaccf (specim, "VELSET") == YES) {
            dispem = TRUE
	    dispabs = TRUE
	    }

#  Set wavelength limts
	call wcs_set (specsh)
	wblue = wcs_p2w (double (ipix01))
	wred = wcs_p2w (double (ipix02))
	reversed = FALSE
	if (wblue > wred) {
	    wred = wcs_p2w (double (ipix01))
	    wblue = wcs_p2w (double (ipix02))
	    reversed = TRUE
	    }

	# Set wavelength shift factor
	switch (vplot) {
	    case VEMISS:
		if (spevel != dindef)
		    velp = spevel
		else
		    velp = 0.d0
	    case VCORREL:
		if (spxvel != dindef)
		    velp = spxvel
		else
		    velp = 0.d0
	    case VGUESS:
		velp = gvel
	    default:
		if (spvel != dindef)
		    velp = spvel
		else
		    velp = 0.d0
	    }

	bcz = 1.d0 + (spechcv / c0)
	shift = (1.d0 + (velp / c0)) / bcz
 
	npix0 = ipix02 - ipix01 + 1
	ipix1 = ipix01
	ipix2 = ipix02
	npix = npix0

	newgraph = FALSE
	dfirst = TRUE
	zfirst = TRUE
	lfit = FALSE
	vfit = FALSE
	sfit = FALSE
	if (clscan ("dispmode") == EOF)
	    dispmode = 2
	else
	    dispmode = clgeti ("dispmode")
	if (dispmode == 4) {
	    contsub = TRUE
	    dispmode = 2
	    }
	else if (dispmode == 5) {
	    contsub = TRUE
	    dispmode = 3
	    }
	else
	    contsub = FALSE

#  Make hardcopy if requested
	hardcopy = clgetb ("hardcopy")
	if (hardcopy) {

	# Open plotter
	    call clgstr("plotter", plotter, SZ_FNAME)
	    gfp = gopen (plotter, APPEND, STDGRAPH)

	# Plot spectrum with line labels
	    if ((dispmode == 3) || (npix > 0 && npix < npix0))
		npix = -npix
	    if (contsub)
		call emeplot (gfp, specfile, specim, mspec, cpix[ipix1],
			      wpix[ipix1], npix)
	    else
		call emeplot (gfp, specfile, specim, mspec, pix[ipix1],
			      wpix[ipix1], npix)

	# Flush plotter and close graphics structure
	    call gflush (gfp)
	    call gclose (gfp)
	    gfp = NULL
	    if (debug) {
		call printf ("EMPLOT: Hardcopy sent to %s\n")
		    call pargstr (plotter)
		call flush (STDOUT)
		}
	    }

#  If not displaying data, return
	displot = clgetb ("displot")
	if (!displot) {
	    return
	    }

#  Sort templates by R-value
	if (ntemp > 0) {
	    do i = 1, ntemp {
		itr[i] = i
		}
	    call gqsort (itr,ntemp,rcompe,czr)
	    it0 = 1
	    }

	if (debug && plotem) {
	    call printf ("EMPLOT: Always labelling %d / %d emission lines\n")
		call pargi (nfound)
		call pargi (nref)
	    call flush (STDOUT)
	    }

#  Open display and eliminate y-axis minor ticks.
	call clgstr("device", device, SZ_FNAME)
	gfd = gopen (device, NEW_FILE, STDGRAPH)

#  Plot spectrum with line labels
	if ((dispmode == 3) || (npix > 0 && npix < npix0))
	    npix = -npix
	if (contsub)
	    call emeplot (gfd, specfile, specim, mspec, cpix[ipix1],
			  wpix[ipix1], npix)
	else
	    call emeplot (gfd, specfile, specim, mspec, pix[ipix1],
			  wpix[ipix1], npix)

#  If not using cursor, return
	cursor = clgetb ("curmode")
	if (!cursor) {
	    call gclose (gfd)
	    return
	    }

#  Initialize graphics pointer
	gt = gt_init()
	mspec1 = mspec
 
	while (gt_gcur ("cursor", x, y, wc, key, command, SZ_FNAME) != EOF) {
	    if (key > 64 && key < 91)
		key = key + 32
	    wx = x

# If wavelength returned is outside of plot limits, return appropriate limit
	    if (reversed) {
		if (wx > wpix[ipix1]) {
#		    wx = wpix[ipix1]
		    ix = ipix1
		    }
		else if (wx < wpix[ipix2]) {
#		    wx = wpix[ipix2]
		    ix = ipix2
		    }
		}
	    else {
		if (wx < wpix[ipix1]) {
#		    wx = wpix[ipix1]
		    ix = ipix1
		    }
		else if (wx > wpix[ipix2]) {
#		    wx = wpix[ipix2]
		    ix = ipix2
		    }
		}
	    wr = wx / shift

# Assign wavelength returned to a specific spectrum pixel
	    px = wcs_w2p (wx)
	    ix = int (px)
	    if (px - double (ix) > 0.5) {
		ix = ix + 1
		}

	    switch (key) {

#	Replot spectrum
		case 'p':
		    newgraph = TRUE

#	Replot spectrum in other mode
		case '/':
		    if (dispmode == 3)
			dispmode = 2
		    else
			dispmode = 3
		    newgraph = TRUE

#	Replot spectrum without header
		case 'h':
		    if (mspec > 0)
			mspec = -mspec
		    else if (mspec == 0)
			mspec = -1
		    else if (mspec < 0)
			mspec = -mspec
		    newgraph = TRUE

#	Plot next aperture/order
		case ')':
		    if (mspec0 > 0) {
			if (mspec0 < nmax)
			    mspec0 = mspec0 + 1
			else
			    mspec0 = 1
			}
		    lfit = FALSE
		    vfit = FALSE
		    sfit = FALSE
		    newgraph = FALSE
		    break

#	Plot previous aperture/order
		case '(':
		    if (mspec0 > 1)
			mspec0 = mspec0 - 1
		    else
			mspec0 = nmax
		    lfit = FALSE
		    vfit = FALSE
		    sfit = FALSE
		    newgraph = FALSE
		    break

#	Make hard copy on STDPLOT
		case '@':
		    call gclose (gfd)
		    call clgstr("plotter", plotter, SZ_FNAME)
		    call printf ("Making hard copy on %s\n")
		 	call pargstr (plotter)
		    gfp = gopen (plotter, APPEND, STDGRAPH)
		    if ((dispmode == 3) || (npix > 0 && npix < npix0))
			npix = -npix
		    if (contsub)
			call emeplot (gfp, specfile, specim, mspec,
				      cpix[ipix1], wpix[ipix1], npix)
		    else
			call emeplot (gfp, specfile, specim, mspec,
				      pix[ipix1], wpix[ipix1], npix)
		    call gflush (gfp)
		    call gclose (gfp)
		    gfp = NULL
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Add closest line to fit
		case '+':
		    iw = 0
		    diff = 1.d2
		    do i = 1, nfound {
			if (dabs (wr - wlrest[i]) < diff) {
			    diff = dabs (wr - wlrest[i])
			    iw = i
			    }
			}
		    if (iw > 0) {
			override[iw] = 1
			call gclose (gfd)
			call printf ("Adding %s line at %d angstroms\n")
			    call pargstr (nmobs[1,iw])
			    call pargd (wlrest[iw])
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			vfit = TRUE
			}

#	Subtract closest line from fit
		case '-':
		    iw = 0
		    diff = 1.d2
		    do i = 1, nfound {
			if (dabs (wr - wlrest[i]) < diff) {
			    diff = dabs (wr - wlrest[i])
			    iw = i
			    }
			}
		    if (iw > 0) {
			override[iw] = -1
			call gclose (gfd)
			call printf ("Dropping %s line at %d angstroms\n")
			    call pargstr (nmobs[1,iw])
			    call pargd (wlrest[iw])
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			vfit = TRUE
			}

#	Set quality flag in header
		case 'y':
		    spvqual = 4
		    call gclose (gfd)
		    call printf ("Setting quality flag to 4=good\n")
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    qplot = TRUE
		    IM_UPDATE(specim) = YES
		case 'n':
		    spvqual = 3
		    call gclose (gfd)
		    call printf ("Setting quality flag to 3=bad\n")
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    qplot = TRUE
		    IM_UPDATE(specim) = YES
		case 'j':
		    spvqual = 1
		    call gclose (gfd)
		    call printf ("Setting quality flag to 1=questionable\n")
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    qplot = TRUE
		    IM_UPDATE(specim) = YES

#	Set intitial radial velocity from an absorption line
		case 'a':
		    call gclose (gfd)
		    call printf ("Absorption line (")
		    do i = 1, nabs {
			call printf (" %s")
			    call pargstr (nmabs[1,i])
			}
		    call printf (" or angstroms): ")
		    call flush (STDOUT)
		    i = getline (STDIN,lname)
		    if (i > 1) {
			waverest = 0.d0
			llnam = strlen (lname) - 1
			do i = 1, llnam {
			    if (i == 1 && lname[i] > 96)
				lname[i] = lname[i] - 32
			    if (i == 3 && lname[i] > 96)
				lname[i] = lname[i] - 32
			    }
			specref = 0
			do i = 1, nabs {
#			    call printf ("%d: %s %s [%7.2fA]")
#				call pargi (i)
#				call pargstr (nmabs[1,i])
#				call pargstr (lname)
#				call pargd (wlabs[i])
			    if (strncmp (nmabs[1,i], lname, llnam) == 0) {
				llref = strlen (nmabs[1,i])
				if (llref == llnam) {
				    specref = i
				    waverest = wlabs[i]
#				    call printf ("*")
				    }
				}
#			    call printf ("\n")
			    }
#			call flush (STDOUT)
			if (specref == 0) {
			    i1 = 1
			    i = ctod (lname,i1,waverest)
			    if (i == 0) waverest = 0.d0
			    }
			if (waverest > 0.d0) {
			    cvel = wx / waverest
			    splvel = (cvel - 1.d0) * c0
			    spvel = splvel
			    if (specref > 0) {
				call printf("%s absorption: Observed %7.2fA -> ")
				call pargstr (nmabs[1,specref])
				call pargd (wx)
				specref = -specref
				}
			    else {
				call printf("Absorption: Observed %7.2fA -> ")
				call pargd (wx)
				}
			    call printf ("Rest %7.2fA vel= %7.1f  z= %6.4f\n")
				call pargd (waverest)
				call pargd (spvel)
				call pargd (cvel - 1.d0)
			    call flush (STDOUT)
			    shift = (1.d0 + (spvel / c0)) / bcz
			    gvel = spvel
			    lfit = TRUE
			    vplot = VGUESS
			    dispem = TRUE
			    dispabs = TRUE
			    newgraph = TRUE
#			    plotem = TRUE
			    }
			else {
			    call printf("* No absorption line found: %s\n")
				call pargstr (lname)
			    call flush (STDOUT)
			    }
			}
		    gfd = gopen (device, APPEND, STDGRAPH)
 
#	Replace a segment of spectrum
		case 'd':
		    if (dfirst) {
			ipx1 = ix
			call gclose (gfd)
			call printf("d again to delete data from pixel %f\n")
			    call pargi (ipx1)
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			dfirst = FALSE
			}
		    else {
			dfirst = TRUE
			ipx2 = ix
			call gclose (gfd)
			call printf("filling from pixel %f to pixel %f\n")
			    call pargi (ipx1)
			    call pargi (ipx2)
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			call fillpix (npoints, spectrum, ipx1, ipx2)
			call fillpix (npoints, pix, ipx1, ipx2)
			edited = TRUE
			newgraph = TRUE
			sfit = TRUE
			}
 
#	Replace a segment of spectrum
		case ',':
		    ipx1 = ix
		    call gclose (gfd)
		    call printf("replacing pixel %.0f with %f\n")
			call pargi (ipx1)
			call pargr (y)
		    call fillpt (npoints, spectrum, ipx1, y, debug)
		    call fillpt (npoints, pix, ipx1, y, debug)
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    edited = TRUE
		    newgraph = TRUE
		    sfit = TRUE

#	Select line labelling options
		case 'o':
		    call gclose (gfd)
		    call printf("Label a>bsorption e>mission b>oth lines else none: ")
		    call flush (STDOUT)
		    i = clgcur ("cursor",x1,y1,wc,key1,command,SZ_LINE)
		    lchar = char (key1)
		    call printf("%c\n")
			call pargc (lchar)
		    call flush (STDOUT)
		    switch (lchar) {
			case 'a':
			    call printf("label absorption lines\n")
			    dispabs = TRUE
			    dispem = FALSE
			    plotem = FALSE
			case 'e':
			    call printf("label emission lines\n")
			    dispem = TRUE
			    dispabs = FALSE
#			    plotem = TRUE
			case 'b':
			    call printf("label absorption and emission lines\n")
			    dispem = TRUE
			    dispabs = TRUE
#			    plotem = TRUE
			default:
			    call printf("do not label absorption and emission lines\n")
			    dispem = FALSE
			    dispabs = FALSE
			}
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    newgraph = TRUE

#	Select different initial velocity source
		case 'i':
		    call gclose (gfd)
		    call printf("s>earch x>cor e>mission f>ile g>uess velocity: ")
		    call flush (STDOUT)
		    i = clgcur ("cursor",x1,y1,wc,key1,command,SZ_LINE)
		    lchar = char (key1)
		    call printf("%c\n")
			call pargc (lchar)
		    call flush (STDOUT)
		    switch (lchar) {
			case 'e':
			    call printf("restart from file emission velocity\n")
			    vinit = VEMISS
			    sfit = TRUE
			case 'x':
			    call printf("restart from correlation velocity\n")
			    vinit = VCORREL
			    sfit = TRUE
			case 'c':
			    call printf("restart from combination velocity\n")
			    vinit = VCOMB
			    sfit = TRUE
			case 'g':
			    call printf("\nRestart using this velocity: ")
			    call flush (STDOUT)
			    i = getline (STDIN,linebuf)
			    if (i > 1) {
				ip = 1
				i = ctod (linebuf,ip,dtemp)
				spvel = dtemp
				vinit = VCOMB
				sfit = TRUE
				}
			default:
			    call printf("restart and search for velocity\n")
			    vinit = VSEARCH
			    sfit = TRUE
			}
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Toggle continuum subtraction
		case 'k':
		    call gclose (gfd)
		    if (contsub) {
			contsub = FALSE
			call printf("Plot spectrum including continuum\n")
			}
		    else {
			contsub = TRUE
			call printf("Plot continuum-subtracted spectrum\n")
			}
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    newgraph = TRUE

#	Set line fitting parameters
		case 'l':
		    call gclose (gfd)
		    call printf("Line search parameters\n")
		    call printf("Number of sigma above continuum (%4.2f): ")
			call pargd (zsig)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctod (linebuf,ip,dtemp)
			zsig = dtemp
			lfit = TRUE
			}

		    call printf("Wavelength to search around redshifted line center (%4.1f):")
			call pargd (wspan)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctod (linebuf,ip,dtemp)
			wspan = dtemp
			sfit = TRUE
			}

		    call printf("Number of pixels to fit around redshifted line center (%d):")
			call pargi (npfit)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,itemp)
			npfit = itemp
			lfit = TRUE
			}
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Set number of times to smooth fit spectrum
		case 'm':
		    call gclose (gfd)
		    call printf("Number of times to smooth fit spectrum (%d): ")
			call pargi (esmooth)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,itemp)
			esmooth = itemp
			lfit = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Set blue limit of spectrum
		case 'b':
		    call gclose (gfd)
		    call printf("Blue wavelength limit in angstroms (%8.1f): ")
			call pargd (wblue)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctod (linebuf,ip,dtemp)
			wblue = dtemp
			if (reversed)
			    ipix2 = int (wcs_w2p (wblue) + 0.99d0)
			else
			    ipix1 = int (wcs_w2p (wblue) + 0.99d0)
			npix = ipix2 - ipix1 + 1
			sfit = TRUE
			newgraph = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Set red limit of spectrum
		case 'r':
		    call gclose (gfd)
		    call printf("Red wavelength limit in angstroms (%8.1f): ")
			call pargd (wred)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctod (linebuf,ip,dtemp)
			wred = dtemp
			if (reversed)
			    ipix1 = int (wcs_w2p (wred) - 0.99d0)
			else
			    ipix2 = int (wcs_w2p (wred) - 0.99d0)
			npix = ipix2 - ipix1 + 1
			newgraph = TRUE
			sfit = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Set number of times to smooth spectrum
		case 'g':
		    call gclose (gfd)
		    call printf("Number of times to smooth spectrum (%d): ")
			call pargi (nsmooth)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,itemp)
			nsmooth = itemp
			call printf("Smoothing spectrum %d times\n")
			    call pargi (nsmooth)
			call flush (STDOUT)
			call amovr (spectrum,pix,npoints)
			call smooth (pix,npoints,nsmooth)
			call amovr (cspec,cpix,npoints)
			call smooth (cpix,npoints,nsmooth)
			sfit = TRUE
			newgraph = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Select different velocity for final redshift
		case 's':
		    call gclose (gfd)
		    call printf("set VELOCITY to e>m x>corr c>omb l>ine g>uess s>et: ")
		    call flush (STDOUT)
		    i = clgcur ("cursor",x1,y1,wc,key1,command,SZ_LINE)
		    lchar = char (key1)
		    call printf("%c\n")
			call pargc (lchar)
		    call flush (STDOUT)
		    switch (lchar) {
			case 'c':
			    call printf("combined velocity redshift used\n")
			    spvel = vcomb
			    sperr = vcerr
			    specref = 0
			    waverest = 0.d0
			case 'e':
			    call printf("emission velocity redshift used\n")
			    spvel = spevel
			    sperr = speerr
			    specref = 0
			    waverest = 0.d0
			case 'x':
			    call printf("correlation velocity redshift used\n")
			    spvel = spxvel
			    sperr = spxerr
			    specref = 0
			    waverest = 0.d0
			case 'l':
			    if (specref > 0) {
				call printf("velocity %.3f km/sec from %s line used\n")
				    call pargd (splvel)
				    call pargstr (nmref[1,specref])
				spvel = splvel
				sperr = 300.d0
				}
			    else if (specref < 0) {
				call printf("velocity %.3f km/sec from %s line used\n")
				    call pargd (splvel)
				    call pargstr (nmabs[1,-specref])
				spvel = splvel
				sperr = 300.d0
				}
			    else if (waverest > 0) {
				call printf("velocity %.3f km/sec from %.3fA line used\n")
				    call pargd (splvel)
				    call pargd (waverest)
				spvel = splvel
				sperr = 300.d0
				}
			    else
				call printf ("no lines have been identified\n")
			case 'g':
			    call printf("last guessed velocity, %.3f km/sec, used\n")
				call pargd (spvel)
			    call printf("Error in km/sec = ")
			    call flush (STDOUT)
			    i = getline (STDIN,linebuf)
			    if (i > 1) {
				ip = 1
				i = ctod (linebuf,ip,sperr)
				}
			    specref = 0
			    waverest = 0.d0
			default:
			    call printf("Velocity in km/sec = ")
			    call flush (STDOUT)
			    i = getline (STDIN,linebuf)
			    if (i > 1) {
				ip = 1
				i = ctod (linebuf,ip,spvel)
				}
			    call printf("Error in km/sec = ")
			    call flush (STDOUT)
			    i = getline (STDIN,linebuf)
			    if (i > 1) {
				ip = 1
				i = ctod (linebuf,ip,sperr)
				}
			}
		    lfit = FALSE
		    vfit = FALSE
		    sfit = FALSE
		    vplot = VCOMB
		    IM_UPDATE(specim) = YES
		    newgraph = TRUE
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Select different velocity for redshift
		case 'v':
		    call gclose (gfd)
			    call printf("e>m x>corr c>ombined g>uessed velocity: ")
			    call flush (STDOUT)
		    i = clgcur ("cursor",x1,y1,wc,key1,command,SZ_LINE)
		    lchar = char (key1)
		    call printf("%c\n")
			call pargc (lchar)
		    call flush (STDOUT)
		    switch (lchar) {
			case 'e':
			    call printf("emission velocity redshift used\n")
			    vplot = VEMISS
			case 'x':
			    call printf("correlation velocity redshift used\n")
			    vplot = VCORREL
			case 'c':
			    call printf("combination velocity redshift used\n")
			    vplot = VCOMB
			default:
			    call printf("last guessed velocity redshift used\n")
			    vplot = VGUESS
			}
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    newgraph = TRUE
		    lfit = TRUE

#	Adjust continuum parameters
		case 'c':
		    call gclose (gfd)
		    call printf("Continuum order: ")
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,itemp)
			call set_cn_order (itemp)
			lfit = TRUE
			}

	    case 't':
		if (ntemp > 0) {
		    call gclose (gfd)
		    call printf("template for redshift = ")
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,it)
			if (it < 1 || it > ntemp)
			    it = it0
			}
		    else
		        it = it0
		    itemp = itr[it]
		    itmax = itr[it]
		    it0 = it
		    spxvel = zvel[itemp]
		    spxerr = czerr[itemp]
		    spxr = czr[itemp]
		    vplot = VCORREL
		    gfd = gopen (device, APPEND, STDGRAPH)
		    newgraph = TRUE
		    IM_UPDATE(specim) = YES
		    }

#	Set intitial radial velocity from an emission line
		case 'e':
		    call gclose (gfd)
		    call printf ("Emission line (")
		    do i = 1, nref {
			call printf (" %s")
			    call pargstr (nmref[1,i])
			}
		    call printf (" or angstroms): ")
		    call flush (STDOUT)
		    i = getline (STDIN,lname)
		    if (i > 1) {
			llnam = strlen (lname) - 1
			do i = 1, llnam {
			    if (i == 1 && lname[i] > 96)
				lname[i] = lname[i] - 32
			    if (i > 1 && (lname[i] == 105 || lname[i] == 118))
				lname[i] = lname[i] - 32
			    }
			specref = 0
			waverest = 0.d0
			do i = 1, nref {
#			    call printf ("%d: %s %s [%7.2fA]")
#				call pargi (i)
#				call pargstr (nmref[1,i])
#				call pargstr (lname)
#				call pargd (wlref[i])
			    if (strncmp (nmref[1,i], lname, llnam) == 0) {
				llref = strlen (nmref[1,i])
				if (llref == llnam) {
				    specref = i
				    waverest = wlref[i]
#				    call printf ("*")
				    }
				}
#			    call printf ("\n")
			    }
			if (specref == 0) {
			    i1 = 1
			    i = ctod (lname,i1,waverest)
			    if (i == 0) waverest = 0.d0
			    }
			if (waverest > 0.d0) {
			    cvel = wx / waverest
			    splvel = (cvel - 1.d0) * c0
			    spvel = splvel
			    if (specref > 0) {
				call printf("%s emission: Observed %7.2fA -> ")
				call pargstr (nmref[1,specref])
				call pargd (wx)
				}
			    else {
				call printf("Emission: Observed %7.2fA -> ")
				call pargd (wx)
				}
			    call printf("Rest %7.2fA vel= %7.1f  z= %6.4f\n")
				call pargd (waverest)
				call pargd (spvel)
				call pargd (cvel-1.d0)
			    call flush (STDOUT)
			    shift = (1.d0 + (spvel / c0)) / bcz
			    gvel = spvel
			    vplot = VGUESS
			    lfit = TRUE
			    dispem = TRUE
			    dispabs = TRUE
#			    plotem = TRUE
			    newgraph = TRUE
			    }
			else {
			    call printf("* No emission line found: %s\n")
				call pargstr (lname)
			    call flush (STDOUT)
			    }
		    }
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Reset zoom and replacement switches
		case '.':
		    if (!zfirst) {
			call gclose (gfd)
			call printf("zoom cancelled\n")
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			ipix1 = ipix01
			ipix2 = ipix02
			npix = ipix2 - ipix1 + 1
			zfirst = TRUE
			}
		    if (!dfirst) {
			call gclose (gfd)
			call printf("delete cancelled\n")
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			dfirst = TRUE
			}

#	Zoom in on section of graph
		case 'z':
		    if (zfirst) {
			if (ix < 1)
			    ipx1 = ipix1
			else
			    ipx1 = ix
			zfirst = FALSE
			call gclose (gfd)
			call printf("z again to zoom from %8.2f angstroms\n")
			    call pargr (wpix[ipx1])
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			}
		    else {
			if (ix < 1)
			    ipx2 = ipix2
			else
			    ipx2 = ix
			zfirst = TRUE
			if (ipx2 < ipx1) {
			    ipix1 = ipx2
			    ipix2 = ipx1
			    }
			else {
			    ipix1 = ipx1
			    ipix2 = ipx2
			    }
			npix = ipix2 - ipix1 + 1
			call gclose (gfd)
			call printf("zooming from %8.2f to %8.2f angstroms")
			if (wpix[ipix1] > wpix[ipix2]) {
			    call pargr (wpix[ipix2])
			    call pargr (wpix[ipix1])
			    }
			else {
			    call pargr (wpix[ipix1])
			    call pargr (wpix[ipix2])
			    }
			call printf(" (pixels %d - %d)\n")
			    call pargi (ipix1)
			    call pargi (ipix2)
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			newgraph = TRUE
			}

#	Unzoom out from section of graph
		case 'u':
		    ipix1 = ipix01
		    ipix2 = ipix02
		    npix = npix0
		    call gclose (gfd)
		    call printf("unzoom to %8.2f to %8.2f angstroms")
		    if (wpix[ipix1] > wpix[ipix2]) {
			call pargr (wpix[ipix2])
			call pargr (wpix[ipix1])
			}
		    else {
			call pargr (wpix[ipix1])
			call pargr (wpix[ipix2])
			}
		    call printf(" (pixels %d - %d)\n")
			call pargi (ipix1)
			call pargi (ipix2)
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    newgraph = TRUE

#	Print list of cursor commands
		case ' ','?':
		    call gclose (gfd)
		    call emcursor
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
 
#	Print rest and observed wavelengths under cursor
		case 'w':
		    call gclose (gfd)
		    call printf ("Wavelength: Rest: %9.3f -> Observed: %9.3f Pixel %9.3f = %8.2f\n")
		    call pargd (wr)
		    call pargd (wx)
		    call pargd (px)
		    call pargr (y)
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Refit emission lines or velocity
		case 'f':
		    xfit = TRUE
		    break

#	Return without refitting
		case ']':
		    lfit = FALSE
		    vfit = FALSE
		    sfit = FALSE
		    newgraph = FALSE
		    savespec = FALSE
		    break

		case 'q':
		    lfit = FALSE
		    vfit = FALSE
		    sfit = FALSE
		    newgraph = FALSE
		    if (edited)
			savespec = TRUE
		    break

		case 'x':
		    lfit = FALSE
		    vfit = FALSE
		    sfit = FALSE
		    newgraph = FALSE
		    edited = FALSE
		    savespec = FALSE
		    plotcorr = TRUE
		    break
		default:
		}
 
	    if (newgraph) {
		if ((dispmode == 3) || (npix > 0 && npix < npix0))
		    npix = -npix
		if (contsub)
		    call emeplot (gfd, specfile, specim, mspec, cpix[ipix1],
				  wpix[ipix1], npix)
		else
		    call emeplot (gfd, specfile, specim, mspec, pix[ipix1],
				  wpix[ipix1], npix)
		newgraph = FALSE
		}
	    }

# Close up graph window
	call gclose (gfd)
	call gt_free (gt)
	if (debug && (mspec != mspec0)) {
	    call printf ("EMPLOT: Aperture/order %d -> %d out of %d\n")
		call pargi (mspec)
		call pargi (mspec0)
		call pargi (nmax)
	    call flush (STDOUT)
	    }
	return
end

procedure emcursor()
begin
	call printf ("*** Emission line display cursor commands ***\n")
	call printf ("a  Set redshift guess from absorption line\n")
	call printf ("b  Set blue limit of line search\n")
	call printf ("c  Change continuum parameters\n")
	call printf ("d  Delete data between 1st and 2nd positions\n")
	call printf ("e  Set redshift guess from emission line\n")
	call printf ("f  Refit redshift\n")
	call printf ("g  Number of times to smooth plotted spectrum\n")
	call printf ("h  Toggle print of heading with filename and redshift\n")
	call printf ("i  Change initial velocity for search\n")
	call printf ("j  Conditional velocity\n")
	call printf ("k  Plot continuum-subtracted spectrum\n")
	call printf ("l  Line search parameters\n")
	call printf ("m  Number of times to smooth fit spectrum\n")
	call printf ("n  Disapprove velocity\n")
	call printf ("o  Turn line labelling on and off\n")
	call printf ("p  Replot current graph\n")
	call printf ("q  Leave plot, writing out edited spectrum\n")
	call printf ("r  Set red limit of line search\n")
	call printf ("s  Set VELOCITY to a specific value\n")
	call printf ("t  Switch template for correlation velocity\n")
	call printf ("u  Unzoom\n")
	call printf ("v  Plot at e>mission x>correlation c>ombined velocity\n")
	call printf ("w  Show rest and observed wavelengths\n")
	call printf ("x  Plot correlation if available or exit\n")
	call printf ("y  Approve velocity\n")
	call printf ("z  Zoom between 1st and 2nd positions\n")
	call printf (".  Cancel delete or zoom\n")
	call printf (",  Set a single pixel value to the cursor location\n")
	call printf ("/  Toggle plot between full screen and lines\n")
	call printf ("+  Add emission line to fit\n")
	call printf ("-  Subtract emission line from fit\n")
	call printf ("(  Plot previous aperture or order\n")
	call printf (")  Plot next aperture or order\n")
	call printf ("]  Leave plot without writing spectrum\n")
	call printf ("@  Make hard copy of screen\n")
	call printf ("?  Display this menu\n")
	return
end


int procedure rcompe (rvalue, itemp1, itemp2)

double	rvalue[ARB]
int	itemp1
int	itemp2

begin
	if (rvalue[itemp1] < rvalue[itemp2])
	    return (1)
	else if (rvalue[itemp1] > rvalue[itemp2])
	    return (-1)
	else
	    return (0)
end

# Dec 19 1991	Fix wavelength zoom logging bug (pargr instead of pargd)

# Feb 14 1992	Add 's' to set redshift from rest wavelength
#		Reset shift when seting vel with a, s, or e
# Mar 24 1992	Use absorption and emission lines from files
#		Allow velocity setting from all lines in files by name
# Mar 31 1992	Accept line names in upper or lower case or rest wavelength
# Apr 22 1992	Reinitialize ctod index to 1 before each call
#		Fix case-changing for absorbtion and emission line names
# May 22 1992	Don't reset vel to velxc if combined used
# May 28 1992	Fix bug with 'e' and 'a' line position setting
# May 29 1992	Clean up line position setting
# Nov 24 1992	Set plot velocity from guess as option
# Nov 30 1992	Fix n command to do immediate resmoothing

# Jan 15 1993	Set ix to limits if outside limits
# Feb 10 1993	Compare lengths as well as characters for 'e' and 'a'
# Dec  2 1993	Pass multispec spectrum number to EMEPLOT

# Feb  3 1994	Compare rest rather than observed velocities for + and -
# Mar 23 1994	Add option to set line search parameters
# Apr  7 1994	Fix smoothing
# Apr  8 1994	Add option to smooth fit spectrum, not just plotted spectrum
# Apr 11 1994	Return from fewer places; edit both raw and smoothed spectra
# Apr 11 1994	Add new flag to restart edited spectrum from scratch
# Apr 11 1994	Add option to reset initial guess
# Apr 12 1994	Add heading to cursor menu; pass esmooth in labelled common
# Apr 13 1994	Clean up paremter entry
# Apr 14 1994	Pass red and blue wavelength limits as arguments
# Apr 22 1994	Drop image structure from argument list
# Apr 26 1994	Fix bug to use default when changing parameters
# May  2 1994	Set guessed velocity appropriately
# Jun 15 1994	Make all input consistent and add error checking
# Jun 23 1994	Set velocity in gquot, not getim labelled common
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 10 1994	If returning with q, do not refit or replot
# Aug 16 1994	Add image header to EMEPLOT calls
# Nov 16 1994	Add y/n/j quality flag setting; change n command to g

# Jan  6 1995	Add feedback to quality flag setting
# Jul 13 1995	Add option to force combined velocty to a specific value
# Jul 19 1995	Use zero if encountering INDEF velocity
# Sep 19 1995	Add x=exit to correlation plot; change undo from x to period
# Sep 20 1995	Update image header with quality flag immediately on setting
# Sep 21 1995	Set NEWRESULTS and QPLOT
# Sep 25 1995	Set XFIT instead of NEWRESULTS
# Oct  2 1995	Use GVEL instead of reading parameter CZGUESS
# Oct  4 1995	Clean up zoom and unzoom code
# Oct 31 1995	Add full-screen display mode
# Oct 31 1995	Add code to switch template used for xcor velocity

# Jan 24 1996	Fix single key responses

# Jan 30 1997	Use FILLPIX instead of FIXIT to remove bad portions of spectra
# Mar 18 1997	Set LFIT true if changing velocity used
# Apr  7 1997	Pass sorted template indices through labelled common
# Sep 22 1997	Add option to turn line labelling on and off
# Sep 30 1997	Allow change in number of points fit per line
# Nov 13 1997	Add option to plot continuum-subtracted spectrum in summary
# Dec 15 1997	Use different graphics descriptor for hard copies

# Jul 31 1998	Add toggling of heading

# Jun 22 1999	Print actual name of plotter device when using @ and hardcopy
# Jun 22 1999	If display mode is 4 or 5, plot continuum-subtracted spectrum
# Jun 25 1999	Change c command to s for setting VELOCITY; old s in a and e
# Jun 25 1999	Add c command to change continuum parameters
# Sep 24 1999	Set IM_UPDATE if quality flag is set

# Mar 28 2001	Allow spectrum to run from right to left
# Jul 31 2001	Update spectrum header if xcor template changed

# May 12 2004	Drop dlw and dw from common/emp/
# Jul 20 2004	Set IM_UPDATE flag if velocity is changed
# Aug 25 2004	Fix bug by changing a call from ctod() to ctoi()

# Feb  7 2006	Print z instead of 1+z
# Apr 14 2006	Add rcompe() to solve problem on Fedora Core 5

# Jan 30 2007	Print message if no matching line found
# Jan 30 2007	Add 'l' option in 's' to set directly from identified line
# Jan 31 2007	Pass in pixel limits instead of wavelength limits
# Jan 31 2007	Return reference line used for manual velocity as specref
# Mar 30 2007	Fix cursor to handle reversed spectrum

# Sep 11 2008	Add a decimal place to wavelengths and print FP pixel in "w"
# Oct  7 2008	Add ')' and '(' commands to move through apertures like splot

# May 20 2009	Add debugging statement if aperture is changed
# Nov 18 2009	Replot spectrum after setting velocity from line with 'a' or 'e'
# Nov 19 2009	When zooming with 'z', set x to limit if off graph (left then right)
# Nov 20 2009	When switching template using 't', switch vplot to VCORREL

# Mar 15 2013	Fix bug which reversed red and blue wavelength limits

# Apr 30 2014	Add , to set one pixel from cursor position
# Apr 30 2014	Fix round-off error when converting to integer pixel
# May  2 2014	Write edited spectrum when exiting
# May  2 2014	Add ] command to exit without writing edited spectrum
