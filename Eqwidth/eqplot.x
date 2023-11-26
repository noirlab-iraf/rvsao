# File rvsao/Eqwidth/eqplot.x
# September 21, 2005
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 2005 Smithsonian Astrophysical Observatory
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
include	"rvsao.h"
 
procedure eqplot (sh,specfile,specim,npoints,spectrum,pix,cspec,
		  cpix,wpix,bands,nbands)
 
pointer	sh		# Spectrum data structure
char	specfile[ARB]	# Name of spectrum file
pointer	specim		# IRAF image descriptor structure
int	npoints		# Number of points in spectrum
real	spectrum[ARB]	# Points in spectrum (unsmoothed)
real	pix[ARB]	# Points to plot (smoothed)
real	cspec[ARB]	# Continuum-subtracted spectrum (unsmoothed)
real	cpix[ARB]	# Continuum-subtracted points to plot (smoothed)
real	wpix[ARB]	# Vector of corresponding wavelengths
pointer	bands		# Equivalent width band data structure
int	nbands		# Number of equivalent width bands
 
int	ipix1,ipix01	# First pixel to plot
int	ipix2,ipix02	# Last pixel to plot
int	npix		# Lumber of points to plot
int	npix0		# Number of points in full plot
bool	dfirst,zfirst, newgraph
bool	heading		# True to plot filename and velocity in heading
char	command[SZ_FNAME]
int	wc, key, ipx1, ipx2, ix, itemp
int	i,ip
real	x, y
double	wx, wr
double	dindef
pointer	gt, gfd, gfp
char	device[SZ_FNAME] # Display on which to plot data
bool	cursor		# true if waiting for cursor after plotting
bool	hardcopy	# true to make automatic hard copy and plot
char	plotter[SZ_FNAME] # Device on which to make hard copies
int	ctoi()
bool	clgetb()
int	gt_gcur()
pointer	gopen(), gt_init()
char	linebuf[SZ_LINE]
int	getline()
extern	rcomp()

common/emp/ shift
double	shift
 
include	"rvsao.com"
include	"eqw.com"
 
begin
	dindef = INDEFD
	heading = TRUE

	# Initial plot of data array
	ipix01 = 1
	ipix02 = npoints
	npix0 = npoints
	ipix1 = ipix01
	ipix2 = ipix02
	npix = npix0

	newgraph = FALSE
	dfirst = TRUE
	zfirst = TRUE

#  Make hardcopy if requested
	hardcopy = clgetb ("hardcopy")
	if (hardcopy) {

	# Open plotter
	    call clgstr("plotter", plotter, SZ_FNAME)
	    gfp = gopen (plotter, APPEND, STDGRAPH)

	# Plot spectrum with line labels
	    call eqwplot (gfp, specfile, specim, cpix[ipix1],
			  wpix[ipix1], npix, bands, nbands)

	# Flush plotter and close graphics structure
	    call gflush (gfp)
	    call gclose (gfp)
	    gfp = NULL
	    if (debug) {
		call printf ("EQPLOT: Hardcopy sent to %s\n")
		    call pargstr (plotter)
		}
	    }

#  Open display and eliminate y-axis minor ticks.
	call clgstr("device", device, SZ_FNAME)
	gfd = gopen (device, NEW_FILE, STDGRAPH)

#  Plot spectrum with line labels
	call eqwplot (gfd, specfile, specim, cpix[ipix1],
		      wpix[ipix1], npix, bands, nbands)

#  If not using cursor, return
	cursor = clgetb ("curmode")
	if (!cursor) {
	    call gclose (gfd)
	    return
	    }

#  Initialize graphics pointer
	gt = gt_init()
 
	while (gt_gcur ("cursor", x, y, wc, key, command, SZ_FNAME) != EOF) {
	    if (key > 64 && key < 91)
		key = key + 32
	    wx = x
	    wr = wx
	    if (wx < wpix[ipix1])
		ix = ipix1
	    if (wx > wpix[ipix2])
		ix = ipix2
	    ix = ipix1
	    do i = ipix1, ipix2 {
		if (wx >= wpix[i] && wx < wpix[i+1]) {
		    if (wx - wpix[i] > 0.5)
			ix = i
		    else
			ix = i + 1
		    break
		    }
		}

	    switch (key) {

#	Replot spectrum
		case 'p':
		    newgraph = TRUE

#	Make hard copy on STDPLOT
		case '@':
		    call gclose (gfd)
		    call clgstr("plotter", plotter, SZ_FNAME)
		    call printf ("Making hard copy on %s\n")
		 	call pargstr (plotter)
		    gfp = gopen (plotter, APPEND, STDGRAPH)
		    call eqwplot (gfp, specfile, specim, cpix[ipix1],
				  wpix[ipix1], npix, bands, nbands)
		    call gflush (gfp)
		    call gclose (gfp)
		    gfp = NULL
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
			newgraph = TRUE
			}

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
			newgraph = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Reset zoom and replacement switches
		case '.':
		    if (!zfirst) {
			call gclose (gfd)
			call printf("zoom cancelled\n")
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			ipix1 = 1
			ipix2 = npoints
			npix = npoints
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
			ipx1 = ix
			zfirst = FALSE
			call gclose (gfd)
			call printf("z again to zoom from %8.2f angstroms\n")
			    call pargr (wpix[ipix1])
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			}
		    else {
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
		    call eqcursor
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
 
#	Print rest and observed wavelengths under cursor
		case 'w':
		    call gclose (gfd)
		    call printf ("Wavelength: Rest: %8.2f -> Observed: %8.2f Pixel %d = %8.2f\n")
		    call pargd (wr)
		    call pargd (wx)
		    call pargi (ix)
		    call pargr (y)
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Return without refitting
		case 'q':
		    newgraph = FALSE
		    break

		case 'x':
		    newgraph = FALSE
		    break
		default:
		}
 
	    if (newgraph) {
		call eqwplot (gfd, specfile, specim, cpix[ipix1],
			      wpix[ipix1], npix, bands, nbands)
		newgraph = FALSE
		}
	    }

# If leaving with q, do not refit
	if (key == 'q') {
	    newgraph = FALSE
	    }
 
# Close up graph window
 
	call gclose (gfd)
	call gt_free (gt)
	return
end

procedure eqcursor()
begin
	call printf ("*** Emission line display cursor commands ***\n")
	call printf ("b  set blue limit of line search\n")
	call printf ("d  delete data between 1st and 2nd positions\n")
	call printf ("g  number of times to smooth plotted spectrum\n")
	call printf ("h  toggle print of heading with filename and redshift\n")
	call printf ("p  replot current graph\n")
	call printf ("q  leave plot\n")
	call printf ("r  set red limit of line search\n")
	call printf ("u  unzoom\n")
	call printf ("w  show rest and observed wavelengths\n")
	call printf ("z  zoom between 1st and 2nd positions\n")
	call printf (".  cancel delete or zoom\n")
	call printf ("@  make hard copy of screen\n")
	call printf ("?  this menu\n")
	return
end
# Sep 21 2005	New subroutine based on Emsao/emplot.x
