# File rvsao/Util/plotspec.x
# Mach 21, 2001
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University

# Copyright(c) 2001 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  PLOTSPEC uses IRAF plot routines to plot a spectrum on STDGRAPH.
#  FIXIT edits a spectrum by interpolating between marked points.		
#  Adapted from splot in the IRAF onedspec package
 
include	<gset.h>
include	<pkg/gtools.h>
include	"rvsao.h"
 
procedure plotspec (npoints, pix, gtitle, xvec, xlabel, ismooth0)
 
int	npoints		# Number of points in spectrum
real	pix[ARB]	# Unsmoothed vector of points to plot
char	gtitle[ARB]	# Title for the plot
real	xvec[ARB]	# X coordinates for plotted vector
char	xlabel[ARB]	# Label for x axis
int	ismooth0	# Number of times to smooth data
			# (if <0, don't wait for cursor response)

int	ipix1		# first pixel to plot
int	ipix2		# last pixel to plot
int	npix		# number of points to plot
char	command[SZ_FNAME]
int	wc, key, itemp, ip, ix
bool	newgraph, dfirst, zfirst, heading, reversed
int	i, d1, d2, ipx1, ipx2, iwidth, i1, i2, ix0, ismooth
real	wx, wy, xmax, wx0
pointer	gt, gfd, gfp, sp
pointer	smpix		# Smoothed vector of points to plot
char    device[SZ_FNAME] # Display on which to plot data
char    plotter[SZ_FNAME] # Device on which to make hard copies
char    linebuf[SZ_LINE] # Line buffer for input
 
int	gt_gcur(), clgeti(), ctoi(), getline()
pointer	gopen(), gt_init()
 
include	"rvsao.com"
 
begin
	heading = TRUE

# Open display and eliminate y-axis minor ticks.
	call clgstr("device", device, SZ_FNAME)
	gfd = gopen (device, NEW_FILE, STDGRAPH)
	If (gfd == ERR) {
	    call printf ("PLOTSPEC:  cannot open display %s\n")
		call pargstr (device)
	    return
	    }
	call gseti (gfd, G_YNMINOR, NO)
	call gseti (gfd, G_XNMINOR, YES)
	call gsetr (gfd, G_TICKLABELSIZE, 0.8)
 
# Initialize graph format
	gt = gt_init()
	call loghead (taskname,linebuf)
	call gt_seti (gt, GTSYSID, NO)
	call gt_sets (gt, GTTITLE, linebuf)
	call gt_sets (gt, GTSUBTITLE, gtitle)
	call gt_sets (gt, GTXLABEL, xlabel)

# Set up vector to hold smoothed data to be plotted
	if (ismooth0 < 0)
	    ismooth = -(ismooth0 + 1)
	else
	    ismooth = ismooth0
	call smark (sp)
	call salloc (smpix, npoints, TY_REAL)
	call amovr (pix,Memr[smpix],npoints)
	call smooth (Memr[smpix],npoints,ismooth)
 
# Initial plot of data array
	ipix1 = 1
	ipix2 = npoints
	npix = npoints
	call replot (gfd, gt, npix, Memr[smpix+ipix1-1], xvec[ipix1])
	call gflush (gfd)

# Return without interaction
	if (ismooth0 < 0) {
	    call gclose (gfd)
	    call gt_free (gt)
	    call sfree (sp)
	    return
	    }

	newgraph = FALSE
	dfirst = true
	zfirst = true
	if (xvec[ipix1] > xvec[ipix2])
	    reversed = TRUE
	else
	    reversed = FALSE
 
# Handle cursor interaction
	while (gt_gcur ("cursor", wx, wy, wc, key, command, SZ_FNAME) != EOF) {

# Find closest pixel within wavelength limits
	    if (reversed) {
		if (wx >= xvec[ipix1])
		    ix = ipix1
		else if (wx <= xvec[ipix2])
		    ix = ipix2
		else {
		    i1 = ipix1
		    i2 = ipix2 - 1
		    do i = i1, i2 {
			if (wx < xvec[i] && wx >= xvec[i+1]) {
			    if (wx - xvec[i+1] > 0.5)
				ix = i + 1
			    else
				ix = i
			    break
			    }
			}
		    }
		}
	    else {
		if (wx <= xvec[ipix1])
		    ix = ipix1
		else if (wx >= xvec[ipix2])
		    ix = ipix2
		else {
		    i1 = ipix1
		    i2 = ipix2 - 1
		    do i = i1, i2 {
			if (wx >= xvec[i] && wx < xvec[i+1]) {
			    if (wx - xvec[i] > 0.5)
				ix = i
			    else
				ix = i + 1
			    break
			    }
			}
		    }
		}
	    switch (key) {
		case 'r':
		    newgraph = TRUE

#	Make a hardcopy of the current display
		case '@':
		    call gflush (gfd)
		    call gclose (gfd)
		    call clgstr("plotter", plotter, SZ_FNAME)
		    call printf ("Making hard copy on %s\n")
			call pargstr (plotter)
		    call flush (STDOUT)
		    gfp = ERR
		    gfp = gopen (plotter, APPEND, STDPLOT)
		    If (gfp == ERR) {
			call printf ("PLOTSPEC:  cannot open plotter %s\n")
			    call pargstr (plotter)
			next
			}
		    call replot (gfp, gt, npix, Memr[smpix+ipix1-1], xvec[ipix1])
		    call gflush (gfp)
		    call gclose (gfp)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    call gseti (gfd, G_YNMINOR, NO)
		    call gseti (gfd, G_XNMINOR, YES)
		    call gsetr (gfd, G_TICKLABELSIZE, 0.8)

#	Toggle presence of heading
		case 'h':
		    if (heading) {
			call gt_sets (gt, GTTITLE, "")
			call gt_sets (gt, GTSUBTITLE, "")
			heading = FALSE
			}
		    else {
			call gt_sets (gt, GTTITLE, linebuf)
			call gt_sets (gt, GTSUBTITLE, gtitle)
			heading = TRUE
			}
		    newgraph = TRUE
 
#	Print wavelength coordinate at cursor
		case 'w':
		    call gclose (gfd)
		    call printf ("Wavelength %8.2f, Pixel %d: %8.2f\n")
		    call pargr (wx)
		    call pargi (ix)
		    call pargr (wy)
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Set number of times to smooth spectrum
		case 'g':
		    call gclose (gfd)
		    call printf("Number of times to smooth spectrum (%d): ")
			call pargi (ismooth)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,itemp)
			ismooth = itemp
			call printf("Smoothing spectrum %d times\n")
			    call pargi (ismooth)
			call flush (STDOUT)
			call amovr (pix,Memr[smpix],npoints)
			call smooth (Memr[smpix],npoints,ismooth)
			newgraph = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Delete portion of data
		case 'd':
		    if (dfirst) {
			d1 = ix
			call gclose (gfd)
			call printf("d again to delete data from pixel %d\n")
			    call pargi (d1)
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			dfirst = false
			}
		    else {
			dfirst = true
			d2 = ix
			call gclose (gfd)
			call printf("filling from pixel %d to pixel %d\n")
			    call pargi (d1)
			    call pargi (d2)
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			call fillpix (npoints, pix, d1, d2)
			call amovr (pix,Memr[smpix],npoints)
			call smooth (Memr[smpix],npoints,ismooth)
			newgraph = TRUE
			}

#	Reinitialize two-keystroke commands
		case 'x':
		    if (!zfirst) {
			call gclose (gfd)
			call printf("zoom cancelled\n")
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			zfirst = true
			}
		    if (!dfirst) {
			call gclose (gfd)
			call printf("delete cancelled\n")
			call flush (STDOUT)
			gfd = gopen (device, APPEND, STDGRAPH)
			dfirst = true
			}

#	Zoom in on section of graph
		case 'z':
		    if (zfirst) {
			ipx1 = ix
			zfirst = false
			call gclose (gfd)
			call printf("z again to zoom from %5.0f\n")
			    call pargr (xvec[ipix1])
			call flush (STDOUT)
			}
		    else {
			ipx2 = ix
			zfirst = true
			if (ipx2 < ipx1) {
			    ipix2 = ipx1
			    ipix1 = ipx2
			    }
			else {
			    ipix1 = ipx1
			    ipix2 = ipx2
			    }
			npix = ipix2 - ipix1 + 1
			call gclose (gfd)
			call printf("zooming from %.0f to %.0f (pix %d - %d)\n")
			if (xvec[ipix1] > xvec[ipix2]) {
			    call pargr (xvec[ipix2])
			    call pargr (xvec[ipix1])
			    }
			else {
			    call pargr (xvec[ipix1])
			    call pargr (xvec[ipix2])
			    }
			    call pargi (ipix1)
			    call pargi (ipix2)
			call flush (STDOUT)
			newgraph = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Unzoom
		case 'u':
		    ipix1 = 1
		    ipix2 = npoints
		    npix = npoints
		    call gclose (gfd)
		    call printf("unzooming back to %.0f to %.0f (pix %d - %d)\n")
		    if (xvec[ipix1] > xvec[ipix2]) {
			call pargr (xvec[ipix2])
			call pargr (xvec[ipix1])
			}
		    else {
			call pargr (xvec[ipix1])
			call pargr (xvec[ipix2])
			}
			call pargi (ipix1)
			call pargi (ipix2)
		    call flush (STDOUT)
		    newgraph = TRUE
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Print menu
		case '?',' ':
		    call gclose (gfd)
		    call printf ("@  make hard copy of screen\n")
		    call printf ("c  show x and y coordinates\n")
		    call printf ("d  delete data between 1st and 2nd positions\n")
		    call printf ("g  smooth plotted data\n")
		    call printf ("h  toggle headings\n")
		    call printf ("p  select correlation peak\n")
		    call printf ("q  leave plot\n")
		    call printf ("r  replot current graph\n")
		    call printf ("s  reset peak search half-width\n")
		    call printf ("u  unzoom\n")
		    call printf ("x  cancel zoom or delete after 1st char\n")
		    call printf ("z  zoom between 1st and 2nd positions\n")
		    call printf (">  overwrite spectrum with edited one\n")
		    call printf ("?  this menu (space works, too)\n")
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Set peak search width
		case 'S','s':
		    call gclose (gfd)
		    iwidth = clgeti ("pksrch")
		    call printf ("Half-width in pixels for peak search (%d) = \n")
			call pargi (iwidth)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 0) {
			ip = 1
			i = ctoi (linebuf,ip,iwidth)
			call clputi ("pksrch",iwidth)
			}
		    gfd = gopen (device, APPEND, STDGRAPH)
 
#	Set estimate of cross-correlation peak to the cursor value
		case 'P','p':
		    call gclose (gfd)
		    ix0 = ix
		    wx0 = wx
		    iwidth = 25
		    iwidth = clgeti ("pksrch")
		    i1 = ix - iwidth + 1
		    if (i1 < 1)
			i1 = 1
		    i2 = ix + iwidth
		    if (i2 > ipix2)
			i2 = ipix2
		    xmax = pix[ix]
		    for (i = i1; i <= i2; i = i + 1) {
			if (pix[i] >= xmax) {
			    xmax = pix[i]
			    ix = i
			    }
			}
		    z[1] = ix
		    wx = xvec[ix]
		    call printf ("Max peak at %d = %8.2f (%d = %8.2f)\n")
			call pargi (ix)
			call pargr (wx)
			call pargi (ix0)
			call pargr (wx0)
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)

		case 'Q','q':
		    break

#	Default is to print coordinates of cursor
		default:
		    call gclose (gfd)
		    call printf ("x,y: %10.3f %10.4g\n")
			call pargr (wx)
			call pargr (wy)
		    call flush (STDOUT)
		    gfd = gopen (device, APPEND, STDGRAPH)
		}
 
#	Replot graph if necessary
	    if (newgraph) {
		call replot (gfd, gt, npix, Memr[smpix+ipix1-1], xvec[ipix1])
		call gflush (gfd)
		newgraph = FALSE
		}
	    }
 
# Close graph window
 
	call gclose (gfd)
	call gt_free (gt)
	call sfree (sp)
	return
end
 
# June	1987	Gerard Kriss
# Feb 	1991	Doug Mink	Clean up code
# May 	1991	Doug Mink	Add zoom and hard-copy
# Sept 	1991	Doug Mink	Separate plot devices
# Dec 18 1991	Plot using wavelength vector instead of end-points
# Dec 19 1991	Add graphics pointer argument to gflush calls

# May 28 1992	Search for maximum near peak set with 'p'

# Feb  3 1994	Add s to change width of peak search
# Mar 17 1994	Add n to smooth spectrum
# Apr  8 1994	Fix n smoothing
# Apr 11 1994	Clean up pointers before returning
# Apr 19 1994	Fix bug freeing smoothed spectrum
# Apr 20 1994	Change GETLINE from subroutine to function
# Apr 26 1994	Fix smoothed spectrum with d command, too
# May  3 1994	Add smoothed vector argument to PLOTSPEC argument list
# May  9 1994	Move smoothed vector out of argument list
# May  9 1994	Add number of times to smooth spectrum to argument list
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 10 1994	Replace system ID with task name and system ID

# Jan 30 1997	Rename FIXIT to FILLPIX and move to file fill.x
# Feb  4 1997	Use NSMOOTH0 instead of NSMOOTH which is now in rvsao.com
# Feb  4 1997	Add option to plot without interacting
# Apr 21 1997	Flush but do not clear graphics buffer after plotting data
# Dec 11 1997	Add option to drop headings and smooth data

# Jan 12 1998	Match menu to actual commands

# Mar 28 2001	Handle right-to-left spectra as well as left-to-right ones
