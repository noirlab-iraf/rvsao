# File rvsao/Xcor/xcplot.x
# July 22, 2010
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1992-2010 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  XCPLOT uses IRAF plot routines to plot a spectrum and its cross-correlation.
#  The spectrum can also be plotted with labelled absorption lines.

###########################################################################

include <gset.h>
include <gio.h>
include <pkg/gtools.h>
include <ttyset.h>
include <fset.h>
include <imio.h>
include <smw.h>
include "rvsao.h"
include "emv.h"

procedure xcplot (wlspec, spectrum, npix, xvel, xcor, nvel, ncor2,
		  specfile, mspec, specim, newfit)

real	wlspec[ARB]	# Wavelengths for object spectrum
real	spectrum[ARB]	# Object spectrum
int	npix		# Number of pixels in object spectrum
real	xvel[ARB]	# Velocities for cross-correlation
real	xcor[ARB]	# Cross-correlation
int	nvel[ARB]	# Number of points in cross-correlation
int	ncor2		# Maximum number of points in cross-correlation
char	specfile[ARB]	# Spectrum file name
int	mspec		# Number of spectrum to read from multispec file
pointer	specim		# Spectrum image header structure
bool	newfit		# Rerun correlation if true

double	strtwav		# Wavelength of first correlation pixel
double	finwav		# Wavelength of last correlation pixel
char    command[SZ_FNAME]
double	dtemp
int	jtemp, i, j, ip, it0, it
char	linebuf[SZ_LINE]
int     wc, key
double	pix1,pix2
real	wx, wy
int	dispmod0
pointer gt, gfd, gfp
char    device[SZ_FNAME] # Display on which to plot data
char    plotter[SZ_FNAME] # Device on which to make hardcopy
bool    cursor          # true if waiting for cursor after plotting
bool    hardcopy        # true to make automatic hard copy and plot
bool	displot
int	dispmode
pointer	openplot()
int	ixi		# Index into xvel and xcor for a given template
int	ipix1, ipix2, ix, npix0
int	itemp		# Template for which to plot cross-correlation peak
double	wl
double	wav1, wav2	# Wavelength limits for spectrum
bool	zfirst
char	xlab[SZ_LINE+1]
int	strlen()

int	getline()
int	gt_gcur()
pointer	gt_init(), gopen()
bool	clgetb()
int	clgeti(), ctoi(), clscan()
double	ctod()
define	endplot_	90
int	rcompx()
extern	rcompx()
double	wcs_p2w()

include "rvsao.com"
include "results.com"
include "emv.com"
include "xplt.com"

begin

#  Set wavelength limts for line-labelled plot
	pix1 = 1.d0
	pix2 = double (npix)
	call wcs_set (specsh)
	wav1 = wcs_p2w (pix1)
	wav2 = wcs_p2w (pix2)

#  Sort templates by R-value
	itemp = itmax
	do i = 1, ntemp {
	    itr[i] = i
	    }
	call gqsort (itr, ntemp, rcompx, czr)
	it0 = 1

# Allocate memory for continuum-substracted spectrum
	if (npix > maxpix) {
	    if (scont != NULL) {
		call mfree (scont, TY_REAL)
		scont = NULL
		}
	    if (cspec != NULL) {
		call mfree (cspec, TY_REAL)
		cspec = NULL
		}
	    if (smspec != NULL) {
		call mfree (smspec, TY_REAL)
		smspec = NULL
		}
	    if (smcspec != NULL) {
		call mfree (smcspec, TY_REAL)
		smcspec = NULL
		}
	    maxpix = npix
	    call malloc (scont, npix, TY_REAL)
	    call malloc (cspec, npix, TY_REAL)
	    call malloc (smspec, npix, TY_REAL)
	    call malloc (smcspec, npix, TY_REAL)
	    }

# Subtract continuum and chop unwanted lines from spectrum using ICFIT
	call amovr (spectrum, Memr[cspec], npix)
        call icsubcon (npix,Memr[cspec],wlspec,specfile,1,nsmooth,Memr[scont])

# Smooth spectrum for plotting
	call amovr (spectrum,Memr[smspec],npix)
	call smooth (Memr[smspec],npix,nsmooth)
	call amovr (Memr[cspec], Memr[smcspec], npix)
	call smooth (Memr[smcspec],npix,nsmooth)
	strtwav = twl1[itemp]
	finwav = twl2[itemp]
#	strtwav = wav1
#	finwav = wav2

	if (clscan ("dispmode") == EOF)
	    dispmode = 2
	else
	    dispmode = clgeti ("dispmode")

	dispabs = TRUE
	zfirst = TRUE
	if (spnl > 0 || plotem)
	    dispem = TRUE
	else {
	    dispem = FALSE
	    do i = 1, ntemp {
		if (!tschop[i] && czr[i] > 4.0) {
		    dispem = TRUE
		    plotem = TRUE
		    }
		}
	    }
	ixi = (itemp-1) * ncor2 + 1

# Print hard copy, if requested

	hardcopy = clgetb ("hardcopy")
	if (hardcopy) {

	# Open plotter
	    call clgstr ("plotter",plotter,SZ_FNAME)
	    gfp = openplot (plotter)

	    if (dispmode == 2)
		call emeplot(gfp,specfile,specim,mspec,Memr[smspec],wlspec,npix)
	    else if (dispmode == 3)
		call emeplot(gfp,specfile,specim,mspec,Memr[smspec],wlspec,-npix)

	# Plot spectrum, cross-correlation, and template information
	    else {
		call xcsplot (gfp,npix,wlspec,Memr[smspec],strtwav,finwav,dispmode)
		call xciplot (gfp,specfile,specim,mspec,strtwav,finwav)
		call xcorplot (gfp,nvel[itemp],xvel[ixi],xcor[ixi],it0)
		}
	    if (debug) {
		call printf ("XCPLOT: Making hardcopy on %s\n")
		    call pargstr (plotter)
		}

	# Close plotter
	    call closeplot (gfp)
	    }

#  If not displaying data, return
	displot = clgetb ("displot")
	if (!displot) {
	    return
	    }

# Plot spectrum with absorption and emission lines labelled
	ipix1 = 1
	ipix2 = npix
	if (dispmode > 1) {
	    xfit = FALSE
	    call emplot (specfile,specim,mspec,npix,spectrum,
			 Memr[smspec],Memr[cspec],Memr[smcspec],
			 wlspec,ipix1,ipix2)
	    if (!plotcorr) {
		if (xfit) {
		    call printf ("Rerun cross-correlation\n")
		    newfit = TRUE
		    }
		return
		}
	    call printf ("Plotting correlation with summary\n")
	    }

#  Open display
	call clgstr ("device",device,SZ_FNAME)
	gfd = openplot (device)

#  Plot spectrum, cross-correlation and template information
	call xcsplot (gfd,npix,wlspec,Memr[smspec],strtwav,finwav,dispmode)
	call xciplot (gfd,specfile,specim,mspec,strtwav,finwav)
	call xcorplot (gfd,nvel[itemp],xvel[ixi],xcor[ixi],it0)
	call gflush (gfd)

#  If not using cursor, return
	cursor = clgetb ("curmode")
	if (!cursor) {
	    call closeplot (gfd)
	    return
	    }

#  Accept cursor commands
	newfit = FALSE
	gt = gt_init()
	npix0 = npix
	while (gt_gcur ("cursor", wx, wy, wc, key, command, SZ_FNAME) != EOF) {
	    wl = wx
	    ix = ipix1
	    do i = ipix1, ipix2 {
		if (wl >= wlspec[i] && wx < wlspec[i+1]) {
		    if (wl - wlspec[i] > 0.5)
			ix = i
		    else
			ix = i + 1
		    break
		    }
		}

	    switch (key) {

	    # Toggle debug flag
		case 'd':
		    call gclose (gfd)
		    if (debug) {
			call printf ("debugging off\n")
			debug = FALSE
			}
		    else {
			debug = TRUE
			call printf ("debugging on\n")
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

	    # Display spectrum so that it can be edited
		case 'e':
		    call gclose (gfd)
		    call printf ("edit spectrum\n")
		    call strcpy (LABEL(specsh),xlab,SZ_LINE)
		    if (strlen (UNITS(specsh)) > 0) {
			call strcat (" in ",xlab,SZ_LINE)
			call strcat (UNITS(specsh),xlab,SZ_LINE)
			}
		    call plotspec (npix,spectrum,specname,wlspec,xlab,nsmooth)
		    gfd = gopen (device, APPEND, STDGRAPH)
		    newfit = TRUE
		    break

#	Set number of times to smooth spectrum
		case 'g':
		    call gclose (gfd)
		    call printf("Number of times to smooth spectrum (%d): ")
			call pargi (nsmooth)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,jtemp)
			nsmooth = jtemp
			call amovr (spectrum,Memr[smspec],npix)
			call smooth (Memr[smspec],npix,nsmooth)
			call amovr (Memr[cspec],Memr[smcspec],npix)
			call smooth (Memr[smcspec],npix,nsmooth)
			}
		    gfd = gopen (device, APPEND, STDGRAPH)

#	Change display mode to specific one
		case 'm':
		    call gclose (gfd)
		    dispmod0 = dispmode
		    call printf("Display mode (0-3) (%d): ")
			call pargi (dispmode)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 0) {
			ip = 1
			j = ctoi (linebuf,ip,jtemp)
			if (jtemp > -4 && jtemp < 4)
			    dispmode = jtemp
			}
		    call printf ("Display mode %d -> %d (%d)\n")
			call pargi (dispmod0)
			call pargi (dispmode)
			call pargi (jtemp)
		    call flush (STDOUT)
		    if (dispmode != dispmod0) {
			if (dispmode > 1) {
			    xfit = FALSE
			    call emplot (specfile,specim,mspec,npix,spectrum,
					 Memr[smspec],Memr[cspec],Memr[smcspec],
					 wlspec,ipix1,ipix2)
			    if (!plotcorr) {
				if (xfit) {
				    newfit = TRUE
				    call printf ("Rerun cross-correlation\n")
				    }
				call gt_free(gt)
				return
				}
			    call printf ("Plotting correlation with summary\n")
			    }
			gfd = openplot (device)
			call xcsplot (gfd,npix,wlspec,Memr[smspec],
				      strtwav,finwav,dispmode)
			call xciplot (gfd,specfile,specim,mspec,strtwav,finwav)
			call xcorplot (gfd,nvel[itemp],xvel[ixi],xcor[ixi],it0)
			call gflush (gfd)
			}
		    else
			gfd = gopen (device, APPEND, STDGRAPH)

#	Select different velocity for final redshift

	    # Repeat correlation plot to select peak using cursor
		case 'p':
		    call gclose (gfd)
		    call printf ("select correlation peak(s)\n")
		    pltcor  = TRUE
		    newfit = TRUE
		    gfd = gopen (device, APPEND, STDGRAPH)
		    break

	    # Plot spectrum with absorption lines marked
		case 'l':
		    call gclose (gfd)
		    call printf ("Plotting absorption lines with summary\n")
#		    do i = 1, npix {
#			call printf ("%4d: %.3fA %9.2f -> %9.2f\n")
#			    call pargi (i)
#			    call pargr (wlspec[i])
#			    call pargr (Memr[smspec+i-1])
#			    call pargr (spectrum[i])
#			}
#		    call printf ("XCPLOT: %d pixels from %.4fA to %.4fA\n")
#			call pargi (npix)
#			call pargd (wav1)
#			call pargd (wav2)

		    xfit = FALSE
		    call emplot (specfile,specim,mspec,npix,spectrum,
				 Memr[smspec],Memr[cspec],Memr[smcspec],
				 wlspec,ipix1,ipix2)
		    if (!plotcorr) {
			if (xfit) {
			    newfit = TRUE
			    call printf ("Rerun cross-correlation\n")
			    }
			call gt_free(gt)
			return
			}
		    call printf ("Plotting correlation with summary\n")
		    gfd = openplot (device)
		    call xcsplot (gfd,npix,wlspec,Memr[smspec],strtwav,finwav,dispmode)
		    ixi = (itemp-1) * ncor2 + 1
		    call xciplot (gfd,specfile,specim,mspec,strtwav,finwav)
		    call xcorplot (gfd,nvel[itemp],xvel[ixi],xcor[ixi],it0)
	      	    call gflush (gfd)

	    # Change peak fitting parameters
		case 'c':
		    call gclose (gfd)
	    	    call printf ("Correlation peak fitting parameters\n")
		    call printf ("Peak-fitting mode (1-parabola 2-quartic 3-cosx/1+x^2) (%d) = ")
			call pargi (pkmode0)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,jtemp)
			pkmode0 = jtemp
			newfit = TRUE
			}

		    call printf("Fraction of peak or number of points for peak fitting (%7.3f) = ")
			call pargd (pkfrac)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctod (linebuf,ip,dtemp)
			pkfrac = dtemp
			newfit = TRUE
			}
		    call printf ("Filter parameters for transform (%d %d %d %d) = ")
			call pargi (lo)
			call pargi (toplo)
			call pargi (nrun)
			call pargi (topnrn)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctoi (linebuf,ip,jtemp)
			lo = jtemp
			if (i > 0) {
			    i = ctoi (linebuf,ip,jtemp)
			    toplo = jtemp
			    }
			if (i > 0) {
			    i = ctoi (linebuf,ip,jtemp)
			    nrun = jtemp
			    }
			if (i > 0) {
			    i = ctoi (linebuf,ip,jtemp)
			    topnrn = jtemp
			    }
			newfit = TRUE
			}
		    gfd = gopen (device, APPEND, STDGRAPH)
		    if (newfit) break

		case 'f':
		    call gclose (gfd)
		    call printf ("Rerunning correlations\n")
		    newfit = TRUE
		    gfd = gopen (device, APPEND, STDGRAPH)
		    break

	    # Change template to use
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
		    ixi = (itemp-1) * ncor2 + 1
		    it0 = it
		    spxvel = zvel[itemp]
		    spxerr = czerr[itemp]
		    spxr = czr[itemp]
		    strtwav = twl1[itemp]
		    finwav = twl2[itemp]
#		    strtwav = wav1
#		    finwav = wav2
		    gfd = openplot (device)
		    call xcsplot (gfd,npix,wlspec,Memr[smspec],strtwav,finwav,dispmode)
		    call xciplot (gfd,specfile,specim,mspec,strtwav,finwav)
		    call xcorplot (gfd,nvel[itemp],xvel[ixi],xcor[ixi],it0)
	      	    call gflush (gfd)
		    IM_UPDATE(specim) = YES
		    }

	    # Change velocity range over which to search for peak
		case 'v':
		    call gclose (gfd)
		    call printf ("Minimum allowable velocity in km/sec (%8.1f) = ")
			call pargd (minvel)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctod (linebuf,ip,dtemp)
			minvel = dtemp
			newfit = TRUE
			}
		    call printf ("Maximum allowable velocity in km/sec (%8.1f) = ")
			call pargd (maxvel)
		    call flush (STDOUT)
		    i = getline (STDIN,linebuf)
		    if (i > 1) {
			ip = 1
			i = ctod (linebuf,ip,dtemp)
			maxvel = dtemp
			newfit = TRUE
			}
		    newfit = TRUE
		    gfd = gopen (device, APPEND, STDGRAPH)
		    break

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

	    # Make a hard copy of the correlation summary display
		case '@':
		    call gclose (gfd)
		    call clgstr ("plotter",plotter,SZ_FNAME)
		    call printf ("Making hard copy on %s\n")
			call pargstr (plotter)
		    gfp = openplot (plotter)
		    call xcsplot (gfp,npix,wlspec,Memr[smspec],strtwav,finwav,dispmode)
		    call xcorplot (gfp,nvel[itemp],xvel[ixi],xcor[ixi],it0)
		    call xciplot (gfp,specfile,specim,mspec,strtwav,finwav)
	      	    call closeplot (gfp)
		    gfd = gopen (device, APPEND, STDGRAPH)

	    # Replot correlation summary
		case 'r':
		    call gclose (gfd)
		    call printf ("Replot XCSAO summary\n")
		    gfd = openplot (device)
		    call xcsplot (gfd,npix,wlspec,Memr[smspec],strtwav,finwav,dispmode)
		    call xciplot (gfd,specfile,specim,mspec,strtwav,finwav)
		    call xcorplot (gfd,nvel[itemp],xvel[ixi],xcor[ixi],it0)
	      	    call gflush (gfd)

	    # Exit without updating the image header
		case 'x':
		    qplot = FALSE
		    break

	    # List available cursor commands
		case ' ','h','?':
		    call gclose (gfd)
		    call printf("XCSAO commands:\n")
		    call printf(" e  edit spectrum\n")
		    call printf(" c  change correlation peak fit parameters\n")
		    call printf(" d  toggle debug flag\n")
		    call printf(" f  rerun cross-correlation\n")
		    call printf(" j  conditional velocity\n")
		    call printf(" l  plot spectrum with absorption lines\n")
		    call printf(" m  change display mode to -1, 0, 1, or 2\n")
		    call printf(" n  disapprove velocity\n")
		    call printf(" p  refit correlation peak interactively\n")
		    call printf(" q  exit from this spectrum\n")
		    call printf(" r  replot summary\n")
		    call printf(" t  change template to use for xcor velocity\n")
		    call printf(" v  change velocity limits for result\n")
		    call printf(" x  Exit without updating header\n")
		    call printf(" y  approve velocity\n")
		    call printf(" @  make hard copy\n")
		    gfd = gopen (device, APPEND, STDGRAPH)

		default:
		    break
		}
	    }

#  Close graphics window
endplot_
	call gt_free(gt)
	call closeplot (gfd)
	return
end


int procedure rcompx (rvalue, itemp1, itemp2)

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

# Jun  2 1992	Fix peak selection mode bug--reopen gfd
# Nov 30 1992	Smooth spectrum here so it can be resmoothed by EMPLOT

# Jan 14 1993	Fix misspellings
# Dec  3 1993	Pass mspec to EMPLOT

# Mar 16 1994	Add option of editing spectrum
# Apr 15 1994	Add wavelength limits to EMPLOT call
# Apr 22 1994	Drop image structure from EMPLOT argument list
# May  3 1994	Add smoothed vector argument to PLOTSPEC call
# May  5 1994	Add correct smoothed vector argument to SPECPLOT call
# May  9 1994	Drop smoothed vector argument from SPECPLOT call
# May  9 1994	Add number of times to smooth spectrum to PLOTSPEC call
# May 23 1994	Locally allocate and free smoothed spectrum
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug  5 1994	Drop nvel as argument to INFOPLOT
# Aug 10 1994	Plot spectrum over wavelength overlap with best template
# Aug 17 1994	Change names of SPECPLOT and INFOPLOT to XCSPLOT and XCIPLOT
# Dec  7 1994	Add multiple display modes; add x command

# Feb 27 1995	Fix hard copy command
# Mar 22 1995	Add options to change fitting parameters
# Aug 18 1995	Set MINVEL and MAXVEL in common, not parameter file
# Aug 25 1995	Add qplot velocity accepatance flags; remove unused variables
# Aug 25 1995	Turn on emission line marking if R>4 for emission template
# Sep 19 1995	Add g to change spectrum smoothing for display
# Sep 25 1995	Use EMPLOT for interaction, not EMEPLOT

# Jan 10 1996	Pass DISPMODE to XCSPLOT to allow more options

# Feb  4 1997	Fix label on spectrum being edited
# Mar 14 1997	Get dispersion axis label from spectrum header
# Apr  7 1997	Pass all correlation vectors so any can be plotted
# Apr  7 1997	Sort template indices by R-value
# Apr  8 1997	Print template number as part of correlation header
# Apr 14 1997	Pass IT0 to XCORPLOT instead of title
# Nov 14 1997	Add arguments to emplot to fake continuum-subtracted spectra

# Jan 13 1998	Add command to change to specific display mode
# Feb 13 1998	Pass correlations and velocities as pointers
# Apr 22 1998	Make boolean assignments consistent (suggested by Bryan Miller)

# Jun 22 1999	Add debug statement showing plotter name
# Jul 29 1999	Recompute continuum-subtraction to pass to emplot
# Aug 19 1999	Set vplot in t_xcsao from vel_plot parameter
# Sep 24 1999	Set IM_UPDATE to yes if quality flag is set

# Jul 31 2001	Update spectrum header if xcor template changed

# May 12 2004	If plotem is set by eminit(), always label emission lines

# Jan 24 2007	Reset WCS before each plot subroutine call
# Jan 30 2007	Always plot all of spectrum
# Jan 31 2007	Send only limiting pixels to emplot(), not wavelength limits
# Apr 27 2007	Send limiting wavelengths for template to xcsplot and xciplot
# Aug 13 2007	Allocate all buffers only when needed
# Oct 15 2007	Change tschop[i] == FALSE test to !tschop[i]

# Jul 22 2010	Check individual pointers before freeing them
