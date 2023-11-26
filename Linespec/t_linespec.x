# File rvsao/Makespec/t_linespec.x
# August 25, 2004
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 2004 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# LINESPEC is an IRAF task for creating a spectrum from a list of emission 
# and/or absorption lines.  Arguments are in the parameter file "linespec.par".
# Information is shared in common blocks defined in "rvsao.com" and "sum.com".
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include "rvsao.h"

define	MAXLINES	100

procedure t_linespec ()

int	i
char	specfile[SZ_PATHNAME]	# Template spectrum file name
char	specdir[SZ_PATHNAME]	# Template spectrum directory
char	specpath[SZ_PATHNAME]	# Template spectrum pathname
char	linefile[SZ_PATHNAME]	# Line list file name
char	dirpath[SZ_PATHNAME]	# Line list directory
char	linepath[SZ_PATHNAME]	# Line list pathname
int	logfiles	# List of log files
char	logfile[SZ_PATHNAME] # Log file name
char	wtitle[20]	# Title for wavelength plots of spectrum

pointer	outspec		# Created spectrum
pointer	specwl		# Wavelength vector for output spectrum
char	str[SZ_LINE]
int	fd
int	ldir
char	lbracket[3]	# "[({"
char	rbracket[3]	# "])}"
double	wl[MAXLINES]
double	ht[MAXLINES]
double	wd[MAXLINES]
char	nm[16,MAXLINES]
double	wdline
bool	maxwidth	# yes to use maximum of line and instrument widths
double	z1		# 1 + z
real	continuum	# Continuum level for final spectrum
 
pointer specim		# Output spectrum image header structure
bool	specint
double	center		# Center of emission line in Angstroms
double	width		# Width of emission line in Angstroms
double	zspec
double	minwav, maxwav, dwav
double	dindef, pi, s2pi, area, eqw
int	npix, lpath, iline, nlines
char	dstr[SZ_LINE]	# Date string
char	hstr[SZ_LINE]	# New history string
char	keyword[SZ_HKWORD]
char	line[SZ_HSTRING]
bool	verbose

bool	clgetb()
int	clpopnu(), clgfil(), open()
int	clgeti()
double	clgetd()
real	clgetr()
int	stridx()
int	strlen()
int	fscan()

include	"rvsao.com"
include	"sum.com"
 
begin
	dindef = INDEFD
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	call sprintf (wtitle,20,"Wavelength")
	c0 = 299792.5d0
	pi = 3.14159265358979323846d0
	s2pi = sqrt (2.d0 * pi)

# Get task parameters.

# Print processing information
	debug  = clgetb ("debug")

# Emission line list for created spectrum
	call clgstr ("linedir",dirpath,SZ_PATHNAME)
	lpath = strlen (dirpath)
	if (lpath > 0 && dirpath[lpath] != '/') {
	    dirpath[lpath+1] = '/'
	    dirpath[lpath+2] = EOS
	    }
	call clgstr ("linefile",linefile,SZ_PATHNAME)
	if (stridx ("/",linefile) > 0 || stridx ("$",linefile) > 0)
	    call strcpy (linefile,linepath,SZ_PATHNAME)
	else {
	    call strcpy (dirpath,linepath,SZ_PATHNAME)
	    call strcat (linefile,linepath,SZ_PATHNAME)
	    }

# Open emission line file; drop out if it cannot be read
	fd = open (linepath, READ_ONLY, TEXT_FILE)
	nlines = 0
	if (fd <= 0) {
	    call printf ("LINESPEC: Cannot read line list in %s\n")
		call pargstr (linepath)
	    return
	    }

#  Read emission line file
	while (fscan(fd) != EOF && nlines < MAXLINES) {
	    nlines = nlines + 1
	    call gargd (wl[nlines])
	    call gargd (wd[nlines])
	    call gargd (ht[nlines])
	    call gargwrd (nm[1,nlines],16)
	    if (debug) {
		call printf ("LINESPEC: %2d: %5s %7.2f +-%6.2f %8.4f\n")
		    call pargi (nlines)
		    call pargstr (nm[1,nlines])
		    call pargd (wl[nlines])
		    call pargd (wd[nlines])
		    call pargd (ht[nlines])
		}
	    }
	call close (fd)

# Emission line half-widths in Angstroms
	wdline = clgetd ("linewidth")

# How to impose instrument resolution
	maxwidth = clgetb ("maxwidth")

# Redshift velocity (cz) for emission lines
	velocity = clgetd ("velspec")
	zspec = clgetd ("zspec")
	if (zspec != 0.d0 && zspec != dindef)
	    velocity = c0 * zspec
	else if (velocity == 0.d0 || velocity == dindef)
	    velocity = 0.d0
	z1 = 1.d0 + (velocity / c0)

# Template file to which to write created spectrum
	call clgstr ("specfile",specfile,SZ_PATHNAME)
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (ldir > 0) {
	    if (specdir[ldir] != '/') {
		specdir[ldir+1] = '/'
		specdir[ldir+2] = EOS
		}
	    call strcpy (specdir,specpath,SZ_PATHNAME)
	    call strcat (specfile,specpath,SZ_PATHNAME)
	    }
	else
	    call strcpy (specfile,specpath,SZ_PATHNAME)
 
# Wavelength limits for output spectrum
        minwav = clgetd ("st_lambda")
        maxwav  = clgetd ("end_lambda")
	dwav = clgetd ("pix_lambda")
	npix = clgeti ("npts")

# If start, increment, and number of pixels are set, compute ending wavelength
	if (npix > 2 && minwav != dindef && dwav != dindef)
	    maxwav = minwav + (dwav * double (npix - 1))

# If start, increment, and last wavelength are set, compute number of pixels
	else if (npix < 3 && minwav != dindef && dwav != dindef && maxwav != dindef)
	    npix = 1 + idnint (((maxwav - minwav) / dwav) + 0.5d0)

# If increment, last, and number of pixels are set, compute starting wavelength
	else if (npix > 3 && maxwav != dindef && dwav != dindef)
	    minwav = maxwav - (dwav * double (npix - 1))

# If first, last, and number of pixels are set, compute wavelength increment
	else if (npix > 3 && minwav != dindef && maxwav != dindef)
	    dwav = (maxwav - minwav) / double (npix - 1)

# Plot and interact with output spectrum
	pltspec = clgetb ("spec_plot")
	specint = clgetb ("spec_int")
	if (specint)
	    tsmooth = 0
	else
	    tsmooth = -1
 
# Open log files and write a header.
	verbose = clgetb ("verbose")
	if (verbose) {
	    logfiles = clpopnu ("logfiles")
	    call fseti (STDOUT, F_FLUSHNL, YES)
	    i = 0
	    call strcpy ("rvsao.linespec",taskname,SZ_LINE)
	    while (clgfil (logfiles, logfile, SZ_PATHNAME) != EOF) {
		fd = open (logfile, APPEND, TEXT_FILE)
		if (fd == ERR) break
		call loghead (taskname,str)
		call fprintf (fd, "%s\n")
		    call pargstr (str)
		i = i + 1
		logfd[i] = fd
		}
	    nlogfd = i
	    call clpcls (logfiles)

	    call sprintf (line,SZ_LINE,"%.3fA - %.3fA by %.3fA: %d points in spectrum")
		call pargd (minwav)
		call pargd (maxwav)
		call pargd (dwav)
		call pargi (npix)
	    do i = 1, nlogfd {
		call fprintf (logfd[i], "%s\n")
		    call pargstr (line)
		}

	    if (velocity != 0.d0) {
		call sprintf (line,SZ_LINE,"Spread by %.2fA, shifted to z = %.4f = %.2f km/sec")
		    call pargd (wdline)
		    call pargd (zspec)
		    call pargd (velocity)
		}
	    else {
		call sprintf (line,SZ_LINE,"Spread by %.2fA\n")
		    call pargd (wdline)
		}
	    do i = 1, nlogfd {
		call fprintf (logfd[i], "%s\n")
		    call pargstr (line)
		}
	    }

#  Open the output file
	call tmp_open (specim,specpath,outspec,npix, 1, NULL)
	if (specim == NULL) {
	    call printf ("LINESPEC: Cannot write spectrum file %s\n")
		call pargstr (specpath)
	    return
	    }

#  Write dispersion world coordinate system information to spectrum header    
	call clgstr ("specname",specname,SZ_LINE)
	call strcpy (specname, IM_TITLE(specim),SZ_IMTITLE)
	call imaddd (specim, "EXPTIME", 1.d0)
	call imaddi (specim, "DISPAXIS", 1)
	call imaddi (specim, "DC-FLAG", 0)
	call imaddr (specim, "CRPIX1", 1.)
	call imaddd (specim, "CRVAL1", minwav)
	call imaddd (specim, "CDELT1", dwav)
	call imaddd (specim, "W0", minwav)
	call imaddd (specim, "WPC", dwav)
	call imaddi (specim, "NP1", 1)
	call imaddi (specim, "NP2", npix)
	call imaddd (specim, "VELOCITY", velocity)
	call imaddd (specim, "BCV", 0.d0)
	call imaddb (specim, "CHOPEM", FALSE)
	call imaddb (specim, "EMCHOP", FALSE)
	call imaddb (specim, "SUBCONT", FALSE)
	call imaddb (specim, "OVERLAP", TRUE)
	call imaddi (specim, "FI-FLAG", 3)
	call imaddi (specim, "DISPAXIS", 1)

# Intitialize spectrum to all zeroes
	call aclrr (Memr[outspec], npix)

# Add lines to spectrum, one at a time
	do iline = 1, nlines {

	# If velocity is nonzero, shift line appropriately
	    center = wl[iline] * z1

	# If line width is positive, it is in Angstroms (physics)
	    if (wd[iline] > 0)
		width = wd[iline]

	# If line width is negative, it is in km/sec (dispersion)
	    else if (wd[iline] < 0)
		width = center * (-wd[iline] / c0)

	# If line width is zero, set it to +- one pixel
	    else
		width = dwav

	# Area in line
	    area = s2pi * width * ht[iline]

	# Equivalent width; area if no continuum
	    if (continuum > 0.0)
		eqw = area / continuum
	    else
		eqw = area

	# If using maximum of instrument and line widths, pick one
	    if (maxwidth && wdline > width)
		width = wdline

	# Add this line to the spectrum
	    call mkgaus (Memr[outspec], center, width, ht[iline],
			 minwav, dwav, npix, debug)

	# Write the position of this line into the header
	    call sprintf (keyword,SZ_HKWORD,"EMLINE%d")
		call pargi (iline)
	    call sprintf(line,67,"%5s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.3f")
		call pargstr (nm[1,iline])
		call pargd (wl[iline])
		call pargd (center)
		call pargd (ht[iline])
		call pargd (width)
		call pargd (velocity)
		call pargd (0.d0)
		call pargd (eqw)
		call pargd (1.d0)
	    call imastr (specim,keyword,line)
	    if (verbose) {
		do i = 1, nlogfd {
		    call fprintf (logfd[i], "%s\n")
			call pargstr (line)
		    }
		}
	    }

# Add continuum level to spectrum if it is non-zero
	continuum = clgetr ("continuum")
	if (continuum != 0.0)
	    call aaddkr (Memr[outspec], continuum, Memr[outspec], npix)

# Plot new spectrum, if requested
	if (pltspec) {
	    call malloc (specwl, npix, TY_REAL)
	    do i = 1, npix {
		Memr[specwl+i-1] = real (minwav + (dwav * double (i - 1)))
		}
	    call plotspec (npix,Memr[outspec],specfile,
			   Memr[specwl],wtitle,tsmooth)
	    }

# Convolve new spectrum with instrument resolution
	if (!maxwidth && wdline > 0.d0) {
	    call wtgaus (Memr[outspec],wdline,minwav,dwav,npix, debug)

# Plot revised spectrum, if requested
	    if (pltspec) {
		call plotspec (npix,Memr[outspec],specfile,
				Memr[specwl],wtitle,tsmooth)
		call mfree (specwl, TY_REAL)
		}
	    }

# Write version information to the header
	call logtime (dstr,SZ_LINE)
	call sprintf (hstr,SZ_LINE,"rvsao.linespec %s %s, %d lines")
                call pargstr (VERSION)
                call pargstr (dstr)
                call pargi (nlines)
	call imputh (specim, "HISTORY", hstr)

# Close the log files
	if (verbose) {
	    do i = 1, nlogfd {
		call close (logfd[i])
		}
	    }                                              

#  Close the output spectrum file
	call tmp_close (specim,outspec,debug)

	return

end

# Mar 31 1997	New task
# Apr 14 1997	Add XCSAO flags to avoid eliminating lines
# Apr 22 1997   Change parameter template to tempfile
# Apr 25 1997	Add tempobj parameter for template title
# Apr 25 1997	Add continuum from parameter continuum
# Apr 29 1997	Change name from LINTEMP to LINESPEC and variable names, too
# May  2 1997	Change specobj parameter to specname
# May  2 1997	Always test against dindef, not INDEFD
# May  8 1997	Fix directory parameter name bug
# May 19 1997	Write all line parameters to the image header, including shift
# Jun 18 1997	Use local double ZSPEC instead of real Z0 in common

# Jun 15 1999	Add nspec argument to tmp_open()
# Jul 27 1999	Add NULL argument to tmp_open()

# Mar  6 2000	Fix major bug which caused improper velocity if Z input
# Jul  6 2000	Add option to specify 3 of 4 dispersion parameters
# Nov  8 2000	Fix bug so dispersion is computed if wave/pix is not specified

# Aug 25 2004	Fix declaration of nm[]
