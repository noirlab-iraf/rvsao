# File rvsao/Emsao/t_emsao.x
# May 2, 2014
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1991-2014 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# EMSAO is an IRAF task for obtaining redshifts using the locations of
# known emission lines.  Arguments are in the IRAF parameter file emsao.par.
# Information is shared in common blocks defined in "rvsao.com" and "emv.com".

include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include <smw.h>
include "rvsao.h"
include "contin.h"
include "emv.h"

procedure t_emsao ()

char	specfile[SZ_PATHNAME]	# Input spectrum file name
char	skynums[SZ_LINE]	# List of multispec sky spectra to read
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	specdir[SZ_PATHNAME]	# Directory for object spectra

char	vel_corr[SZ_LINE]	# type of velocity correction for spectrum
				# none | file | heliocentric | barycentric
char	vel_init[SZ_LINE]	# type of velocity for initial guess
				# correlation|emission|combination|search|guess
char	vel_plot[SZ_LINE]	# type of velocity for redshifting plot
				# correlation|emission|combination|search
char	logfile[SZ_PATHNAME]	# Log file name
pointer	speclist	# List of spectrum files
int	logfiles	# IRAF list pointer for log files
int	rmode		# report format (1=normal,2=one-line)
double 	minwav0		# Starting wavelength from parameter file
double 	maxwav0		# Ending wavelength from parameter file
double	minwav		# Actual wavelength of first pixel used
double	maxwav		# Actual wavelength of last pixel used
bool	arcflag		# archive record writing flag
int	mspec		# Aperture of spectrum to read from multispec file
int	mband		# Band of spectrum to read from multispec file
int	sspec		# Aperture of sky spectrum to read from multispec file
int	sband		# Band of sky spectrum to read from multispec file
pointer spectrum	# Object spectrum
pointer specim		# Object image header structure
pointer specsky		# Sky spectrum
pointer skyim		# Sky image header structure
pointer skysh		# Sky spectrum header structure
pointer	pxspec		# Pixel vector for object spectrum
pointer wlspec		# Wavelength vector for object spectrum
pointer	pxsky		# Pixel vector for sky spectrum
int	nmspec0		# Number of object multispec spectra
int	nspec		# Number of object multispec spectra
int	mspec0		# Number of spectrum in file to read
int	nsspec		# Number of sky multispec spectra
int	mspec_range[3,MAX_RANGES]
int	sspec_range[3,MAX_RANGES]
int	pix1,pix2	# Pixel range for revised wavelength range
int	ip,jp,lfile	# Limits for multispec aperture decoding
char	lbracket[3]	# "[({"
char	rbracket[3]	# "])}"
int	npix, ipix
double	dindef
bool	linefit
bool	wenable

char	str[SZ_LINE]
int	i,j,iline,ldir
int	fd
int	nsp, npix0
double	sptot, spmax, spmin, speci
pointer	ispec
 
pointer	imtopenp()
int	imaccess(), imtgetim()
bool	clgetb()
int	clgeti()
int	clpopnu(), clgfil(), open(), clscan()
int	strdic(),strncmp()
double	clgetd()
int	decode_ranges(),get_next_number(),stridx(),strlen(),stridxs()
double	wcs_w2p(), wcs_p2w()

define	endspec_ 90
define	nextspec_ 80

include	"rvsao.com"
include	"emv.com"
include "contin.com"
include "results.com"
 
begin
	c0 = 299792.5d0
	npix0 = 0
	dindef = INDEFD
	schop = FALSE
	wlspec = NULL
	specsh = NULL
	skysh = NULL
	specref = 0
	waverest = 0.d0
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")

# Initialize emission line information
	do iline = 1, MAXREF {
	    do j = 1, 2 {
		do i = 1, 10 {
		    emparams[i,j,iline] = 0.d0
		    }
		}
	    wlrest[iline] = 0.d0
	    pxobs[iline] = 0.d0
	    wlobs[iline] = 0.d0
	    wtobs[iline] = 0.d0
	    }

# Get task parameters.

# Get list of images to process
	speclist = imtopenp ("spectra")

# Multspec spectrum numbers
	call clgstr ("specnum",specnums,SZ_LINE)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR) {
	    call sprintf (str, SZ_LINE, "T_EMSAO: Illegal multispec object list <%s>")
		call pargstr (specnums)
	    call error (1, str)
	    }

# Multispec band number
	mband = clgeti ("specband")

# Multispec sky spectrum numbers
	call clgstr ("skynum",skynums,SZ_LINE)
	if (skynums[1] == EOS)
	    call strcpy ("0",skynums,SZ_LINE)
	skyspec = FALSE
	if (strncmp (skynums,"0",1) != 0) {
	    if (decode_ranges (skynums, sspec_range, MAX_RANGES, nsspec) == ERR) {
		call sprintf (str, SZ_LINE, "T_EMSAO: Illegal multispec sky list <%s>")
		    call pargstr (specnums)
		call error (1, str)
		}
	    else
		skyspec = TRUE
	    }

# Multispec band number for sky
	sband = clgeti ("skyband")

# Spectrum directory 
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }
 
# Wavelength limits if not from spectrum
	minwav0 = clgetd ("st_lambda")
	maxwav0 = clgetd ("end_lambda")

# Data renormalization
	renorm = clgetb ("renormalize")

	pltspec = clgetb("obj_plot")
	if (clscan ("contsub_plot") != EOF)
	    pltcon  = clgetb ("contsub_plot")
	else
	    pltcon = FALSE
	debug  = clgetb("debug")

# Type of heliocentric velocity correction to make
	call clgstr ("vel_corr",vel_corr,SZ_LINE)
	svcor = strdic (vel_corr,vel_corr,SZ_LINE,HC_VTYPES)

# Type of velocity for plotting emission and absorption lines
	call clgstr ("vel_plot",vel_plot,SZ_LINE)
	vplot = strdic (vel_plot,vel_plot,SZ_LINE,PL_VTYPES)

# Type of velocity for initial redshift for search
	call clgstr ("vel_init",vel_init,SZ_LINE)
	vinit = strdic (vel_init,vel_init,SZ_LINE,IEM_VTYPES)
	call clgstr ("cortemp",cortemp,SZ_FNAME)

# Report mode for log file
	rmode = 1
	rmode = clgeti ("report_mode")

# SAO TDC archive record writing flag
        arcflag = FALSE
        arcflag = clgetb ("archive")

# IRAF image header result writing flag
	savevel = FALSE
	savevel = clgetb("save_vel")
	wenable = savevel
	linefit = TRUE
	linefit = clgetb ("linefit")
        if (!linefit) {
	    if (savevel) {
		savevel = FALSE
		wenable = TRUE
		}
	    else
		wenable = FALSE
	    }


# Emission line on data plot flag
	dispem = false
	dispem = clgetb ("dispem")

# Absorption lines on data plot flag
	dispabs = false
	dispabs = clgetb ("dispabs")

# Open log files

	logfiles = clpopnu("logfiles")
	call fseti (STDOUT, F_FLUSHNL, YES)
	nlogfd = 0
	while (clgfil (logfiles, logfile, SZ_PATHNAME) != EOF) {
	    fd = open (logfile, APPEND, TEXT_FILE)
	    if (fd == ERR) next
	    nlogfd = nlogfd + 1
	    logfd[nlogfd] = fd
	    }
	call clpcls (logfiles)

# Write a header for multiline report.

	call strcpy ("rvsao.emsao",taskname,SZ_LINE)
	if (rmode == 1) {
	    do i = 1, nlogfd {
		call loghead (taskname,str)
		fd = logfd[i]
		call fprintf (fd, "%s\n")
		    call pargstr (str)
		}
	    }

# Initialize wavelengths from files
	call eminit (TRUE)

# Print tab table heading, if that form of output is used
	if (rmode < 0)
	    call emrshead (rmode)

# Get next object spectrum file name from the list
	while (imtgetim (speclist, specfile, SZ_PATHNAME) != EOF) {

# Decode specified apertures in multispec spectrum file
	ip = stridxs (lbracket,specfile)
	if (ip > 0) {
	    lfile = strlen (specfile)
	    specfile[ip] = EOS
	    jp = 0
	    ip = ip + 1
	    while (stridx (specfile[ip],rbracket) == 0 && ip <= lfile) {
		jp = jp + 1
		specnums[jp] = specfile[ip]
		specfile[ip] = EOS
		ip = ip + 1
		}
	    if (jp > 0)
		specnums[jp+1] = EOS
	    else
		call strcpy ("0",specnums,SZ_LINE)
	    if (decode_ranges (specnums,mspec_range,MAX_RANGES,nmspec) == ERR) {
		call sprintf (str, SZ_LINE, "T_EMSAO: Illegal multispec list <%s>")
		    call pargstr (specnums)
		call error (1, str)
		}
	    }
	else
	    nmspec = nmspec0
	nspec = nmspec
	if (debug) {
	    call printf ("EMSAO: Next file is %s [%s] = %d aps\n")
		call pargstr (specfile)
		call pargstr (specnums)
		call pargi (nspec)
	    call flush (STDOUT)
	    }

	call strcpy (specdir,specpath,SZ_PATHNAME)
	call strcat (specfile,specpath,SZ_PATHNAME)

# Check for readability of object spectrum
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("EMSAO: cannot read spectrum path %s \n")
		call pargstr (specpath)
	    next
	    }

	mspec = -1
	sspec = -1

# Get next multispec number from list
	while (nspec > 0 && get_next_number (mspec_range, mspec) != EOF) {

nextspec_
	if (debug) {
	    call printf ("EMSAO: Next aperture is %s [%d] = %d aps\n")
		call pargstr (specfile)
		call pargi (mspec)
		call pargi (nspec)
	    }

# Load spectrum
	call getspec (specpath, mspec, mband, spectrum, specim, wenable)
	if (specim == ERR)
	    go to endspec_
	call xcrhead (mspec, specim)
	if (!linefit)
	    call emrhead (mspec, specim)

# Turn off spectrum continuum subtraction, if not needed
	tcont = TRUE
	call imgbpar (specim,"SUBCONT",tcont)
	if (tcont)
	    conproc = SUBCONT
	else
	    conproc = NOCONT
	if (debug) {
	    if (tcont)
		call printf ("T_EMSAO: Spectrum continuum will be subtracted\n")
	    else
		call printf ("T_EMSAO: Spectrum continuum already subtracted\n")
	    call flush (STDOUT)
	    }

# If RENORM is not set, set it if maximum counts are less than 1
        if (!renorm) {
            spmax = Memr[spectrum]
            spmin = Memr[spectrum]
            do i = 2, specpix {
                speci = Memr[spectrum+i-1]
#		if (debug) {
#		    call printf ("T_EMSAO: %4d/%4d: %f\n")
#			call pargi (i)
#			call pargi (specpix)
#			call pargd (speci)
#		    call flush (STDOUT)
#		    }
                if (speci > spmax)
                    spmax = speci
                if (speci < spmin)
                    spmin = speci
                }
            if (spmax < 1.0 && spmin > -1.0)
                renorm = TRUE
            else
                renorm = FALSE
            }
	if (debug) {
	    call printf ("T_EMSAO: maximum count in spectrum is %f\n")
		call pargd (spmax)
	    call flush (STDOUT)
	    }

# Renormalize spectrum if requested
	if (renorm) {
	    sptot = 0.d0
	    nsp = 0
	    do i = 1, specpix {
		ispec = spectrum+i-1
		if (Memr[ispec] != 0.) {
		    sptot = sptot + double (Memr[ispec])
		    nsp = nsp + 1
		    }
		}
	    if (nsp > 0) {
		spmean = 0.001d0 * sptot / double (nsp)
		do i = 1, specpix {
		    ispec = spectrum+i-1
		    if (Memr[ispec] != 0.) {
			Memr[ispec] = Memr[ispec] / spmean
			}
		    }
		}
	    else {
		call eprintf ("*** Spectrum is all zeroes\n")
		call close_image (specim, specsh)
		next
		}
	    }

# Set up wavelength vector for intermediate graphs
	npix = NP2(specsh) - NP1(specsh)
	if (NP2(specsh) < NP1(specsh)) {
	    npix = NP1(specsh) - NP2(specsh) + 1
	    }
	else {
	    npix = NP2(specsh) - NP1(specsh) + 1
	    }
	if (npix > npix0) {
	    if (npix0 > 0) {
		call mfree (wlspec, TY_REAL)
		}
	    call malloc (wlspec, npix, TY_REAL)
	    npix0 = npix
	    }
	do ipix = NP1(specsh), NP2(specsh) {
	    Memr[wlspec+ipix-1] = wcs_p2w (double (ipix))
	    }
#	wlspec = SX(specsh) + NP1(specsh) - 1
	pxspec = spectrum + NP1(specsh) - 1
	if (debug) {
	    call printf ("T_EMSAO: %d-pixel spectrum set up\n")
		call pargi (specpix)
	    call flush (STDOUT)
	    }

# Load sky spectrum
	if (skyspec) {
	    if (nsspec > 0) {
		if (get_next_number (sspec_range, sspec) != EOF) {
		    if (sspec > 0) 
			call getsky (specpath,sspec,sband,specsky,skyim,skysh)
		    else {
			specsky = NULL
			skyspec = FALSE
			}
		    }
		}
	    else
		call getsky (specpath, sspec, sband, specsky, skyim, skysh)
	    if (debug) {
		call printf ("EMSAO: Sky spectrum %d read %s\n")
		    call pargi (sspec)
		    call pargstr (skyname)
		}
	    pxsky = specsky + NP1(skysh) - 1
	    }
	else
	    pxsky = NULL

# Starting wavelength from image header or parameter file
	minwav = W0(specsh)
	maxwav = W1(specsh)
	if (maxwav < minwav) {
	    maxwav = W0(specsh)
	    minwav = W1(specsh)
	    }

# Use limit parameter if greater than minimum wavelength of spectrum
	if (minwav0 != dindef) {
	    if (minwav0 > 0.0) {
		if (minwav0 > minwav)
		    minwav = minwav0
		}

# If minimum parameter is <0, increase lower limit by -minimum parameter
	    else if (minwav0 < 0.0)
		minwav = minwav - minwav0
	    }
	pix1 = idnint (wcs_w2p (minwav))
	minwav = wcs_p2w (double(pix1))

# Ending wavelength from image header or parameter file
# Use limit parameter if less than maximum wavelength of spectrum
	if (minwav0 != dindef) {
	    if (maxwav0 > 0.0) {
		if (maxwav0 < maxwav)
		    maxwav = maxwav0
		}

# If maximum parameter is <0, decrease upper limit by -maximum parameter
	    else if (maxwav0 < 0.0)
		maxwav = maxwav + maxwav0
	    }
	pix2 = idnint (wcs_w2p (maxwav))
	maxwav = wcs_p2w (double(pix2))

	if (pix1 > pix2) {
	    pix2 = idnint (wcs_w2p (minwav))
	    pix1 = idnint (wcs_w2p (maxwav))
	    }

	if (debug) {
	    npix = pix2 - pix1 + 1
	    call printf("EMSAO: from %10.3fa(%d) to %10.3fa(%d), npts = %d \n")
	    call pargd (minwav)
	    call pargi (pix1)
	    call pargd (maxwav)
	    call pargi (pix2)
	    call pargi (npix)
	    call flush (STDOUT)
	    }

# Find emission line velocity for this spectrum
	IM_UPDATE(specim) = NO
	newresults = FALSE
	mspec0 = mspec
	call emfit (Memr[pxspec],Memr[wlspec],Memr[pxsky],
		    specfile,mspec,specim,pix1,pix2,rmode)
	if (debug) {
	    call printf ("T_EMSAO: back from EMFIT\n")
	    call flush (STDOUT)
	    }
	if (mspec != mspec0) {
	    call close_image (specim, specsh)
	    if (skyspec)
		call close_image (skyim, skysh)
	    go to nextspec_
	    }

# Save results to TDC archive file, if requested
	if (arcflag)
	    call emarch (specim, specfile)

# Save Cz and error to image header, if requested
	if (savespec || savevel || qplot) {
	    if (imaccess (specpath, READ_WRITE) == NO) {
		call eprintf ("EMSAO: cannot write to %s; not saving results\n")
		    call pargstr (specpath)
		}
	    else if (savevel) {
		call vwhead (mspec,specim)
		call emwhead (mspec,specim)
		}
	    else if (qplot && newresults) {
		call qwhead (mspec,specim,"emsao")
		call emwhead (mspec,specim)
		}
	    else if (qplot) {
		call qwhead (mspec,specim,"emsao")
		}
	    if (savespec) {
		call putspec (specim, mspec, spectrum, TY_REAL, debug)
		savespec = FALSE
		}
	    }
	if (debug) {
	    call printf ("T_EMSAO: about to close image\n")
	    call flush (STDOUT)
	    }

# Close the object spectrum image and headers
endspec_
	call close_image (specim, specsh)

# Close the sky spectrum image
	if (skyspec)
	    call close_image (skyim, skysh)

# Move on to next aperture or next image
	nspec = nspec - 1

# End of multispec loop within single image
	}

# End of image loop
	}

	if (debug) {
	    call printf ("T_EMSAO: about to close logs and files\n")
	    call flush (STDOUT)
	    }
 
# Close the log files
	do i = 1, nlogfd {
	    call close (logfd[i])
	    }                                              
	call imtclose (speclist)
	if (npix0 > 0)
	    call mfree (wlspec, TY_REAL)

end
# Nov 20 1991	Free buffers before malloc'ing them
# Dec 18 1991	Set up and pass wavelength vector

# Feb 18 1992	Rewrite to add multispec file reading and writing
# Mar 24 1992	Read line information using eminit and save using emarch
# Apr 22 1992	Fix multispec file handling
#		Implement wavelength limits
# May 22 1992	Initialize mspec to -1
# Aug 10 1992	Add specim to emarch arguments
# Aug 11 1992	Write HISTORY line to image header
# Sep 25 1992	Add option to renormalize data (useful if in flux)
# Oct 22 1992	Make hstr and dstr strings, not pointers
# Nov 30 1992	Initialize emission lines after log files opened

# Feb  1 1993	Allow INDEF for wavelength limits
# Feb  2 1993	Exit if a zero spectrum is found when renormalizing
# Jun  4 1993	Add mwcs for wavelengths
# Jun 14 1993	Get wavelength mapping from spectrum header
# Jun 16 1993	Set version to 1.1
# Jul  7 1993	Add spectrum header to getspec and getsky
# Jul  8 1993	Close sh properly
# Aug 20 1993	Deal with double wavelength vector
# Dec  1 1993	Write history only after last multispec line is done
# Dec  3 1993	Pass mspec to EMFIT
# Dec  3 1993	Update to version 1.2

# Feb  3 1994	Proceed gracefully if files are not useable
# Mar 22 1994	Add multispec aperture specification for each image
# Mar 22 1994	Update to version 1.3
# Apr 11 1994	Pass full spectrum to emfit
# Apr 12 1994	Allocate rather than reallocate wlspec
# Apr 12 1994	Return MWCS header pointer from getspec and getsky
# Apr 14 1994	Fix bug by not changing specpix
# Apr 21 1994	Initialize sspec
# Apr 22 1994	Let SHDR_CLOSE free data vectors
# Apr 26 1994	Don't malloc anything
# Apr 28 1994	Use while instead of goto
# May  3 1994	Clean up code
# May  3 1994	Add smoothed vector argument to PLOTSPEC call
# May  5 1994	Make smoothed vector array, not pointer
# May 16 1994	Add report mode so 1-line reports can be generated
# May 23 1994	Clean up write-protect check
# Jun 23 1994	Do not use getim labelled common
# Jun 23 1994	Keep MWCS pointer in SHDR structure
# Jun 23 1994	Update to version 1.4
# Jul 29 1994	Read version from fquot.h
# Aug  1 1994	Move header writing to emhead
# Aug  2 1994	Fix logging bug
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug  3 1994	Read dispersion error as parameter
# Aug  4 1994	Move specnums and nmspec to labelled common
# Aug  8 1994   Add specdir as parameter
# Aug 10 1994	Write log heading in one line instead of two
# Aug 12 1994	Change name of header-writing subroutine
# Aug 17 1994	Set nmspec correctly
# Aug 24 1994	Fix bug so SPECDIR can be null
# Sep  9 1994	Keep from crashing if log-lambda file
# Sep 21 1994	Print emsao header correctly
# Sep 23 1994	Reset NTEMP to 0 each time a new spectrum is read
# Dec  7 1994	Do not set NTEMP
# Dec  8 1994	Move mwcs initialization to GETSPEC

# Jan 31 1995	Change lengths of file and directory names to sz_pathname
# Mar 13 1995	Use CLOSE_IMAGE to free allocated memory
# Mar 29 1995	Improve error handling
# Apr  5 1995	Change selection choice vector names
# May 15 1995	Change all sz_fname to sz_line, which is 100 chars longer
# Jul 19 1995	Call VWHEAD as well as EMWHEAD to write velocity to header
# Sep 21 1995	Set im_update to no before EMPLOT so qplotting works
# Sep 21 1995	Move check for writeability to GETSPEC
# Oct  4 1995	Default SPECNUMS to 0 if null string is input

# Aug  7 1996	Use smw.h

# Feb 12 1997	Add switch to plot continuum-subtracted spectrum
# Mar 14 1997	Drop SPECSH from GETSPEC call; pass it in common
# Apr  9 1997	Deal with null aperture list in DECODE_RANGES, not here
# May  2 1997	Always test against dindef, not INDEFD
# May 19 1997	Add parameter specifying template for initial velocity or report
# Aug 27 1997	Add bands for sky and spectrum for multispec spectra

# Apr 22 1998	Fix bug in sky spectrum retrieval
# Jun 12 1998	Use only portion of spectrum covered by WCS

# Mar 11 1999	Add write argument to getspec() call
# Mar 17 1999	Multiply by 1000 when renormalizing
# Jul 16 1999	Fix bug so QPLOT updates header

# Jan 25 2000	Read spectrum header flag to decide whether to remove contin.
# Jul  5 2000	Renormalize if maximum value is less than one

# Mar 22 2001	Fix limit-setting to work when wavelength runs right to left

# Aug 25 2004	Fix bug involving misuse for spectrum pointer variable
# Aug 25 2004	Add third argument to first qwhead() call

# May 25 2005	Compute wavelength vector so shift is included
# Jul 14 2005	Fix bug of calling mfree without data type
# Nov  3 2005	Force renormalization only if -1 < spectrum < 1

# Jan 31 2007	Drop wavelength limits from call to emfit()

# Mar 14 2008	Use minwav0 and maxwav0 as spacing from edge if negative
# Oct  7 2008	If mspec changed by emfit, run with new value

# May 20 2009	Set mspec to closest limit if exceeded

# May  2 2014	Write out spectrum if it has been changed
