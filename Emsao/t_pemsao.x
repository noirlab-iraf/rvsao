# File rvsao/Emsao/t_pemsao.x
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

procedure t_pemsao ()

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
int	nsum
int	nsp, npix0, nflines
double	sptot, spmax, spmin, speci
pointer	ispec
double	vl, dwdp
 
pointer	imtopenp()
int	imaccess(), imtgetim()
bool	clgetb()
int	clgeti()
int	clpopnu(), clgfil(), open()
int	strdic(),strncmp()
double	clgetd()
int	decode_ranges(),get_next_number(),stridx(),strlen(),stridxs()
double	wcs_w2p(), wcs_p2w()

define	endspec_ 90

include	"rvsao.com"
include	"emv.com"
include "contin.com"
include "results.com"
include "emsums.com"
 
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
	    call sprintf (str, SZ_LINE, "PEMSAO: Illegal multispec object list <%s>")
		call pargstr (specnums)
	    call error (1, str)
	    }

# Multispec band number
	mband = clgeti ("specband")

# FITS Extension number (passed through rvsao.com)
	mext = clgeti ("specext")

# Multispec sky spectrum numbers
	call clgstr ("skynum",skynums,SZ_LINE)
	if (skynums[1] == EOS)
	    call strcpy ("0",skynums,SZ_LINE)
	skyspec = FALSE
	if (strncmp (skynums,"0",1) != 0) {
	    if (decode_ranges (skynums, sspec_range, MAX_RANGES, nsspec) == ERR) {
		call sprintf (str, SZ_LINE, "PEMSAO: Illegal multispec sky list <%s>")
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
	pltcon  = clgetb ("contsub_plot")
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
	    nsum = 0

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
		call sprintf (str, SZ_LINE, "PEMSAO: Illegal multispec list <%s>")
		    call pargstr (specnums)
		call error (1, str)
		}
	    }
	else
	    nmspec = nmspec0
	nspec = nmspec
	if (debug) {
	    call printf ("PEMSAO: Next file is %s [%s] = %d aps\n")
		call pargstr (specfile)
		call pargstr (specnums)
		call pargi (nspec)
	    call flush (STDOUT)
	    }

	call strcpy (specdir,specpath,SZ_PATHNAME)
	call strcat (specfile,specpath,SZ_PATHNAME)

# Allocate vector for median shift
	call malloc (evel, nmspec, TY_DOUBLE)
	call malloc (ewl, nmspec, TY_DOUBLE)
	call malloc (epix, nmspec, TY_DOUBLE)

# Check for readability of object spectrum
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("PEMSAO: cannot read spectrum path %s \n")
		call pargstr (specpath)
	    next
	    }

# Get next multispec number from list
	mspec = -1
	sspec = -1
	sumvel = 0.d0
	sumerr = 0.d0
	sumdwl = 0.d0
	sumpix = 0.d0
	nap = 0
	nflines = 0

	while (nspec > 0 && get_next_number (mspec_range, mspec) != EOF) {

	if (debug) {
	    call printf ("PEMSAO: Next aperture is %s [%d] = %d aps\n")
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
		call printf ("PEMSAO: Spectrum continuum will be subtracted\n")
	    else
		call printf ("PEMSAO: Spectrum continuum already subtracted\n")
	    call flush (STDOUT)
	    }

# If RENORM is not set, set it if maximum counts are less than 1
        if (!renorm) {
            spmax = Memr[spectrum]
            spmin = Memr[spectrum]
            do i = 2, specpix {
                speci = Memr[spectrum+i-1]
#		if (debug) {
#		    call printf ("PEMSAO: %4d/%4d: %f\n")
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
	    call printf ("PEMSAO: maximum count in spectrum is %f\n")
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
	    call printf ("PEMSAO: %d-pixel spectrum set up\n")
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
		call printf ("PEMSAO: Sky spectrum %d read %s\n")
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
	if (minwav0 > 0.0 && minwav0 != dindef && minwav0 > minwav)
	    minwav = minwav0
	pix1 = idnint (wcs_w2p (minwav))
	minwav = wcs_p2w (double(pix1))

# Ending wavelength from image header or parameter file
	if (maxwav0 > 0.0 && maxwav0 != dindef && maxwav0 < maxwav)
	    maxwav = maxwav0
	pix2 = idnint (wcs_w2p (maxwav))
	maxwav = wcs_p2w (double(pix2))
	if (pix1 > pix2) {
	    pix2 = idnint (wcs_w2p (minwav))
	    pix1 = idnint (wcs_w2p (maxwav))
	    }

	if (debug) {
	    npix = pix2 - pix1 + 1
	    call printf("PEMSAO: from %10.3fa(%d) to %10.3fa(%d), npts = %d \n")
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
	call emfit (Memr[pxspec],Memr[wlspec],Memr[pxsky],
		    specfile,mspec,specim,pix1,pix2,rmode)

	if (spnlf > 0) {
	    spx = wcs_w2p (wlobs[1]) - wcs_w2p (wlref[1])
	    swl = wlobs[1] - wlref[1]
	    vl = wcs_w2p (wlobs[1])
	    dwdp = wcs_p2w (vl + 1.d0) - wcs_p2w (vl)
	    }
	else {
	    spevel = 0.d0
	    speerr = 0.d0
	    spx = 0.d0
	    swl = 0.d0
	    vl = wcs_w2p (wlref[1])
	    dwdp = wcs_p2w (vl + 1.d0) - wcs_p2w (vl)
	    }
	if (nmspec > 0) {
	    if (spnlf > 0) {
		nsum = nsum + 1
		if (debug) {
		    call printf ("PEMSAO: %3d / %3d: %.3f (%.3f) %.3fA %.3fp %d lines\n")
			call pargi (nap+1)
			call pargi (nmspec)
			call pargd (spevel)
			call pargd (speerr)
			call pargd (swl)
			call pargd (spx)
			call pargi (spnlf)
		    }
		Memd[evel+nsum-1] = spevel
		sumvel = sumvel + spevel	
		sumerr = sumerr + (speerr * speerr)
		Memd[epix+nsum-1] = spx
		sumpix = sumpix + spx
		Memd[ewl+nsum-1] = swl
		sumdwl = sumdwl + swl
		nflines = nflines + spnlf
		}
	    else if (debug) {
		call printf ("PEMSAO: %3d / %3d: %.3f (%.3f) %d lines\n")
		    call pargi (nap+1)
		    call pargi (nmspec)
		    call pargd (spevel)
		    call pargd (speerr)
		    call pargi (spnlf)
		}
	    nap = nap + 1
	    }

# Save results to TDC archive file, if requested
	if (arcflag)
	    call emarch (specim, specfile)

# Save Cz and error to image header, if requested
	if ((savespec || savevel || qplot) && nmspec < 10) {
	    if (imaccess (specpath, READ_WRITE) == NO) {
		call eprintf ("PEMSAO: cannot write to %s; not saving results\n")
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
#	if (debug) {
#	    call printf ("PEMSAO: %3d / %3d / %3d completed\n")
#		call pargi (nsum)
#		call pargi (nap)
#		call pargi (nmspec)
#	    call flush (STDOUT)
#	    }

# Close the object spectrum image and headers
endspec_

	if (nmspec < 1 || mspec < nmspec)
	    call close_image (specim, specsh)

# Close the sky spectrum image
	if (skyspec)
	    call close_image (skyim, skysh)

#	if (debug) {
#	    call printf ("PEMSAO: (%3d / %3d) %3d / %3d images closed\n")
#		call pargi (nsum)
#		call pargi (nap)
#		call pargi (mspec)
#		call pargi (nmspec)
#	    call flush (STDOUT)
#	    }

# Move on to next aperture or next image
	nspec = nspec - 1

# End of multispec loop within single image
	}

# End of image loop
	}

# Save global results in the parameter file
	if (nmspec > 0) {
	    call printsum (specim, nsum)
	    }

# Save some results for first line in the parameter file
	if (!savevel)
	    call printline (nflines, dwdp)

# Close the log files
	do i = 1, nlogfd {
	    call close (logfd[i])
	    }                                              
	call imtclose (speclist)
	if (npix0 > 0)
	    call mfree (wlspec, TY_REAL)

end


# Save velocity results in image header if savevel is true
# Save velocity results in parameter file if savevel is false

procedure printsum (specim, nsum)

pointer specim	# Object image header structure
int	nsum	# Number of apertures for which lines were found

int	ip,jp,lfile	# Limits for multispec aperture decoding

char	str[SZ_LINE], sstr[SZ_LINE]
int	i,j,iline,ldir
int	fd
double	avgvel, avgerr, medvel, avgpix, medpix, avgdwl, meddwl, q1, q2
double	vl
 
include	"rvsao.com"
include	"emv.com"
include "results.com"
include "emsums.com"
 
begin

	if (savevel) {
	    IM_UPDATE(specim) = YES
	    }

# Save global results in the parameter or image file
	if (debug) {
	    if (savevel) {
		call printf("PRINTSUMS: About to save results in header\n")
		call strcpy ("would be", sstr, SZ_LINE)
		}
	    else {
		call printf("PRINTSUMS: About to save results as parameters\n")
		call strcpy ("saved as", sstr, SZ_LINE)
		}
	    call flush (STDOUT)
	    }

# Mean velocity and error
	if (nsum > 0)
	    avgvel = sumvel / double (nsum)
	else
	    avgvel = 0.d0
	if (debug) {
	    call printf("PEMSAO: Mean velocity = %.3f %s meanvel\n")
		call pargd (avgvel)
		call pargstr (sstr)
	    call flush (STDOUT)
	    }
	if (nsum > 0)
	    avgerr = dsqrt (sumerr) / double (nsum)
	else
	    avgerr = 0.d0
	if (debug) {
	    call printf("PEMSAO: Mean velocity error = %.3f %s meanerr\n")
		call pargd (avgerr)
		call pargstr (sstr)
	    call flush (STDOUT)
	    }

# Median velocity and error
	call median (Memd[evel], nsum, medvel, q1, q2)
	if (debug) {
	    call printf("PEMSAO: Median velocity = %.3f %s medvel\n")
		call pargd (medvel)
		call pargstr (sstr)
	    call printf("PEMSAO: First quartile = %.3f %s medq1\n")
		call pargd (q1)
		call pargstr (sstr)
	    call printf("PEMSAO: Second quartile = %.3f %s medq2\n")
		call pargd (q2)
		call pargstr (sstr)
	    call flush (STDOUT)
	    }

# If not saving in header, save in parameter file
	if (!savevel) {
	    call clputd ("meanvel", avgvel)
	    call clputd ("meanerr", avgerr)
	    call clputd ("medvel",medvel)
	    call clputd ("medq1",q1)
	    call clputd ("medq2",q2)
	    }

# Mean wavelength and pixel shift
	if (nsum > 0) {
	    avgdwl = sumdwl / double (nsum)
	    avgpix = sumpix / double (nsum)
	    }
	else {
	    avgdwl = 0.d0
	    avgpix = 0.d0
	    }
	if (debug) {
	    call printf("PEMSAO: Mean pixel shift = %.3f %s meanpix\n")
		call pargd (avgpix)
		call pargstr (sstr)
	    call printf("PEMSAO: Mean wavelength shift = %.3f %s meanwl\n")
		call pargd (avgdwl)
		call pargstr (sstr)
	    call flush (STDOUT)
	    }

# Median wavelength and pixel shift
	call median (Memd[ewl], nsum, meddwl, q1, q2)
	call median (Memd[epix], nsum, medpix, q1, q2)
	if (debug) {
	    call printf("PEMSAO: Median pixel shift = %.3f %s medpix\n")
		call pargd (medpix)
		call pargstr (sstr)
	    call printf("PEMSAO: Median wavelength shift = %.3f %s medpix\n")
		call pargd (meddwl)
		call pargstr (sstr)
	    call flush (STDOUT)
	    }

	if (savevel) {
	    call imaddd (specim, "MEANPIX", avgpix, "Mean pixel shift")
	    call imaddd (specim, "MEDPIX", medpix, "Median pixel shift")
	    call imaddd (specim, "MEANDWL", avgdwl, "Mean wavelength shift")
	    call imaddd (specim, "MEDDWL", meddwl, "Median wavelength shift")
	    }
	else {
	    call clputd ("meanpix", avgpix)
	    call clputd ("medpix",medpix)
	    call clputd ("meandwl", avgdwl)
	    call clputd ("meddwl",meddwl)
	    }

	call close_image (specim, specsh)
	return
end


procedure printline (nflines, dwdp)

int	nflines
double	dwdp	# Wavelength per pixel
 
double linvel, linerr, linwidth

include	"rvsao.com"
include	"emv.com"
include "results.com"

begin

# Save some results for first line in the parameter file
	if (debug) {
	    call printf("PEMSAO: Number of lines fit = %d saved\n")
		call pargi (nflines)
	    call printf("PEMSAO: First line fit = %s saved as emline\n")
		call pargstr (nmobs[1,1])
	    call printf("PEMSAO: Rest wavelength = %.3f saved as wlrest\n")
		call pargd (wlrest[1])
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputi ("nlfit", nflines)
            call clpstr ("emline", nmobs[1,1])
            call clputd ("wlrest", wlrest[1])
	    }

# Line velocity
	linvel = (c0 * emparams[9,1,1]) + spechcv
	if (debug) {
	    call printf("PEMSAO: Velocity = %.4f saved as velocity\n")
		call pargd (linvel)
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputd ("velocity", linvel)
	    }

# Line velocity error
	linerr = c0 * emparams[9,2,1]
	if (debug) {
	    call printf("PEMSAO: Velocity error = %.4f saved as velerr\n")
		call pargd (linerr)
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputd ("velerr", linerr)
	    }

# Line height
	if (debug) {
	    call printf("PEMSAO: Line height = %.2f saved as lineheight\n")
		call pargd (emparams[5,1,1])
	    call flush (STDOUT)
	    }
	if (!savevel) {
	    call clputd ("lineheight", emparams[5,1,1])
	    }

# Line width
	linwidth = emparams[6,1,1] * abs (dwdp)
	if (debug) {
	    call printf("PEMSAO: Line width = %.4f saved as linewidth\n")
		call pargd (linwidth)
	    call printf("PEMSAO: Line equivalent width = %.4f saved as lineeqw\n")
		call pargd (emparams[7,1,1])
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputd ("linewidth", linwidth)
	    call clputd ("lineeqw", emparams[7,1,1])
	    }

	return
end
# Aug 23 2007	New task based on EMSAO to save parameters of one line

# Mar 19 2008	Save mean velocity and error as parameters, too
# Mar 19 2008	Save total number of lines fit across apertures
# Mar 25 2008	Do not compute averages if no lines fit
# Aug  5 2008	Add median and quartiles as output parameters
# Aug  6 2008	Decrease number of variables so it will compile
# Aug  6 2008	Add mean and median pixel and wavelength shifts
# Sep 15 2008	Add debugging printout for mean and median shifts
# Sep 15 2008	Save linewidth as positive when wavelength scale is reversed

# Jun  4 2009	Compute wavelength per pixel by pixel offset for observed center

# Mar  1 2010	Write mean and median pixel shifts to image header if save_vel
# Mar  2 2010	Do not write to parameter file if save_vel=yes
# Mar  8 2010	Skip apertures where no lines were found
# Jul 08 2010	Read specext FITS extension number into mext in rvsao.com
# Aug 13 2010	In debug mode, print shift in pixels and wavelength too
# Aug 13 2010	Add mean and median wavelength shifts to FITS file header
# Aug 25 2010	Write mean and median wavelength shift parameters

# May 25 2012	Drop unused declarations to get back below limit

# May  2 2014	Write out edited spectrum
