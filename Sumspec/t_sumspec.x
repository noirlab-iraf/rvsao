# File rvsao/Makespec/t_sumspec.x
# June 10, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1997-2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# SUMSPEC is an IRAF task for combining spectra with different velocities
# For arguments, see parameter file sumspec.par.
# Information is shared in common blocks defined in "rvsao.com".
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include "rvsao.h"

procedure t_sumspec ()

int	i
char	spectra[SZ_PATHNAME]	# List of input spectra
char	specfile[SZ_PATHNAME]	# Input spectrum file name
char	specdir[SZ_PATHNAME]	# Directory for input spectra
char	specpath[SZ_PATHNAME]	# Input spectrum path name
char	compfiles[SZ_PATHNAME]	# List of input spectra
char	compfile[SZ_PATHNAME]	# Output spectrum file name
char	compdir[SZ_PATHNAME]	# Output spectrum directory
char	comppath[SZ_PATHNAME]	# Output spectrum pathname
#char	filt_type[SZ_LINE]	# Filter for transform
				# (ramp | welch | hanning | cos-bell)
char	svel_corr[SZ_LINE]	# Type of velocity correction for spectrum
				# (none | file | heliocentric | barycentric)
int	logfiles		# List of log files
char	logfile[SZ_PATHNAME]	# Log file name
char	wtitle[20]		# Title for wavelength plots of spectrum

int	mspec		# Object aperture to read from multispec file
int	mband		# Object band to read from multispec file

pointer	speclist	# List of spectrum files
pointer	complist	# List of output files
pointer	compspec	# Summed spectrum
char	str[SZ_LINE]
int	fd
int	ldir
bool	fproc

int	nmspec0		# Number of object multispec spectra
int	mspec_range[3,MAX_RANGES]
int	ip,jp,lfile	# Limits for multispec aperture decoding
double	zcomp
char	lbracket[3]	# "[({"
char	rbracket[3]	# "])}"
pointer compim		# Output spectrum header structure
int	ispec, jspec
int	lcomp
int	nspec		# Number of spectra in output file
int	npix
int	indefi
int	nspecap
bool	specint, compint, perspec
double	minwav, maxwav, pixwav
double	dindef
 
bool	clgetb()
int	clpopnu(), clgeti(), clgfil(), open()
double	clgetd()
int	strdic(), stridx(), stridxs()
int	decode_ranges(),get_next_number()
int	imtgetim(), imaccess(), strlen()
pointer	imtopen()
real	clgetr()

define	newspec_ 10
define	newap_	 20
define	endcomp_	 90

include	"rvsao.com"
include	"sum.com"
 
begin
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	call sprintf (wtitle,20,"Wavelength")
	c0 = 299792.5
	dindef = INDEFD
	indefi = INDEFI
	nlogfd = 0
	compim = NULL
	velhc = dindef
	npix = 0
	fproc = FALSE

# Get task parameters.

# Spectra to combine
	call clgstr ("spectra",spectra,SZ_PATHNAME)
	speclist = imtopen (spectra)

# Multispec spectrum apertures (use only first if multiple files)
	call clgstr ("specnum",specnums,SZ_LINE)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR)
	    call error (1, "SUMSPEC: Illegal multispec list")
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }
	mband = clgeti ("specband")

# Number of output spectra in output file
	nspec = clgeti ("nspec")

# File to which to write summed spectrum
	call clgstr ("compfile",compfiles,SZ_PATHNAME)
	lcomp = strlen (compfiles)
	complist = imtopen (compfiles)

	save_names = FALSE
	save_names = clgetb ("save_names")

# Optional intermediate data plot switches
	pltspec = clgetb ("spec_plot")
	plttemp = clgetb ("comp_plot")

# Print processing information
	debug  = clgetb ("debug")

# Continuum fit parameter pset
	call csum_get_pars()
	pltcon  = clgetb ("cont_plot")

# Number of times to smooth (1-2-1) final data plot
	nsmooth = clgeti ("nsmooth")

# Minimum and maximum values for data in graphs
	ymin = clgetr ("ymin")
	ymax = clgetr ("ymax")

# Keep header from first spectrum
	copyhead = TRUE
	copyhead = clgetb ("copy_header")
	if (nspec == 1)
	    copyhead = TRUE

# Interact with display of final composite spectrum
	compint = clgetb ("comp_int")
	if (compint)
	    tsmooth = nsmooth
	else
	    tsmooth = -(nsmooth + 1)

# Interact with display of input spectra
	specint = clgetb ("spec_int")
	if (!specint)
	    nsmooth = -(nsmooth + 1)

# Type of heliocentric velocity correction to be used
	call clgstr ("svel_corr",svel_corr,SZ_LINE)
	svcor = strdic (svel_corr,svel_corr,SZ_LINE, HC_VTYPES)

# Optional wavelength limits for output spectrum
	minwav0 = clgetd ("st_lambda")
	maxwav0 = clgetd ("end_lambda")
	if (minwav0 == dindef && maxwav0 == dindef && lcomp == 0)
	    perspec = TRUE
	else
	    perspec = FALSE
	pixwav = clgetd ("pix_lambda")
	npts = clgeti ("npts")

# If start, increment, and number of pixels are set, compute ending wavelength
	if (npts != indefi && minwav0 != dindef && pixwav != dindef &&
	    maxwav0 == dindef)
	    maxwav0 = minwav0 + (pixwav * double (npts - 1))

# If start, increment, and last wavelength are set, compute number of pixels
	else if (npts == indefi && minwav0 != dindef && pixwav != dindef &&
		 maxwav0 != dindef)
	    npts = 1 + idnint (((maxwav0 - minwav0) / pixwav) + 0.0001d0)

# If increment, last, and number of pixels are set, compute starting wavelength
	else if (npts != indefi && maxwav0 != dindef && pixwav != dindef &&
		 minwav0 == dindef)
	    minwav0 = maxwav0 - (pixwav * double (npts - 1))

# If start, last, and number of pixels are set, compute wavelength increment
	else if (npts != indefi && maxwav0 != dindef && minwav0 != dindef &&
		 pixwav == dindef)
	    pixwav = (maxwav0 - minwav0) / double (npts - 1)

# If not completely specified, get wavelength limits from input spectra
	else if (minwav0 == dindef || maxwav0 == dindef && !perspec) {
	    call wlrange (spectra, minwav, maxwav, npix)
	    if (minwav0 == dindef)
		minwav0 = minwav
	    if (maxwav0 == dindef)
		maxwav0 = maxwav
	    if (npts == indefi)
		npts = npix
	    pixwav = (maxwav0 - minwav0) / double (npts - 1)
	    speclist = imtopen (spectra)
	    }
	if (debug || npts < 3) {
	    call printf ("SUMSPEC: %d-point spectrum from %.3fA to %.3fA by %.3fA\n")
		call pargi (npts)
		call pargd (minwav0)
		call pargd (maxwav0)
		call pargd (pixwav)
	    call flush (STDOUT)
	    }
	if (npts < 3)
	    go to endcomp_

# Redshift velocity for composite spectrum
	velocity = clgetd ("velcomp")
	zcomp = clgetd ("zcomp")
	if (zcomp != dindef)
	    velocity = c0 * zcomp
	if (debug) {
	    call printf ("SUMSPEC: Z = %.5f, CZ = %.2f\n")
		call pargd (zcomp)
		call pargd (velocity)
	    call flush (STDOUT)
	    }
 
# Open log files and write a header.
	logfiles = clpopnu ("logfiles")
	call fseti (STDOUT, F_FLUSHNL, YES)
	i = 0
	call strcpy ("rvsao.sumspec",taskname,SZ_LINE)
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
	ispec = 0

# Get next object spectrum file name from the list
newspec_
	if (imtgetim (speclist, specfile, SZ_PATHNAME) == EOF) {
	    if (debug) {
		call printf ("SUMSPEC: No more files in input list\n")
		call flush (STDOUT)
		}
	    go to endcomp_
	    }
	if (debug) {
	    call printf ("SUMSPEC: next file is %s\n")
		call pargstr (specfile)
	    call flush (STDOUT)
	    }

# Get next output file from list if first spectrum or first output filled
	if (ispec == 0 || (ispec >= nspec && fproc)) {
	    if (imtgetim (complist,compfile,SZ_PATHNAME) == EOF) {
		if (debug) {
		    call printf ("SUMSPEC: Cannot read image %s\n")
			call pargstr (compfile)
		    call flush (STDOUT)
		    }
		go to endcomp_
		}
	    lcomp = strlen (compfile)
	    call clgstr ("compdir",compdir,SZ_PATHNAME)
	    ldir = strlen (compdir)
	    if (ldir > 0) {
		if (compdir[ldir] != '/') {
		    compdir[ldir+1] = '/'
		    compdir[ldir+2] = EOS
		    }
		call strcpy (compdir,comppath,SZ_PATHNAME)
		if (lcomp > 0)
		    call strcat (compfile,comppath,SZ_PATHNAME)
		}
	    else if (lcomp > 0)
		call strcpy (compfile,comppath,SZ_PATHNAME)
	    else
		comppath[0] = EOS
	    if (lcomp == 0)
		fproc = TRUE
	    else
		fproc = FALSE
	    if (debug) {
		call printf ("SUMSPEC: next output file is %s\n")
		call pargstr (comppath)
		call flush (STDOUT)
		}
	    }

# Check for specified apertures in multispec spectrum file
	ip = stridxs (lbracket,specfile)
	if (ip > 0) {
	    lfile = strlen (specfile)
	    specfile[ip] = EOS
	    ip = ip + 1
	    jp = 1
	    while (stridx (specfile[ip],rbracket) == 0 && ip <= lfile) {
		specnums[jp] = specfile[ip]
		specfile[ip] = EOS
		ip = ip + 1
		jp = jp + 1
		}
	    specnums[jp] = EOS
	    if (decode_ranges (specnums,mspec_range,MAX_RANGES,nmspec) == ERR)
		call error (1, "SUMSPEC: Illegal multispec list")
	    }
	else {
	    nmspec = nmspec0
	    }
	nspecap = nmspec
	if (debug) {
	    call printf ("SUMSPEC: next file is %s [%s] = %d aps\n")
		call pargstr (specfile)
		call pargstr (specnums)
		call pargi (nmspec)
	    call flush (STDOUT)
	    }

# Check for readability of object spectrum
	call strcpy (specdir,specpath,SZ_PATHNAME)
	call strcat (specfile,specpath,SZ_PATHNAME)
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("SUMSPEC: cannot read spectrum file %s \n")
		call pargstr (specpath)
	    go to newspec_
	    }

# Get next multispec number from list
	mspec = -1
newap_
	if (nmspec <= 0)
	    go to newspec_
	if (get_next_number (mspec_range, mspec) == EOF)
	    go to newspec_

	ispec = ispec + 1
	call addspec (specfile, specdir, mspec, mband, nspecap, fproc,
		      compim ,comppath, compspec, ispec, nspec, perspec)

	if (nspec > 1)
	    call tmp_write_iraf (compim, ispec, compspec, TY_REAL, debug)

	if (fproc) {
	    call tmp_close (compim,compspec,debug)
	    compim = NULL
	    }

# Move on to next aperture or next image
	nmspec = nmspec - 1
	if (nmspec > 0)
	    go to newap_
	if (debug) {
	    call printf ("SUMSPEC: Reading file %d next\n")
		call pargi (ispec+1)
	    call flush (STDOUT)
	    }
	go to newspec_
 
endcomp_

# If fewer than nspec spectra read, fill out the file with zeroes
	if (ispec < nspec) {
	    call aclrr (Memr[compspec], npts)
	    do jspec = ispec+1, nspec {
		call tmp_write_iraf (compim, ispec, compspec, TY_REAL, debug)
		}
	    }

# Close the log files
	if (nlogfd > 0) {
	    do i = 1, nlogfd {
		call close (logfd[i])
		}                                              
	    }

#  Close spectrum list
	call imtclose (speclist)

#  Close the output composite spectrum file
	if (compim != NULL) {
	    if (nspec > 1) {
		call imunmap (compim)
		call mfree (compspec, TY_REAL)
		}
	    else
		call tmp_close (compim,compspec,debug)
	    }

end
# Jan 14 1997	Add sum.com labelled common
# Jan 14 1997	Change parameter tempvel to veltemp
# Feb  3 1997	Call CSUM_GET_PARS instead of CONT_GET_PARS
# Mar 17 1997	Fix log file heading
# Apr 18 1997	Fix final template interaction
# Apr 21 1997	If limits are indef, use overlap between all input spectra
# Apr 22 1997	Change parameter template to tempfile
# Apr 29 1997	Change name from SUMTEMP to SUMSPEC and variable names, too
# May  2 1997	Fix wavelength limit problem
# May  5 1997	Pass SPECLIST to WLRANGE
# May 16 1997	Leave VELOCITY set to parameter file value, even INDEF
# Jun 18 1997	Fix handling of INDEF as VELOCITY or Z0 value
# Jul 22 1997	Add flag to make saving filenames in output file optional
# Jul 27 1997	If composite spectrum filename is null, print range and exit
# Aug 27 1997	Add parameter for mul;tispec spectrum band

# Dec 18 1998	Add parameter nspec to write to 2D output file

# May 11 1999	Fix bug so 2-D output files are written correctly
# Jun 16 1999	If fewer input spectra than output 2-D dimension, zero-pad
# Jul  6 2000	Add option to specify 3 of 4 dispersion parameters
# Jul 21 2000	Initialize barycentric velocity correction to INDEF

# Apr 25 2001	Allow null output filename to write one file per input file
# Apr 25 2001	Allow wavelength increment to be computed
# Apr 26 2001	Fix bug which miscomputed npts if parameter INDEF
# Apr 27 2001	Add one-at-a-time file rebinning

# Mar 29 2002	Add perspec to set limits separately for each spectrum

# Apr 13 2004	Add option to use a list of output files
# Jun  4 2004	Check npts against INDEFI

# Mar 23 2005	Add ymin and ymax to scale all graphs the same
# Aug 30 2005	Add copyhead to copy header information from first spectrum

# Dec 18 2006	Pass number of input apertures to addspec

# Jan 11 2008	Get fproc from t_sumspec.x, where it is set correctly
# Jan 14 2008	Add more debugging; reopen speclist after wlrange
# Sep 16 2008	Pass unchanging number of apertures to addspec()

# Jun 10 2009	Set lcomp as soon as compfiles parameter is read
