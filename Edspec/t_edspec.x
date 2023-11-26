# File rvsao/Makespec/t_edspec.x
# May 13, 2013
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 2013 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# EDSPEC is an IRAF task for editing individual spectra within a file
# For arguments, see parameter file edspec.par.
# Information is shared in common blocks defined in "rvsao.com".
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include "rvsao.h"

procedure t_edspec ()

int	i
char	spectra[SZ_PATHNAME]	# List of input spectra
char	specfile[SZ_PATHNAME]	# Input spectrum file name
char	specdir[SZ_PATHNAME]	# Directory for input spectra
char	specpath[SZ_PATHNAME]	# Input spectrum path name
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

int	nmspec0		# Number of object multispec apertures
int	nmband0		# Number of object multispec bands
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

# Spectra to edit
	call clgstr ("spectra",spectra,SZ_PATHNAME)
	speclist = imtopen (spectra)

# Multispec spectrum apertures (use only first if multiple files)
	call clgstr ("specnum",specnums,SZ_LINE)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR)
	    call error (1, "EDSPEC: Illegal multispec list")
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }

# Multispec spectrum bands (use only first if multiple files)
	call clgstr ("specband",specbands,SZ_LINE)
	if (decode_ranges (specbands, mband_range, MAX_RANGES, nmband0) == ERR)
	    call error (1, "EDSPEC: Illegal band list")
	mband = clgeti ("specband")

# Directory containing input spectra
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }

# Number of output spectra in output file
	nspec = clgeti ("nspec")

# File to which to write summed spectrum
	call clgstr ("compfile",compfiles,SZ_PATHNAME)
	lcomp = strlen (compfiles)
	complist = imtopen (compfiles)

	save_names = FALSE
	save_names = clgetb ("save_names")

# Optional intermediate data plot switches
	pltspec = YES
	plttemp = NO

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
	specint = YES

# Type of heliocentric velocity correction to be used
	call clgstr ("svel_corr",svel_corr,SZ_LINE)
	svcor = strdic (svel_corr,svel_corr,SZ_LINE, HC_VTYPES)

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
		call printf ("EDSPEC: No more files in input list\n")
		call flush (STDOUT)
		}
	    go to endcomp_
	    }
	if (debug) {
	    call printf ("EDSPEC: next file is %s\n")
		call pargstr (specfile)
	    call flush (STDOUT)
	    }

# Get next output file from list if first spectrum or first output filled
	if (ispec == 0 || (ispec >= nspec && fproc)) {
	    if (imtgetim (complist,compfile,SZ_PATHNAME) == EOF) {
		if (debug) {
		    call printf ("EDSPEC: Cannot read image %s\n")
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
		call printf ("EDSPEC: next output file is %s\n")
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
		call error (1, "EDSPEC: Illegal multispec list")
	    }
	else {
	    nmspec = nmspec0
	    }
	nspecband = nmband0
	nspecap = nmspec
	if (debug) {
	    call printf ("EDSPEC: next file is %s [%s] = %d aps, %d bands\n")
		call pargstr (specfile)
		call pargstr (specnums)
		call pargi (nmspec)
		call pargi (nmband)
	    call flush (STDOUT)
	    }

# Check for readability of object spectrum
	call strcpy (specdir,specpath,SZ_PATHNAME)
	call strcat (specfile,specpath,SZ_PATHNAME)
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("EDSPEC: cannot read spectrum file %s \n")
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
	    call printf ("EDSPEC: Reading file %d next\n")
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

end

# May 13 2013	New task based on SUMSPEC
