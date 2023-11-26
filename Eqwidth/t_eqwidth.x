# File rvsao/Eqwidth/t_eqwidth.x
# July 20, 2009
# IRAF noao.onedspec.sbands modified to include known redshift and fix sky flux
# (Modified by Doug Mink, Smithsonian Astrophysical Observatory)

include	<error.h>
include	<smw.h>
include	<fset.h>
include "rvsao.h"
include "emv.h"
include "contin.h"
include "eqw.h"

define	SZ_HSTRING	68
define	SZ_HKWORD	9

# T_EQWIDTH -- Compute band fluxes, indices, and equivalent widths.
# A list of bandpasses with rest wavelengths is supplied in a text file,
# and all of them are applied at redshifted wavelengths to each spectrum
# in the list.  The output is written to an output file in multicolumn format.
# This program is the same as onedspec.sbands, but shifts the line list
# to match the observed redshift of each spectrum.

procedure t_eqwidth ()

pointer	speclist			# Input list of spectra
char	dbands[SZ_PATHNAME]
char	fbands[SZ_FNAME]
bool	verbose			# Verbose header?

int	nbands, nsubbands, nimages, fd
char	specfile[SZ_PATHNAME]	# Input spectrum file name
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	specdir[SZ_PATHNAME]	# Directory for object spectra
pointer spectrum		# Object spectrum
pointer specim			# Object image header structure
pointer sspec			# Sky spectrum
pointer skyim			# Sky image header structure
int	nmspec0			# Number of object multispec spectra
int	nspec			# Number of object multispec spectra
int	mspec			# Aperture of spectrum from multispec file
int	mband			# Band of object spectrum from multispec file
int	sband			# Band of sky spectrum from multispec file
int	ip,jp,lfile		# Limits for multispec aperture decoding
int	mspec_range[3,MAX_RANGES]
int	row_range[3,MAX_RANGES]
char	lbracket[3]		# "[({"
char	rbracket[3]		# "])}"
char	title_line[SZ_LINE]
char	imsec[SZ_LINE]
char	rowsec[SZ_LINE]
pointer	bands
pointer	smspec
pointer	smcspec
bool	wenable
double	dindef
double	wl
char	str[SZ_LINE]
char	rows[SZ_LINE]
int	logfiles	# List of log files
char	logfile[SZ_PATHNAME] # Log file name
char	comma
int	ldir
int	im10
int	stype
pointer	cont		# Vector for fit continuum
pointer	work		# Temporary working storage for continuum removal
pointer	wlspec		# Wavelength vector for spectrum overlap
int	i, j, nch
char	wlab[SZ_LINE]
bool	pltflat, intflat, intspec
int	nsm
int	nrows
int	nrsec
int	bin
int	fitpix
pointer	skysh

int	imaccess(), imtgetim()
bool	clgetb()
int	clgeti()
int	open()
double	clgetd()
pointer	imtopenp()
int	decode_ranges(),get_next_number(),stridx(),strlen(),stridxs()
double	wcs_p2w()
int	clpopnu(), clgfil()

include "rvsao.com"
include "results.com"
include "contin.com"
#include "emv.com"
include "eqw.com"

define  endspec_ 90

begin
	c0 = 299792.5d0
	comma = ','
	dindef = INDEFD
	specsh = NULL
	work = NULL
	cont = NULL
	fitpix = 0
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	wenable = FALSE
	conproc = CONTFIT
	rmode = 1

# Get task parameters.

# Get list of images to process
	speclist = imtopenp ("spectra")

# Multspec spectrum numbers
	call clgstr ("specnum",specnums,SZ_LINE)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR) {
	    call sprintf (str, SZ_LINE, "T_EQWIDTH: Illegal multispec object list <%s>")
		call pargstr (specnums)
	    call error (1, str)
	    }

# Multispec object spectrum band number
	mband = clgeti ("specband")

# Multispec sky spectrum band number (0 if none)
	sband = clgeti ("skyband")

# Spectrum directory 
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }

# Open log files and write a header.
	logfiles = clpopnu ("output")
	call fseti (STDOUT, F_FLUSHNL, YES)
	nlogfd = 0
	call strcpy ("rvsao.eqwidth",taskname,SZ_LINE)
	while (clgfil (logfiles, logfile, SZ_PATHNAME) != EOF) {
	    fd = open (logfile, APPEND, TEXT_FILE)
	    if (fd == ERR) next
	    nlogfd = nlogfd + 1
	    logfd[nlogfd] = fd
	    }
	if (nlogfd < 1) {
	    nlogfd = 1
	    logfd[nlogfd] = STDOUT
	    }
	call clpcls (logfiles)

	call clgstr ("bands", fbands, SZ_FNAME)
	call clgstr ("banddir", dbands, SZ_FNAME)
	readfilt = clgetb ("bandfilt")
	contname = clgetb ("bandcont")
	norm = clgetb ("normalize")
	mag = clgetb ("mag")
	magzero = clgetd ("magzero")
	bindim = clgeti ("bindim")
	byexp = clgetb ("byexp")
	bypix = clgetb ("bypix")
	net = clgetb ("netflux")
	torest = clgetb ("torest")
	if (torest)
	    call strcpy ("Rest Wavelength in Angstroms", wlab, SZ_LINE)
	else
	    call strcpy ("Wavelength in Angstroms", wlab, SZ_LINE)
	fitcont = clgetb ("fitcont")
	nsmooth = clgeti ("nsmooth")

	rmode = clgeti ("report_mode")
	pltcon = clgetb ("plot_fitcont")
	pltflat = clgetb ("plot_contsub")
	intflat = clgetb ("int_contsub")
	pltspec = clgetb ("plot_obj")
	intspec = clgetb ("int_obj")
	verbose = clgetb ("verbose")
	debug = clgetb ("debug")


# Continuum fit parameter pset
	call cont_get_pars()

	if (debug) {
	    call printf ("T_EQWIDTH: Ready to read in bands\n")
	    call flush (STDOUT)
	    }

	# Read bands from the band file.
	call eq_bands (fbands, dbands, bands, nbands, nsubbands)

	if (debug) {
	    call printf ("T_EQWIDTH: %d bands read\n")
		call pargi (nbands)
	    call flush (STDOUT)
	    }

	# Loop over the input spectra.
	nimages = 0
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
		    call sprintf (str, SZ_LINE, "EQWIDTH: Illegal multispec list <%s>")
			call pargstr (specnums)
		    call error (1, str)
		    }
		}
	    else
		nmspec = nmspec0
	    nspec = nmspec
	    if (debug) {
		call printf ("EQWIDTH: Next file is %s [%s] : %d aps left\n")
		    call pargstr (specfile)
		    call pargstr (specnums)
		    call pargi (nspec-1)
		}

	# Make spectrum pathname
	    call strcpy (specdir,specpath,SZ_PATHNAME)
	    call strcat (specfile,specpath,SZ_PATHNAME)

	# Check for readability of object spectrum
	    if (imaccess (specpath, READ_ONLY) == NO) {
		call eprintf ("EQWIDTH: cannot read spectrum path %s \n")
		    call pargstr (specpath)
		next
		}

	    mspec = -1

	# Get next multispec number from list
	    while (nspec > 0 && get_next_number (mspec_range, mspec) != EOF) {
		if (debug) {
		    call printf ("EQWIDTH: Next aperture is %s [%d] : %d aps left\n")
			call pargstr (specfile)
			call pargi (mspec)
			call pargi (nspec - mspec)
		    }

	    # Load sky spectrum for error computation
		if (sband > 0) {
		    call getspec (specpath,mspec,sband,sspec,skyim,wenable)
		    skysh = specsh
		    }
		else {
		    skysh = NULL
		    }

	    # Load spectrum
		spvel = 0.d0
		spechcv = 0.d0
		call getspec (specpath,mspec,mband,spectrum,specim,wenable)
#		if (specim == ERR)
#		    go to endspec_

	    # Compute the wavelength shift observed/rest, removing helio corr.
		if (torest) {
		    if (spvel == dindef) {
			spvel = 0.d0
			spechcv = 0.d0
			}
		    if (spechcv == dindef)
			spechcv = 0.d0
		    z1 = (1.d0 + (spvel / c0)) / (1.d0 + (spechcv / c0))
		    if (debug) {
			call printf ("EQWIDTH: Redshift z = %9.6f\n")
			    call pargd (z1)
			call flush (STDOUT)
			}
		    }
		else
		    z1 = 1.d0

	    # Allocate buffers
		if (specpix > fitpix) {
		    if (work != NULL)
			call mfree (work, TY_REAL)
		    call malloc (work, specpix, TY_REAL)
		    if (cont != NULL)
			call mfree (cont, TY_REAL)
		    call malloc (cont, specpix, TY_REAL)
		    if (smspec != NULL)
			call mfree (smspec, TY_REAL)
		    call malloc (smspec, specpix, TY_REAL)
		    if (smcspec != NULL)
			call mfree (smcspec, TY_REAL)
	 	    call malloc (smcspec, specpix, TY_REAL)
		    if (wlspec != NULL)
			call mfree (wlspec, TY_REAL)
	 	    call malloc (wlspec, specpix, TY_REAL)
		    fitpix = specpix
		    }

	    # Set up wavelength vector
		call wcs_set (specsh)
		do j = 1, specpix {
		    wl = wcs_p2w (double (j))
		    Memr[wlspec+j-1] = real (wl)
#		    call printf ("%4d: %8.3f\n")
#			call pargi (j)
#			call pargr (Memr[wlspec+j-1])
		    }

		nimages = nimages + 1
		if (debug) {
		    im10 = mod (nimages, 10)
		    if (im10 == 1)
			call printf ("EQWIDTH: %dst spectrum read\n")
		    else if (im10 == 2)
			call printf ("EQWIDTH: %dnd spectrum read\n")
		    else if (im10 == 3)
			call printf ("EQWIDTH: %drd spectrum read\n")
		    else
			call printf ("EQWIDTH: %dth spectrum read\n")
			call pargi (nimages)
		    call flush (STDOUT)
		    }

	    # Get exposure time if it is needed for normalization
		exptime = 1.d0
		if (byexp) {
		    exptime = 0.0
		    call imgdpar (specim, "EXPTIME", exptime)
		    if (exptime == 0.0)
			call imgdpar (specim, "EXPOSURE", exptime)
		    if (exptime == 0.0)
		    	exptime = 1.d0
		    }
		if (debug) {
		    call printf ("EQWIDTH: Exposure time is %6.2f sec\n")
			call pargd (exptime)
		    call flush (STDOUT)
		    }

	    # Get number of rows if it is needed for normalization
		dnrows = 1.d0
		rows[1] = EOS
		bin = 1
		if (bypix) {
		    call strcpy ("1", rows, SZ_LINE)
		    call imgspar (specim, "FINDOBJ", rows, SZ_LINE)
		    if (strlen (rows) == 0)
			nrows = 1
		    else if (decode_ranges (rows,row_range,MAX_RANGES,nrows) == ERR)
			nrows = 1
		    if (debug) {
			if (nrows > 1)
			    call printf ("EQWIDTH: rows %s = %d rows\n")
			else
			    call printf ("EQWIDTH: rows %s = %d row\n")
			    call pargstr (rows)
			    call pargi (nrows)
			}
		    if (bindim > 0) {
			call imgspar (specim, "CCDSEC", imsec, SZ_LINE)
			if (strlen (imsec) > 0) {
			    i = stridx (comma, imsec)
			    j = stridx (rbracket[1], imsec)
			    nch = j - i - 1
			    call strcpy (imsec[i+1], rowsec, nch)
		    	    if (decode_ranges (rowsec,row_range,MAX_RANGES,nrsec) != ERR) {
				if (nrsec < (bindim / 2)) {
				    dnrows = double (nrows) * 4.d0
				    bin = 4
				    }
				else if (nrsec < bindim) {
				    dnrows = double (nrows) * 2.d0
				    bin = 2
				    }
				else
				    dnrows = double (nrows)
				}
			    }
			else
			    dnrows = double (nrows)
			}
		    else
			dnrows = double (nrows)
		    }
		if (debug) {
		    if (dnrows > 1.d0)
			call printf ("EQWIDTH: %.1f (bin %d) rows added for flux\n")
		    else
			call printf ("EQWIDTH: %.1f (%d) row added for flux\n")
			call pargd (dnrows)
			call pargi (bin)
		    call flush (STDOUT)
		    }

	    # Open output file and write a verbose header if desired.
	    # It is delayed until now to avoid output if an error occurs
	    # such as image not found.
		if (nimages == 1 && (verbose || debug)) {
		    call eq_header (specid,fbands,bands,nbands,nsubbands)
		    }

	    # Set up smoothed spectrum for plotting
		call amovr (Memr[spectrum], Memr[smspec], specpix)
		call smooth (Memr[smspec], specpix, nsmooth)

		if (pltspec) {
		    if (intspec)
			nsm = nsmooth
		    else
			nsm = -nsmooth
		    call plotspec (specpix,Memr[spectrum],specid,Memr[wlspec],wlab,nsm)
		    }

	    # Fit continuum if requested
		if (fitcont) {
		    npts = specpix
		    if (debug) {
			call printf ("EQWIDTH: Continuum fit %d spectrum pixels\n")
			    call pargi (npts)
			call flush (STDOUT)
			}
		    schop = TRUE
		    stype = 1
		    do j = 1, npts {
			Memr[cont+j-1] = Memr[spectrum+j-1]
			}
		    call icsubcon (npts,Memr[cont],Memr[wlspec],specfile,stype,
				   nsmooth,Memr[work])

	    # Set up smoothed continuum-subtracted spectrum for plotting
		    call amovr (Memr[work], Memr[smcspec], specpix)
		    call smooth (Memr[smcspec], specpix, nsmooth)

		    if (pltflat) {
			call strcpy (specid, title_line, SZ_LINE)
			call strcat (" - continuum", title_line, SZ_LINE)
			do j = 1, npts {
			    Memr[work+j-1] = Memr[spectrum+j-1] - Memr[cont+j-1]
			    }
			if (intflat)
			    nsm = nsmooth
			else
			    nsm = -nsmooth
			call plotspec (specpix,Memr[work],title_line,
					Memr[wlspec],wlab,nsm)
#			call eqplot (specsh, specid, specim, specpix,
#				     Memr[spectrum], Memr[smspec], Memr[work],
#				     Memr[smcspec], Memr[wlspec], bands, nbands)
			}

		    if (debug) {
			call printf ("EQWIDTH: Continuum fit %d spectrum pixels\n")
			        call pargi (npts)
			call flush (STDOUT)
			}
		    }

	    # Compute equivalent widths for this spectrum
		call eq_proc (specsh, Memr[spectrum], skysh, Memr[sspec],
			      Memr[wlspec], specpix, mspec, nmspec, 
			      Memr[cont], bands, nbands)

	    # Close the object spectrum image and headers
endspec_
		call close_image (specim, specsh)
		if (sband > 0)
		    call close_image (skyim, skysh)

	    # move on to next aperture or next image
		nspec = nspec - 1

	    # End of multispec loop within single image
		}

	# End of image loop
	    }

# Finish up.
	call eq_free (bands, nbands)
	do i = 1, nlogfd {
	    call close (logfd[i])
	    }
#	call imtclose (speclist)
	if (work != NULL)
	    call mfree (work, TY_REAL)
	if (cont != NULL)
	    call mfree (work, TY_REAL)
	if (smspec != NULL)
	    call mfree (smspec, TY_REAL)
	if (smcspec != NULL)
	    call mfree (smcspec, TY_REAL)
	if (wlspec != NULL)
	    call mfree (wlspec, TY_REAL)
end


# EQ_BANDS - Read bands from the band file and put them into an array
# of band pointers.

procedure eq_bands (fbands, dbands, bands, nbands, nsubbands)

char	fbands[ARB]		#I File containing band information
char	dbands[ARB]		#I Directory containing band information
pointer	bands			#O Bandpass table descriptor
int	nbands			#O Number of bandpasses
int	nsubbands		#O Number of individual bands

bool	bandok, cdone
int	ip, lb
int	fd			#I Bandpass file descriptor
double	center, width, w1, w2
char	line[SZ_LINE]
char	ffile[SZ_FNAME]
char	bandid[SZ_FNAME]
char	bandpath[SZ_PATHNAME]
char	id[SZ_LINE]

int	getline(), ctowrd(), ctod()
int	open(), strlen()

include "eqw.com"
include "rvsao.com"

begin
	call strcpy ("none", ffile, SZ_FNAME)
	bandid[1] = EOS
	call strcpy (dbands, bandpath, SZ_PATHNAME)
	lb = strlen (bandpath)
	if (lb < 1)
	    call strcpy (fbands, bandpath, SZ_PATHNAME)
	else if (dbands[lb] != '/') {
	    call strcat ("/", bandpath, SZ_PATHNAME)
	    call strcat (fbands, bandpath, SZ_PATHNAME)
	    }
	else
	    call strcat (fbands, bandpath, SZ_PATHNAME)

	if (debug) {
	    call printf ("EQ_BANDS: Reading band information from %s\n")
		call pargstr (bandpath)
	    }

	fd = open (bandpath, READ_ONLY, TEXT_FILE)

	# Read the bands.  If the first band is not seen
	# skip the line.  Check for 1, 2, or 3 bandpasses.
	# Can't use fscan() because fscan() will be called later to
	# read any filter file.

	bands = NULL
	nbands = 0
	nsubbands = 0
	while (getline (fd, line) != EOF) {
	    if (debug) {
		call printf ("EQ_BANDS: %s\n")
		    call pargstr (line)
		}
	    if (line[1] == '#')
		next

# Read first entry in the table
	    ip = 1
	    bandok = (ctowrd (line, ip, id, SZ_FNAME) > 0)
	    bandok = (bandok && ctod (line, ip, w1) > 0)
	    bandok = (bandok && ctod (line, ip, w2) > 0)
	    if (readfilt)
		bandok = (bandok && ctowrd (line,ip,ffile,SZ_FNAME)>0)
	    if (!bandok || id[1] == '#')
		next

# Allocate and reallocate the array of band pointers.
	    if (nbands == 0)
		call malloc (bands, 10 * NBANDS, TY_POINTER)
	    else if (mod (nbands, 10) == 0)
		call realloc (bands, (nbands + 10) * NBANDS, TY_POINTER)
	    nbands = nbands + 1
	    call strcpy (id, bandid, SZ_FNAME)

# Convert from center and width to limits
	    cdone = FALSE
	    if (w2 < w1) {
		cdone = TRUE
		center = w1
		width = w2
		w1 = center - (width / 2.d0)
		w2 = center + (width / 2.d0)
		call strcpy (bandid, id, SZ_FNAME)
		call eq_alloc (BAND(bands,nbands,1), id, ffile, w1, w2)
		}
	    else {
		call strcpy ("c1", id, SZ_FNAME)
		call strcat (bandid, id, SZ_FNAME)
		call eq_alloc (BAND(bands,nbands,2), id, ffile, w1, w2)
		}
	    nsubbands = nsubbands + 1

	    if (debug) {
		call printf ("EQ_BANDS: %d %d %s %6.1f - %6.1f A\n");
		    call pargi (nbands)
		    call pargi (nsubbands)
		    call pargstr (bandid)
		    call pargd (w1)
		    call pargd (w2)
		call flush (STDOUT)
		}

	    if (contname)
		bandok = (ctowrd (line, ip, id, SZ_FNAME) > 0)
	    else {
		if (cdone) {
		    call strcpy ("c1", id, SZ_FNAME)
		    call strcat (bandid, id, SZ_FNAME)
		    }
		else {
		    call strcpy (bandid, id, SZ_FNAME)
		    }
		bandok = TRUE
		}
	    bandok = (bandok && ctod (line, ip, w1) > 0)
	    bandok = (bandok && ctod (line, ip, w2) > 0)
	    if (readfilt)
		bandok = (bandok && ctowrd (line,ip,ffile,SZ_FNAME)>0)
	    if (bandok) {
		if (w2 < w1) {
		    center = w1
		    width = w2
		    w1 = center - (width / 2.0d0)
		    w2 = center + (width / 2.0d0)
		    }
		if (cdone)
		    call eq_alloc (BAND(bands,nbands,2), id, ffile, w1, w2)
		else
		    call eq_alloc (BAND(bands,nbands,1), id, ffile, w1, w2)
		nsubbands = nsubbands + 1

		if (debug) {
		    call printf ("EQ_BANDS: %d %d %6.1f - %6.1f A\n");
			call pargi (nbands)
			call pargi (nsubbands)
			call pargd (w1)
			call pargd (w2)
		    call flush (STDOUT)
		    }
		}
	    else
		BAND(bands,nbands,2) = NULL

	    if (contname)
		bandok = (ctowrd (line, ip, id, SZ_FNAME) > 0)
	    else {
		call strcpy ("c2", id, SZ_FNAME)
		call strcat (bandid, id, SZ_FNAME)
		bandok = TRUE
		}
	    bandok = (bandok && ctod (line, ip, w1) > 0)
	    bandok = (bandok && ctod (line, ip, w2) > 0)
	    if (readfilt)
		bandok = (bandok && ctowrd (line,ip,ffile,SZ_FNAME)>0)
	    if (bandok) {
		if (w2 < w1) {
		    center = w1
		    width = w2
		    w1 = center - (width / 2.0d0)
		    w2 = center + (width / 2.0d0)
		    }
		call eq_alloc (BAND(bands,nbands,3),id,ffile,w1, w2)
		nsubbands = nsubbands + 1

		if (debug) {
		    call printf ("EQ_BANDS: %d %d %6.1f - %6.1f A\n");
			call pargi (nbands)
			call pargi (nsubbands)
			call pargd (w1)
			call pargd (w2)
		    call flush (STDOUT)
		    }
		}
	    else
		BAND(bands,nbands,3) = NULL

	    if (BAND(bands,nbands,3) == NULL) {
		if (BAND(bands,nbands,2) != NULL) {
		    if (debug) {
			call printf ("EQ_BANDS: %d No second continuum region specified\n")
			    call pargi (nbands)
			}
		    }
		else {
		    fitcont = TRUE
		    if (debug && nbands == 1)
			call printf ("EQ_BANDS: no continuum region specified\n")
		    }
		}
	    }

	call close (fd)
end


# EQ_ALLOC -- Allocate a band structure.

procedure eq_alloc (band, id, ffile, w1, w2)

pointer	band			#O Band pointer
char	id[ARB]			#I Band id
char	ffile[ARB]		#I Band filter file name
double	w1			#I Band lower wavelength limit
double	w2			#I Band upper wavelength limit

int	fn, fd, strlen(), open(), fscan(), nscan()
double	w, r
pointer	fw, fr
bool	streq()
errchk	open()

begin
	call calloc (band, LEN_BAND, TY_STRUCT)
	call malloc (BAND_ID(band), strlen(id), TY_CHAR)
	call malloc (BAND_FILTER(band), strlen(ffile), TY_CHAR)
	call strcpy (id, Memc[BAND_ID(band)], ARB)
	call strcpy (ffile, Memc[BAND_FILTER(band)], ARB)
	BAND_W1(band) = w1
	BAND_W2(band) = w2
	BAND_FN(band) = 0
	BAND_FW(band) = NULL
	BAND_FR(band) = NULL

	if (streq (ffile, "none"))
	    return

	# Read the filter file.
	fd = open (ffile, READ_ONLY, TEXT_FILE)
	fn = 0
	while (fscan (fd) != EOF) {
	    call gargd (w)
	    call gargd (r)
	    if (nscan() != 2)
		next
	    if (fn == 0) {
		call malloc (fw, 100, TY_DOUBLE)
		call malloc (fr, 100, TY_DOUBLE)
	    } else if (mod (fn, 100) == 0) {
		call realloc (fw, fn+100, TY_DOUBLE)
		call realloc (fr, fn+100, TY_DOUBLE)
	    }
	    Memd[fw+fn] = w
	    Memd[fr+fn] = r
	    fn = fn + 1
	}
	call close (fd)

	BAND_FN(band) = fn
	BAND_FW(band) = fw
	BAND_FR(band) = fr
	return
end


# EQ_FREE -- Free band structures.

procedure eq_free (bands, nbands)

pointer	bands			#I bands descriptor
int	nbands			#I number of bands

int	i, j
pointer	band

begin
	do i = 1, nbands {
	    do j = 1, NBANDS {
		band = BAND(bands,i,j)
		if (band != NULL) {
		    call mfree (BAND_ID(band), TY_CHAR)
		    call mfree (BAND_FILTER(band), TY_CHAR)
		    if (BAND_FN(band) > 0) {
			call mfree (BAND_FW(band), TY_DOUBLE)
			call mfree (BAND_FR(band), TY_DOUBLE)
			}
		    call mfree (band, TY_STRUCT)
		}
	    }
	}
	call mfree (bands, TY_POINTER)
	return
end


# EQ_HEADER -- Print output header.

procedure eq_header (spec, fbands, bands, nbands, nsubbands)

char	spec[ARB]		#I Typical spectrum name
char	fbands[ARB]		#I Band file
pointer	bands			#I Pointer to array of bands
int	nbands			#I Number of bands
int	nsubbands		#I Number of subbands

int	i, j, k, lspec, iband
pointer	band
char	str[SZ_LINE]

int	strlen()

include "eqw.com"
include "rvsao.com"

begin
	if (debug) {
	    call printf ("EQ_HEADER: printing headers for %d bands to")
		call pargi (nbands)
	    do k = 1, nlogfd {
		call printf (" %d")
		    call pargi (logfd[k])
		}
	    call printf ("\n")
	    call flush (STDOUT)
	    }

	# Output a banner and task parameters.
	call sysid (str, SZ_LINE)
	if (debug) {
	    call printf ("EQ_HEADER: %s\n")
		call pargi (str)
	    call flush (STDOUT)
	    }

	do k = 1, nlogfd {
	    call fprintf (logfd[k], "\n# EQWIDTH: %s\n#  ")
		call pargstr (str)
	    if (fbands[1] != EOS) {
		call fprintf (logfd[k], " bands = %s,")
		    call pargstr (fbands)
		}
	    call fprintf (logfd[k], " norm = %b, mag = %b")
		call pargb (norm)
		call pargb (mag)
	    if (mag) {
		call fprintf (logfd[k], ", magzero = %.2f")
		    call pargd (magzero)
		call strcpy ("mag", str, SZ_LINE)
		}
	    else
		call strcpy ("flux", str, SZ_LINE)
	    if (torest)
		call fprintf (logfd[k], " torest = yes")
	    if (byexp)
		call fprintf (logfd[k], " byexp = yes")
	    if (bypix)
		call fprintf (logfd[k], " bypix = yes")

	    # Output the bands.
	    call fprintf (logfd[k], "\n# bno           band     filter wavelength1  wavelength2\n")
	    do iband = 1, nbands {
		do j = 1, NBANDS {
		    band = BAND(bands,iband,j)
		    if (band == NULL)
			next
		    call fprintf (logfd[k], "# %3d %14s %10s %10.4f - %10.4f\n")
			call pargi (iband)
			call pargstr (Memc[BAND_ID(band)])
			call pargstr (Memc[BAND_FILTER(band)])
			call pargd (BAND_W1(band))
			call pargd (BAND_W2(band))
		    }
		}
	    call fprintf (logfd[k], "#\n")

	    # Output column headings.

	    # Print indices and equivalent widths
	    lspec = strlen (spec)
	    call sprintf (fspec, 16, "%%%ds")
		call pargi (lspec)
	    call fprintf (logfd[k], fspec)
		call pargstr ("spectrum")
	    do iband = 1, nbands {
		if (rmode == 1) {
		    call fprintf (logfd[k],"\tband%d\tflux%d\tcont%d\tcflux%d\tindex%d\teqwidth%d")
			call pargi (iband)
			call pargi (iband)
			call pargi (iband)
			call pargi (iband)
			call pargi (iband)
			call pargi (iband)
		    }

	    # Print indices and errors
		else if (rmode > 1) {
		    band = BAND(bands,iband,1)
		    if (rmode < 4) {
			call fprintf (logfd[k],"\tflux%s\tcflux%s")
			    call pargstr (Memc[BAND_ID(band)])
			    call pargstr (Memc[BAND_ID(band)])
			}
		    if (rmode == 2 || rmode == 4) {
			call fprintf (logfd[k],"\tindex%s\terr%s")
			    call pargstr (Memc[BAND_ID(band)])
			    call pargstr (Memc[BAND_ID(band)])
			}
		    else if (rmode == 3 || rmode == 5) {
			call fprintf (logfd[k],"\teqw%s\terr%s")
			    call pargstr (Memc[BAND_ID(band)])
			    call pargstr (Memc[BAND_ID(band)])
			}
		    }
		}
	    if (torest) {
		call fprintf (logfd[k], "\t%7s\n")
		    call pargstr ("1+z")
		}
	    else
		call fprintf (logfd[k], "\n")

	    # Print underlines for Starbase table
	    do i = 1, lspec {
		call fprintf (logfd[k], "-")
		}
	    do iband = 1, nbands {
		if (rmode == 1)
		    call fprintf (logfd[k],"\t-------\t--------")
		if (rmode < 3)
		    call fprintf (logfd[k],"\t-------\t--------")
		call fprintf (logfd[k],"\t--------\t--------")
		}
	    if (torest) {
		call fprintf (logfd[k], "\t-------\n")
		}
	    else
		call fprintf (logfd[k], "\n")
	    }
	return
end


# EQ_PROC -- Measure the band fluxes and possibly a band index and eq. width.

procedure eq_proc (sh,spec,skysh,sspec,wl,npix,ispec,nspec,cont,bands,nbands)

pointer	sh			#I Object spectrum descriptor
real	spec[ARB]		#I Object spectrum
pointer	skysh			#I Sky spectrum descriptor
real	sspec[ARB]		#I Sky spectrum
real	wl[ARB]			#I Wavelengths for spectrum
int	npix			#I Number of pixels in spectrum
int	ispec			#I number of spectrum in image
int	nspec			#I number of spectra in image
real	cont[ARB]		#I Object continuum vector, if one has been fit
pointer	bands			#I Bandpass table pointer
int	nbands			#I Number of bandpasses

char	imname[SZ_FNAME]
int	j, api, iband
double	flux, contval, index, eqwidth, w1, w2, dw1
double	flux1, norm1, flux2, norm2, flux3, norm3, a, b, fluxnet, zero
double	sflux1, snorm1, sflux2, snorm2, sflux3, snorm3, index1, index2
double	sigma1, sigma2, sigma3, sigma, inderr, eqwerr
#double	noise
double	wc1, wc2, wc3, fluxratio, fluxrat1, fluxrat2
pointer	band1, band2, band3
real	rindef, wl1, wl2
bool	skystat

include "eqw.com"
include "rvsao.com"

begin
	rindef = INDEF
	api = AP(sh)
	wl1 = wl[1]
	wl2 = wl[2]

	if (nspec < 2)
	    call strcpy (IMNAME(sh), imname, SZ_FNAME)
	else {
	    call sprintf (imname, SZ_FNAME, "%s[%d]")
		call pargstr (IMNAME(sh))
		call pargd (ispec)
	    }

	if (debug) {
	    call printf ("EQ_PROC: About to process %s(%d)\n")
		call pargstr (imname)
		call pargi (api)
	    call flush (STDOUT)
	    }

	zero = 0.0d0
	if (skysh != NULL)
	    skystat = TRUE
	else
	    skystat = FALSE

	if (debug) {
	    call printf ("EQ_PROC: About to compute %d bands in %s\n")
		call pargi (nbands)
		call pargstr (imname)
	    call flush (STDOUT)
	    }
	do j = 1, nlogfd {
	    call fprintf (logfd[j], fspec)
		call pargstr (imname)
	    }

	# Loop over all bandpasses
	do iband = 1, nbands {

	    # Measure primary band flux, normalize, and print result.
	    band1 = BAND(bands,iband,1)
	    w1 = BAND_W1(band1)
	    w2 = BAND_W2(band1)
	    dw1 = w2 - w1
	    wc1 = (w1 + w2) * 0.5d0
	    if (debug) {
		call printf ("\nEQPROC: Band %d_1 %7.7s %11.4f-%11.4f (%.4f)\n")
		    call pargi (iband)
		    call pargstr (Memc[BAND_ID(band1)])
		    call pargd (w1)
		    call pargd (w2)
		    call pargd (wc1)
		call flush (STDOUT)
		}

	    # Measure primary band flux, normalize, and print result.
	    call eq_flux (sh,spec,wl,npix,w1,w2,band1,flux1,norm1,debug)
	    if (IS_INDEFD(flux1))
		flux1 = 0.d0

	    # Divide out exposure time
	    if (byexp && exptime > 0.d0)
		flux1 = flux1 / exptime

	    # Measure primary band sky flux, normalize, and print result.
	    if (skysh != NULL) {
		call eq_flux (skysh,sspec,wl,npix,w1,w2,band1,sflux1,snorm1,debug)
		if (IS_INDEFD(sflux1))
		    sflux1 = 0.d0

		# Divide out exposure time
		if (byexp && exptime > 0.d0)
		    sflux1 = sflux1 / exptime

		if (flux1 > 0.d0)
		    sigma1 = flux1 - sflux1
		else
		    sigma1 = 0.d0
		}
	    else {
		sflux1 = 0.d0
		if (flux1 > 0.d0)
		    sigma1 = flux1
		else 
		    sigma1 = 0.d0
		}

	    # Divide out number of rows for skies
	    if (bypix && dnrows > 0.0) 
		flux1 = flux1 / dnrows

	    # Divide out filter-weighted number of pixels
	    if (norm) {
		if (norm1 != 0.0) {
		    flux1 = flux1 / norm1
		    sigma1 = sigma1 / norm1
		    }
		norm1 = 1.d0
		}

	    # Convert to magnitude if requested
	    if (mag && flux1 > 0.)
		flux = magzero - 2.5 * log10 (flux1)
	    else
		flux = flux1

	    do j = 1, nlogfd {
		if (rmode == 1) {
		    call fprintf (logfd[j], "	%7.7s")
			call pargstr (Memc[BAND_ID(band1)])
		    }
		if (rmode < 4) {
		    call fprintf (logfd[j], "	%11.6g")
			call pargd (flux)
		    }
		}
	    if (debug)
		call printf ("\n")

	    # Measure continuum in same place if it has been fit
	    if (fitcont) {
		call eq_flux (sh,cont,wl,npix,w1,w2,band1,flux2,norm2,debug)

		# Print zeroes and go to next band if no continuum is specified
		if (IS_INDEFD(flux2)) {
		    do j = 1, nlogfd {
			if (rmode ==1 ) {
			    call fprintf (logfd[j], "	%7.7s")
				call pargstr ("cont")
			    }
			if (rmode < 4) {
			    call fprintf (logfd[j], " %11.6g")
				call pargd (zero)
			    call fprintf (logfd[j], " %9.6g %9.6g")
				call pargd (zero)
				call pargd (zero)
			    }
			else {
			    call fprintf (logfd[j], " %9.6g %9.6g")
				call pargd (zero)
				call pargd (zero)
			    }
			}
		    next
		    }
		else if (flux2 > 0.d0) {

		    # Divide out exposure time
		    if (byexp && exptime > 0.d0)
			flux2 = flux2 / exptime
		    sigma2 = flux2 + sflux1
		    if (norm && norm2 > 0)
			flux2 = flux2 / norm2
		    sigma3 = 1.d0
		    }
		else {
		    sigma2 = 1.d0
		    sigma3 = 1.d0
		    }
		}

	    # Otherwise measure it in the indicated places
	    else {

		# Measure the first continuum band object flux
		band2 = BAND(bands,iband,2)
		w1 = BAND_W1(band2)
		w2 = BAND_W2(band2)
		wc2 = (w1 + w2) * 0.5d0
		if (debug) {
		    call printf ("EQPROC: Band %d_2 %.4f - %.4f (%.4f)\n")
			call pargi (iband)
			call pargd (w1)
			call pargd (w2)
			call pargd (wc2)
			}
		call eq_flux (sh,spec,wl,npix,w1,w2,band2,flux2,norm2,debug)

		# Divide out exposure time
		if (byexp && exptime > 0.d0)
		    flux2 = flux2 / exptime

		# Measure the first continuum band sky flux
		if (skystat) {
			call eq_flux (skysh,sspec,wl,npix,w1,w2,band2,sflux2,snorm2,debug)

		    # Divide out exposure time
		    if (byexp && exptime > 0.d0)
			sflux2 = sflux2 / exptime

		    if (IS_INDEFD(sflux2))
			sflux2 = zero
		    if (flux2 > 0)
			sigma2 = flux2 + sflux2
		    else
			sigma2 = 1.d0
		    }
		else {
		    sflux2 = zero
		    if (flux2 > 0.d0)
			sigma2 = flux2
		    else
			sigma2 = 1.d0
		    }
		if (norm && norm2 > 0)
		    sigma2 = sigma2 / norm2

		# Measure the second continuum band flux
		band3 = BAND(bands,iband,3)
		if (band3 == NULL) {
		    flux3 = INDEFD
		    sigma3 = 1.0
		    }
		else {
		    w1 = BAND_W1(band3)
		    w2 = BAND_W2(band3)
		    wc3 = (w1 + w2) * 0.5d0
		    if (debug) {
			call printf ("\nEQPROC: Band %d_3 %.4f - %.4f (%.4f)\n")
			    call pargi (iband)
			    call pargd (w1)
			    call pargd (w2)
			    call pargd (wc3)
			call flush (STDOUT)
			}
		    call eq_flux (sh,spec,wl,npix,w1,w2,band3,flux3,norm3,debug)

		    # Divide out exposure time
		    if (byexp && exptime > 0.d0)
			flux3 = flux3 / exptime

		    # Measure the second continuum band sky flux
		    if (skystat) {
			call eq_flux (skysh,sspec,wl,npix,w1,w2,band3,sflux3,snorm3,debug)

			# Divide out exposure time
			if (byexp && exptime > 0.d0)
			    sflux3 = sflux3 / exptime

			if (IS_INDEFD(sflux3))
			    sflux3 = zero
			if (flux3 > 0.d0)
			    sigma3 = flux3 + sflux3
			else
			    sigma3 = 1.d0
			}
		    else {
			sflux3 = zero
			if (flux3 > 0.d0)
			    sigma3 = flux3
			else
			    sigma3 = 1.d0
			}
		    if (norm && norm3 > 0)
			sigma3 = sigma3 / norm3
		    }

	    if (debug) {
		call printf ("\nEQPROC: Band %d %.4f, %.4f + %.4f (%.4f, %.4f, %.4f)\n")
		    call pargi (iband)
		    call pargd (wc1)
		    call pargd (wc2)
		    call pargd (wc3)
		    call pargd (norm1)
		    call pargd (norm2)
		    call pargd (norm3)
		call printf ("EQPROC: Band %d sigmas= %.4f, %.4f, %.4f\n")
		    call pargi (iband)
		    call pargd (sigma1)
		    call pargd (sigma2)
		    call pargd (sigma3)
		call flush (STDOUT)
		}

		# Print zeroes and go to next band if no continuum is specified
		if ((IS_INDEFD(flux2)) && (IS_INDEFD(flux3))) {
		    do j = 1, nlogfd {
			if (rmode ==1 ) {
			    call fprintf (logfd[j], "	%7.7s")
				call pargstr ("cont")
			    }
			if (rmode < 4) {
			    call fprintf (logfd[j], " %11.6g")
				call pargd (zero)
			    call fprintf (logfd[j], " %9.6g %9.6g")
				call pargd (zero)
				call pargd (zero)
			    }
			else {
			    call fprintf (logfd[j], "	%9.6g	%9.6g")
				call pargd (zero)
				call pargd (zero)
			    }
			}
		    next
		    }

		# Compute and output the band index and equivalent width.
		if (net || norm) {
		    if (!IS_INDEFD(flux2)) {
			if (norm2 != 0.0)
			    flux2 = flux2 / norm2
			norm2 = 1.d0
			}
		    if (!IS_INDEFD(flux3)) {
			if (norm3 != 0.0)
			    flux3 = flux3 / norm3
			norm3 = 1.d0
			}
		    }

		else if (IS_INDEFD(flux2))
		    flux2 = 0.0
		else if (IS_INDEFD(flux3))
		    flux3 = 0.0
		}

	    contval = INDEFD
	    index = INDEFD
	    eqwidth = INDEFD

	    if (debug) {
		call printf ("\nEQPROC: Band %d %.4f, %.4f + %.4f (%.4f, %.4f, %.4f)\n")
		    call pargi (iband)
		    call pargd (wc1)
		    call pargd (wc2)
		    call pargd (wc3)
		    call pargd (norm1)
		    call pargd (norm2)
		    call pargd (norm3)
		call flush (STDOUT)
		}

	    # Use underlying fit continuum
	    if (fitcont) {
		contval = flux2
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr ("cont")
			}
		    }
		sigma = 1.d0 / dsqrt (sigma1 + sigma2)
		}

	    # Interpolate to the center of the primary band.
	    else if (!IS_INDEFD(flux2) && !IS_INDEFD(flux3)) {
		a = ((flux2 / norm2) - (flux3 / norm3)) /
		    (wc2 - wc3)
		b = (flux2 / norm2) - (a * wc2)
		contval = (a * wc1 + b) * norm1
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr ("cont")
			}
		    }
		if (sigma2 > 0.d0 && sigma3 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1 + sigma2 + sigma3)
		else if (sigma2 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1 + sigma2)
		else if (sigma3 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1 + sigma3)
		else if (sigma1 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1)
		else
		    sigma = 0.d0
		}

	    # Use first continuum if second is not computable
	    else if (!IS_INDEFD(flux2)) {
		contval = 0.d0
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr (Memc[BAND_ID(band2)])
			}
		    }
		if (sigma2 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1 + sigma2)
		else if (sigma1 > 0.d0)
		    sigma = 1.d0 / dsqrt (dsqrt (sigma1))
		else
		    sigma = 0.d0
		}

	    # Use second continuum if first is not computable
	    else if (!IS_INDEFD(flux3)) {
		contval = 0.d0
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr (Memc[BAND_ID(band3)])
			}
		    }
		if (sigma3 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1 + sigma3)
		else if (sigma1 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1)
		else
		    sigma = 0.d0
		}

	    else {
		contval = 0.d0
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr ("cont")
			}
		    }
		if (sigma1 > 0.d0)
		    sigma = 1.d0 / dsqrt (sigma1)
		else
		    sigma = 0.d0
		}

	    if (mag && contval > 0.)
		flux = magzero - 2.5 * log10 (contval)
	    else
		flux = contval
	    if (rmode < 4) {
		do j = 1, nlogfd {
		    call fprintf (logfd[j], "	%11.6g")
			call pargd (flux)
		    }
		}

	    if (debug) {
		call printf ("EQPROC: flux1=%.2f flux2=%.2f flux3=%.2f\n")
		    call pargd (flux1)
		    call pargd (flux2)
		    call pargd (flux3)
		call flush (STDOUT)
		}

# If line or either continuum flux is INDEF, skip this line
	    if (IS_INDEFD (flux1) || IS_INDEFD (flux1) || IS_INDEFD (flux1)) {
		fluxratio = INDEFD
		eqwidth = 0.d0
		eqwerr = 0.d0
		fluxnet = flux1
		}

# Skip this if any of the fluxes is less than zero
	    else if (flux1 < 0.d0 || flux2 < 0.d0 || flux3 < 0.d0) {
		fluxratio = INDEFD
		eqwidth = 0.d0
		eqwerr = 0.d0
		fluxnet = flux1
		}
	    else if (flux1>0. && contval>0. && norm1>0.d0 && norm2>0.d0 && norm3>0.d0) {
		if (norm || net)
		    fluxratio = flux1 / (0.5d0 * (flux2 + flux3))
		else
		    fluxratio = (flux1/norm1)/(0.5d0*((flux2/norm2)+(flux3/norm3)))
		eqwidth = (1.d0 - fluxratio) * dw1
		eqwerr = sigma * eqwidth
		fluxnet = flux1 - contval
		if (eqwerr < 0.d0)
		    eqwerr = -eqwerr
		}
	    else if (flux1 > 0.) {
		fluxratio = INDEFD
		eqwidth = 0.d0
		eqwerr = 0.d0
		fluxnet = flux1
		}
	    else {
		fluxratio = INDEFD
		eqwidth = 0.0
		eqwerr = 0.d0
		fluxnet = 0.0
		}
	    if (mag) {
		if (!IS_INDEFD(contval) && contval > 0.)
		    contval = magzero - 2.5 * log10 (contval)
		}

	    if (debug) {
		call printf ("EQPROC: net flux=%.2f fluxratio=%.2f sigma=%.4f eqwidth=%.3f (%.3f)\n")
		    call pargd (fluxnet)
		    call pargd (fluxratio)
		    call pargd (sigma)
		    call pargd (eqwidth)
		    call pargd (eqwerr)
		call flush (STDOUT)
		}

# Compute 
	    if (!IS_INDEFD(fluxratio)) {
		index = -2.5 * log10 (fluxratio)
		fluxrat1 = fluxratio * (1.d0 - sigma)
		if (fluxrat1 < 1.d0)
		    fluxrat1 = 1.d0
		fluxrat2 = fluxratio * (1.d0 + sigma)
		if (fluxrat2 < 1.d0)
		    fluxrat2 = 1.d0
		index1 = -2.5 * log10 (fluxrat1)
		index2 = -2.5 * log10 (fluxrat2)
		inderr = 0.5 * (index2 - index1)
		if (inderr < 0.d0)
		    inderr = -inderr
		}
	    else {
		index = 0.0
		inderr = 0.0
		}

	    if (net) {
		do j = 1, nlogfd {
		    call fprintf (logfd[j], "\t%8.6g")
			call pargd (fluxnet)
		    }
		}
	    else {
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "\t%8.6g\t%8.6g")
			    call pargd (index)
			    call pargd (eqwidth)
			}
		    else if (rmode == 2 || rmode == 4) {
			call fprintf (logfd[j], "\t%8.6g\t%8.6g")
			    call pargd (index)
			    call pargd (inderr)
			}
		    else {
			call fprintf (logfd[j], "\t%8.6g\t%8.6g")
			    call pargd (eqwidth)
			    call pargd (eqwerr)
			}
		    }
		}
	    }

# Print redshift as 1+z, if shifted, and flush output buffers
	do j = 1, nlogfd {
	    if (torest) {
		call fprintf (logfd[j], "\t%7.5f")
		    call pargd (z1)
		}
	    call fprintf (logfd[j], "\n")
	    call flush (logfd[j])
	    }

	return
end


# EQ_FLUX - Compute the flux and total response in a given band.
# Return INDEF if the band is outside of the spectrum.

procedure eq_flux (sh, spec, wl, npix, w1, w2, band, flux, normb, debug)

pointer	sh		#I Spectrum descriptor
real	spec[ARB]	#I spectrum
real	wl[ARB]		#I wavelengths for spectrum
int	npix		#I Number of pixels in spectrum
double	w1		#I lower wavelength for band
double	w2		#I upper wavelength for band
pointer	band		#I band descriptor
double	flux		#O flux
double	normb		#O normalization
bool	debug

int	i, i1, i2
bool	noflux
double	a, b, r1, r2, di1, di2, wt, wl1, wl2
double	eq_filter(), shdr_wl()

include "eqw.com"

begin
	flux = INDEFD
	noflux = FALSE
	normb = 1.0
	if (band == NULL)
	    return

	wl1 = wl[1]
	wl2 = wl[npix]
	if (wl1 > wl2) {
	    wl2 = wl[1]
	    wl1 = wl[npix]
	    }

	# Redshift wavelength limits to get correct pixel limits
	r1 = w1 * z1
	r2 = w2 * z1
	if (r2 < r1) {
	    r2 = w1 * z1
	    r1 = w2 * z1
	    }
	if (debug) {
#	    call printf ("\n")
            call printf ("EQFLUX: Rest %.3f - %.3f -> Obs %.3f-%.3f at %.3f\n")
                call pargd (w1)
                call pargd (w2)
                call pargd (r1)
                call pargd (r2)
                call pargd (z1)
            call printf ("EQFLUX: Spectrum: %.3f - %.3f\n")
                call pargd (wl1)
                call pargd (wl2)
	    call flush (STDOUT)
	    }

	if (r1 < wl1) {
	    if (debug) {
		call printf ("EQFLUX: Obs %.3f < Spectrum %.3f\n")
		    call pargd (r1)
		    call pargd (wl1)
		}
	    return
	    }

	if (r2 > wl2) {
	    if (debug) {
		call printf ("EQFLUX: Obs %.3f > Spectrum %.3f\n")
		    call pargd (r2)
		    call pargd (wl2)
		}
	    return
	    }

	a = shdr_wl (sh, r1)
	b = shdr_wl (sh, r2)
	di1 = min (a, b)
	di2 = max (a, b)
	i1 = nint (di1)
	i2 = nint (di2)
	if (debug) {
	    call printf ("EQFLUX: Rest %.3f - %.3f -> Obs %.3f-%.3f = pixels %.3f-%.3f\n")
		call pargd (w1)
		call pargd (w2)
		call pargd (r1)
		call pargd (r2)
		call pargd (di1)
		call pargd (di2)
	    if (i1 > npix) {
		call printf (" off spectrum\n")
		noflux = TRUE
		}
	    else if (i2 < 1) {
		call printf (" off spectrum\n")
		noflux = TRUE
		}
	    else if (i2 > npix) {
		call printf (" > %d\n")
		    call pargi (npix)
		}
	    else if (i1 < 1) {
		call printf (" < 1")
		}
	    call flush (STDOUT)
	    }
	if (r1 > wl2) {
	    flux = INDEFD
	    noflux = TRUE
	    }
	if (r2 < wl1) {
	    flux = INDEFD
	    noflux = TRUE
	    }
	if (i1 < 1) {
	    i1 = 1
	    di2 = 1.d0
	    }
	if (i2 > npix) {
	    i2 = npix
	    di2 = double (npix)
	    }
	if (noflux) {
#	    call printf ("\n")
#	    if (debug) {
#		call printf (": 0.000000 (%9.6g)\n")
#		    call pargd (normb)
#		call flush (STDOUT)
#		}
	    return
	    }

# If only one pixel is involved
	if (i1 == i2) {
	    wt = eq_filter (wl[i1], band) * (di2 - di1)
	    flux = wt * spec[i1]
	    normb = wt
	    }

# Otherwise, use appropriate fractions of first and last pixels
	else {

# First pixel
	    wt = eq_filter (wl[i1], band) * (double(i1) + 0.5d0 - di1)
	    flux = wt * spec[i1]
	    normb = wt

# Middle pixels
	    do i = i1+1, i2-1, 1 {
		wt = eq_filter (wl[i], band)
		flux = flux + wt * spec[i]
		normb = normb + wt
		}

# Last pixel
	    wt = eq_filter (double (wl[i2]), band) * (di2 - double (i2) + 0.5d0)
	    flux = flux + wt * spec[i2]
	    normb = normb + wt
	    }

	if (debug) {
#	    call printf ("\nEQFLUX: Rest %7.2f-%7.2f -> Obs %7.2f-%7.2f -> Pix %d-%d: %9.6g (%9.6g)\n")
#		call pargd (w1)
#		call pargd (w2)
#		call pargd (r1)
#		call pargd (r2)
#		call pargi (i1)
#		call pargi (i2)
	    call printf (": %9.6g (%9.6g)\n")
		call pargd (flux)
		call pargd (normb)
	    call flush (STDOUT)
	    }
	return
end


# EQ_FILTER -- Given a filter array interpolate to the specified wavelength.

double procedure eq_filter (w, band)

double	w		# Wavelength desired
pointer	band		# Band pointer

int	i, n
double	x1, x2
pointer	x, y

begin
	n = BAND_FN(band)
	if (n == 0)
	    return (1.)

	x = BAND_FW(band)
	y = BAND_FR(band)
	x1 = Memd[x]
	x2 = Memd[x+n-1]

	if (w <= x1)
	    return (Memd[y])
	else if (w >= x2)
	    return (Memd[y+n-1])
	
	if ((w - x1) < (x2 - w))
	    for (i = 1; w > Memd[x+i]; i=i+1)
		;
	else
	    for (i = n - 1; w < Memd[x+i-1]; i=i-1)
		;
		
	x1 = Memd[x+i-1]
	x2 = Memd[x+i]
	return ((w - x1) / (x2 - x1) * (Memd[y+i] - Memd[y+i-1]) + Memd[y+i-1])
end
# Sep 30 2002	Add redshift to noao.onedspec.sbands
# Oct  1 2002	Print redshift
# Oct  8 2002	Rewrite spectrum input to match other RVSAO packages

# Sep 28 2004	Add options to normalize by exposure (byexp) and pixels (bypix)

# Jan 14 2005	Add debugging
# Jan 19 2005	Add /eqw/ common; fix width bug
# Jan 19 2005	Print output even if fluxes are zero
# Jan 28 2005	Print flux - continuum
# Feb  9 2005	Add parameters and optional continuum fitting
# Feb 23 2005	Deal with binned data and fix bug reading number of rows
# Mar  1 2005	Clean up code; read EXPOSUIRE if EXPTIME not in header
# Mar  3 2005	Add string length to calls to strcpy() and strcat()
# Jul 27 2005	Add option to log to multiple files
# Jul 28 2005	If no continuum region found, print "none" and zeroes
# Sep 21 2005	Fix flux computation when fit continuum subtracted
# Sep 22 2005	Clean up image name for logging

# May 28 2008	Add output option 2, including index and error
# May 28 2008	Add skyband and read sky spectrum if not zero
# May 29 2008	Allow either limiting wavelengths or centers and widths
# Jul  2 2008	Add bandcont and bandfilt so filter and cont names are optional
# Jul 11 2008	If limits, cont1-center-cont2; if centers, center-cont1-cont2
# Jul 11 2008	Print one line from EQ_FLUX
# Jul 11 2008	Improve error handling and print 1+z heading only once
# Aug 13 2008	Drop unused variable i from eq_proc
# Sep  5 2008	If all of spectral region is off spectrum set flux to zero
# Dec  9 2008	Compute Eqw index error using John Huchra's method

# Jan 26 2009	Divide line flux by n-1 in noise computation
# Jan 29 2009	Set velocity to zero if not in header
# Jan 30 2009	Normalize *all* fluxes by exposure, if requested
# Apr 15 2009	Rename variables for clarity in index and eqw computation
# Apr 15 2009	Use bandpass sigmas for index and eqw error computation
# May 14 2009	Use wavelength limits instead of centers and widths
# May 21 2009	Check shifted instead of rest wavelengths agains spectrum limits
# Jun  4 2009	Drop unused variable twidth from eq_bands()
# Jul 16 2009	Update sigma computation after John derives it again
# Jul 17 2009	Check sigmas for >0 before dividing
# Jul 20 2009	Fix log of negative number when fluxratio near 1.0
# Jul 20 2009	Fix equivalent width error computation again
# Jul 20 2009	Make index error always positive
