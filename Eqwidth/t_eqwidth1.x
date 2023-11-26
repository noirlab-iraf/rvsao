# File rvsao/Eqwidth/t_eqwidth.x
# May 29, 2008
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
pointer	sp
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
include "emv.com"
include "eqw.com"

define  endspec_ 90

begin
	c0 = 299792.5d0
	comma = ','
	dindef = INDEFD
	specsh = NULL
	work = NULL
	cont = NULL
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	wenable = FALSE
	conproc = CONTFIT
	rmode = 1

	call smark (sp)

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
	i = 0
	call strcpy ("rvsao.eqwidth",taskname,SZ_LINE)
	while (clgfil (logfiles, logfile, SZ_PATHNAME) != EOF) {
	    fd = open (logfile, APPEND, TEXT_FILE)
	    if (fd == ERR) next
	    i = i + 1
	    logfd[i] = fd
	    }
	nlogfd = i
	if (nlogfd < 1) {
	    nlogfd = 1
	    logfd[1] = STDOUT
	    }
	call clpcls (logfiles)

	call clgstr ("bands", fbands, SZ_FNAME)
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
	call eq_bands (fbands, bands, nbands, nsubbands)

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
			call pargi (nspec)
		    }

	    # Load sky spectrum for error computation
		if (sband > 0) {
		    call getspec (specpath,mspec,sband,sspec,skyim,wenable)
		    skysh = specsh
		    }

	    # Load spectrum
		spvel = 0.d0
		spechcv = 0.d0
		call getspec (specpath,mspec,mband,spectrum,specim,wenable)
#		if (specim == ERR)
#		    go to endspec_

	    # Compute the wavelength shift observed/rest, removing helio corr.
		if (torest) {
		    z1 = (1.d0 + (spvel / c0)) / (1.d0 + (spechcv / c0))
		    if (debug) {
			call printf ("EQWIDTH: Redshift z = %9.6f\n")
			    call pargd (z1)
			call flush (STDOUT)
			}
		    }
		else
		    z1 = 1.d0

	    # Set up wavelength vector for plotting
		call wcs_set (specsh)
		if (wlspec == NULL)
	 	    call salloc (wlspec, specpix, TY_REAL)
		do j = 1, specpix {
		    wl = wcs_p2w (double (j))
		    Memr[wlspec+j-1] = real (wl / z1)
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
		    call imgspar (specim, "FINDOBJ", rows, SZ_LINE)
		    if (strlen (rows) == 0)
			nrows = 1
		    else if (decode_ranges (rows,row_range,MAX_RANGES,nrows) == ERR)
			nrows = 1
		    if (debug) {
			call printf ("EQWIDTH: rows %s = %d rows\n")
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
		if (nimages == 1) {
		    if (verbose)
			call eq_header (rmode,specid,fbands,bands,nbands,nsubbands)
		    }

		if (work == NULL)
		    call salloc (work, specpix, TY_REAL)
		if (cont == NULL)
		    call salloc (cont, specpix, TY_REAL)

	    # Set up smoothed spectrum for plotting
		if (smspec == NULL)
	 	    call salloc (smspec, specpix, TY_REAL)
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
		    }

	    # Set up smoothed continuum-subtracted spectrum for plotting
		if (smcspec == NULL)
	 	    call salloc (smcspec, specpix, TY_REAL)
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
#		    call eqplot (specsh, specid, specim, specpix,
#				 Memr[spectrum], Memr[smspec], Memr[work],
#				 Memr[smcspec], Memr[wlspec], bands, nbands)
		    }

	    # Compute equivalent widths for this spectrum
		call eq_proc (specsh, Memr[spectrum], Memr[sspec],
			      Memr[wlspec], specpix, Memr[cont], bands, nbands)

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
	call imtclose (speclist)
	call sfree (sp)
end


# EQ_BANDS - Read bands from the band file and put them into an array
# of band pointers.

procedure eq_bands (fbands, bands, nbands, nsubbands)

char	fbands[ARB]		#I File containing band information
pointer	bands			#O Bandpass table descriptor
int	nbands			#O Number of bandpasses
int	nsubbands		#O Number of individual bands

bool	bandok
int	ip
int	fd			#I Bandpass file descriptor
double	center, width, twidth
pointer	sp, line, id, filter

int	getline(), ctowrd(), ctod()
int	open()

begin
	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)
	call salloc (id, SZ_FNAME, TY_CHAR)
	call salloc (filter, SZ_FNAME, TY_CHAR)

	fd = open (fbands, READ_ONLY, TEXT_FILE)

	# Read the bands.  If the first band is not seen
	# skip the line.  Check for 1, 2, or 3 bandpasses.
	# Can't use fscan() because fscan() will be called later to
	# read any filter file.

	bands = NULL
	nbands = 0
	nsubbands = 0
	while (getline (fd, Memc[line]) != EOF) {

# Read line center
	    ip = 1
	    bandok = (ctowrd (Memc[line], ip, Memc[id], SZ_FNAME) > 0)
	    bandok = (bandok && ctod (Memc[line], ip, center) > 0)
	    bandok = (bandok && ctod (Memc[line], ip, width) > 0)
	    bandok = (bandok && ctowrd (Memc[line],ip,Memc[filter],SZ_FNAME)>0)
	    if (!bandok || Memc[id] == '#')
		next

# Convert from limits to center and width
	    if (width > center) {
		twidth = width - center
		center = (center + width) / 2.d0
		width = twidth
		}

	    # Allocate and reallocate the array of band pointers.
	    if (nbands == 0)
		call malloc (bands, 10 * NBANDS, TY_POINTER)
	    else if (mod (nbands, 10) == 0)
		call realloc (bands, (nbands + 10) * NBANDS, TY_POINTER)
	    nbands = nbands + 1

	    call eq_alloc (BAND(bands,nbands,BAND1),
		Memc[id], Memc[filter], center, width)
	    nsubbands = nsubbands + 1

	    bandok = (ctowrd (Memc[line], ip, Memc[id], SZ_FNAME) > 0)
	    bandok = (bandok && ctod (Memc[line], ip, center) > 0)
	    bandok = (bandok && ctod (Memc[line], ip, width) > 0)
	    bandok = (bandok && ctowrd (Memc[line],ip,Memc[filter],SZ_FNAME)>0)
	    if (bandok) {
		if (width > center) {
		    twidth = width - center
		    center = (center + width) / 2.0
		    width = twidth
		    }
		call eq_alloc (BAND(bands,nbands,BAND2),
		    Memc[id], Memc[filter], center, width)
		nsubbands = nsubbands + 1
		}
	    else
		BAND(bands,nbands,BAND2) = NULL

	    bandok = (ctowrd (Memc[line], ip, Memc[id], SZ_FNAME) > 0)
	    bandok = (bandok && ctod (Memc[line], ip, center) > 0)
	    bandok = (bandok && ctod (Memc[line], ip, width) > 0)
	    bandok = (bandok && ctowrd (Memc[line],ip,Memc[filter],SZ_FNAME)>0)
	    if (bandok) {
		if (width > center) {
		    twidth = width - center
		    center = (center + width) / 2.0
		    width = twidth
		    }
		call eq_alloc (BAND(bands,nbands,BAND3),
		    Memc[id], Memc[filter], center, width)
		nsubbands = nsubbands + 1
		}
	    else
		BAND(bands,nbands,BAND3) = NULL
	    }

	call close (fd)
	call sfree (sp)
end


# EQ_ALLOC -- Allocate a band structure.

procedure eq_alloc (band, id, filter, center, width)

pointer	band			#O Band pointer
char	id[ARB]			#I Band id
char	filter[ARB]		#I Band filter
double	center			#I Band wavelength
double	width			#I Band width

int	fn, fd, strlen(), open(), fscan(), nscan()
double	w, r
pointer	fw, fr
bool	streq()
errchk	open()

begin
	call calloc (band, LEN_BAND, TY_STRUCT)
	call malloc (BAND_ID(band), strlen(id), TY_CHAR)
	call malloc (BAND_FILTER(band), strlen(filter), TY_CHAR)
	call strcpy (id, Memc[BAND_ID(band)], ARB)
	call strcpy (filter, Memc[BAND_FILTER(band)], ARB)
	BAND_WC(band) = center
	BAND_DW(band) = width
	BAND_FN(band) = 0
	BAND_FW(band) = NULL
	BAND_FR(band) = NULL

	if (streq (filter, "none"))
	    return

	# Read the filter file.
	fd = open (filter, READ_ONLY, TEXT_FILE)
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
end


# EQ_HEADER -- Print output header.

procedure eq_header (spec, fbands, bands, nbands, nsubbands)

char	spec[ARB]		#I Typical spectrum name
char	fbands[ARB]		#I Band file
pointer	bands			#I Pointer to array of bands
int	nbands			#I Number of bands
int	nsubbands		#I Number of subbands

int	i, j, k, lspec
pointer	sp, str, band
char	fspec[16]

int	strlen()

include "eqw.com"
include "rvsao.com"

begin
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)

	# Output a banner and task parameters.
	call sysid (Memc[str], SZ_LINE)

	do k = 1, nlogfd {
	    call fprintf (logfd[k], "\n# EQWIDTH: %s\n#  ")
		call pargstr (Memc[str])
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
		call strcpy ("mag", Memc[str], SZ_LINE)
		}
	    else
		call strcpy ("flux", Memc[str], SZ_LINE)
	    if (torest)
		call fprintf (logfd[k], " torest = yes")
	    if (byexp)
		call fprintf (logfd[k], " byexp = yes")
	    if (bypix)
		call fprintf (logfd[k], " bypix = yes")

	    # Output the bands.
	    call fprintf (logfd[k], "\n# %14s %10s %10s %10s\n")
		call pargstr ("band")
		call pargstr ("filter")
		call pargstr ("wavelength")
		call pargstr ("width")
	    do i = 1, nbands {
		do j = 1, NBANDS {
		    band = BAND(bands,i,j)
		    if (band == NULL)
			next
		    call fprintf (logfd[k], "# %14s %10s %10g %10g\n")
			call pargstr (Memc[BAND_ID(band)])
			call pargstr (Memc[BAND_FILTER(band)])
			call pargd (BAND_WC(band))
			call pargd (BAND_DW(band))
		    }
		}
	    call fprintf (logfd[k], "#\n")

	    # Output column headings.

	    # Print indices and equivalent widths
	    lspec = strlen (spec)
	    if (rmode == 1) {
		call sprintf (fspec, 16, "%ds	")
		    call pargi (lspec)
		call fprintf (logfd[k], fspec)
		    call pargstr ("spectrum")
		do i = 1, nbands {
		    call fprintf (logfd[k],"	band%d	flux%d	cont%d	cflux%d	index%d	eqwidth%d")
			call pargi (i)
			call pargi (i)
			call pargi (i)
			call pargi (i)
			call pargi (i)
			call pargi (i)
		    }
		if (torest) {
		    call fprintf (logfd[k], " %7s")
			call pargstr ("  1+z  ")
		    }
		call fprintf (logfd[k], "\n")

		# Print underlines for Starbase table
		call fprintf (logfd[k], "%24s	")
		call pargstr ("------------------------")
		do i = 1, nbands {
		    call fprintf (logfd[k],"	----	-------	----	--------	-------	-------")
		    }
		}

	    # Print indices and equivalent widths
	    else if (rmode == 2) {
		call fprintf (logfd[k], "%24s	")
		call pargstr ("spectrum")
		band = BAND(bands,i,1)
		do i = 1, nbands {
		    call fprintf (logfd[k],"	flux%s	cflux%s		index%s	error%s")
			call pargstr (Memc[BAND_ID(band)])
			call pargstr (Memc[BAND_ID(band)])
			call pargstr (Memc[BAND_ID(band)])
			call pargstr (Memc[BAND_ID(band)])
		    }
		}
	    }

	call sfree (sp)
end


# EQ_PROC -- Measure the band fluxes and possibly a band index and eq. width.

procedure eq_proc (sh, spec, skysh, sspec, wl, npix, cont, bands, nbands)

pointer	sh			#I Object spectrum descriptor
real	spec[ARB]		#I Object spectrum
pointer	skysh			#I Sky spectrum descriptor
real	sspec[ARB]		#I Sky spectrum
real	wl[ARB]			#I Wavelengths for spectrum
int	npix			#I Number of pixels in spectrum
real	cont[ARB]		#I Object continuum vector, if one has been fit
pointer	bands			#I Bandpass table pointer
int	nbands			#I Number of bandpasses

int	i, j
double	flux, contval, index, eqwidth, wc1, wc2, wc3, dw1, dw2, dw3
double	flux1, norm1, flux2, norm2, flux3, norm3, a, b, fluxnet, zero
double	sflux1, snorm1, sflux2, snorm2, sflux3, snorm3
double	sigma1, sigma2, sigma3, sigma
pointer	sp, imname, band1, band2, band3
bool	skystat

int	strcmp()

include "eqw.com"
include "rvsao.com"

begin
	call smark (sp)
	call salloc (imname, SZ_FNAME, TY_CHAR)
	zero = 0.0d0
	if (skysh != NULL)
	    skystat = TRUE
	else
	    skystat = FALSE

	if (APHIGH(sh) > 1) {
	    if (strcmp (IMSEC(sh), "[\*,1,1]") != 0) {
		call sprintf (Memc[imname], SZ_FNAME, "%s%s(%d)")
		    call pargstr (IMNAME(sh))
		    call pargstr (IMSEC(sh))
		    call pargi (AP(sh))
		}
	    else {
		call sprintf (Memc[imname], SZ_FNAME, "%s(%d)")
		    call pargstr (IMNAME(sh))
		    call pargi (AP(sh))
		}
	    }
	else {
	    if (strcmp (IMSEC(sh), "[\*,1,1]") != 0) {
		call sprintf (Memc[imname], SZ_FNAME, "%s%s")
		    call pargstr (IMNAME(sh))
		    call pargstr (IMSEC(sh))
		}
	    else {
		call sprintf (Memc[imname], SZ_FNAME, "%s")
		    call pargstr (IMNAME(sh))
		}
	    }

	do j = 1, nlogfd {
	    call fprintf (logfd[j], "%26s")
		call pargstr (Memc[imname])
	    }

	# Loop over all bandpasses
	do i = 1, nbands {

	    # Measure primary band flux, normalize, and print result.
	    band1 = BAND(bands,i,BAND1)
	    wc1 = BAND_WC(band1)
	    dw1 = BAND_DW(band1)
	    if (debug) {
		call printf ("\nBand %d %7.7s %11.6g %11.6g\n")
		    call pargi (i)
		    call pargstr (Memc[BAND_ID(band1)])
		    call pargd (wc1)
		    call pargd (dw1)
		call flush (STDOUT)
		}

	    # Measure primary band sky flux, normalize, and print result.
	    call eq_flux (sh,spec,wl,npix,wc1,dw1,band1,flux1,norm1,debug)
	    if (IS_INDEFD(flux1))
		flux1 = 0.d0
	    if (skysh != NULL) {
		call eq_flux (skysh,sspec,wl,npix,wc1,dw1,band1,sflux1,snorm1,debug)
		if (IS_INDEFD(sflux1))
		    sflux1 = 0.d0
		sigma1 = dsqrt (flux1 + (2.d0 * sflux1)) / flux1
		}
	    else {
		sflux1 = 0.d0
		sigma1 = dsqrt (flux1) / flux1
		}

	    # Divide out exposure time
	    if (byexp)
		flux1 = flux1 / exptime

	    # Divide out number of rows for skies
	    if (bypix)
		flux1 = flux1 / dnrows

	    if (norm) {
		flux1 = flux1 / norm1
		norm1 = 1.d0
		}
	    if (mag && flux1 > 0.)
		flux = magzero - 2.5 * log10 (flux1)
	    else
		    flux = flux1

	    do j = 1, nlogfd {
		if (rmode == 1) {
		    call fprintf (logfd[j], "	%7.7s	%11.6g")
			call pargstr (Memc[BAND_ID(band1)])
		    }
		    call pargd (flux)
		}
	    if (debug)
		call printf ("\n")

	    # Measure continuum in same place if it has been fit
	    if (fitcont) {
		call eq_flux (sh,cont,wl,npix,wc1,dw1,band1,flux2,norm2,debug)

		# Print zeroes and go to next band if no continuum is specified
		if (IS_INDEFD(flux2)) {
		    do j = 1, nlogfd {
			if (rmode ==1 ) {
			    call fprintf (logfd[j], "	%7.7s")
				call pargstr ("cont")
			    }
			call fprintf (logfd[j], " %11.6g")
			    call pargd (zero)
			call fprintf (logfd[j], " %9.6g %9.6g")
			    call pargd (zero)
			    call pargd (zero)
			}
		    next
		    }
		else {
		    sigma2 = dsqrt (flux2 + (2.d0 * sflux1)) / flux2
		    sigma3 = 0.d0
		    }
		}

	    # Otherwise measure it in the indicated places
	    else {

		# Measure the first continuum band object flux
		band2 = BAND(bands,i,BAND2)
		wc2 = BAND_WC(band2)
		dw2 = BAND_DW(band2)
		call eq_flux (sh,spec,wl,npix,wc2,dw2,band2,flux2,norm2,debug)

		# Measure the first continuum band sky flux
		if (skystat) {
		    call eq_flux (skysh,sspec,wl,npix,wc2,dw2,band2,sflux2,snorm2,debug)
		    if (IS_INDEFD(sflux2))
			sflux2 = 0.d0
		    if (IS_INDEFD(sflux2))
			sflux2 = 0.d0
		    sigma2 = dsqrt (flux2 + (2.d0 * sflux2)) / flux2
		    }
		else {
		    sflux2 = 0.d0
		    sigma2 = dsqrt (flux2) / flux2
		    }

		# Measure the second continuum band flux
		band3 = BAND(bands,i,BAND3)
		wc3 = BAND_WC(band3)
		dw3 = BAND_DW(band3)
		call eq_flux (sh,spec,wl,npix,wc3,dw3,band3,flux3,norm3,debug)

		# Measure the second continuum band sky flux
		if (skystat) {
		    call eq_flux (skysh,sspec,wl,npix,wc3,dw3,band3,sflux3,snorm3,debug)
		    if (IS_INDEFD(sflux2))
			sflux3 = 0.d0
		    if (IS_INDEFD(sflux2))
			sflux3 = 0.d0
		    sigma3 = dsqrt (flux3 + (2.d0 * sflux3)) / flux3
		    }
		else {
		    sflux3 = 0.d0
		    sigma3 = dsqrt (flux3) / flux3
		    }

		# Print zeroes and go to next band if no continuum is specified
		if ((IS_INDEFD(flux2)) && (IS_INDEFD(flux3))) {
		    do j = 1, nlogfd {
			if (rmode ==1 ) {
			    call fprintf (logfd[j], "	%7.7s")
				call pargstr ("cont")
			    }
			call fprintf (logfd[j], " %11.6g")
			    call pargd (zero)
			call fprintf (logfd[j], " %9.6g %9.6g")
			    call pargd (zero)
			    call pargd (zero)
			}
		    next
		    }

		# Compute and output the band index and equivalent width.
		if (net || norm) {
		    if (!IS_INDEFD(flux2)) {
			flux2 = flux2 / norm2
			norm2 = 1.d0
			}
		    if (!IS_INDEFD(flux3)) {
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

	    # Use underlying fit continuum
	    if (fitcont) {
		contval = flux2
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr ("cont")
			}
		    }
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
		sigma = dsqrt (sigma1*sigma1 + sigma2*sigma2 + sigma3*sigma3)
		}

	    # Use first continuum if second is not computable
	    else if (!IS_INDEFD(flux2)) {
		contval = flux2
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr (Memc[BAND_ID(band2)])
			}
		    }
		sigma = dsqrt (sigma1*sigma1 + sigma2*sigma2)
		}

	    # Use second continuum if first is not computable
	    else if (!IS_INDEFD(flux3)) {
		contval = flux3
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr (Memc[BAND_ID(band3)])
			}
		    }
		sigma = dsqrt (sigma1*sigma1 + sigma3*sigma3)
		}

	    else {
		contval = 0.0
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%7.7s")
			    call pargstr ("cont")
			}
		    }
		sigma = sigma1
		}

	    if (mag && contval > 0.)
		flux = magzero - 2.5 * log10 (contval)
	    else
		flux = contval
	    do j = 1, nlogfd {
		call fprintf (logfd[j], "	%11.6g")
		    call pargd (flux)
		}

	    if (debug) {
		call printf ("\n")
		call flush (STDOUT)
		}

	    if (flux1 > 0. && contval > 0.) {
		index = flux1 / contval
		eqwidth = (1.d0 - index) * dw1
		fluxnet = flux1 - contval
	    }
	    else if (flux1 > 0.)
		fluxnet = flux1
	    else {
		index = 0.0
		eqwidth = 0.0
		fluxnet = 0.0
		}
	    if (mag) {
		if (!IS_INDEFD(contval) && contval > 0.)
		    contval = magzero - 2.5 * log10 (contval)
		if (!IS_INDEFD(index))
		    index = -2.5 * log10 (index)
		}

	    if (net) {
		do j = 1, nlogfd {
		    call fprintf (logfd[j], "	%9.6g")
			call pargd (fluxnet)
		    }
		}
	    else {
		do j = 1, nlogfd {
		    if (rmode == 1) {
			call fprintf (logfd[j], "	%9.6g	%9.6g")
			    call pargd (index)
			    call pargd (eqwidth)
			}
		    else if (rmode == 2) {
			call fprintf (logfd[j], "	%9.6g	%9.6g")
			    call pargd (index)
			    call pargd (sigma)
			}
		    }
		}
	    }

# Print redshift as 1+z if not zero
	do j = 1, nlogfd {
	    if (z1 != 1.0) {
		call fprintf (logfd[j], " %7.5f")
		    call pargd (z1)
		}
	    call fprintf (logfd[j], "\n")
	    }

	# Flush output and finish up.
	do j = 1, nlogfd {
	    call flush (logfd[j])
	    }
	call sfree (sp)
end


# EQ_FLUX - Compute the flux and total response in a given band.
# Return INDEF if the band is outside of the spectrum.

procedure eq_flux (sh, spec, wl, npix, cw, dw, band, flux, normb, debug)

pointer	sh			#I Spectrum descriptor
real	spec[ARB]		#I spectrum
real	wl[ARB]			#I wavelengths for spectrum
int	npix			#I Number of pixels in spectrum
double	cw			#I center wavelength for band
double	dw			#I wavelength width for band
pointer	band			#I band descriptor
double	flux			#O flux
double	normb			#O normalization
bool	debug

int	i, i1, i2, fluxi, flux2
double	a, b, r1, r2, w1, w2, di1, di2, wt, sigma
double	eq_filter(), shdr_wl()

include "eqw.com"

begin
	flux = INDEFD
	normb = 1.0
	if (band == NULL)
	    return

	r1 = cw - (dw / 2.)
	r2 = cw + (dw / 2.)

	# Redshift wavelength limits to get correct pixel limits
	w1 = r1 * z1
	w2 = r2 * z1
	a = shdr_wl (sh, w1)
	b = shdr_wl (sh, w2)
	di1 = min (a, b)
	di2 = max (a, b)
	i1 = nint (di1)
	i2 = nint (di2)
	if (debug) {
	    call printf ("EQFLUX: %.3f-%.3f = pixels %.3f-%.3f")
		call pargd (w1)
		call pargd (w2)
		call pargd (di1)
		call pargd (di2)
	    if (i1 > npix)
		call printf (" off spectrum\n")
	    else if (i2 < 1)
		call printf (" off spectrum\n")
	    else if (i2 > npix) {
		call printf (" > %d\n")
		    call pargi (npix)
		}
	    else if (i1 < 1)
		call printf (" < 1\n")
	    else
		call printf ("\n")
	    }
	if (i2 > npix) {
	    i2 = npix
	    di2 = double (npix)
	    }
	if (i1 < 1) {
	    i1 = 1
	    di2 = 1.d0
	    }
	if (di1 == di2 || i2 < 1 || i1 > npix)
	    return

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
	    fluxi = wt * spec[i1]
	    flux = fluxi
	    flux2 = fluxi * fluxi
	    normb = wt

# Middle pixels
	    do i = i1+1, i2-1, 1 {
		wt = eq_filter (wl[i], band)
		fluxi = wt * spec[i]
		flux = flux + fluxi
		flux2 = flux2 + (fluxi * fluxi)
		normb = normb + wt
		}

# Last pixel
	    wt = eq_filter (double (wl[i2]), band) * (di2 - double (i2) + 0.5d0)
	    fluxi = wt * spec[i2]
	    flux = flux + fluxi
	    flux2 = flux2 + (fluxi * fluxi)
	    normb = normb + wt

	    mean = flux / normb
	    sigma = (flux2 - (flux * flux)) / normb
	    }

	if (debug) {
	    call printf ("\nEQFLUX: Band %7.2f %7.2f-%7.2f (%7.2f-%7.2f)(%d-%d): %9.6g (%9.6g)\n")
		call pargd (cw)
		call pargd (r1)
		call pargd (r2)
		call pargd (w1)
		call pargd (w2)
		call pargi (i1)
		call pargi (i2)
		call pargd (flux)
		call pargd (normb)
	    }
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
