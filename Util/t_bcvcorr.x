# File rvsao/Util/t_bcvcorr.x
# April 21, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1995-2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# BCVCORR is an IRAF task which computes heliocentric and solar system
# barycentric corrections for spectra


include	<imhdr.h>
include	<fset.h>
include <imhdr.h>
include <imio.h>

define	LEN_USER_AREA	100000
define	MAX_RANGES	50
define	SZ_HKWORD	9
define	SZ_MONTH	3
define	UT_WHEN		"|start|mid|end|"
define	UT_START	1
define	UT_MID		2
define	UT_END		3


procedure t_bcvcorr ()

pointer	specim		# image structure for spectrum
char	spectra[SZ_LINE]	# Object spectra list
char	specfile[SZ_PATHNAME]	# Object spectrum file name
char	specdir[SZ_PATHNAME]	# Directory for object spectra
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	obsname[SZ_LINE]	# Observatory name
char	keyword[SZ_LINE]	# Various keywords
char	keyobs[SZ_LINE]		# Observatory name keyword
char	keyhjd[SZ_LINE]		# Header keyword for Heliocentric Julian Day
pointer	speclist	# List of spectrum files
bool	savebcv,savebcv0	# Save velocity correction in data file header
bool	savejd,savejd0	# Save Julian Date and midtime in data file header
int	mode		# Mode (READ_ONLY, READ_WRITE)
int	printmode	# Output mode (1=all, 2=bcv only...)
pointer	obs		# pointer to observatory structure
bool	gottime		# True if observation time obtained
bool	gotdate		# True if observation date obtained
bool	gotjd		# True if observation JD is from parameter list
bool	keydate		# True if observation date from keywords
bool	gotpos		# True if position of observed object obtained
bool	keypos		# True if position of observed object from keywords
bool	specsky		# True if spectrum is of twilight sky
bool	subgrav		# True if solar gravitaional redshift is removed
#bool	debug
bool	verbose, debug, clgetb()
int	mm, dd, yyyy, ldir, ndim, clgeti()
double	ra, dec, eq, exp, dindef
double	ut1, ut2, utm, ut, dlong, dlat, dalt
double	hjd		# Heliocentric Julian Date
double	gjd		# Geocentric Julian Date
double	dpi,dtr, lt, dt
double	dgcv, dbcv, dhcv, fbcv, fhcv, tpi, sgrv
real	rbcv, rhcv
int	lspec
int	reltime
char	colon
char	specnums[SZ_LINE]
int	mspec_range[3,MAX_RANGES]
int	nmspec0,mspec
char	str[SZ_LINE]
char	uts[32]
char	when[SZ_LINE]
char	rstr[32], dstr[32]
string  month   "JanFebMarAprMayJunJulAugSepOctNovDec"

int	strlen(), imtgetim(), imaccess(), strncmp()
double	juldat(), clgetd()
int	ctod(), ip
int	strdic()
pointer	imtopen()
pointer	immap()
errchk	immap()
double	obsgetd()
pointer	obsvopen()
errchk	obsvopen
int	decode_ranges(),get_next_number()
int	imaccf()
int	stridx(), strcmp(), strmatch()

define	newspec_	10
define	bcvcomp_	50
define	endvc_	90

begin
	colon = ':'
	dpi = 3.1415926535897932d0
	tpi = 2.d0 * dpi
	dtr = dpi / 180.0d0
	dindef = INDEFD
	gjd = dindef
	hjd = dindef
	obs = NULL
	gotdate = FALSE
	gottime = FALSE

	verbose = clgetb ("verbose")
	debug = clgetb ("debug")
	printmode = clgeti ("printmode")

# Figure out which values come from paramter list and which from header

# Position (convert RA to degrees if sexigesimal)
	ra = dindef
	call strcpy ("INDEF", str, SZ_LINE)
	call clgstr ("ra", str, SZ_LINE)
	ip = 1
	if (strcmp (str, "INDEF") != 0) {
	    if (stridx (colon, str) > 0) {
		if (ctod (str, ip, ra) <= 0)
		    ra = dindef
		if (debug) {
		    call printf ("BCVCORR: RA %s -> %0h\n")
			call pargstr (str)
			call pargd (ra)
		    call flush (STDOUT)
		    }
		}
	    else {
		if (ctod (str, ip, ra) <= 0)
		    ra = dindef
		else
		    ra = ra / 15.d0
		if (debug) {
		    call printf ("BCVCORR: RA %s -> %0h\n")
			call pargstr (str)
			call pargd (ra)
		    call flush (STDOUT)
		    }
		}
	    }
	dec = clgetd ("dec")
	eq = clgetd ("equinox")
	if (ra != dindef && dec != dindef)
	    gotpos = TRUE
	else
	    gotpos = FALSE
	if (gotpos)
	    keypos = FALSE
	else
	    keypos = TRUE

# Date
	hjd = clgetd ("hjd")
	gjd = clgetd ("gjd")
	yyyy = clgeti ("year")
	mm = clgeti ("month")
	dd = clgeti ("day")
	if (yyyy != INDEFI && mm != INDEFI && dd != INDEFI)
	    gotdate = TRUE
	else if (gjd != 0.d0 && gjd != dindef)
	    gotdate = TRUE
	else
	    gottime = FALSE
	ut = clgetd ("ut")
	gotjd = FALSE
	if (ut != dindef)
	    gottime = TRUE
	else if (gjd != 0.d0 && gjd != dindef) {
	    gotjd = TRUE
	    gotdate = TRUE
	    }
	else
	    gottime = FALSE
	if (gottime)
	    keydate = FALSE
	else
	    keydate = TRUE

# Spectra for which to compute solar system center velocity corrections
	call clgstr ("spectra",spectra,SZ_PATHNAME)
	lspec = strlen (spectra)
	if (lspec == 0) {
	    specim = NULL
	    savebcv = FALSE
	    savejd = FALSE
	    go to bcvcomp_
	    }
	speclist = imtopen (spectra)
	call clgstr ("specdir",specpath,SZ_PATHNAME)

# Multispec spectrum numbers (use only first if multiple files)
	call clgstr ("specnum",specnums,SZ_LINE)
	if (specnums[1] == EOS)
	    call strcpy ("0",specnums,SZ_LINE)
	call flush (STDOUT)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR){
	    call sprintf (str, SZ_LINE, "T_XCSAO: Illegal multispec list <%s>")
	    call pargstr (specnums)
	    call error (1, str)
	    }

	call clgstr ("specdir",specdir,SZ_PATHNAME)
	savebcv0 = clgetb ("savebcv")
	savejd0 = clgetb ("savejd")

	specsky = FALSE
	specsky = clgetb ("specsky")

# Get next object spectrum file name from the list
newspec_
	savebcv = savebcv0
	savejd = savejd0
	dindef = INDEFD
	if (!keydate) {
	    gjd = dindef
	    hjd = dindef
	    }
	if (imtgetim (speclist, specfile, SZ_PATHNAME) == EOF)
	   go to endvc_

	if (get_next_number (mspec_range, mspec) == EOF)
	    mspec = 0

# Check for readability of object spectrum
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

	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("BCVCORR: cannot read spectrum file %s \n")
		call pargstr (specpath)
	    go to newspec_
	    }
	else if (debug) {
	    call printf ("BCVCORR: Reading spectrum file %s \n")
		call pargstr (specpath)
	    call flush (STDOUT)
	    }

# Load spectrum header
	if (savebcv0 || savejd0)
	    mode = READ_WRITE
	else
	    mode = READ_ONLY
        iferr (specim = immap (specpath, mode, LEN_USER_AREA)) {
            mode = READ_ONLY
	    if (debug) {
		call printf ("Opening %s read-only\n")
		    call pargstr (specpath)
		call flush (STDOUT)
		}
	    iferr (specim = immap (specpath, mode, LEN_USER_AREA)) {
		if (specim != NULL) call imunmap (specim)
		call eprintf ("BCVCORR: Cannot read image %s\n")
		call pargstr (specpath)
		specim = ERR
		go to endvc_
		}
	    call eprintf ("BCVCORR: Cannot write to image %s\n")
	    call pargstr (specpath)
	    savebcv = FALSE
	    savejd = FALSE
	    }
        ndim = IM_NDIM(specim)
	if (debug || (verbose && printmode == 1)) {
	    call printf ("%s: %d x %d x %d %d-D image\n")
		call pargstr (specpath)
		call pargi (IM_LEN(specim,1))
		call pargi (IM_LEN(specim,2))
		call pargi (IM_LEN(specim,3))
		call pargi (ndim)
	    call flush (STDOUT)
	    }

# Direction for velocity correction (convert RA to degrees if sexigesimal)
	if (keypos) {
	    call clgstr ("keyra",keyword,SZ_LINE)
	    ra = dindef
 	    call imgdpar (specim, keyword, ra)
	    if (ra != dindef) {
		call imgspar (specim, "ra", str, SZ_LINE)
		if (stridx (colon, str) == 0)
		    ra = ra / 15.d0
		}
	    call clgstr ("keydec",keyword,SZ_LINE)
	    dec = dindef
	    call imgdpar (specim,keyword, dec)
	    call clgstr ("keyeqnx",keyword,SZ_LINE)
	    eq = 1950.0d0
	    call imgdpar (specim,keyword,eq)
	    if (ra != dindef && dec != dindef)
		gotpos = TRUE
	    else
		gotpos = FALSE
	    }

# Time of observation from image header
	if (keydate) {
	    call clgstr ("keydate",keyword,SZ_LINE)

# Julian date (JD or jd)
	    if (strmatch (keyword, "JD") > 0 ||
		strmatch (keyword, "jd") > 0 ) {
		call flush (STDOUT)
		call imgdpar (specim, keyword, gjd)

# Modified Julian Date (MJD or mjd)
		if (strmatch (keyword, "MJD") > 0 ||
		    strmatch (keyword, "mjd") > 0 ) {
		    if (debug) {
			call printf ("MJD %.5f -> ")
			    call pargd (gjd)
			}
		    gjd = gjd + 2400000.5d0
		    }
		call jd2ut (gjd, yyyy, mm, dd, ut)
		if (debug) {
		    call printf ("JD %.5f -> %4d-%02d-%02dT%13.4h\n")
			call pargd (gjd)
			call pargi (yyyy)
			call pargi (mm)
			call pargi (dd)
			call pargd (ut)
		    }
		call flush (STDOUT)
		gottime = TRUE
		gotdate = TRUE
		gotjd = TRUE
		}

# FITS date (DATE-OBS or date-obs)
	    else if (strmatch (keyword, "DATE") > 0 ||
		     strmatch (keyword, "date") > 0 ) {
		call imgdate (specim, keyword, mm, dd, yyyy)
		if (debug) {
		    call printf ("%s -> %4d-%02d-%02d\n")
			call pargstr (keyword)
			call pargi (yyyy)
			call pargi (mm)
			call pargi (dd)
		    call flush (STDOUT)
		    }
		gotdate = TRUE
		ut = dindef
		call imgdtim (specim, keyword, ut)
		if (debug) {
		    call printf ("%s -> %h UT\n")
			call pargstr (keyword)
			call pargd (ut)
		    call flush (STDOUT)
		    }
		if (ut != dindef) {
		    gottime = TRUE
		    }
		else {
		    call clgstr ("keytime",keyword,SZ_LINE)
		    call imgdpar (specim, keyword, ut)
		    if (ut != dindef) {
			gottime = TRUE
			}
		    }
		if (gotdate && gottime) {
		    gjd = juldat (mm, dd, yyyy, ut)
		    gotjd = TRUE
		    }
		}

# HJD direct from header keyword
	    call clgstr ("keyhjd",keyhjd,SZ_LINE)
	    if (strlen (keyhjd) > 0) {
		call imgdpar (specim, keyhjd, hjd)
		}

# UT from separate keyword
	    if (!gottime) {
		call clgstr ("keytime",keyword,SZ_LINE)
		if (strlen (keyword) > 0) {
		    call imgdpar (specim,keyword, ut)
			call pargi (dd)
			call pargi (dd)
		    if (ut == dindef)
			ut = 0.d0
		    }
		else
		    ut = 0.d0
		gottime = TRUE
		}

# Exposure time in seconds
	    exp = 0.
	    call clgstr ("keyexp",keyword,SZ_LINE)
	    call imgdpar (specim,keyword, exp)
	    dt = exp * 0.5d0 / 3600.d0

	    call clgstr ("keywhen",when,SZ_LINE)
	    reltime = strdic (when,when,SZ_LINE,UT_WHEN)
	    if (reltime == UT_START) {
		ut1 = ut
		ut = ut + dt
		ut2 = ut + dt
		gjd = gjd + dt / 24.d0
		}
	    else if (reltime == UT_END) {
		ut2 = ut
		ut = ut - dt
		ut1 = ut - dt
		gjd = gjd - dt / 24.d0
		}
	    else {
		ut1 = ut - dt
		ut2 = ut + dt
		}
	    if (debug ) {
		call printf("UT start: %0h, mid: %0h, end: %0h, exp: %d")
		    call pargd (ut1)
		    call pargd (ut)
		    call pargd (ut2)
		    call pargd (exp)
		if (strlen (keyword) > 0) {
		    call printf (" from %s\n")
			call pargstr (keyword)
		    }
		else
		    call printf ("\n")
		}

# Print date, time, and pointing direction
	    if (debug) {
		call printf("%04d-%3.3s-%02d %0h UT\n")
		    call pargi (yyyy)
		    call pargstr (month[(mm-1) * SZ_MONTH + 1])
		    call pargi (dd)
		    call pargd (ut)
		}
	    if (mm != dindef && dd != dindef && yyyy != INDEFI) {
		gotdate = TRUE
		}
	    else {
		gotdate = FALSE
		}
	    }
	else {
	    if (debug) {
		call printf("Julian date: %.4f\n")
		    call pargd (gjd)
		}
	    gotdate = TRUE
	    }

# If no Julian Date yet, compute it from time of observation
	if (!gotjd)
       	    gjd = juldat (mm, dd, yyyy, ut)

# If solar sky spectrum, use direction of sun as pointing direction
	if (specsky) {
	    subgrav = FALSE
	    subgrav = clgetb ("subgrav")
	    call suncoords (gjd, ra, dec)

# If solar sky spectrum, coordinates are equinox of date
	    eq =  (gjd - 2415020.0 - 0.313)/365.24219572 + 1900.0;
	    if (debug) {
		call printf("RA: %0.3h, Dec: %0.2h  %9.4f (sun)\n")
		    call pargd (ra)
		    call pargd (dec)
		    call pargd (eq)
		}

# Save solar coordinates to image header
	    if (savebcv) {
		call sprintf (rstr, 32, "%.3h")
		    call pargd (ra)
		call sprintf (dstr, 32, "%.3h")
		    call pargd (dec)
		call imastr (specim,"RA      ",rstr)
		call imastr (specim,"DEC     ",dstr)
		call imaddd (specim,"EPOCH   ",eq)
		if (debug) {
		    call eprintf ("BCVCORR: writing SGRV = %.4f to %s\n")
			call pargd (sgrv)
			call pargstr (specpath)
		    }
		}
	    }

	else if (verbose && printmode == 1) {
	    call printf("RA: %0h, Dec: %0h  %6.1f\n")
		call pargd (ra)
		call pargd (dec)
		call pargd (eq)
	    }

# Location of observatory
bcvcomp_
	obsname[1] = EOS
	dlong = dindef
	dlat = dindef
	dalt = 0.d0
	call clgstr ("obsname",obsname,SZ_LINE)
	call clgstr ("keyobs",keyobs,SZ_LINE)

# Read position from image header
	if (specim!=NULL && specim!=ERR && strncmp (obsname,"file",4)==0) {
 	    call imgspar (specim,keyobs,obsname,SZ_LINE)
	    if (obsname[1] != EOS) {
		call clgstr ("keylong",keyword,SZ_LINE)
		if (strlen (keyword) > 0) {
 		    call imgdpar (specim,keyword,dlong)
		    call clgstr ("keylat",keyword,SZ_LINE)
 		    call imgdpar (specim,keyword,dlat)
		    if (dlat == dindef || dlong == dindef)
			obs = obsvopen (obsname, debug)
		    else {
			call clgstr ("keyobs",keyword,SZ_LINE)
 			call imgspar (specim,keyword,obsname,SZ_LINE)
			call clgstr ("keyalt",keyword,SZ_LINE)
 			call imgdpar (specim,keyword,dalt)
			}
		    }
		else
		    obs = obsvopen (obsname, debug)
		if (obs != NULL) {
		    if (debug)
			call obslog (obs,"BCVCORR","latitude longitude altitude",STDOUT)
		    dlat = obsgetd (obs, "latitude")
		    dlong = obsgetd (obs, "longitude")
		    dalt = obsgetd (obs, "altitude")
		    call obsclose (obs)
		    }
		}
	    }

# Read position from BCVCORR parameters
	if (dlong == dindef || dlat == dindef) {
	    dlong = clgetd ("obslong")
	    dlat = clgetd ("obslat")
	    dalt = clgetd ("obsalt")
	    }

# Read position from IRAF observatory database using name from file
	if ((dlong == dindef || dlat == dindef) && keyobs[1] != EOS) {
 	    call imgspar (specim,keyobs,obsname,SZ_LINE)
	    if (obsname[1] != EOS) {
        	obs = obsvopen (obsname, debug)
		if (obs != NULL) {
		    if (debug)
			call obslog (obs,"BCVCORR","latitude longitude altitude",STDOUT)
		    dlat = obsgetd (obs, "latitude")
		    dlong = obsgetd (obs, "longitude")
		    dalt = obsgetd (obs, "altitude")
		    call obsclose (obs)
		    }
		}
	    }

# Read position from IRAF observatory database using name from parameter list
	if ((dlong == dindef || dlat == dindef) && obsname[1] != EOS) {
            obs = obsvopen (obsname, debug)
	    if (obs != NULL) {
		if (debug)
		    call obslog (obs, "BCVCORR", "latitude longitude altitude", STDOUT)
		dlat = obsgetd (obs, "latitude")
		dlong = obsgetd (obs, "longitude")
		dalt = obsgetd (obs, "altitude")
		call obsclose (obs)
		}
	    }

	if (verbose && printmode == 1) {
	    call printf("%s lat %0h , long %0h, alt %.1f\n")
		call pargstr (obsname)
		call pargd (dlat)
		call pargd (dlong)
		call pargd (dalt)
	    }
	call flush (STDOUT)

# Reformat so juldat and heliocvel can use the info
	rhcv = 0.
	rbcv = 0.
	if (gotdate && gottime) {

	# Calculate Julian Date
	    if (yyyy == INDEFI || mm == INDEFI || dd == INDEFI || ut == dindef)
		call jd2ut (gjd, yyyy, mm, dd, ut)
	    else if (gjd == 0.0 || gjd == dindef) {
        	gjd = juldat (mm, dd, yyyy, ut)
		}
	    if (verbose && printmode == 1) {
		call printf("Julian date is %.5f at %4d-%3.3s-%02d %0h UT\n")
		    call pargd (gjd)
		    call pargi (yyyy)
		    call pargstr (month[(mm-1) * SZ_MONTH + 1])
		    call pargi (dd)
		    call pargd (ut)
		call flush (STDOUT)
		}

	# Compute the Julian Date at the time when the light reached the Sun
	    if ((hjd == 0.0 || hjd == dindef) &&
		ra != dindef && dec != dindef) {
		call jd2hjd (ra, dec, gjd, lt, hjd)
		}
	    if (verbose && printmode == 1) {
		if (specsky)
		    call printf("Sun at ra %0h dec %0h eq %.4d\n")
		else
		    call printf("Object at ra %0.3h dec %0.2h eq %.1d\n")
		    call pargd (ra)
		    call pargd (dec)
		    call pargd (eq)
		call printf("Heliocentric Julian date: %.5f\n")
		    call pargd (hjd)
		call flush (STDOUT)
		}

	# Compute radial velocity corrections
	    call bcv (gjd, dlong,dlat,dalt, ra,dec,eq, dbcv, dhcv, dgcv)

	# Solar gravitional redshift if observing twilight sky
	    if (specsky)
		sgrv = 0.636d0
	    else
		sgrv = 0.d0

	    if (verbose && printmode == 1) {
		call printf ("gbcvel = %.4f  ghcvel = %.4f  geovel = %.4f\n")
		    call pargd (dbcv)
		    call pargd (dhcv)
		    call pargd (dgcv)
		if (specsky) {
		    if (subgrav) {
			call printf ("solar gravitational redshift removed = %.4f\n")
			    call pargd (sgrv)
			}
		    else {
			call printf ("solar gravitational redshift is = %.4f\n")
			    call pargd (sgrv)
			}
		    }
		call flush (STDOUT)
		}

	# Add up both corrections
	    rbcv = dbcv + dgcv
	    rhcv = dhcv + dgcv
	    if (subgrav) {
		rbcv = rbcv - sgrv
		rhcv = rhcv - sgrv
		}
	    if (verbose && printmode == 1) {
		call printf ("bcv = %.4f  hcv = %.4f computed\n")
		    call pargr (rbcv)
		    call pargr (rhcv)
		fbcv = dindef
		fhcv = dindef
		if (specim != NULL && specim != ERR) {
		    call imgdpar (specim, "BCV", fbcv)
		    call imgdpar (specim, "HCV", fhcv)
		    if (fbcv == dindef && fhcv == dindef)
			call printf ("No BCV or HCV in file header\n")
		    else {
			if (fbcv != dindef) {
			    call printf ("bcv = %.4f ")
				call pargd (fbcv)
			    }
			if (fhcv != dindef) {
			    call printf (" hcv = %.4f ")
				call pargd (fhcv)
			    }
			call printf ("from file\n")
			}
		    }
		call flush (STDOUT)
		}
	    if (verbose && printmode == 2) {
		call printf ("%.4f\n")
		    call pargr (rbcv)
		call flush (STDOUT)
		}
	    call clputr ("bcv", rbcv)
	    call clputr ("hcv", rhcv)

# Check for writability if results are being saved in image header
	    IM_UPDATE(specim) = NO
	    if (savejd) {
		if (imaccess (specpath, READ_WRITE) == NO) {
		call eprintf ("BCVCORR: cannot write to %s; not saving results\n")
		    call pargstr (specpath)
		    IM_UPDATE(specim) = NO
		    }
		else {

#		Write midtime of observation in UT
		    utm = dindef
		    call imgdpar (specim,"UTMID   ", utm)
		    if (utm == dindef) {
			call sprintf (uts, 32, "%0h")
			    call pargd (ut)
			call imastr (specim, "UTMID", uts)
			if (debug) {
			    call eprintf ("BCVCORR: writing UTMID = %.1h to %s\n")
				call pargd (ut)
				call pargstr (specpath)
			    }
			}
		    else {
			if (debug) {
			    call eprintf ("BCVCORR: writing UTMID = %.1h to %s\n")
				call pargd (utm)
				call pargstr (specpath)
			    }
			}

#		Write midtime of observation as geocentric Julian Date
		    if (imaccf (specim, "GJDN") == NO) {
			call imaddd (specim,"GJDN",gjd)
			if (debug) {
			    call eprintf ("BCVCORR: writing GJDN = %.5f to %s\n")
				call pargd (gjd)
				call pargstr (specpath)
			    }
			}

#		Write midtime of observation as Heliocentric Julian Date
		    if (imaccf (specim, "HJDN") == NO) {
			call imaddd (specim,"HJDN",hjd)
			if (debug) {
			    call eprintf ("BCVCORR: writing HJDN = %.5f to %s\n")
				call pargd (hjd)
				call pargstr (specpath)
			    }
			}
		    IM_UPDATE(specim) = YES
		    }
		}

	    if (savebcv) {
		if (imaccess (specpath, READ_WRITE) == NO) {
		call eprintf ("BCVCORR: cannot write to %s; not saving results\n")
		    call pargstr (specpath)
		    IM_UPDATE(specim) = NO
		    }
		else {

#		Write barycentric velocity correction
		    if (mspec == 0)
			call imaddr (specim,"BCV     ",rbcv)
		    else {
			call sprintf (keyword,SZ_HKWORD,"APVEL%d")
			    call pargi (mspec)
			call sprintf (str,SZ_LINE," 0.0 0.0 _ %.3f %.3f")
			    call pargr (rbcv)
			    call pargr (rhcv)
			call imastr (specim,keyword,str)
			}
		    if (debug) {
			call eprintf ("BCVCORR: writing BCV = %.4f to %s\n")
			    call pargr (rbcv)
			    call pargstr (specpath)
			}

#		Write heliocentric velocity correction
		    if (mspec == 0)
			call imaddr (specim,"HCV     ",rhcv)
		    if (debug) {
			call eprintf ("BCVCORR: writing HCV = %.4f to %s\n")
			    call pargr (rhcv)
			    call pargstr (specpath)
			}

#		Write solar gravitational redshift correction if twilight sky
		    if (specsky)
			call imaddd (specim,"SGRV    ",sgrv)
		    if (debug) {
			call eprintf ("BCVCORR: writing SGRV = %.4f to %s\n")
			    call pargd (sgrv)
			    call pargstr (specpath)
			}
		    IM_UPDATE(specim) = YES
		    }
		}
	    }

# Close current image and move on to next image
	if (specim != ERR && specim != NULL)
	    call imunmap (specim)
	if (lspec > 0)
	    go to newspec_

#  Close spectrum list
endvc_	if (lspec > 0)
	    call imtclose (speclist)

end

double procedure juldat (mm, dd, yyyy, ut)

int	mm	# Month
int	dd	# Day of month
int	yyyy	# Year
double	ut	# Universal Time

double	julday
int	jday
#int	igreg, jy, jm, ja

begin

	jday = ( 1461 * ( yyyy + 4800 + ( mm - 14 ) / 12 ) ) / 4 +
	       ( 367 * ( mm - 2 - 12 * ( ( mm - 14 ) / 12 ) ) ) / 12 -
	       ( 3 * ( ( yyyy + 4900 + ( mm - 14 ) / 12 ) / 100 ) ) / 4 +
	         dd - 32075

#        #igreg = 15 + 31 * (10 + 12 * 1582)
#
#	if (yyyy < 0)
#	    return 0.d0
#
#	if (mm > 2) {
#	    jy = yyyy
#	    jm = mm + 1
#	    }
#	else {
#	    jy = yyyy - 1
#	    jm = mm + 13
#	    }
#
#	julday = int (365.25*jy) + int (30.6001*jm) + dd + 1720995

#	if (dd+31*(mm+12*yyyy) >= igreg) {
#	    ja = int (0.01 * jy)
#	    julday = julday + 2 - ja + int (0.25 * ja)
#	    }

	julday = double (jday) - 0.5d0 + (ut / 24.d0)
	return julday
end


procedure jd2ut (jd, yyyy, mm, dd, ut)

double	jd	#Julian Date
int	yyyy	# Year	(returned)
int	mm	# Month	(returned)
int	dd	# Day	(returned)
double	ut	# UT	(returned)

double	tsec	# Seconds since 1/1/1950 0:00
int	hr	# Hours
int	mn	# Minutes
double	sec	# Seconds

double	t, days
int	ihms,nc,nc4,nly,ny,m,im
define	setmonth_ 10

begin

# Convert Julian Date to seconds since 1950 Jan 1 0:00
	tsec = (jd - 2433282.5d0) * 8.64d4

#  Calculate time of day (hours, minutes, seconds, .1 msec)
	t = dint ((tsec + 61530883200.d0) * 1.d4 + .5d0)
	hr = idint (dmod (t/36.d6,24.d0))
	mn = idint (dmod (t/60.d4,60.d0))
	if (tsec >= 0) {
	    ihms = idint (dmod (tsec+1.d-6,1.d0) * 1.d4)
	    sec = dmod (tsec+1.d-6,60.d0)
	    }
	else {
	    ihms = idint (dmod (tsec-1.d-6,1.d0) * 1.d4)
	    sec = dmod (tsec-1.d-6,60.d0)
	    }
	ut = double (hr) + (double (mn) / 60.d0) + (sec / 3.6d3)

#  Calculate number of days since 0 hr 0/0/0000
	days = dint ((t / 864.0d6) + 1.0d-6)

#  Number of leap centuries (400 years)
	nc4 = idint ((days / 146097.0d0) + 1.0d-5)

#  Number of centuries since last /400
	days = days - (146097.d0 * double (nc4))
	nc = idint ((days / 36524.d0) + 1.0d-4)
	If (nc > 3) nc = 3

#  Number of leap years since last century
	days = days - (36524.d0 * nc)
	nly = idint ((days / 1461.d0) + 1.0d-10)

#  Number of years since last leap year
	days = days - (1461.d0 * double (nly))
	ny = idint ((days / 365.d0) + 1.d-8)
	If (ny > 3) ny = 3

#  Calculate day of month and month
	days = days - (365.d0 * double (ny))
	dd = idint (days + 1.0d-8) + 1
	mm = 0
	do m = 1, 12 {
	    im = mod ((m + ((m - 1) / 5)), 2)
	    if (dd - 1 < im + 30)
		go to setmonth_
	    dd = dd - im - 30
	    }
setmonth_
	mm  =  mod (m+1, 12) + 1

# Year
	yyyy = nc4*400 + nc*100 + nly*4 + ny + m/11

	return
end


# Compute RA and Dec of sun in apparent coordinates
procedure suncoords (jd, ra, dec)

double jd	# Julian Date
double ra	# Right ascension of sun (hours returned)
double dec	# Declination of sun (degrees returned)

double lsun, mlsun, masun
double secliptic	# sine and cosine of obliquity
double cecliptic	# of ecliptic (23 deg. 27 min.)
double pi, pi2, pi32, d2r, r2d, tpi
double du, du1

begin

	pi = 3.14159265358979323846d0
	pi2 = pi * 0.5d0
	pi32 = pi * 1.5d0
	d2r = pi / 180.d0
	r2d = 180.d0 / pi
	tpi = 2.d0 * pi
	cecliptic = 0.9174650d0
	secliptic = 0.3978166d0

# Days from 0 Jan 1900 12:00 noon
	du = jd - 2415020.d0
	du1 = du / 1.d4

# Mean longitude of sun in degrees
	mlsun = dmod (279.696678d0 + 0.9856473354d0*du + 2.267d-5*du1*du1, 3.6d2)

# Compute mean anomaly of sun in degrees and convert to radians
	masun = dmod (358.475833d0 + 0.9856002670d0*du -
		   du1 * (1.12d-5 * du1 + 7.d-8 * du1*du1), 3.6d2)
	masun = d2r * masun

# Longitude of sun in degrees
	lsun = (0.020d0 * sin (2.0 * masun)) +
	       (1.916d0 * sin (masun)) + mlsun

# Convert from degrees to radians
	mlsun = d2r * mlsun
	lsun = d2r * lsun

# Right ascension of sun
	ra = atan (tan (lsun) * cecliptic)

# Restore to 4 quad, based on lsun
	if (lsun > pi32)
	    ra = ra + tpi
	else if (lsun > pi2)
	    ra = ra + pi

# Declination of sun
	dec = asin (sin (lsun) * secliptic)

# Convert right ascension from radians to hours
	ra = (ra * r2d) / 15.d0

# Convert declination from radians to degrees
	dec = dec * r2d

	return
end

# Jul 15 1995	New task
# Oct 18 1995	If observatory position not indef, use parameters
# Oct 18 1995	Always print RA, Dec, Equinox if verbose
# Dec  4 1995	Initialize altitude to 0 if reading from keyword

# Jan 12 1996	Allow position and time entry from parameter file, not image
# May 28 1996	If SPECNUM > 0, write to APVELnnn instead of BCV keyword

# Jan  8 1997	Compute heliocentric julian date
# Jan 17 1997	Compute BCV at geocentric Julian Date
# Jan 30 1997	Minimize repition in output
# Oct 23 1997	Fix egregious error which reversed BCV and HCV in call to BCV
# Dec 17 1997   Change order of date from dd-mm-yyyy to yyyy-mmm-dd
# Dec 17 1997   Print names of Julian Date keywords, if written to header
# Dec 17 1997	Always write keywords expected by other RVSAO programs

# Feb 11 1998	Separate date and time existence checks
# Feb 11 1998	Use end or start of observation if center not computable
# Feb 18 1998	Print values when verifying writing to header
# Mar  4 1998	Compute UT if only Julian date is given

# Mar 29 1999	Read IRAF observatory database if name keyword but not lat,long
# Dec  8 1999	Fix bug which failed to use image observatory name

# Jun  2 2003	Add keypos and keydate so task parameter input works

# Aug 24 2004	Fix minor bug in error message 

# Feb  1 2005	Always figure out source of values, not just if no image
# Feb  1 2005	Read observatory name from header along with position
# Feb  2 2005	Set dindef to INDEFD
# Mar 16 2005	Convert RA to hours from degrees if not sexigesimal

# Jan 27 2006	Add specsky parameter use sun position if twilight sky
# Feb 16 2006	Change parameters to accomodate JD, MJD, or FITS date/time
# Aug  1 2006	Fix bug which printed hcv instead of bcv in debug mode
# Aug  2 2006	Write HCV as well as BCV to image header
# Aug  2 2006	Subtract solar gravitational redshift for twilight skies
# Aug  2 2006	Write sun RA, Dec, equinox to image header for twilight skies

# Jun 20 2007	Print only BCV if not in verbose mode
# Jun 20 2007	Fix bug when reading observatory code from file keyword
# Jun 21 2007	Return BCV and HCV as parameters; add printmode parameter
# Jun 21 2007	Add savejd parameter to separate saving of times and velocities

# Apr 21 2009	Make gravitational redshift subtraction optional
