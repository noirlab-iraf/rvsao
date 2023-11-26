# July 7, 1995
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1995 Smithsonian Astrophysical Observatory
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

procedure t_bcvcorr ()

pointer	specim		# image structure for spectrum
char	specfile[SZ_PATHNAME]	# Object spectrum file name
char	specdir[SZ_PATHNAME]	# Directory for object spectra
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	obsname[SZ_LINE]	# Observatory name
char	keyword[SZ_LINE]	# Observatory name
pointer	speclist	# List of spectrum files
bool	savebcv,savebcv0	# Save velocity correction in data file header
int	mode		# Mode (READ_ONLY, READ_WRITE)
pointer	obs		# pointer to observatory structure
bool	gottime		# True if observation time obtained
bool	gotpos		# True if position of observed object obtained
#bool	debug
bool	verbose, debug, clgetb()
int	mm, dd, yyyy, ldir, ndim, k
double	ra, dec, eq, exp, dindef
double	ut1,ut2,ut, dlong, dlat, dalt, gjd
double	dc[3],dcc[3],dprema[3,3],dvelh[3],dvelb[3]
double	djd,dpi,daukm,dctrop,dcbes,dc1900,dct0,dtr,dst
double	dlongs, dlats, deqt, dras, decs, dha, dra2, dec2
double	dgcvel, dbcvel, dhcvel, fbcv, fhcv, tpi, darg
real	rbcv, rhcv
int	strlen(), imtgetim(), imaccess(), strncmp()
double	juldat(), clgetd()
pointer	imtopenp()
pointer	immap()
errchk	immap()
double	obsgetd()
pointer	obsopen()
errchk	obsopen

define	newspec_	10
define	endvc_	90

begin
	dpi = 3.1415926535897932d0
	tpi = 2.d0 * dpi
	daukm  =1.4959787d08
	dctrop = 365.24219572d0
	dcbes = 0.313d0
	dc1900 = 1900.0d0
	dctrop = 365.24219572d0
	dcbes = 0.313d0
	dc1900 = 1900.0d0
	dtr = dpi / 180.0d0

# Spectra for which to compute solar system center velocity corrections
	speclist = imtopenp ("spectra")
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	verbose = clgetb ("verbose")
	debug = clgetb ("debug")
	savebcv0 = clgetb ("savebcv")

# Get next object spectrum file name from the list
newspec_
	savebcv = savebcv0
	if (imtgetim (speclist, specfile, SZ_PATHNAME) == EOF)
	   go to endvc_

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

# Load spectrum header
	mode = READ_WRITE
        iferr (specim = immap (specpath, mode, LEN_USER_AREA)) {
            mode = READ_ONLY
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
	    }
        ndim = IM_NDIM(specim)
	if (debug) {
	    call printf ("%s: %d x %d x %d %d-D image\n")
		call pargstr (specpath)
		call pargi (IM_LEN(specim,1))
		call pargi (IM_LEN(specim,2))
		call pargi (IM_LEN(specim,3))
		call pargi (ndim)
	    }

	dindef = INDEFD

# Direction for velocity correction
	call clgstr ("keyra",keyword,SZ_LINE)
	ra = dindef
 	call imgdpar (specim, keyword, ra)
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

# Midtime of observation
	call clgstr ("keyjd",keyword,SZ_LINE)
	gjd = dindef
	call imgdpar (specim, keyword, gjd)
	if (gjd == dindef) {
	    call clgstr ("keymid",keyword,SZ_LINE)
	    ut = dindef
	    call imgdpar (specim,keyword, ut)
	    if (ut == dindef) {
		call clgstr ("keyend",keyword,SZ_LINE)
		ut2 = dindef
		call imgdpar (specim, keyword, ut2)
		ut1 = dindef
		call clgstr ("keystart",keyword,SZ_LINE)
		call imgdpar (specim,keyword, ut1)
		exp = 0.
		call clgstr ("keyexp",keyword,SZ_LINE)
		call imgdpar (specim,keyword, exp)
		if (debug ) {
		    call printf("UT: %h, UTOPEN: %h, UTMID: %h, EXP: %d\n")
			call pargd (ut2)
			call pargd (ut1)
			call pargd (ut)
			call pargd (exp)
		    }
		if (ut1 != dindef && ut2 != dindef) {
		    ut = (ut1 + ut2) * 0.5
		    gottime = TRUE
		    }
		else if (ut2 != dindef && exp > 0.d0) {
		    ut = ut2 - (exp * 0.5 / 3600.)
		    gottime = TRUE
		    }
		else if (ut1 != dindef && exp > 0.d0) {
		    ut = ut1 + (exp * 0.5 / 3600.)
		    gottime = TRUE
		    }
		else {
		    ut = dindef
		    gottime = FALSE
		    }
		}
	    else {
		ut = dindef
		gottime = FALSE
		}

# Date of observation, if midtime is not Julian Date
	    mm = INDEFI
	    dd = INDEFI
	    yyyy = INDEFI
	    call clgstr ("keydate",keyword,SZ_LINE)
	    call imgdate (specim, keyword, mm, dd, yyyy)

# Print date, time, and pointing direction
	    if (verbose ) {
		call printf("%d/%d/%d %h UT  ra: %h, dec: %h  %6.1f\n")
		    call pargi (mm)
		    call pargi (dd)
		    call pargi (yyyy)
		    call pargd (ut)
		    call pargd (ra)
		    call pargd (dec)
		    call pargd (eq)
		}
	    if (mm != dindef && dd != dindef && yyyy != INDEFI)
		gottime = TRUE
	    else
		gottime = FALSE
	    }
	else {
	    if (verbose ) {
		call printf("Julian date: %.4f\n")
		    call pargd (gjd)
		}
	    gottime = TRUE
	    }

# Location of observatory
	obsname[1] = EOS
	dlong = dindef
	dlat = dindef
	dalt = dindef
	call clgstr ("obsname",obsname,SZ_LINE)
	if (strncmp (obsname,"file",4) == 0) {
	    call clgstr ("keyobs",keyword,SZ_LINE)
 	    call imgspar (specim,keyword,obsname,SZ_LINE)
	    if (obsname[1] != EOS) {
		call clgstr ("keylong",keyword,SZ_LINE)
 		call imgdpar (specim,keyword,dlong)
		call clgstr ("keylat",keyword,SZ_LINE)
 		call imgdpar (specim,keyword,dlat)
		call clgstr ("keyalt",keyword,SZ_LINE)
 		call imgdpar (specim,keyword,dalt)
		}
	    }
	else {
            obs = obsopen (obsname)
            call obslog (obs, "BCVCORR", "latitude longitude altitude", STDOUT)
            dlat = obsgetd (obs, "latitude")
            dlong = obsgetd (obs, "longitude")
            dalt = obsgetd (obs, "altitude")
            call obsclose (obs)
	    }
	if (dlong == dindef && dlat == dindef &&dalt == dindef) {
	    dlong = clgetd ("obslong")
	    dlat = clgetd ("obslat")
	    dalt = clgetd ("obsalt")
	    }
	if (verbose) {
	    call printf("%s lat %h , long %h, alt %.1fm\n")
		call pargstr (obsname)
		call pargd (dlat)
		call pargd (dlong)
		call pargd (dalt)
	    }
         
# Reformat so juldat and heliocvel can use the info
	rhcv = 0.
	rbcv = 0.
	if (gottime && gotpos) {

	# Calculate julian date
	    if (gjd == dindef) {
        	djd = juldat (mm,dd,yyyy) - 0.5d0
        	djd = djd + ut / 24.0d0
		}
	    else
		djd = gjd
	    if (verbose ) {
		call printf("Julian date: %.4f\n")
		    call pargd (djd)
		}

	# Calculate local sidereal time
	    dlongs = dlong*dtr
	    call sidtim (djd,dlongs,dst)

	# Precess r.a. and dec. to mean equator and equinox of date (deqt)
	    deqt = (djd - dct0 - dcbes) / dctrop + dc1900
	    dras = ra * 15.0d0 * dtr
	    decs = dec * dtr
	    dc(1) = dcos(dras) * dcos(decs)
	    dc(2) = dsin(dras) * dcos(decs)
	    dc(3) =              dsin(decs)
	    call pre (eq,deqt,dprema)
	    do k = 1, 3 {
		dcc[k]=dc[1]*dprema[k,1]+dc[2]*dprema[k,2]+dc[3]*dprema[k,3]
		}

	    dra2 = datan2 (dcc[1],dcc[2])
	    if (dcc[1] != 0.d0) {
		darg = dcc[2] / dcc[1]
		dra2 = datan (darg)
		if (dcc[1] < 0.d0)
		    dra2 = dra2 + dpi
		else if (dcc[2] < 0.d0)
		    dra2 = dra2 + tpi
		}
	    else {
		if (dcc[2] > 0.d0)
		    dra2 = 0.5d0 * dpi
		else
		    dra2 = 1.5d0 * dpi
		}
	    dec2 = dasin(dcc(3])

	# Calculate hour angle = local sidereal time - r.a.
	    dha = dst - dra2
	    dha = dmod (dha + tpi , tpi)

	# Calculate observer's geocentric velocity
	# (altitude assumed to be zero)
	    dlats = dlat*dtr
	    call geovel (dlats,dalt,dec2,-dha,dgcvel)
	    if (debug) {
		call printf ("BCVCORR: RA = %.3h  Dec = %.2h\n")
		    call pargd (ra)
		    call pargd (dec)
		call printf ("BCVCORR: current RA = %.3h  Dec = %.2h\n")
		    call pargd (dra2)
		    call pargd (dec2)
		call printf ("BCVCORR: current ST = %.3h  HA = %.3h\n")
		    call pargd (dst)
		    call pargd (dha)
		call printf ("BCVCORR: Obs lat = %.3h  Alt = %.2f\n")
		    call pargd (dlat)
		    call pargd (dalt)
		}

	# Calculate components of earth's barycentric veolcity,
	# dvelb(i), i=1,2,3  in units of a.u./s
	    call barvel (djd,deqt,dvelh,dvelb)

	# Project barycentric velocity to the direction of the star, and
	# convert to km/s
	    dbcvel = 0.0d0
	    dhcvel = 0.0d0
	    do k = 1, 3 {
		dbcvel = dbcvel + dvelb[k]*dcc[k]*daukm
		dhcvel = dhcvel + dvelh[k]*dcc[k]*daukm
		}
	    if (verbose) {
		call printf ("bcvel = %.4f  hcvel = %.4f  geovel = %.4f\n")
		    call pargd (dbcvel)
		    call pargd (dhcvel)
		    call pargd (dgcvel)
		}

	# Add up both corrections
	    rbcv = dbcvel + dgcvel
	    rhcv = dhcvel + dgcvel
	    if (verbose) {
		call printf ("bcv = %.4f  hcv = %.4f computed\n")
		    call pargr (rbcv)
		    call pargr (rhcv)
		fbcv = dindef
		call imgdpar (specim, "BCV     ", fbcv)
		fbcv = dindef
		call imgdpar (specim, "HCV     ", fhcv)
		call printf ("bcv = %.4f  hcv = %.4f from file\n")
		    call pargd (fbcv)
		    call pargd (fhcv)
		}

# Check for writability if results are being saved in image header
	    if (savebcv) {
		if (imaccess (specpath, READ_WRITE) == NO) {
		call eprintf ("BCVCORR: cannot write to %s; not saving results\n")
		    call pargstr (specpath)
		    IM_UPDATE(specim) = NO
		    }
		else {
		    call imaddr (specim,"BCV     ",rbcv)
		    IM_UPDATE(specim) = YES
		    if (debug) {
			call eprintf ("BCVCORR: writing BCV to %s\n")
			    call pargstr (specpath)
			}
		    }
		}
	    else
		IM_UPDATE(specim) = NO
	    }

# Close current image and move on to next image
	if (specim != ERR && specim != NULL)
	    call imunmap (specim)
	go to newspec_

#  Close spectrum list
endvc_	call imtclose (speclist)

end

double procedure juldat (mm,dd,yyyy)

int	mm
int	dd
int	yyyy

double	julday
int	igreg, jy, jm, ja

begin

        igreg = 15 + 31 * (10 + 12 * 1582)

	if (yyyy < 0)
	    return 0.d0

	if (mm > 2) {
	    jy = yyyy
	    jm = mm + 1
	    }
	else {
	    jy = yyyy - 1
	    jm = mm + 13
	    }

	julday = int (365.25*jy) + int (30.6001*jm) + dd + 1720995

	if (dd+31*(mm+12*yyyy) >= igreg) {
            ja = int (0.01 * jy)
            julday = julday + 2 - ja + int (0.25 * ja)
            }
	return julday
end
