# File rvsao/Emvel/emrslts.x
# May 20, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Log radial velocity from emission line shift

include <imhdr.h>
include "rvsao.h"
include "emv.h"
define	SZ_MONTH	3
define	SZ_YEAR		40

procedure emrslts (specfile, mspec, image, rmode0)

char	specfile[ARB]	# Name of spectrum file
int	mspec		# Spectrum number in multispec file (else 0)
pointer	image		# Header structure for spectrum
int	rmode0		# Report format (1=normal,2=one-line...)

int	i, j
char	wch[8]		# Line weight indicator
double	emz		# z from emission lines
double	xcz		# z from cross-correlation
double	spz		# z from cross-correlation and emission
double	spexvel		# Emission line template cross-correlation velocity
double	spexerr		# Emission line template cross-correlation error
int	iczexvel	# Emission line template cross-correlation velocity
int	iczexerr	# Emission line template cross-correlation error
double	spexr		# Emission line template cross-correlation R-value
int	rmode		# Report format (1=normal,2=one-line)
char	instrument[10]
char	lbr, rbr
int	ibr, ibl

bool	verbose, clgetb()
double	ra, dec, ut, dj, hdj
double	twt, wt, exptime, alt, az, ha, lat, sha,cha,slat,clat,sdec,cdec,epoch
double	ac,dpi,tpi,rcon, rdec, rha, dindef, dwdp, wl1, pix1, pix2
int	mm,dd,yyyy, itemp, lspec
char	month[SZ_YEAR]
char	nolines[SZ_LINE]
char	nodate[SZ_LINE]
char	kra[8]
char	kha[8]
char	kdec[8]
int	icz,iczxc,iczem,icze,iczxce,iczeme
int	strcmp(), strlen()
char	ssvcor[16]
double	wcs_p2w(), wcs_w2p()
int	imaccf(), stridx()
char	spchar[8]

include "emv.com"
include "rvsao.com"
include "results.com"
include "lineinfo.com"

begin

# Return if no log file, including STDOUT
	if (nlogfd < 1)
	    return

	dindef = INDEFD
        lbr = char (91)
        rbr = char (93)
	ibr = 0
	ibl = 0
	verbose = clgetb ("verbose")
        if (!verbose) return
	call strcpy ("JanFebMarAprMayJunJulAugSepOctNovDec",month, SZ_YEAR)
	call strcpy ("No emission lines found\n", nolines, SZ_LINE)
	call strcpy (" (No observation date)", nodate, SZ_LINE)
	call strcpy ("RA", kra, 8)
	call strcpy ("DEC", kdec, 8)
	call strcpy ("HA", kha, 8)

	call juldate (image, ut, dj, hdj, debug)

	mm = INDEFI
	dd = INDEFI
	yyyy = INDEFI
	if (imaccf (image, "DATE-OBS") == YES)
	    call imgdate (image, "DATE-OBS", mm, dd, yyyy)
	else if (imaccf (image, "DATE") == YES)
	    call imgdate (image, "DATE", mm, dd, yyyy)

	if (spvel != dindef) {
	    icz = idnint (spvel)
	    icze = idnint (sperr)
	    }
	else {
	    icz = 0
	    icze = 0
	    }

	if (spxvel != dindef) {
	    iczxc = idnint (spxvel)
	    iczxce = idnint (spxerr)
	    }
	else {
	    iczxc = 0
	    iczxce = 0
	    }

	if (spevel != dindef) {
	    iczem = idnint (spevel)
	    iczeme = idnint (speerr)
	    }
	else {
	    iczem = 0
	    iczeme = 0
	    }

# Set emission line template velocity, if there is one
	spexvel = 0.0
	spexerr = 0.0
	spexr = 0.0
	iczexvel = 0
	iczexerr = 0
	if (strlen (cortemp) > 0 && ntemp > 0) {
	    do itemp = 1, ntemp {
		if (strcmp (cortemp,tempid[1,itemp]) == 0) {
		    spexvel = zvel[itemp]
		    spexerr = czerr[itemp]
		    iczexvel = int (zvel[itemp])
		    iczexerr = int (czerr[itemp])
		    spexr = czr[itemp]
		    }
		}
	    }
	else {
	    spexvel = spxvel
	    spexerr = spxerr
	    if (spxvel != dindef) {
		iczexvel = int (spxvel)
		iczexerr = int (spxerr)
		}
	    else {
		iczexvel = 0
		iczexerr = 0
		}
	    spexr = spxr
	    }

	lspec = strlen (specid)

# Set up spectrum number in specid
#	if (mspec > 0) {
#	    call sprintf (specid,SZ_PATHNAME,"%s[%d]")
#		call pargstr (specfile)
#		call pargi (mspec)
#	    }
#	else
#	    call strcpy (specfile,specid,SZ_PATHNAME)

	if (debug) {
	    call printf ("EMRSLTS: ready to print results\n")
	    call flush (STDOUT)
	    }

# Set separating character to tab if report mode is less than 0
	if (rmode0 >= 0) {
	    call strcpy (" ", spchar, 8)
	    rmode = rmode0
	    }
	else {
	    call strcpy ("\t", spchar, 8)
	    rmode = -rmode0
	    }

# Convert brackets to spaces in all one-line report modes
	if (rmode > 0) {
            ibr = stridx (rbr, specid)
            if (ibr > 0)
                specid[ibr] = ' '
            ibl = stridx (lbr, specid)
            if (ibl > 0)
                specid[ibl] = ' '
	    }

	if (rmode < 2 || rmode == 10) {

#	Read position of object from spectrum header
	    ra = dindef
	    call imgdpar (image,kra,ra)
	    dec = dindef
	    call imgdpar (image,kdec,dec)
	    epoch = 1950.d0
	    call imgdpar (image,"EPOCH",epoch)
	    call imgdpar (image,"EQUINOX",epoch)

	    call printf ("\n")
	    call flush (STDOUT)

#	Set type and source of velocity correction to center of solar system
	    switch (svcor) {
		case NONE:
		    call strcpy ("BCV",ssvcor,16)
		case HCV:
		    call strcpy ("HCV",ssvcor,16)
		case BCV:
		    call strcpy ("BCV",ssvcor,16)
		case FHCV:
		    call strcpy ("File HCV",ssvcor,16)
		case FBCV:
		    if (specvb) 
			call strcpy ("File BCV",ssvcor,16)
		    else
			call strcpy ("File HCV",ssvcor,16)
		default:
		    call strcpy ("BCV",ssvcor,16)
		}

#	Log object information
	    do j = 1, nlogfd {
		if (rmode < 2) {
		call fprintf(logfd[j], "File: %s\n")
		    call pargstr (specid)
		if (ra != dindef && dec != dindef) {
		    call fprintf(logfd[j], "Object: %s  RA: %011.2h Dec: %010.1h %.1f\n")
			call pargstr (specname)
			call pargd (ra)
			call pargd (dec)
			call pargd (epoch)
		    }
		else {
		    call fprintf(logfd[j], "Object: %s  (No position)\n")
			call pargstr (specname)
		    }

		if (debug) {
		    call flush (STDOUT)
		    }
		if (dd != INDEFI && mm != INDEFI && yyyy != INDEFI) {
		    call fprintf(logfd[j], "Observed %4d-%3.3s-%02d %011.2h = JD %.4f %s: %6.2f\n")
			call pargi (yyyy)
			call pargstr (month[(mm - 1) * SZ_MONTH + 1])
			call pargi (dd)
			call pargd (ut)
			call pargd (dj)
			call pargstr (ssvcor)
			call pargd (spechcv)
		    }
		else {
		    call fprintf(logfd[j], "%s %s: %6.2f\n")
			call pargstr (nodate)
			call pargstr (ssvcor)
			call pargd (spechcv)
		    }
		call flush (logfd[j])

	#  Print combined velocity
		if (spvel != dindef) {
		    spz = spvel / c0
		    call fprintf(logfd[j],
		    "Combined vel = %8.2f +- %7.2f km/sec, z= %6.4f \n")
			call pargd (spvel)
			call pargd (sperr)
			call pargd (spz)
		    call flush (logfd[j])
		    }

	# Print cross-correlation velocity
		if (spxvel != dindef) {
		    xcz = spxvel / c0
		    call fprintf(logfd[j],
		    "Correlation vel = %8.2f +- %7.2f km/sec, z= %6.4f R= %6.1f \n")
			call pargd (spxvel)
			call pargd (spxerr)
			call pargd (xcz)
			call pargd (spxr)
		    call flush (logfd[j])
		    }

	# Print emission line velocity
		if (spevel != dindef) {
		    emz = spevel / c0
		    call fprintf(logfd[j],
		    "Emission vel = %8.2f +- %7.2f km/sec, z= %6.4f for %d/%d lines\n")
			call pargd (spevel)
			call pargd (speerr)
			call pargd (emz)
			call pargi (nfit)
			call pargi (nfound)
		    call flush (logfd[j])
		    }
		}

	# Print emission line information if any were found
		if (nfound > 0) {
		    twt = 0.d0
		    do i = 1, nfound {
			wt = emparams[10,1,i]
			if (wt != 0.d0)
			    twt = twt + 1.d0 / (wt * wt)
			}
		    call fprintf(logfd[j], "\nLine  Rest lam  Obs. lam    ")
		    call fprintf(logfd[j], "Pixel     z         vel    dvel")
		    if (rmode == 1)
			call fprintf(logfd[j], "    eqw    wt\n")
		    else
			call fprintf(logfd[j], "    eqw    eqe    wt\n")
		    do i = 1, nfound {
			wtobs[i] = emparams[10,1,i]
			if (wtobs[i] > 0) {
			    if (override[i] < 0)
				call strcpy ("-",wch,8)
			    else if (override[i] > 0)
				call strcpy ("+",wch,8)
			    else
				call strcpy (" ",wch,8)
			    }
			else if (linedrop[i] > 0) {
			    call sprintf(wch,8,"X%d")
				call pargi (linedrop[i])
			    }
			else
			    call strcpy ("X",wch,8)
			call fprintf(logfd[j],"%4s   %7.2f   %7.2f  %7.2f")
			    call pargstr (nmobs[1,i])
			    call pargd (wlrest[i])
			    call pargd (wlobs[i])
			    call pargd (emparams[4,1,i])
			wt = emparams[10,1,i]
			if (wt != 0.d0)
			    wt = (1.d0 / (wt * wt)) / twt
			if (rmode == 1 || rmode == 10)
			    call fprintf (logfd[j],"  %6.4f  %8.2f %7.2f %6.2f %6.3f %s\n")
			else
			    call fprintf (logfd[j],"  %6.4f  %8.2f %7.2f %6.2f %6.2f %6.3f %s\n")
			    call pargd (emparams[9,1,i])
			    call pargd ((c0*emparams[9,1,i]) + spechcv)
			    call pargd (c0*emparams[9,2,i])
			    call pargd (emparams[7,1,i])
			if (rmode == 0)
			    call pargd (emparams[7,2,i])
			    call pargd (wt)
			    call pargstr (wch)
			}
		    }
		else {
		    call fprintf(logfd[j], nolines)
		    }
		}
	    }

# Single line report
	else if (rmode == 2 || rmode == 5) {
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strlen (instrument) == 0)
		call strcpy ("________", instrument, 8)
	    do j = 1, nlogfd {
		if (lspec > 32) {
		    call fprintf(logfd[j],"%s")
			call pargstr (specid)
		    }
		else {
		    call fprintf(logfd[j],"%-32.32s")
			call pargstr (specid)
		    }
		call fprintf(logfd[j],"%s%-8.8s%s%-16.16s%s%13.5f")
		    call pargstr (spchar)
		    call pargstr (instrument)
		    call pargstr (spchar)
		    call pargstr (IM_TITLE(image))
		    call pargstr (spchar)
		    call pargd (dj)
		call fprintf(logfd[j],"%s%6d%s%3d")
		    call pargstr (spchar)
		    call pargi (iczem)
		    call pargstr (spchar)
		    call pargi (iczeme)
		call fprintf(logfd[j],"%s%6d%s%3d%s%5.1f")
		    call pargstr (spchar)
		    if (rmode == 2) {
			call pargi (iczxc)
			call pargstr (spchar)
			call pargi (iczxce)
			call pargstr (spchar)
			call pargd (spxr)
			}
		    else {
			call pargi (iczexvel)
			call pargstr (spchar)
			call pargi (iczexerr)
			call pargstr (spchar)
			call pargd (spexr)
			}
		call fprintf(logfd[j],"%s%6d%s%3d%sE%s%2d%s%2d")
		    call pargstr (spchar)
		    call pargi (icz)
		    call pargstr (spchar)
		    call pargi (iczxce)
		    call pargstr (spchar)
		    call pargstr (spchar)
		    call pargi (nfound)
		    call pargstr (spchar)
		    call pargi (nfit)
		if (nfit > 0) {
		    do i = 1, nfound {
			wtobs[i] = emparams[10,1,i]
			if (wtobs[i] > 0 && override[i] >= 0) {
			    call fprintf(logfd[j], "%s%s(%7.2f)")
				call pargstr (spchar)
				call pargstr (nmobs[1,i])
				call pargd (wlrest[i])
			    }
			}
		    }
		call fprintf(logfd[j], "\n")
		}
	    }

# Single line report with one velocity per possible line
	else if (rmode == 3 || rmode == 6 || rmode == 9 || rmode == 11) {
	    call imgspar (image,"INSTRUMEN",instrument, 9)
	    if (strlen (instrument) == 0)
		call strcpy ("_________", instrument, 9)
	    call aclrd (velline,MAXREF)

#	Compute individual line velocities from shift
	    if (nfound > 0) {
		do i = 1, nfound {
		    wtobs[i] = emparams[10,1,i]
		    if (wtobs[i] > 0 && override[i] >= 0) {
			do j = 1, nref {
			    if (wlref[j] == wlrest[i]) {
				if (rmode == 9)
				    velline[j] = wlobs[i] - wlref[j]
				else if (rmode == 11) {
				    pix1 = wcs_w2p (wlref[j])
				    pix2 = wcs_w2p (wlobs[i])
				    velline[j] = pix2 - pix1
				    }
				else
				    velline[j] = (c0 * emparams[9,1,i]) + spechcv
				}
			    }
			}
		    }
		}

	# Print one-line report
	    do j = 1, nlogfd {
		if (lspec > 32) {
		    call fprintf(logfd[j],"%s")
			call pargstr (specid)
		    }
		else {
		    call fprintf(logfd[j],"%-32.32s")
			call pargstr (specid)
		    }
		call fprintf(logfd[j],"%s%-9.9s%s%-16.16s%s%13.5f")
		    call pargstr (spchar)
		    call pargstr (instrument)
		    call pargstr (spchar)
		    call pargstr (IM_TITLE(image))
		    call pargstr (spchar)
		    call pargd (dj)
		if (rmode < 9) {
		    call fprintf(logfd[j],"%s%8.2f%s%6.2f")
			call pargstr (spchar)
			call pargd (spevel)
			call pargstr (spchar)
			call pargd (speerr)
		    call fprintf(logfd[j],"%s%8.2f%s%6.2f%s%5.1f")
			call pargstr (spchar)
		    if (rmode == 3) {
			call pargd (spxvel)
			call pargstr (spchar)
			call pargd (spxerr)
			call pargstr (spchar)
			call pargd (spxr)
			}
		    else {
			call pargd (spexvel)
			call pargstr (spchar)
			call pargd (spexerr)
			call pargstr (spchar)
			call pargd (spexr)
			}
		    call fprintf(logfd[j],"%s%8.2f%s%6.2f%sE")
			call pargstr (spchar)
			call pargd (spvel)
			call pargstr (spchar)
			call pargd (sperr)
			call pargstr (spchar)
		    }
		call fprintf(logfd[j],"%s%2d%s%2d")
		    call pargstr (spchar)
		    call pargi (spnl)
		    call pargstr (spchar)
		    call pargi (spnlf)
		do i = 1, nref {
		    call fprintf (logfd[j],"%s%9.3f")
			call pargstr (spchar)
			call pargd (velline[i])
		    }
		call fprintf(logfd[j], "\n")
		}
	    }

# Single line report with velocity, error, height, width, and equivalent width
# for each emission line
	else if (rmode == 4 || rmode == 7) {
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strlen (instrument) == 0)
		call strcpy ("________", instrument, 8)
	    call aclrd (velline,MAXREF)
	    call aclrd (errline,MAXREF)
	    call aclrd (htline,MAXREF)
	    call aclrd (widline,MAXREF)
	    call aclrd (eqwline,MAXREF)
	    call aclrd (eqeline,MAXREF)

#	Compute individual line velocities from shift
	    if (nfound > 0) {
		do i = 1, nfound {
		    wtobs[i] = emparams[10,1,i]
		    if (wtobs[i] > 0 && override[i] >= 0) {
			do j = 1, nref {
			    if (wlref[j] == wlrest[i]) {
				velline[j] = (c0 * emparams[9,1,i]) + spechcv
				errline[j] = c0 * emparams[9,2,i]
				htline[j] = emparams[5,1,i]
				wl1 = wcs_p2w (emparams[4,1,i] - 1.d0)
				dwdp = wlobs[i] - wl1
				widline[j] = emparams[6,1,i] * dwdp
				eqwline[j] = emparams[7,1,i]
				eqeline[j] = emparams[7,2,i]
				}
			    }
			}
		    }
		}

	# Get sky position
	    ra = 99.d0
	    call imgdpar (image,kra,ra)
	    dec = 99.d0
	    call imgdpar (image,kdec,dec)
	    ha = 99.d0
	    call imgdpar (image,kha,ha)
	    lat = 99.d0
	    call imgdpar (image,"SITELAT",lat)

	    if (debug) {
		call printf ("EMRSLTS: RA: %011.2h Dec: %010.1h HA: %011.2h Lat: %010.1h\n")
		    call pargd (ra)
		    call pargd (dec)
		    call pargd (ha)
		    call pargd (lat)
		}

	# Compute altitude and azimuth
	    if (ha < 12.d0 && ha > -12.d0 && dec <= 90.d0 && dec >= -90.d0 &&
 		lat >= -90.d0 && lat <= 90.d0) {
		dpi = 3.1415926535897932d0
		tpi = 6.283185307179588d0
		rcon = 1.74532925199433d-2
		rha = ha * 15.d0 * rcon
		if (rha == 0)
		    rha = 1.d-5
		sha = sin (rha)
		cha = cos (rha)
		rdec = dec * rcon
		sdec = sin (rdec)
		cdec = cos (rdec)
		slat = sin (lat)
		clat = cos (lat)
		ac = (clat * sdec) - (slat * cdec * cha)
		az = atan (-cdec * sha / abs (ac))
		if (ac < 0) az = dpi - az
		if (az < 0) az = tpi + az
		ac = -cdec * sha / sin (az)
		alt = atan ((slat * sdec) + (clat * cdec * cha) / abs (ac))
		if (ac < 0) alt = dpi - alt
		az = az / rcon
		alt = alt / rcon
		}
	    else {
		alt = 0.d0
		az = 0.d0
		}

	    exptime = 0.d0
	    call imgdpar (image,"EXPTIME",exptime)

	# Print one-line report
	    do j = 1, nlogfd {
		if (lspec > 32) {
		    call fprintf(logfd[j],"%s")
			call pargstr (specid)
		    }
		else {
		    call fprintf(logfd[j],"%-32.32s")
			call pargstr (specid)
		    }
		call fprintf(logfd[j],"%s%-8.8s%s%-16.16s")
		    call pargstr (spchar)
		    call pargstr (instrument)
		    call pargstr (spchar)
		    call pargstr (IM_TITLE(image))
		call fprintf(logfd[j],"%s%011.2h%s%010.1h%s%5.1f%s%5.1f")
		    call pargstr (spchar)
		    call pargd (ra)
		    call pargstr (spchar)
		    call pargd (dec)
		    call pargstr (spchar)
		    call pargd (az)
		    call pargstr (spchar)
		    call pargd (alt)
		call fprintf(logfd[j],"%s%13.5f%s%6.1f")
		    call pargstr (spchar)
		    call pargd (dj)
		    call pargstr (spchar)
		    call pargd (exptime)
		call fprintf(logfd[j],"%s%8.2f%s%6.2f%s%8.2f%s%6.2f%s%5.1f")
		    call pargstr (spchar)
		    call pargd (spevel)
		    call pargstr (spchar)
		    call pargd (speerr)
		    call pargstr (spchar)
		    call pargd (spexvel)
		    call pargstr (spchar)
		    call pargd (spexerr)
		    call pargstr (spchar)
		    call pargd (spexr)
		call fprintf(logfd[j],"%s%8.2f%s%6.2f%s%2d%s%2d")
		    call pargstr (spchar)
		    call pargd (spvel)
		    call pargstr (spchar)
		    call pargd (sperr)
		    call pargstr (spchar)
		    call pargi (spnl)
		    call pargstr (spchar)
		    call pargi (spnlf)
		do i = 1, nref {
		    if (rmode == 4)
			call fprintf (logfd[j],"%s%8.2f%s%6.2f%s%9.2f%s%5.2f%s%6.2f")
		    else
			call fprintf (logfd[j],"%s%8.2f%s%6.2f%s%9.2f%s%5.2f%s%6.2f%s%6.2f")
			call pargstr (spchar)
			call pargd (velline[i])
			call pargstr (spchar)
			call pargd (errline[i])
			call pargstr (spchar)
			call pargd (htline[i])
			call pargstr (spchar)
			call pargd (widline[i])
			call pargstr (spchar)
			call pargd (eqwline[i])
		    if (rmode == 7)
			call pargd (eqeline[i])
		    }
		call fprintf(logfd[j], "\n")
		}
	    }

# Single line report with line offset, error, height, width, and
# equivalent width for each emission line
	else if (rmode == 8) {
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strlen (instrument) == 0)
		call strcpy ("________", instrument, 8)
	    call aclrd (velline,MAXREF)
	    call aclrd (errline,MAXREF)
	    call aclrd (htline,MAXREF)
	    call aclrd (widline,MAXREF)
	    call aclrd (eqwline,MAXREF)
	    call aclrd (eqeline,MAXREF)

#	Compute individual line velocities from shift
	    if (nfound > 0) {
		do j = 1, nref {
		    velline[j] = 0.0
		    htline[j] = 0.0
		    widline[j] = 0.0
		    errline[j] = 0.0
		    eqwline[j] = 0.0
		    eqeline[j] = 0.0
		    }
		do i = 1, nfound {
		    wtobs[i] = emparams[10,1,i]
		    if (wtobs[i] > 0 && override[i] >= 0) {
			do j = 1, nref {
			    if (wlref[j] == wlrest[i]) {
				velline[j] = wlobs[i] - wlrest[i]
				htline[j] = emparams[5,1,i]
				wl1 = wcs_p2w (emparams[4,1,i] - 1.d0)
				dwdp = wlobs[i] - wl1
				widline[j] = emparams[6,1,i] * dwdp
				errline[j] = emparams[4,2,i] * dwdp
				eqwline[j] = emparams[7,1,i]
				eqeline[j] = emparams[7,2,i]
				}
			    }
			}
		    }
		}

	# Get sky position
	    ra = 99.d0
	    call imgdpar (image,kra,ra)
	    dec = 99.d0
	    call imgdpar (image,kdec,dec)
	    ha = 99.d0
	    call imgdpar (image,kha,ha)
	    lat = 99.d0
	    call imgdpar (image,"SITELAT",lat)

	    if (debug) {
		call printf ("EMRSLTS: RA: %011.2h Dec: %010.1h HA: %011.2h Lat: %010.1h\n")
		    call pargd (ra)
		    call pargd (dec)
		    call pargd (ha)
		    call pargd (lat)
		}

	# Compute altitude and azimuth
	    if (ha < 12.d0 && ha > -12.d0 && dec <= 90.d0 && dec >= -90.d0 &&
 		lat >= -90.d0 && lat <= 90.d0) {
		dpi = 3.1415926535897932d0
		tpi = 6.283185307179588d0
		rcon = 1.74532925199433d-2
		rha = ha * 15.d0 * rcon
		if (rha == 0)
		    rha = 1.d-5
		sha = sin (rha)
		cha = cos (rha)
		rdec = dec * rcon
		sdec = sin (rdec)
		cdec = cos (rdec)
		slat = sin (lat)
		clat = cos (lat)
		ac = (clat * sdec) - (slat * cdec * cha)
		az = atan (-cdec * sha / abs (ac))
		if (ac < 0) az = dpi - az
		if (az < 0) az = tpi + az
		ac = -cdec * sha / sin (az)
		alt = atan ((slat * sdec) + (clat * cdec * cha) / abs (ac))
		if (ac < 0) alt = dpi - alt
		az = az / rcon
		alt = alt / rcon
		}
	    else {
		alt = 0.d0
		az = 0.d0
		}

	    exptime = 0.d0
	    call imgdpar (image,"EXPTIME",exptime)

	# Print one-line report
	    do j = 1, nlogfd {
		if (lspec > 32) {
		    call fprintf(logfd[j],"%s")
			call pargstr (specid)
		    }
		else {
		    call fprintf(logfd[j],"%-32.32s")
			call pargstr (specid)
		    }
		call fprintf(logfd[j],"%s%-8.8s%s%-16.16s")
		    call pargstr (spchar)
		    call pargstr (instrument)
		    call pargstr (spchar)
		    call pargstr (IM_TITLE(image))
		call fprintf(logfd[j],"%s%011.2h%s%010.1h%s%5.1f%s%5.1f")
		    call pargstr (spchar)
		    call pargd (ra)
		    call pargstr (spchar)
		    call pargd (dec)
		    call pargstr (spchar)
		    call pargd (az)
		    call pargstr (spchar)
		    call pargd (alt)
		call fprintf(logfd[j],"%s%13.5f%s%6.1f")
		    call pargstr (spchar)
		    call pargd (dj)
		    call pargstr (spchar)
		    call pargd (exptime)
		call fprintf(logfd[j],"%s%2d%s%2d")
		    call pargstr (spchar)
		    call pargi (spnl)
		    call pargstr (spchar)
		    call pargi (spnlf)
		do i = 1, nref {
		    call fprintf (logfd[j],"%s%8.3f%s%8.4f%s%6.4f%s%9.2f%s%6.3f%s%6.2f%s%6.2f")
			call pargstr (spchar)
			call pargd (wlref[i])
			call pargstr (spchar)
			call pargd (velline[i])
			call pargstr (spchar)
			call pargd (errline[i])
			call pargstr (spchar)
			call pargd (htline[i])
			call pargstr (spchar)
			call pargd (widline[i])
			call pargstr (spchar)
			call pargd (eqwline[i])
			call pargstr (spchar)
			call pargd (eqeline[i])
		    }
		call fprintf(logfd[j], "\n")
		}
	    }

# For echelle with more precision, fewer lines
	else if (rmode == 12) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j], "# %s")
		    call pargstr (specid)
		if (dd != INDEFI && mm != INDEFI && yyyy != INDEFI) {
		    call fprintf(logfd[j], "%4d-%3.3s-%02d %011.2h %.4f")
			call pargi (yyyy)
			call pargstr (month[(mm - 1) * SZ_MONTH + 1])
			call pargi (dd)
			call pargd (ut)
			call pargd (dj)
		    }
		else {
		    call fprintf(logfd[j], nodate)
		    }
		call fprintf(logfd[j], " %8.4f %7.4f %d/%d\n")
			call pargd (spevel)
			call pargd (speerr)
			call pargi (nfit)
			call pargi (nfound)

	# Print emission line information if any were found
		if (nfound > 0) {
		    twt = 0.d0
		    do i = 1, nfound {
			wt = emparams[10,1,i]
			if (wt != 0.d0)
			    twt = twt + 1.d0 / (wt * wt)
			}
		    call fprintf(logfd[j], "#Rest wl  Obs. wl   Line Order Pixel       Vel    Dvel     Eqw    Wt\n")
		    do i = 1, nfound {
			wtobs[i] = emparams[10,1,i]
			if (wtobs[i] > 0) {
			    if (override[i] < 0)
				call strcpy ("-",wch,8)
			    else if (override[i] > 0)
				call strcpy ("+",wch,8)
			    else
				call strcpy (" ",wch,8)
			    }
			else if (linedrop[i] > 0) {
			    call sprintf(wch,8,"X%d")
				call pargi (linedrop[i])
			    }
			else
			    call strcpy ("X",wch,8)
			call fprintf(logfd[j],"%9.4f%s%9.4f%s%4s%s%3d%s%9.4f")
			    call pargd (wlrest[i])
			    call pargstr (spchar)
			    call pargd (wlobs[i])
			    call pargstr (spchar)
			    call pargstr (nmobs[1,i])
			    call pargstr (spchar)
			    call pargi (mspec)
			    call pargstr (spchar)
			    call pargd (emparams[4,1,i])
			wt = emparams[10,1,i]
			if (wt != 0.d0)
			    wt = (1.d0 / (wt * wt)) / twt
			call fprintf (logfd[j],"%s%8.4f%s%6.4f%s%8.2f%s%5.3f%s%s\n")
			    call pargstr (spchar)
			    call pargd ((c0*emparams[9,1,i]) + spechcv)
			    call pargstr (spchar)
			    call pargd (c0*emparams[9,2,i])
			    call pargstr (spchar)
			    call pargd (emparams[7,1,i])
			    call pargstr (spchar)
			    call pargd (wt)
			    call pargstr (spchar)
			    call pargstr (wch)
			}
		    }
		else {
		    call fprintf(logfd[j], nolines)
		    }
		}
	    }
	if (ibr > 0)
	    specid[ibr] = ']'
	if (ibl > 0)
	    specid[ibl] = '['
end

procedure emrshead (rmode0)

int	rmode0		# Report format (1=normal,2=one-line...)

int	i, j
int	rmode
char	spchar[8]
char	uscore[40]
char	tstring[16]

include "emv.com"
include "rvsao.com"
include "results.com"

begin
	if (rmode0 >= 0) {
	    call strcpy (" ", spchar, 8)
	    rmode = rmode0
	    }
	else {
	    call strcpy ("\t", spchar, 8)
	    rmode = -rmode0
	    }
	call strcpy ("--------------------------------", uscore, 40)

	# Print one-line report
	if (rmode == 2 || rmode == 5) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%-8.8s%s%-16.16s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("instrument")
		    call pargstr (spchar)
		    call pargstr ("title")
		call fprintf(logfd[j],"%s%-13.13s%s%-5.5s%s%3.3s")
		    call pargstr (spchar)
		    call pargstr ("Julian Date")
		    call pargstr (spchar)
		    call pargstr ("emvel")
		    call pargstr (spchar)
		    call pargstr ("emerr")
		if (rmode == 3) {
		    call fprintf(logfd[j],"%s%-6.6s%s%3.3s%s%-5.5s")
			call pargstr (spchar)
			call pargstr ("xcvel")
			call pargstr (spchar)
			call pargstr ("xcerr")
			call pargstr (spchar)
			call pargstr ("xcr")
		    }
		else {
		    call fprintf(logfd[j],"%s%-6.6s%s%3.3s%s%-5.5s")
			call pargstr (spchar)
			call pargstr ("xcevel")
			call pargstr (spchar)
			call pargstr ("xcever")
			call pargstr (spchar)
			call pargstr ("xcer")
		    }
		call fprintf(logfd[j],"%s%6.6s%s%3.3s%s%1.1s%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr ("velocity")
		    call pargstr (spchar)
		    call pargstr ("verr")
		    call pargstr (spchar)
		    call pargstr ("E")
		    call pargstr (spchar)
		    call pargstr ("nl")
		    call pargstr (spchar)
		    call pargstr ("nf")
		do i = 1, nref {
		    call sprintf (tstring,16,"line%02d(rest)")
			call pargi (i)
		    call fprintf (logfd[j],"%s%12s")
			call pargstr (spchar)
			call pargstr (tstring)
		    }
		call fprintf(logfd[j], "\n")

		call fprintf(logfd[j],"%32.32s%s%8.8s%s%16.16s")
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%13.13s%s%6.6s%s%3.3s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%6.6s%s%3.3s%s%5.5s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%6.6s%s%3.3s%s%1.1s%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		do i = 1, nref {
		    call fprintf (logfd[j],"%s%12.12s")
			call pargstr (spchar)
			call pargstr (uscore)
		    }
		call fprintf(logfd[j], "\n")
		}
	    }
	else if (rmode == 3 || rmode == 6) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%-8.8s%s%-16.16s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("instrument")
		    call pargstr (spchar)
		    call pargstr ("title")
		call fprintf(logfd[j],"%s%-13.13s%s%-8.8s%s%6.6s")
		    call pargstr (spchar)
		    call pargstr ("Julian Date")
		    call pargstr (spchar)
		    call pargstr ("emvel")
		    call pargstr (spchar)
		    call pargstr ("emverr")
		    if (rmode == 3) {
			call fprintf(logfd[j],"%s%-8.8s%s%6.6s%s%-5.5s")
			    call pargstr (spchar)
			    call pargstr ("xcvel")
			    call pargstr (spchar)
			    call pargstr ("xcverr")
			    call pargstr (spchar)
			    call pargstr ("xcr")
			}
		    else {
			call fprintf(logfd[j],"%s%-8.8s%s%6.6s%s%-5.5s")
			    call pargstr (spchar)
			    call pargstr ("xcevel")
			    call pargstr (spchar)
			    call pargstr ("xcever")
			    call pargstr (spchar)
			    call pargstr ("xcer")
			}
		call fprintf(logfd[j],"%s%8.8s%s%6.6s%s%1.1s%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr ("velocity")
		    call pargstr (spchar)
		    call pargstr ("velerr")
		    call pargstr (spchar)
		    call pargstr ("E")
		    call pargstr (spchar)
		    call pargstr ("nl")
		    call pargstr (spchar)
		    call pargstr ("nf")
		do i = 1, nref {
		    call sprintf (tstring,16,"lin%02dvel")
			call pargi (i)
		    call fprintf (logfd[j],"%s%8s")
			call pargstr (spchar)
			call pargstr (tstring)
		    }
		call fprintf(logfd[j], "\n")

		call fprintf(logfd[j],"%32.32s%s%8.8s%s%16.16s")
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%13.13s%s%8.8s%s%6.6s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%8.8s%s%6.6s%s%5.5s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%8.8s%s%6.6s%s%1.1s%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		do i = 1, nref {
		    call fprintf (logfd[j],"%s%8.8s")
			call pargstr (spchar)
			call pargstr (uscore)
		    }
		call fprintf(logfd[j], "\n")
		}
	    }
	else if (rmode == 4 || rmode == 8) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%8s%s%-16.16s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("instrume")
		    call pargstr (spchar)
		    call pargstr ("title")
		call fprintf(logfd[j],"%s%-11.11s%s%-10.10s%s%-5.5s%s%-5.5s")
		    call pargstr (spchar)
		    call pargstr ("ra")
		    call pargstr (spchar)
		    call pargstr ("dec")
		    call pargstr (spchar)
		    call pargstr ("az")
		    call pargstr (spchar)
		    call pargstr ("alt")
		call fprintf(logfd[j],"%s%-13.13s%s%6.6s")
		    call pargstr (spchar)
		    call pargstr ("Julian Date")
		    call pargstr (spchar)
		    call pargstr ("exposure")
		if (rmode == 4) {
		    call fprintf(logfd[j],"%s%-8.8s%s%-6.6s%s%-8.8s%s%-6.6s%s%-5.5s")
			call pargstr (spchar)
			call pargstr ("emvel")
			call pargstr (spchar)
			call pargstr ("emerr")
			call pargstr (spchar)
			call pargstr ("xcvel")
			call pargstr (spchar)
			call pargstr ("xcerr")
			call pargstr (spchar)
			call pargstr ("xcr")
		    call fprintf(logfd[j],"%s%8.8s%s%6.6s")
			call pargstr (spchar)
			call pargstr ("velocity")
			call pargstr (spchar)
			call pargstr ("velerr")
		    }
		call fprintf(logfd[j],"%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr ("nl")
		    call pargstr (spchar)
		    call pargstr ("nf")
		do i = 1, nref {
		    if (rmode == 8) {
			call fprintf (logfd[j],"%swlrest%02d%swlshft%02d")
			    call pargstr (spchar)
			    call pargi (i)
			    call pargstr (spchar)
			    call pargi (i)
			}
		    else if (rmode == 9) {
			call fprintf (logfd[j],"%swlshft%02d ")
			    call pargstr (spchar)
			    call pargi (i)
			}
		    else if (rmode == 11) {
			call fprintf (logfd[j],"%spxshft%02d")
			    call pargstr (spchar)
			    call pargi (i)
			}
		    else {
			call fprintf (logfd[j],"%slinvel%02d")
			    call pargstr (spchar)
			    call pargi (i)
			}
		    call fprintf (logfd[j],"%slerr%02d%slht%02d    %slwd%02d%sleqw%02d")
			call pargstr (spchar)
			call pargi (i)
			call pargstr (spchar)
			call pargi (i)
			call pargstr (spchar)
			call pargi (i)
			call pargstr (spchar)
			call pargi (i)
		    }
		call fprintf(logfd[j], "\n")
		}
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%32.32s%s%8.8s%s%16.16s")
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%11.11s%s%10.10s%s%5.5s%s%5.5s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		call fprintf(logfd[j],"%s%13.13s%s%6.6s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		if (rmode == 4) {
		    call fprintf(logfd[j],"%s%8.8s%s%6.6s%s%8.8s%s%6.6s%s%5.5s")
			call pargstr (spchar)
			call pargstr (uscore)
			call pargstr (spchar)
			call pargstr (uscore)
			call pargstr (spchar)
			call pargstr (uscore)
			call pargstr (spchar)
			call pargstr (uscore)
			call pargstr (spchar)
			call pargstr (uscore)
		    call fprintf(logfd[j],"%s%8.8s%s%6.6s")
			call pargstr (spchar)
			call pargstr (uscore)
			call pargstr (spchar)
			call pargstr (uscore)
		    }
		call fprintf(logfd[j],"%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr (uscore)
		    call pargstr (spchar)
		    call pargstr (uscore)
		do i = 1, nref {
		    if (rmode == 9 || rmode == 11) {
			call fprintf (logfd[j],"%s%8.8s")
			    call pargstr (spchar)
			    call pargstr (uscore)
			}
		    else {
			call fprintf (logfd[j],"%s%8.8s%s%6.6s%s%9.9s%s%5.5s%s%6.6s")
			    call pargstr (spchar)
			    call pargstr (uscore)
			    call pargstr (spchar)
			    call pargstr (uscore)
			    call pargstr (spchar)
			    call pargstr (uscore)
			    call pargstr (spchar)
			    call pargstr (uscore)
			    call pargstr (spchar)
			    call pargstr (uscore)
			}
		    }
		call fprintf(logfd[j], "\n")
		}
	    }


end

# Aug 12 1992	Correctly print z, not 1+z for fit velocities
# Dec  2 1993	Print spectrum number if multispec file

# Apr  8 1993	Set speed of light locally

# Apr 13 1994	Remove unused variables k0, k1 and rlim
# Apr 22 1994	Remove unused variable gwidth
# Apr 25 1994	Declare month as char
# May 17 1994	Add second option for single-line output format
# Jun 15 1994	Add third output option for per-line velocities
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 10 1994	Fix ut time
# Aug 17 1994	Add weight per line for mode 1 report
# Aug 22 1994	Don't add weight if zero

# Feb  9 1995	Give 2 decimal places for all velocities
# Feb 15 1995	Add heliocentric vel. correction to individual line velocities
# Apr  3 1995	Add fourth output option for single-line per-line velocities
# Apr  5 1995	Add sky coordinates, alt/az and exposure time to 4th option
# Apr 10 1995	Fix alt/az
# May 15 1995	Change all sz_fname to sz_pathname, which is longer
# Jul 13 1995	Get Julian Date from JULDATE instead of JULDAY
# Jul 13 1995	Get UT from JULDATE
# Jul 19 1995	Handle INDEF values for SPVEL and SPXVEL
# Sep 21 1995	Handle INDEF value for SPEVEL

# Jan 15 1997	Add HDJ to JULDATE arguments
# Feb 12 1997	Add modes 5 and 6 to plot emission line template velocity
# May  2 1997	Always test against dindef, not INDEFD
# May  6 1997	Print messages if no position or no observation date
# May 19 1997	Add height and width to mode 4 report
# May 19 1997	Use CORTEMP parameter instead of inline template name
# May 21 1997	Reinitialize all line parameter arrays for each spectrum
# May 22 1997	If CORTEMP is null string, use best template
# Sep 25 1997	Add null print line so this works in Digital Unix
# Sep 25 1997	In mode 1, put file name on separate line
# Oct  1 1997	Print rejection code in mode 1 if line is dropped from fit
# Oct  6 1997	Fix bug in mode 1 BCV flag code
# Dec 17 1997   Change order of date from dd-mmm-yyyy to yyyy-mmm-dd
# Dec 17 1997   Use EQUINOX if it is present instead of EPOCH

# May 21 1999	Read DATE if DATE-OBS is no found
# Jul 22 1999	Put _ in null strings
# Jul 23 1999	Add report mode 5, same as 4, but with tabs
# Jul 23 1999	Add emrshead() to print report heading

# Jun  8 2000	Add equivalent width error to new report modes 0 and 7
# Oct  3 2000	Add Mode 8 to print wavelength shifts

# Jan 30 2001	Fix bug in Mode 8 printed wavelength 
# Aug  7 2001	Add Mode 10 to print only line information

# May 11 2005	Print 3 decimal places in output mode 9 shifts
# May 19 2005	Add mode 11 to print line shifts in pixels
# Jul 14 2005	Initialize ibr and ibl
# Mar 18 2008	Add mode 12 for high precision echelle velocities
# May 30 2008	Use strings more efficiently for IRAF Solaris limits

# May 20 2009	Return immediately if not logging results
