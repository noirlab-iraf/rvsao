# File rvsao/Xcsao/xcrslts.x
# May 25, 2012
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1990-2012 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Produce a statistical summary of all the cross-correlations of templates
#  and the object spectrum.  Most input to summary is in common in "rvsao.com"

include "rvsao.h"
include <imhdr.h>
include "contin.h"
define  SZ_MONTH	3

procedure xcrslts (specfile, mspec, image, newpts, rmode0, filexcor, txfile)

char	specfile[ARB]		# spectrum file name
int	mspec			# spectrum number within file
pointer	image			# spectrum data structure
int	newpts			# number of points in cross-correlation
int	rmode0			# report format (0=archive,1=normal,2=one-line)
bool	filexcor		# true if correlation vector file is written
char	txfile[SZ_PATHNAME,ARB]	# names of correlation vector files

char	filename[SZ_PATHNAME]
char	instrument[8]
char	incode
char	bcvform[16], hcvform[16]
char	jdot, lbr, rbr, colon
int	mm, dd, yyyy
int	i, j, k, nc, imax, itemp
int	strldx()
int	icz, icze, iczxc, iczxce, iczem, iczeme
double	dj		# Julian date of observation
double	djh		# Heliocentric Julian date of observation
double	width
double	dindef
int	strncmp(),strcmp(), strlen()
int	inova,jnova,qnova,ltemp,ic
char	pname[20]	# image header parameter name
int	rmode			# report format (0=archive,1=normal,2=one-line)
 
double	novavel[MAXTEMPS]	# Correlation velocity of each template
double	novaerr[MAXTEMPS]	# Correlation velocity error of each template
double	novar[MAXTEMPS]		# Correlation R-value of each template
real	novahght[MAXTEMPS]	# Height of correlation peak for this template
real	novawdth[MAXTEMPS]	# Width of correlation peak for this template
real	novarms[MAXTEMPS]	# Asymetric RMS for this template
int	nnova			# Total number of templates
char	novaname[SZ_PATHNAME,MAXTEMPS]	# Names of the templates
char	tempfile[SZ_PATHNAME]	# Template filename
char	str[SZ_LINE]
char	tempinfo[32]
char	spchar[8]
int	ibl, ibr
int	fiber, aperture, beam
char	nameform[16]

double	ra, dec, ut, vel, fracpeak, z1, epoch
string	month	"JanFebMarAprMayJunJulAugSepOctNovDec"
int	imaccf(), stridx()

include "rvsao.com"
include "results.com"
include "oldres.com"

begin
	dindef = INDEFD
	jdot = '.'
	colon = ':'
	lbr = char (91)
	rbr = char (93)
	mm = INDEFI
	dd = INDEFI
	yyyy = INDEFI
	if (imaccf (image, "DATE-OBS") == YES)
	    call imgdate (image, "DATE-OBS", mm, dd, yyyy)
	else if (imaccf (image, "DATE") == YES)
	    call imgdate (image, "DATE", mm, dd, yyyy)
	djh = 0.d0
	dj = 0.d0
	call juldate (image,ut,dj,djh,debug)
	call strcpy ("HCV: %6.2f", hcvform, 16)
	call strcpy ("BCV: %6.2f", bcvform, 16)

# Set separating character to tab if report mode is less than 0
	if (rmode0 >= 0) {
	    call strcpy (" ", spchar, 8)
	    rmode = rmode0
	    }
	else {
	    call strcpy ("\t", spchar, 8)
	    rmode = -rmode0
	    }

# Mode 1: Print everything
	if (rmode == 1) {
	    if (pkfrac < 0.d0)
		fracpeak = -pkfrac
	    else
		fracpeak = pkfrac

#	Set up heading in logfiles
	    ra = dindef
	    dec = dindef
	    call imgdpar (image,"RA",ra)
	    if (ra != dindef) {
		call imgspar (image, "RA", str, SZ_LINE)
		if (stridx (colon, str) == 0)
		    ra = ra / 15.0
		}
	    call imgdpar (image,"DEC",dec)
	    epoch = 1950.0
	    call imgdpar (image,"EPOCH",epoch)
	    call imgdpar (image,"EQUINOX",epoch)
	    do i = 1, nlogfd {
		call fprintf(logfd[i], "\n%s Object: %s \n")
		    call pargstr (specid)
		    call pargstr (specname)
		if (ra != dindef && dec != dindef) {
		    call fprintf(logfd[i], "RA: %011.2h Dec: %010.1h %.1f\n")
			call pargd (ra)
			call pargd (dec)
			call pargd (epoch)
		    }
		else
		    call fprintf(logfd[i], "(No position)\n")
		if (dd == INDEFI || mm == INDEFI || yyyy == INDEFI)
		    call fprintf(logfd[i], "(No observation date) ")
		else if (djh > 0.d0)  {
		    call fprintf(logfd[i], "%4d-%3.3s-%02d %011.2h =HJD%.4f ")
			call pargi (yyyy)
			call pargstr (month[(mm - 1) * SZ_MONTH + 1])
			call pargi (dd)
			call pargd (ut)
			call pargd (djh)
		    }
		else {
		    call fprintf(logfd[i], "%4d-%3.3s-%02d %011.2h =JD%.4f ")
			call pargi (yyyy)
			call pargstr (month[(mm - 1) * SZ_MONTH + 1])
			call pargi (dd)
			call pargd (ut)
			call pargd (dj)
		    }
		switch (svcor) {
		    case NONE:
			call pargstr ("no BCV")
		    case HCV:
			call fprintf(logfd[i], hcvform)
		    case BCV:
			call fprintf(logfd[i], bcvform)
		    case FHCV:
			if (specvb) 
			    call fprintf(logfd[i], bcvform)
			else
			    call fprintf(logfd[i], hcvform)
		    case FBCV:
			if (specvb) 
			    call fprintf(logfd[i], bcvform)
			else
			    call fprintf(logfd[i], hcvform)
		    default:
		    }
		    call pargd (spechcv)
		call fprintf(logfd[i], "\n")

		call fprintf(logfd[i],"%.1fA- %.1fA ")
		    call pargd (twl1[itmax])
		    call pargd (twl2[itmax])
		if (zpad) {
		    call fprintf(logfd[i],"%dx2 points, ")
			call pargi (npts)
		    }
		else {
		    call fprintf(logfd[i],"%d points, ")
			call pargi (newpts)
		    }
		call fprintf(logfd[i],"filter: %d %d %d %d apodize %4.2f fit %3.2f best %d\n")
		    call pargi (lo)
		    call pargi (toplo)
		    call pargi (topnrn)
		    call pargi (nrun)
		    call pargr (han)
		    call pargd (fracpeak)
		    call pargi (itmax)
		if (tshift != 0.0) {
		    call fprintf(logfd[i],"Tshift: %6.2f \n")
		    call pargd (tshift)
		    }
		}

# Record results in log files (one line per template per peak fit mode)

	    do itemp = 1, ntemp {
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"Temp: %16s vel: %6.2f ")
			call pargstr (tempid[1,itemp])
			call pargd (tempvel[itemp])
		    call fprintf(logfd[j]," tsh: %6.2f HCV: %6.2f Peak: %7.3f")
			call pargd (tempshift[itemp])
			call pargd (temphcv[itemp])
			call pargd (cz[itemp])
		    if (correlate == COR_PIX || correlate == COR_WAV)
			call fprintf(logfd[j]," h:%5.3f R: %6.2f Shift: %7.3f +/- %7.3f %d")
		    else
			call fprintf(logfd[j]," h:%5.3f R: %6.2f CZ: %7.3f +/- %7.3f %d")
			call pargd (thght[itemp])
			call pargd (czr[itemp])
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
			call pargi (pkmode0)
		    if (tschop[itemp])
			call fprintf(logfd[j]," -el")
		    if (tachop[itemp])
			call fprintf(logfd[j]," -al")
		    if (tempfilt[itemp] == 1)
			call fprintf(logfd[j]," -tf")
		    else if (tempfilt[itemp] == 2)
			call fprintf(logfd[j]," +hi")
		    else if (tempfilt[itemp] == 3)
			call fprintf(logfd[j]," +hi -tf")
		    if (!tscont[itemp])
			call fprintf(logfd[j]," +tc")
		    else if (tconproc[itemp] == ZEROCONT)
			call fprintf(logfd[j]," divc")
		    call fprintf(logfd[j],"\n")
		    }
		}
	    }

# Mode 2:  One line per template
	else if (rmode == 2) {
	    do itemp = 1, ntemp {
		if (taa[itemp] != 0) {
		    dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		    width = c0 * ((10.d0 ** (twdth[itemp] * dlogw)) - 1.d0)
		    }
		else
		    width = tvw[itemp]
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%-16.16s %-16.16s %6.2f %9.3f %8.3f %5.3f %7.3f\n")
			call pargstr (specid)
			call pargstr (tempid[1,itemp])
			call pargd (czr[itemp])
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
			call pargd (thght[itemp])
			call pargd (width)
		    }
		}
	    }

# Mode 3:  One line per object spectrum with past results and best template
	else if (rmode == 3) {
	    nc = 16
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)
	    imax = 1
	    i = imax
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-16.16s %6.2f %6.2f %9.3f %9.3f %8.3f %8.3f %6.3f %6.3f %12.2f %c %s\n")
		call pargstr (filename)
		call pargd (spr0)
		call pargd (czr[itmax])
		call pargd (spvel0)
		call pargd (zvel[itmax])
		call pargd (sperr0)
		call pargd (czerr[itmax])
		call pargd (spechcv0)
		call pargd (spechcv)
		call pargd (djh)
		call pargc (spqual0)
		call pargstr (tempid[1,itmax])
		}
	    }

	else if (rmode == 4) {
	    nc = strldx (jdot,specfile) - 1
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)
	    do itemp = 1, ntemp {
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%-16.16s %s %6.2f %9.3f %8.3f\n")
			call pargstr (filename)
			call pargstr (tempid[1,itemp])
			call pargd (czr[itemp])
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
		    }
		}
	    }

	else if (rmode == 5) {
	    nc = strldx (jdot,specid) - 1
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)
	    do itemp = 1, ntemp {
		if (taa[itemp] != 0) {
		    dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		    width = c0 * ((10.d0 ** (twdth[itemp] * dlogw)) - 1.d0)
		    }
		else
		    width = tvw[itemp]
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%2d %2d %3d %3d %-16.16s %5.2f %-14.14s")
			call pargi (lo)
			call pargi (toplo)
			call pargi (topnrn)
			call pargi (nrun)
			call pargstr (tempid[1,itemp])
			call pargd (tempshift[itemp])
			call pargstr (specid)
		    call fprintf(logfd[j]," %6.2f %4.1f %4.2f %5.1f %5.2f\n")
			call pargd (zvel[itemp])
			call pargd (czr[itemp])
			call pargd (thght[itemp])
			call pargd (width)
			call pargd (spechcv)
		    }
		}
	    }

	else if (rmode == 6) {
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strncmp (instrument,"echelle",7) == 0)
		incode = 'T'
	    else if (strncmp (instrument,"mmtech",6) == 0)
		incode = 'M'
	    else if (strncmp (instrument,"oroech",6) == 0)
		incode = 'W'
	    else
		incode = '_'
	    nc = strldx (jdot,specid) - 1
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)
	    do itemp = 1, ntemp {
		if (taa[itemp] != 0) {
		    dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		    width = c0 * ((10.d0 ** (twdth[itemp] * dlogw)) - 1.d0)
		    }
		else
		    width = tvw[itemp]
		if (width < 0.d0)
		    width = 0.d0
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%c %-14.14s %-16.16s %14.5f %7.2f %6.2f %5.3g %4.2f %5.1f %6.2f")
			call pargc (incode)
			call pargstr (specid)
			call pargstr (tempid[1,itemp])
			call pargd (djh)
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
			call pargd (czr[itemp])
			call pargd (thght[itemp])
			call pargd (width)
			call pargd (spechcv)
		    if (filexcor) {
			call fprintf(logfd[j]," %7.2f %7.2f %s")
			    call pargd (tsig1[itemp])
			    call pargd (tsig2[itemp])
			    call pargstr (txfile[1,itemp])
			}
		    call fprintf(logfd[j],"\n")
		    }
		}
	    }

	else if (rmode == 7) {
#	    nc = strldx (jdot,specid) - 1
	    nc = 16
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)

	# Read previous cross-correlation information from image file
	    call imgipar (image,"NTEMP",nnova)
	    call imgipar (image,"QCSTAT",qnova)
	    do inova= 1, nnova {
		call sprintf (pname,6,"TEMPL%d")
		    call pargi (inova)
		call imgspar (image,pname,novaname[1,inova],SZ_PATHNAME)
		ltemp = strlen (novaname[1,inova])
		if (ltemp > 0) {
		    do itemp = 1, ltemp {
			ic = novaname[itemp,inova]
			if (ic == 32) ic = EOS
			if (ic > 64 & ic < 96) ic = ic + 32
			novaname[itemp,inova] = ic
			}
		    }
		else if (inova == 1)
		call imgspar (image,pname,novaname[1,inova],SZ_PATHNAME)
		    call strcpy ("ztemp",novaname[1,inova],SZ_PATHNAME)
		call sprintf (pname,8,"TCZXC%d")
		    call pargi (inova)
		call imgdpar (image,pname,novavel[inova])
		call sprintf (pname,8,"TCZXCR%d")
		    call pargi (inova)
		call imgdpar (image,pname,novar[inova])
		call sprintf (pname,8,"TCZXCER%d")
		    call pargi (inova)
		call imgdpar (image,pname,novaerr[inova])
		call sprintf (pname,8,"TXCHGHT%d")
		    call pargi (inova)
		call imgrpar (image,pname,novahght[inova])
		call sprintf (pname,8,"TXCWDTH%d")
		    call pargi (inova)
		call imgrpar (image,pname,novawdth[inova])
		call sprintf (pname,8,"TXCARMS%d")
		    call pargi (inova)
		call imgrpar (image,pname,novarms[inova])
		call sprintf (pname,8,"QCSTATS%d")
		    call pargi (inova)
		}
	    call strcpy ("none",novaname[1,nnova+1],SZ_PATHNAME)
	    novavel[nnova+1] = 0.d0
	    novar[nnova+1] = 0.d0
	    novaerr[nnova+1] = 0.d0
	    novahght[nnova+1] = 0.d0
	    novarms[nnova+1] = 0.d0

	    do itemp = 1, ntemp {
		jnova = nnova+1
		do inova = 1, nnova {
		    if (strcmp(novaname[1,inova],tempid[1,itemp]) == 0)
			jnova = inova
		    }
		do j = 1, nlogfd {
		call fprintf(logfd[j],"%16.16s %6.2f %6.2f %9.3f %9.3f %8.3f %8.3f")
		    call pargstr (filename)
		    call pargd (novar[jnova])
		    call pargd (czr[itemp])
		    call pargd (novavel[jnova])
		    call pargd (zvel[itemp])
		    call pargd (novaerr[jnova])
		    call pargd (czerr[itemp])
		call fprintf(logfd[j]," %7.5f %7.5f %7.5f %7.5f %12.2f %s %s %s\n")
		    call pargr (novahght[jnova])
		    call pargd (thght[itemp])
		    call pargr (novarms[jnova])
		    call pargd (tarms[itemp])
		    call pargd (dj)
		    call pargstr (tempid[1,itemp])
		    call pargstr (novaname[1,jnova])
		    if (qnova == 1)
			call pargc ("?")
		    else if (qnova == 2)
			call pargc ("?")
		    else if (qnova == 3)
			call pargc ("X")
		    else if (qnova == 4)
			call pargc ("Q")
		    else
			call pargc ("_")
		}
		}
	    }

# Single line report, including emission line information
	else if (rmode == 8) {
	    spxvel = zvel[itmax]
	    spxerr = czerr[itmax]
	    spxr = czr[itmax]
	    call vcombine (spxvel,spxerr,spxr,spevel,speerr,spnlf,spvel,sperr,debug)
	    if (spvel != dindef) {
		icz = nint (spvel)
		icze = nint (sperr)
		}
	    else {
		icz = 0
		icze = 0
		}
	    if (spxvel != dindef) {
		iczxc = nint (spxvel)
		iczxce = nint (spxerr)
		}
	    else {
		iczxc = 0
		iczxce = 0
		}
	    if (spevel != dindef) {
		iczem = nint (spevel)
		iczeme = nint (speerr)
		}
	    else {
		iczem = 0
		iczeme = 0
		}
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strlen (instrument) == 0)
		call strcpy ("________", instrument, 8)
	    if (strlen (IM_TITLE(image)) == 0)
		call strcpy ("________________", IM_TITLE(image), SZ_IMTITLE)
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%6s%s%-8s%s%-16s%s%15.5f")
		    call pargstr (specid)
		    call pargstr (spchar)
		    call pargstr (instrument)
		    call pargstr (spchar)
		    call pargstr (IM_TITLE(image))
		    call pargstr (spchar)
		    call pargd (dj)
		call fprintf(logfd[j],"%s%5d%s%3d%s%5d%s%3d%s%5.1f")
		    call pargstr (spchar)
		    call pargi (iczem)
		    call pargstr (spchar)
		    call pargi (iczeme)
		    call pargstr (spchar)
		    call pargi (iczxc)   
		    call pargstr (spchar)
		    call pargi (iczxce)
		    call pargstr (spchar)
		    call pargd (spxr)
		call fprintf(logfd[j],"%s%5d%s%3d%sC%s%2d%s%2d%s%s\n")
		    call pargstr (spchar)
		    call pargi (icz)
		    call pargstr (spchar)
		    call pargi (icze)
		    call pargstr (spchar)
		    call pargstr (spchar)
		    call pargi (spnl)
		    call pargstr (spchar)
		    call pargi (spnlf)
		    call pargstr (spchar)
		    call pargstr (tempid[1,itmax])
		}
	    }

# Print one line per template with correlation output file name, if used
	else if (rmode == 9) {
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strncmp (instrument,"echelle",7) == 0)
		incode = 'T'
	    else if (strncmp (instrument,"mmtech",6) == 0)
		incode = 'M'
	    else if (strncmp (instrument,"oroech",6) == 0)
		incode = 'W'
	    else
		incode = '_'
	    nc = strldx (jdot,specid) - 1
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)
	    do itemp = 1, ntemp {
		if (taa[itemp] != 0) {
		    dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		    width = c0 * ((10.d0 ** (twdth[itemp] * dlogw)) - 1.d0)
		    }
		else
		    width = tvw[itemp]
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%c %-14.14s %-16.16s %14.5f %8.3f")
			call pargc (incode)
			call pargstr (specid)
			call pargstr (tempid[1,itemp])
			call pargd (djh)
			call pargd (zvel[itemp])
		    call fprintf(logfd[j]," %7.3f %4.1f %4.2f %5.1f %6.2f")
			call pargd (czerr[itemp])
			call pargd (czr[itemp])
			call pargd (thght[itemp])
			call pargd (width)
			call pargd (spechcv)
		    if (filexcor) {
			call fprintf(logfd[j]," %7.2f %7.2f %s")
			    call pargd (tsig1[itemp])
			    call pargd (tsig2[itemp])
			    call pargstr (txfile[1,itemp])
			}
		    call fprintf(logfd[j],"\n")
		    }
		}
	    }

# Print all templates on one line
	else if (rmode == 10) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-14.14s%s%14.5f%s%d")
		    call pargstr (specid)
		    call pargstr (spchar)
		    call pargd (dj)
		    call pargstr (spchar)
		    call pargi (itmax)
		if (ntemp > 1) {
		    call fprintf(logfd[j],"%s%-12.12s%s%7.1f%s%5.1f%s%4.2f")
			call pargstr (spchar)
			call pargstr (tempid[1,itmax])
			call pargstr (spchar)
			call pargd (zvel[itmax])
			call pargstr (spchar)
			call pargd (czerr[itmax])
			call pargstr (spchar)
			call pargd (czr[itmax])
		    }
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%-12.12s%s%7.1f%s%5.1f%s%4.2f")
			call pargstr (spchar)
			call pargstr (tempid[1,itemp])
			call pargstr (spchar)
			call pargd (zvel[itemp])
			call pargstr (spchar)
			call pargd (czerr[itemp])
			call pargstr (spchar)
			call pargd (czr[itemp])
		    }
		call fprintf(logfd[j],"\n")
		}
	    }

# Print all templates on one line with long file names
	else if (rmode == 11) {
	    if (maxlength > 0) {
		call sprintf (nameform, 16, "%%-%d.%ds")
		    call pargi (maxlength)
		    call pargi (maxlength)
		}
	    else {
		call strcpy ("%-32.32s", nameform, 16)
		}
	    do j = 1, nlogfd {
		call fprintf(logfd[j],nameform)
		    call pargstr (specid)
		call fprintf(logfd[j],"%s%14.5f%s%d")
		    call pargstr (spchar)
		    call pargd (dj)
		    call pargstr (spchar)
		    if (czr[itmax] < 1.d0)
			call pargi (0)
		    else
			call pargi (itmax)
		if (ntemp > 1) {
		    if (czr[itmax] < 1.d0) {
		    call fprintf(logfd[j],"%s%-14.14s%s%7.1f%s%5.1f%s%4.2f")
			call pargstr (spchar)
			call pargstr ("none")
			call pargstr (spchar)
			call pargd (0.d0)
			call pargstr (spchar)
			call pargd (0.d0)
			call pargstr (spchar)
			call pargd (czr[itmax])
			}
		    else {
		    call fprintf(logfd[j],"%s%-14.14s%s%7.1f%s%5.1f%s%4.2f")
			call pargstr (spchar)
			call pargstr (tempid[1,itmax])
			call pargstr (spchar)
			call pargd (zvel[itmax])
			call pargstr (spchar)
			call pargd (czerr[itmax])
			call pargstr (spchar)
			call pargd (czr[itmax])
			}
		    }
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%-14.14s%s%7.7g%s%5.5g%s%4.2f")
			call pargstr (spchar)
			call pargstr (tempid[1,itemp])
			call pargstr (spchar)
			call pargd (zvel[itemp])
			call pargstr (spchar)
			call pargd (czerr[itemp])
			call pargstr (spchar)
			call pargd (czr[itemp])
		    }
		call fprintf(logfd[j],"\n")
		}
	    }

# Print all templates on one line with long file names - no best velocity
	else if (rmode == 12) {
	    ibr = stridx (rbr, specid)
	    if (ibr > 0)
		specid[ibr] = ' '
	    ibl = stridx (lbr, specid)
	    if (ibl > 0)
		specid[ibl] = ' '
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%14.5f%s%d")
		    call pargstr (specid)
		    call pargstr (spchar)
		    call pargd (dj)
		    call pargstr (spchar)
		    call pargi (itmax)
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%-14.14s%s%7.7g%s%5.5g%s%4.2f")
			call pargstr (spchar)
			call pargstr (tempid[1,itemp])
			call pargstr (spchar)
			call pargd (zvel[itemp])
			call pargstr (spchar)
			call pargd (czerr[itemp])
			call pargstr (spchar)
			call pargd (czr[itemp])
		    }
		call fprintf(logfd[j],"\n")
		}
	    if (ibl > 0)
		specid[ibl] = '['
	    if (ibr > 0)
		specid[ibl] = ']'
	    }

# Print measured rather than fit center and height
	else if (rmode == 13) {
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strncmp (instrument,"echelle",7) == 0)
		incode = 'T'
	    else if (strncmp (instrument,"mmtech",6) == 0)
		incode = 'M'
	    else if (strncmp (instrument,"oroech",6) == 0)
		incode = 'W'
	    else
		incode = '_'
	    nc = strldx (jdot,specid) - 1
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)
	    do itemp = 1, ntemp {
		if (taa[itemp] != 0) {
		    dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		    width = c0 * ((10.d0 ** (tpwdth[itemp] * dlogw)) - 1.d0)
		    z1 = (10.d0 ** (tpcent[itemp]*dlogw)) * (1.d0+(spechcv/c0))
		    vel = c0 * (z1 - 1.d0)
		    }
		else {
		    width = tvw[itemp]
		    vel =  cz[itemp]
		    }
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%c %-14.14s %-16.16s %14.5f %7.2f %5.3f %6.2f %6.2f")
			call pargc (incode)
			call pargstr (specid)
			call pargstr (tempid[1,itemp])
			call pargd (djh)
			call pargd (vel)
			call pargd (tphght[itemp])
			call pargd (width)
			call pargd (spechcv)
		    if (filexcor) {
			call fprintf(logfd[j]," %7.2f %7.2f %s")
			    call pargd (tsig1[itemp])
			    call pargd (tsig2[itemp])
			    call pargstr (txfile[1,itemp])
			}
		    call fprintf(logfd[j],"\n")
		    }
		}
	    }

# Print emission line velocity plus multiple templates
	else if (rmode == 14) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%14.5f")
		    call pargstr (specid)
		    call pargstr (spchar)
		    call pargd (dj)
		call fprintf(logfd[j],"%s%7.1f%s%5.1f%s%d%s%d")
		    call pargstr (spchar)
		    call pargd (spevel)
		    call pargstr (spchar)
		    call pargd (speerr)
		    call pargstr (spchar)
		    call pargi (spnl)
		    call pargstr (spchar)
		    call pargi (spnlf)
		call fprintf(logfd[j],"%s%d%s%-14.14s")
		    call pargstr (spchar)
		    call pargi (itmax)
		    call pargstr (spchar)
		    call pargstr (tempid[1,itmax])
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%-14.14s%s%7.1f%s%5.1f%s%4.2f")
			call pargstr (spchar)
			call pargstr (tempid[1,itemp])
			call pargstr (spchar)
			call pargd (zvel[itemp])
			call pargstr (spchar)
			call pargd (czerr[itemp])
			call pargstr (spchar)
			call pargd (czr[itemp])
		    }
		call fprintf(logfd[j],"\n")
		}
	    }
# Like mode 6 but with longer template and more digits
	else if (rmode == 15) {
	    call imgspar (image,"INSTRUME",instrument,8)
	    if (strncmp (instrument,"echelle",7) == 0)
		incode = 'T'
	    else if (strncmp (instrument,"mmtech",6) == 0)
		incode = 'M'
	    else if (strncmp (instrument,"oroech",6) == 0)
		incode = 'W'
	    else
		incode = '_'
	    nc = strldx (jdot,specid) - 1
	    call strcpy (specid,filename,nc)
#	    call printf("%16s -> %16s\n")
#		call pargstr (specid)
#		call pargstr (filename)
	    do itemp = 1, ntemp {
		if (taa[itemp] != 0) {
		    dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		    width = c0 * ((10.d0 ** (twdth[itemp] * dlogw)) - 1.d0)
		    }
		else
		    width = tvw[itemp]
		if (width < 0.d0)
		    width = 0.d0
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%c %-14.14s %-18.18s %14.5f %9.3f %7.3f %5.3g %5.3f %7.3f %7.3f")
			call pargc (incode)
			call pargstr (specid)
			call pargstr (tempid[1,itemp])
			call pargd (djh)
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
			call pargd (czr[itemp])
			call pargd (thght[itemp])
			call pargd (width)
			call pargd (spechcv)
		    if (filexcor) {
			call fprintf(logfd[j]," %7.2f %7.2f %s")
			    call pargd (tsig1[itemp])
			    call pargd (tsig2[itemp])
			    call pargstr (txfile[1,itemp])
			}
		    call fprintf(logfd[j],"\n")
		    }
		}
	    }

# Mode 16:  One line per template
	else if (rmode == 16) {
	    do itemp = 1, ntemp {
		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%-24.24s %-16.16s %6.2f %9.3f %8.3f %5.3f")
			call pargstr (specid)
			call pargstr (tempid[1,itemp])
			call pargd (czr[itemp])
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
			call pargd (thght[itemp])
		    call fprintf (logfd[j]," %7.1f %7.1f %7.1f\n")
			call pargd (tempwl1[itemp])
			call pargd (tempwl2[itemp])
			call pargd ((tempwl1[itemp]+tempwl2[itemp])*0.5d0)
		    }
		}
	    }

# Mode 17:  One line per template for Hectochelle
	else if (rmode == 17) {
	    aperture = 1
	    if (mspec > 0)
		aperture = mspec
	    else
		call imgipar (image,"APERTURE",aperture)
	    fiber = 0
	    call imgipar (image,"FIBER",fiber)
	    beam = 1
	    call imgipar (image,"BEAM",beam)
	    do itemp = 1, ntemp {
		if (taa[itemp] != 0) {
		    dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		    width = c0 * ((10.d0 ** (twdth[itemp] * dlogw)) - 1.d0)
		    }
		else
		    width = tvw[itemp]
		if (width < 0.d0)
		    width = 0.d0

	    # Keep only last 35 characters of template file name - .fits
		tempinfo[28] = EOS
		ltemp = strlen (tempid[1,itemp])
		i = ltemp - 35
		if (i < 0) {
		    i = 0
		    nc = ltemp
		    }
		else {
		    nc = 35
		    }
		k = strldx ('/',tempid[1,itemp])
		if (k > 0) {
		    i = k
		    nc = nc - (k - i)
		    }
		do j = 1, nc {
		    i = i + 1
		    tempinfo[j] = tempid[i,itemp]
		    }

		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%03d %03d %d %-24.24s %-35.35s %14.5f %9.3f %8.3f %6.3f %6.4f %6.2f %7.3f\n")
			call pargi (aperture)
			call pargi (fiber)
			call pargi (beam)
			call pargstr (specname)
			call pargstr (tempinfo)
			call pargd (djh)
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
			call pargd (czr[itemp])
			call pargd (thght[itemp])
			call pargd (width)
			call pargd (spechcv)
		    }
		}
	    }

# Mode 18:  One line per template for TRES
	else if (rmode == 18) {
	    aperture = 1
	    if (mspec > 0)
		aperture = mspec
	    else
		call imgipar (image,"APERTURE",aperture)
	    fiber = 0
	    call imgipar (image,"FIBER",fiber)
	    beam = 1
	    call imgipar (image,"BEAM",beam)
	    do itemp = 1, ntemp {
		ibl = stridx (lbr, tempid[1,itemp])
		if (ibl > 0)
		    call strcpy (tempid[1,itemp], tempfile, ibl-1)
		else
		    call strcpy (tempid[1,itemp], tempfile, SZ_PATHNAME)

		do j = 1, nlogfd {
		    call fprintf(logfd[j],"%02d %s %s %14.5f %9.5f %8.5f %6.3f %7.3f\n")
			call pargi (aperture)
			call pargstr (specfile)
			call pargstr (tempfile)
			call pargd (djh)
			call pargd (zvel[itemp])
			call pargd (czerr[itemp])
			call pargd (czr[itemp])
			call pargd (spechcv)
		    }
		}
	    }
end

procedure xcrshead (rmode0)

int	rmode0		# Report format (1=normal,2=one-line...)

int	itemp, j
int	rmode
char	spchar[8]
char	dashes[40]
char	tstring[16]

include "rvsao.com"

begin
	if (rmode0 >= 0) {
	    call strcpy (" ", spchar, 8)
	    rmode = rmode0
	    }
	else {
	    call strcpy ("\t", spchar, 8)
	    rmode = -rmode0
	    }
	call strcpy ("--------------------------------", dashes, 40)

	# Print one-line report
	if (rmode == 8) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%-8.8s%s%-16.16s%s%-15.15s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("instrument")
		    call pargstr (spchar)
		    call pargstr ("title")
		    call pargstr (spchar)
		    call pargstr ("julian_date")
		call fprintf(logfd[j],"%s%-5.5s%s%5.5s%s%-5.5s%s%5.5s%s%-5.5s")
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
		call fprintf(logfd[j],"%s%5.5s%s%3.3s%sC%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr ("vel")
		    call pargstr (spchar)
		    call pargstr ("err")
		    call pargstr (spchar)
		    call pargstr (spchar)
		    call pargstr ("nl")
		    call pargstr (spchar)
		    call pargstr ("nf")
		call fprintf(logfd[j],"%s%-8.8s\n")
		    call pargstr (spchar)
		    call pargstr ("best")

		call fprintf(logfd[j],"%32.32s%s%8.8s%s%16.16s%s%15.15s")
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		call fprintf(logfd[j],"%s%5.5s%s%5.5s%s%5.5s%s%5.5s%s%5.5s")
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		call fprintf(logfd[j],"%s%5.5s%s%3.3s%s%1.1s%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		call fprintf(logfd[j],"%s%8.8s\n")
		    call pargstr (spchar)
		    call pargstr (dashes)
		}
	    }
	else if (rmode == 10) {
	
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-14.14s%s%-14.14s%s%-4.4s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("julian_date")
		    call pargstr (spchar)
		    call pargstr ("best")
		if (ntemp > 1) {
		    call fprintf(logfd[j],"%s%-12.12s%s%-7.7s%s%-7.7s%s%-5.5s")
			call pargstr (spchar)
			call pargstr ("besttemp")
			call pargstr (spchar)
			call pargstr ("bestvel")
			call pargstr (spchar)
			call pargstr ("besterr")
			call pargstr (spchar)
			call pargstr ("bestr")
		    }
		do itemp = 1, ntemp {
		    call sprintf (tstring,16,"temp%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-12.12s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"vel%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-7.7s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"err%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-5.5s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"xcr%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-4.4s")
			call pargstr (spchar)
			call pargstr (tstring)
		    }
		call fprintf(logfd[j],"\n")

		call fprintf(logfd[j],"%14.14s%s%14.14s%s%2.2s")
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		if (ntemp > 1) {
		    call fprintf(logfd[j],"%s%12.12s%s%7.7s%s%5.5s%s%4.4s")
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
		    }
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%12.12s%s%7.7s%s%5.5s%s%4.4s")
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
		    }
		call fprintf(logfd[j],"\n")
		}
	    }

# Print all templates on one line with long file names
	else if (rmode == 11) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%-14.14s%s%s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("julian_date")
		    call pargstr (spchar)
		    call pargstr ("best")
		if (ntemp > 1) {
		    call fprintf(logfd[j],"%s%-14.14s%s%-7.7s%s%-5.5s%s%4.4s")
			call pargstr (spchar)
			call pargstr ("besttemp")
			call pargstr (spchar)
			call pargstr ("xcvel")
			call pargstr (spchar)
			call pargstr ("xcerr")
			call pargstr (spchar)
			call pargstr ("xcr")
		    }
		do itemp = 1, ntemp {
		    call sprintf (tstring,16,"temp%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-14.14s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"vel%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-7.7s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"err%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-5.5s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"xcr%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-4.4s")
			call pargstr (spchar)
			call pargstr (tstring)
		    }
		call fprintf(logfd[j],"\n")

		call fprintf(logfd[j],"%32.32s%s%14.14s%s%2.2s")
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		if (ntemp > 1) {
		    call fprintf(logfd[j],"%s%14.14s%s%7.7s%s%5.5s%s%4.4s")
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
		    }
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%14.14s%s%7.7s%s%5.5s%s%4.4s")
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
		    }
		call fprintf(logfd[j],"\n")
		}
	    }

# Print all templates on one line with long file names - no best velocity
	else if (rmode == 12) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%-14.14s%s%-4.4s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("julian_date")
		    call pargstr (spchar)
		    call pargstr ("best")
		do itemp = 1, ntemp {
		    call sprintf (tstring,16,"temp%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-14.14s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"vel%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%7s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"err%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-5s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"xcr%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%4s")
			call pargstr (spchar)
			call pargstr (tstring)
		    }
		call fprintf(logfd[j],"\n")

		call fprintf(logfd[j],"%32.32s%s%14.14s%s%4.4s")
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%-14.14s%s%-7.7s%s%-5.5s%s%4.4s")
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
		    }
		call fprintf(logfd[j],"\n")
		}
	    }

# Print emission line velocity plus multiple templates
	else if (rmode == 14) {
	    do j = 1, nlogfd {
		call fprintf(logfd[j],"%-32.32s%s%-14.14s")
		    call pargstr ("specid")
		    call pargstr (spchar)
		    call pargstr ("julian_date")
		call fprintf(logfd[j],"%s%-7.7s%s%-5s%s%s%s%s")
		    call pargstr (spchar)
		    call pargstr ("emvel")
		    call pargstr (spchar)
		    call pargstr ("emerr")
		    call pargstr (spchar)
		    call pargstr ("nl")
		    call pargstr (spchar)
		    call pargstr ("nf")
		call fprintf(logfd[j],"%s%s%s%-14.14s")
		    call pargstr (spchar)
		    call pargstr ("best")
		    call pargstr (spchar)
		    call pargstr ("besttemp")
		do itemp = 1, ntemp {
		    call sprintf (tstring,16,"temp%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-14.14s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"vel%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%7s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"err%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%-5s")
			call pargstr (spchar)
			call pargstr (tstring)
		    call sprintf (tstring,16,"xcr%d")
			call pargi (itemp)
		    call fprintf(logfd[j],"%s%4s")
			call pargstr (spchar)
			call pargstr (tstring)
		    }
		call fprintf(logfd[j],"\n")

		call fprintf(logfd[j],"%32.32s%s%14.14s")
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		call fprintf(logfd[j],"%s%7.7s%s%5.5s%s%2.2s%s%2.2s")
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		call fprintf(logfd[j],"%s%2.2s%s%14.14s")
		    call pargstr (spchar)
		    call pargstr (dashes)
		    call pargstr (spchar)
		    call pargstr (dashes)
		do itemp = 1, ntemp {
		    call fprintf(logfd[j],"%s%14.14s%s%7.7s%s%5.5s%s%4.4s")
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
			call pargstr (spchar)
			call pargstr (dashes)
		    }
		call fprintf(logfd[j],"\n")
		}
	    }
end
# Oct    1988	Stephen Levine wrote original version

# Oct 	 1990	Add hcv or bcv from file
# Nov	 1990	Add 1 line/template output
# Jan	 1991	Use zvel computed by main program
# Feb	 1991	Print 16-character template names
# May	 1991	Make all velocities double
# June	 1991	Use best template for mode 3
# June	 1991	Print Julian date for mode 3
# July	 1991	Add mode 5 for echelle testing
# Aug	 1991	Add mode 6 for echelle 
# Aug	 1991	Add mode 7 for error testing

# Jun 29 1992	Left-justify mode 6 names and add one digit to velocity

# Feb 12 1993	Print template filename in mode 4 correctly, fix mode 5 bug
# Dec  2 1993	Print mspec if non-zero

# Jan 20 1994	Print correct peak fitting mode
# Feb 11 1994	Add extra information to type 6 report for echelles
# Mar 23 1994	Pass template file id's through labelled common
# Apr 13 1994	Drop unused variables it and rmax
# Jun 24 1994	Add mode 8 to print both emission line and correlation results
# Jun 28 1994	Add mode 9 to print correlation velocity and error to meters
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 19 1994	Use wavelength range of best template
# Nov 17 1994	Add mode 10 to print multiple templates on one line
# Dec 19 1994	Print filter and emission line chopping info in mode 1

# Jan 11 1995	Add filter mode 3
# Jan 25 1995	Add mode 10 to allow long file names

# Feb 24 1995	Add quality flag to mode 7
# Mar 13 1995	Add absorption line removal note
# May  9 1995	Print all of modes 10 and 11 at once
# May 15 1995	Change all sz_fname to sz_line, which is 100 chars longer
# May 15 1995	Add template-driven option to divide, not subtract, continuum
# May 18 1995	Print 14-character templates in mode 11
# Jun  9 1995	Add mode 12 to avoid printing any template data twice
# Jun 19 1995	Get SPECID from rvsao common
# Jul  3 1995	Add mode 13 to print measured not fit numbers for echelle group
# Jul 13 1995	Add debugging argument to vcombine
# Jul 13 1995	Get Julian Date from JULDATE, not JULDAY
# Jul 13 1995	(except echelle modes which get HJDN from header)
# Jul 13 1995	Get UT from JULDATE
# Jul 17 1995	Fix regular Julian Date computation
# Aug  7 1995	Print heliocentric Julian Date if present, else gecoentric JD
# Aug 22 1995	Use absolute value of PKFRAC
# Sep 19 1995	Use width from header, not computed if TAA is zero
# Sep 21 1995	Handle indef values of spvel, spxvel and spevel
# Sep 25 1995	Fix mode 13 to print center correctly

# Feb  2 1996	Fix mode 6 to avoid right shift for large R values

# Jan 15 1997	Add JDH to JULDATE arguments
# Feb  4 1997	Drop declaration of C0; it is now in rvsao.com
# Apr 15 1997	Add mode 14 for template comparison to emission velocity
# May  2 1997	Always test against dindef, not INDEFD
# May  6 1997	Print message if no position or observation date
# Oct 23 1997	Print position and observation time on separate lines in mode 1
# Dec 17 1997   Change order of date from dd-mmm-yyyy to yyyy-mmm-dd
# Dec 17 1997   Use EQUINOX if it is present instead of EPOCH

# Jan 21 1998	Add mode 15 for Dave Latham--longer template name, m/sec
# Jul 14 1998	Increase filename length from 24 to 32 in modes 11, 12, and 14

# May 21 1999	Read DATE if DATE-OBS is no found
# Jul 23 1999	Put _ in place of space in incode and null in instrument
# Jul 29 1999	Add tab table option for single line report modes
# Sep 15 1999	Clean up tab table output

# Sep 20 2000	Print Shift instead of CZ if wavelength or pixel correlation
# Sep 27 2000	Print mode 11 and 12 results in g, not f format

# Feb  9 2001	Print aperture as separate column in report mode 12
# Aug  7 2002	Print template range and center as additional columns

# Jun  1 2006	Add report mode 17 for Soren Meibold for Hectochelle
# Dec  1 2006	Do not print best template info if only one template

# Apr  5 2007	Print last 27 characters of template name in mode 17 (-.fits)
# Apr  6 2007	Drop directory from template name
# Jun 21 2007	Add oldres.com to save previous results for mode 3 output

# Mar 10 2008	Add mode 18 for TRES multi-order echelle spectra
# May 30 2008	Use string constants more efficiently for Solaris limits

# Feb 20 2009	Assume decimal RA keyword value to be in degrees
# Apr  9 2009	Allow mode 17 template names to be up to 35 characters long

# Jul 21 2010	If no good template in mode 11, print zeroes

# May 25 2012	Add maxlength to print long filenames in mode 11
