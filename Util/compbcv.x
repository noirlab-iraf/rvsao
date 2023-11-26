# File rvsao/Util/compbcv.x
# March 16, 2005
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1995-2005 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

include	<imhdr.h>

real procedure compbcv (im,ihcv,debug)

pointer	im		# image structure for spectrum
int	ihcv		# If =1, return barycentric velocity, else heliocentric
bool	debug		# If true, display diagnostics

double	ra, dec, epoch, dindef
double	dj, hdj, dlong, dlat, dalt, ut
double	dhcv, dbcv, dgcv
real	rbcv, rhcv
char	obsname[SZ_LINE]
char	str[SZ_LINE]
char	colon
pointer	obs		# pointer to observatory structure

int	stridx()
double	obsgetd()
pointer obsvopen()
errchk	obsvopen

begin
	dindef = INDEFD
	colon = ':'

# Julian Date of midtime of observation
	call juldate (im, ut, dj, hdj, debug)
	if (dj == 0.d0)
	    dj = hdj

# Direction of observed object (Convert RA to hours if not sexigesimal)
	ra = dindef
 	call imgdpar (im,"RA",ra)
	if (ra != dindef) {
	    call imgspar (im, "RA", str, SZ_LINE)
	    if (stridx (colon, str) == 0)
		ra = ra / 15.0
	    }
	dec = dindef
	call imgdpar (im,"DEC",dec)
	epoch = 1950.00
	call imgdpar (im,"EPOCH",epoch)
	call imgdpar (im,"EQUINOX",epoch)

# Position of observatory from header keywords
	dlong = dindef
 	call imgdpar (im,"SITELONG",dlong)
	if (dlong == dindef)
	    call imgdpar (im,"OBS-LONG",dlong)
	dlat = dindef
 	call imgdpar (im,"SITELAT",dlat)
	if (dlat == dindef)
	    call imgdpar (im,"OBS-LAT",dlat)
	dalt = dindef
 	call imgdpar (im,"SITEELEV",dalt)
	if (dalt == dindef)
	    call imgdpar (im,"OBS-ELEV",dalt)

# Position of observatory from IRAF database
	if (dlong == dindef) {
	    call imgspar (im,"OBSERVAT",obsname,SZ_LINE)
	    if (obsname[1] != EOS) {
		obs = obsvopen (obsname, NO)
		if (obs != NULL) {
		    dlat = obsgetd (obs, "latitude")
		    dlong = obsgetd (obs, "longitude")
		    dalt = obsgetd (obs, "altitude")
		    call obsclose (obs)
		    }
		}
	    }

# Barycentric velocity correction
	dhcv = 0.d0
	dbcv = 0.d0
	dgcv = 0.d0
#	call printf ("COMPBCV: dj= %.5f  ra= %.3h dec= %.2h\n")
#	    call pargd (dj)
#	    call pargd (ra)
#	    call pargd (dec)
#	call printf ("COMPBCV: lat= %.5f  long= %.6f alt= %.6f\n")
#	    call pargd (dlat)
#	    call pargd (dlong)
#	    call pargd (dalt)
	if (ra != dindef && dec != dindef && dj != 0.d0 &&
	    dlong != dindef && dlat != dindef && dalt != dindef) {
	    call bcv (dj, dlong,dlat,dalt, ra,dec,epoch, dbcv, dhcv, dgcv)
#	    call printf ("bcv is %f\n")
#	    call pargr (rbcv)
	    }
	else {
	    if (debug) {
		if (ra == dindef || dec == dindef) {
		    call printf ("COMPBCV: Missing pointing; ra= %f, dec= %f\n")
			call pargd (ra)
			call pargd (dec)
		    }
	    	if (dlong == dindef || dlat == dindef || dalt == dindef) {
		    call printf ("COMPBCV: Missing obs; long= %f, lat= %f, alt= %f\n")
			call pargd (dlong)
			call pargd (dlat)
			call pargd (dalt)
		    }
	    	if (dj == dindef) {
		    call printf ("COMPBCV: Missing time; jd= %f\n")
			call pargd (dj)
		    }
		}
	    dhcv = 0.d0
	    dbcv = 0.d0
	    dgcv = 0.d0
	    }
	if (debug) {
	    call printf ("BCVCORR: JD = %.4f RA = %.3h  Dec = %.2h\n")
		call pargd (hdj)
		call pargd (ra)
		call pargd (dec)
	    call printf ("BCVCORR: Obs lat = %.3h  long = %.3h  alt = %.2f\n")
		call pargd (dlat)
		call pargd (dlong)
		call pargd (dalt)
	    }
	rhcv = dhcv + dgcv
	rbcv = dbcv + dgcv

	if (debug) {
	    call printf ("COMPBCV: %.4f: hcv= %.3f, bcv = %.3f, gcv = %.3f\n")
		call pargd (hdj)
		call pargd (dhcv)
		call pargd (dbcv)
		call pargd (dgcv)
	    call printf ("COMPBCV:       helio.corr.= %.3f, bary.corr. = %.3f\n")
		call pargr (rhcv)
		call pargr (rbcv)
	    }
	if (ihcv == 1)
	    return (rbcv)       
	else
	    return (rhcv)       
end
# Apr 13 1994	Clean up code after running ftnchek
# Aug  3 1994	Change common and header from fquot to rvsao

# Jun 27 1995	Initialize all header parameters to indef
# Jul 13 1995	Compute Julian Date in a separate subroutine
# Jul 17 1995	Pass geocentric back from bcv, too.
# Sep 22 1995	Initialize velocity correction to zero

# Jan 12 1996	Add JD to diagnostic listing

# Jan 15 1997	Add HJD to JULDATE arguments and use for BCV/HCV computation
# Dec 17 1997	Use EQUINOX if it is present instead of EPOCH

# Mar 29 1999	Add code to read observatory from IRAF observatory database
# Jun 30 1999	Add more error reporting if missing parameters

# Feb  5 2002	Check for new standard WCS keywords OBS-LONG, OBS-LAT, and OBS-ELEV

# Mar 16 2005	Convert RA from degrees to hours if not sexigesimal
