# File rvsao/Util/juldate.x
# June 27, 2011
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

include	<imhdr.h>
include <math.h>

# JULDATE -- Compute Julian Date at midtime of observation

procedure juldate (im,ut,jd,hjd,debug)

pointer	im		# image structure for spectrum
double	ut		# UT of observation midtime (returned)
double	jd		# Julian date of observation (returned)
double	hjd		# Heliocentric Julian date of observation (returned)
bool	debug		# True for diagnostic display

int	mm, dd, yyyy, jday
#int	 igreg, ja, jm, jy
double	dindef
double	ut1, ut2, exp
double	ra, dec, lt
int	imaccf()

begin
	dindef = INDEFD
	jd = 0.d0
	hjd = 0.d0

# Get date of observation
	mm = INDEFI
	dd = INDEFI
	yyyy = INDEFI
	if (imaccf (im, "DATE-OBS") == YES)
	    call imgdate (im, "DATE-OBS", mm, dd, yyyy)
	else if (imaccf (im, "DATE") == YES)
	    call imgdate (im, "DATE", mm, dd, yyyy)

# Get UT at middle of observation
	ut = dindef
	call imgdpar (im, "UTMID", ut)
	if (ut == dindef) {
	    ut1 = dindef
	    ut2 = dindef
	    call imgdpar (im,"UTOPEN", ut1)
	    if (ut1 == dindef)
		call imgdpar (im,"UTSTART", ut1)
	    if (ut1 == dindef)
		call imgdpar (im, "UT", ut1)
	    call imgdpar (im,"UTEND", ut2)
	    exp = 0.d0
	    call imgdpar (im,"EXPOSURE", exp)
	    if (exp == 0.d0)
		call imgdpar (im,"EXPTIME", exp)
	    if (ut1 == dindef) {
		if (ut2 != dindef)
		    ut = ut2 - (exp * 0.5 / 3600.)
		else
		    ut = 0.d0
		}
	    else if (ut2 != dindef)
		ut = (ut1 + ut2) * 0.5
	    else if (ut1 != dindef)
		ut = ut1 + (exp * 0.5 / 3600.)
	    else
		ut = 0.d0
	    }

# Find heliocentric Julian Day of observation
	hjd = 0.d0
	call imgdpar (im, "HJDN", hjd)
	if (hjd == 0.d0)
	    call imgdpar (im, "HJD", hjd)
	if (hjd == 0.d0)
	    call imgdpar (im, "BJD", hjd)
	jd = 0.d0
	call imgdpar (im, "GJDN", jd)
	if (jd == 0.d0)
	    call imgdpar (im, "GJD", jd)
	if (jd == 0.d0) {
	    call imgdpar (im, "MJD-OBS", jd)
	    if (jd == 0.d0)
		call imgdpar (im, "MJD", jd)
	    if (jd != 0.d0)
		jd = jd + 2400000.5
	    }

# If both JD and HJD are available in header, return (dropped due to FAST HJDN)
#	if (jd != 0.d0 && hjd != 0.d0)
#	    return

# Set Julian date from month, day, and year if it is not in image header
	if (jd == 0.d0) {

#	If year, month, and day are not all set, JD cannot be set
	    if (yyyy < 0 || mm == INDEFI || dd == INDEFI || yyyy == INDEFI)
		return

#	    igreg = 15 + (31 * (10 + (12 * 1582)))
#	    if (mm > 2) {
#		jy = yyyy
#		jm = mm + 1
#		}
#	    else {
#		jy = yyyy - 1
#		jm = mm + 13
#		}
#	    jday = int (365.25*jy) + int (30.6001*jm) + dd + 1720995
#	    if (dd+31*(mm+12*yyyy) >= igreg) {
#		ja = int (0.01 * jy)
#		jday = jday + 2 - ja + int (0.25 * ja)
#		}
	    jday = ( 1461 * ( yyyy + 4800 + ( mm - 14 ) / 12 ) ) / 4 +
		   ( 367 * ( mm - 2 - 12 * ( ( mm - 14 ) / 12 ) ) ) / 12 -
		   ( 3 * ( ( yyyy + 4900 + ( mm - 14 ) / 12 ) / 100 ) ) / 4 +
		   dd - 32075

# Compute Julian Date of observation
	    jd = double (jday) - 0.5d0 + (ut / 24.d0)
	    if (debug) {
		call printf("JULDATE: date = %d/%d/%d ut = %.6f = %h\n")
		    call pargi(dd)
		    call pargi(mm)
		    call pargi(yyyy)
		    call pargd(ut)
		    call pargd(ut)
	        }
	    }
	if (debug) {
	    call printf("JULDATE: date = %.6f\n")
		call pargd (jd)
	    }

# If Heliocentric Julian Date is already set, do not recompute it
	if (hjd != 0.0)
	    return

# Get pointing direction for light travel time correction to sun 
	ra = dindef
	call imgdpar (im, "RA", ra)
	dec = dindef
	call imgdpar (im,"DEC", dec)

# If there is a pointing direction in the image header, compute the Julian Date
# when the light from the object reached the Sun
	if (ra != dindef && dec != dindef) {
	    call jd2hjd (ra, dec, jd, lt, hjd)
	    if (debug) {
		call printf("JULDATE: RA = %.5f hours  Dec = %.5f degrees\n")
		    call pargd (ra)
		    call pargd (dec)
		call printf("JULDATE: heliocentric date = %.6f\n")
		    call pargd (hjd)
		}
	    }
	else
	    hjd = jd

end


# JD2HJD -- Helocentric Julian Day from UT (Geocentric) Julian Day
# Copied from noao.astutils.asttools.asthjd.x

procedure jd2hjd (ra, dec, jd, lt, hjd)

double	ra		# Right ascension of observation (hours)
double	dec		# Declination of observation (degrees)
double	jd		# Geocentric Julian date of observation
double	lt		# Light travel time in seconds
double	hjd		# Helocentric Julian Day

double	t, manom, lperi, oblq, eccen, tanom, slong, r, d, l, b, rsun

begin
	# JD is the geocentric Julian date.
	# T is the number of Julian centuries since J1900.

	t = (jd - 2415020d0) / 36525d0

	# MANOM is the mean anomaly of the Earth's orbit (degrees)
	# LPERI is the mean longitude of perihelion (degrees)
	# OBLQ is the mean obliquity of the ecliptic (degrees)
	# ECCEN is the eccentricity of the Earth's orbit (dimensionless)

	manom = 358.47583d0 +
	    t * (35999.04975d0 - t * (0.000150d0 + t * 0.000003d0))
	lperi = 101.22083d0 +
	    t * (1.7191733d0 + t * (0.000453d0 + t * 0.000003d0))
	oblq = 23.452294d0 -
	    t * (0.0130125d0 + t * (0.00000164d0 - t * 0.000000503d0))
	eccen = 0.01675104d0 - t * (0.00004180d0 + t * 0.000000126d0)

	# Convert to principle angles
	manom = mod (manom, 360.0D0)
	lperi = mod (lperi, 360.0D0)

	# Convert to radians
	r = DEGTORAD(ra * 15)
	d = DEGTORAD(dec)
	manom = DEGTORAD(manom)
	lperi = DEGTORAD(lperi)
	oblq = DEGTORAD(oblq)

	# TANOM is the true anomaly (approximate formula) (radians)
	tanom = manom + (2 * eccen - 0.25 * eccen**3) * sin (manom) +
	    1.25 * eccen**2 * sin (2 * manom) +
	    13./12. * eccen**3 * sin (3 * manom)

	# SLONG is the true longitude of the Sun seen from the Earth (radians)
	slong = lperi + tanom + PI

	# L and B are the longitude and latitude of the star in the orbital
	# plane of the Earth (radians)

	call concrd (double (0.), double (0.), double (-HALFPI),
	    HALFPI - oblq, r, d, l, b)

	# R is the distance to the Sun.
	rsun = (1. - eccen**2) / (1. + eccen * cos (tanom))

	# LT is the light travel difference to the Sun.
	lt = -0.005770d0 * rsun * cos (b) * cos (l - slong)
	hjd = jd + lt
end


# CONCRD -- Convert spherical coordinates to new system.
# Copied from noao.astutils.asttools.astcoord.x
#
# This procedure converts the longitude-latitude coordinates (a1, b1)
# of a point on a sphere into corresponding coordinates (a2, b2) in a
# different coordinate system that is specified by the coordinates of its
# origin (ao, bo).  The range of a2 will be from -pi to pi.

procedure concrd (ao, bo, ap, bp, a1, b1, a2, b2)

double	ao, bo		# Origin of new coordinates (radians)
double	ap, bp		# Pole of new coordinates (radians)
double	a1, b1		# Coordinates to be converted (radians)
double	a2, b2		# Converted coordinates (radians)

double	sao, cao, sbo, cbo, sbp, cbp
double	x, y, z, xp, yp, zp, temp

begin
	x = cos (a1) * cos (b1)
	y = sin (a1) * cos (b1)
	z = sin (b1)
	xp = cos (ap) * cos (bp)
	yp = sin (ap) * cos (bp)
	zp = sin (bp)

	# Rotate the origin about z.
	sao = sin (ao)
	cao = cos (ao)
	sbo = sin (bo)
	cbo = cos (bo)
	temp = -xp * sao + yp * cao
	xp = xp * cao + yp * sao
	yp = temp
	temp = -x * sao + y * cao
	x = x * cao + y * sao
	y = temp

	# Rotate the origin about y.
	temp = -xp * sbo + zp * cbo
	xp = xp * cbo + zp * sbo
	zp = temp
	temp = -x * sbo + z * cbo
	x = x * cbo + z * sbo
	z = temp

	# Rotate pole around x.
	sbp = zp
	cbp = yp
	temp = y * cbp + z * sbp
	y = y * sbp - z * cbp
	z = temp

	# Final angular coordinates.
	a2 = atan2 (y, x)
	b2 = asin (z)
end

# Jul 13 1995	New program
# Sep 19 1995	Format floating point numbers
# Sep 22 1995	Set UT to 0 if no times are set
# Oct 18 1995	Use either EXPOSURE or EXPTIME
# Jul 13 1995	New program
# Sep 19 1995	Format floating point numbers
# Sep 22 1995	Set UT to 0 if no times are set
# Oct 18 1995	Use either EXPOSURE or EXPTIME

# Jan  3 1997	Compute Heliocentric, not Geocentric Julian Date
# Jan  8 1997	Return Helicoentric Julian Date from HJDN in header, if there
# Jan 15 1997	Return both Heliocentric and Geocentric Julian Date
# Dec 17 1997	Write midtime of observation as UTMID
# Dec 17 1997	Always write keywords expected by other RVSAO programs

# Jan 15 1998	Allow input date to be MJD-OBS

# May 21 1999	Assume UT to be start, not end of exposure, UTEND to be end
# May 21 1999	Try DATE if DATE-OBS is not present

# Jun 12 2000	Do not use HJDN from header as it was not always correct in FAST

# Jul 19 2002	Check for UTSTART as well as UTOPEN; HJD and GJD, too

# Dec 14 2005	Fix but which set hjd from GJD (though GJDN set jd correctly)

# Dec 15 2006	Try new Julian Date computation

# Apr 10 2009	Use MJD if in header; return hjd = jd if no pointing direction

# Jun 27 2011	Do not recompute HJD if already in header
# Jun 27 2011	Check for BJD as well as HJD
