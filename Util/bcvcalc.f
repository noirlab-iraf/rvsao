C SCCS ID = @(#)bcvcalc.f	1.1  3/18/93
c-------------------------------------------------------------------------
c
c       subroutine:     BCVCALC ( DJD, DLONG, DLAT, DELEV,
c                               DRA, DEC, DEQ, BC, HC)
c
c       purpose:        calculates the correction required to reduce
c                       observed (topocentric) radial velocities of a given
c                       star to the barycenter of the solar system.
c                       - includes correction for the effect of the earth's
c                       rotation.
c                       - the maximum error of this routine is not expected
c                       to be larger than 0.6 m/s.
c
c       input:          DJD = julian day number
c
c                       DLONG = geodetic longitude (radians)
c                       DLAT  = geodetic latitide (radians)
c                       DELEV = observer altitude (meters)
c
c                       DRA = right ascension of star (radians)
c                       DEC = declination of star (radians)
c                       DEQ = mean equator and equinox of coordinates, years
c                             (e.g., 1950.0)
c
c       output:         BC = barycentric correction (km/s)
c                       HC = barycentric correction (km/s)
c
c       author:         G. Torres (1989)
c			D. Mink (1990) return HCV as well as BCV
c                       W. Wyatt (1990) compartmentalization and nutation
c                                       for apparent sidereal time, elevation
c------------------------------------------------------------------------------

	SUBROUTINE BCVCALC (DJD, DLONG, DLAT, DELEV, DRA, DEC, DEQ, BC, HC)

c        implicit real*8 (d)

        real*8 djd, dlong, dlat, delev, dra, dec, deq, deqt
	real*8 dcc, dvelh, dvelb, dctrop, dc1900, dcbes,
     $         dct0, dpi, daukm, dra2, dec2, dst, dha,
     $         dgcvel, dbcvel, dhcvel

        real*4 bc, hc

        integer*2 k

        dimension dcc(3), dvelh(3), dvelb(3)

        data dctrop /365.24219572d0/, dcbes /0.313d0/,
     $       dc1900 /1900.0d0/,        dct0 /2415020.0d0/

        data dpi/3.1415926535897932d0/, daukm/1.4959787d08/

c Calculate local sidereal time (1 for apparent time, 0 for mean time)

        call sid1950 (djd, dlong, dst, 1)

c Precess and nutate R.A. and Dec. to mean equator and equinox of date 
c (the fourth argument to the call)
        deqt = (djd - dct0 - dcbes)/dctrop + dc1900
        call nutprec (dra, dec, deq, deqt, dra2, dec2)

c Calculate hour angle = local sidereal time - R.A.

        dha = dst - dra2
        dha = dmod(dha + 2.0d0*dpi , 2.0d0*dpi)

c Calculate observer's geocentric velocity

        call geovel (dlat, delev, dec2, -dha, dgcvel)

c Calculate components of Earth's barycentric velocity,
c   dvelb(i), i=1,2,3  in units of a.u./s

        call barvel (djd, 0.0d0, dvelh, dvelb)

c Project barycentric velocity to the direction of the star, and
c convert to km/s
        dbcvel = 0.0d0
        dhcvel = 0.0d0

        dcc(1) = dcos(dra2) * dcos(dec2)
        dcc(2) = dsin(dra2) * dcos(dec2)
        dcc(3) =              dsin(dec2)

        do 200 k=1,3
            dbcvel = dbcvel + dvelb(k)*dcc(k)*daukm
            dhcvel = dhcvel + dvelh(k)*dcc(k)*daukm
200     continue

c Add up both corrections
        bc = dbcvel + dgcvel
        hc = dhcvel + dgcvel

        return
        end
c------------------------------------------------------------------------------
