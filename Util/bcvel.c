static char SccsID[] = "@(#)bcvel.c	2.1\t04/22/99";

/* Routines to return heliocentric and barycentric velocity
 * 
 */

#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#include "trig.h"
#include "poly.h"

static double secliptic = 0.3978166;  /* sine and cosine of obliquity  */
static double cecliptic = 0.9174650;  /* of ecliptic (23 deg. 27 min.) */

/* Heliocentric velocity correction, km/sec */
/*
 *   *****  NEW CALCULATION BASED ON BARVEL, July 1990  *****
 *   *                                                      *
 *   *  This velocity should be better than 0.6 meters/sec  *
 *   *                                                      *
 *   ********************************************************
 */
void
bcv_vel(jdn, rastar, decstar, epoch, longitude, latitude, elevation, bcv, hcv)
     double jdn;               /* Julian day of observation */
     double rastar, decstar;   /* Star coordinate, radians */
     double epoch;             /* Star coords epoch, in Years */
     double longitude;         /* in radians west of 0 */
     double latitude;          /* in radians */
     double elevation;         /* observatory altitude, meters */
     double *bcv;              /* RETURNED = barycentric velocity */
     double *hcv;              /* RETURNED = heliocentric velocity */
{
    float bcvf, hcvf;
    void bcvcalc_();  /* FORTRAN routine for doing velocities */

    bcvcalc_(&jdn, &longitude, &latitude, &elevation, 
	     &rastar, &decstar, &epoch, &bcvf, &hcvf);
    *bcv = (double)bcvf;
    *hcv = (double)hcvf;

    return;
}

void
sun_coords(jdn, rasun, decsun)
     double jdn;
     double *rasun, *decsun;  /* RETURN - coords in radians */
{
    double lsun, mlsun, masun;
    double lperisun;
    double epoch;
    double sun3long();
    
    lsun = sun3long(jdn, &mlsun, &masun);
    lperisun = mlsun - masun + M_PI;   /* longitude of perihelion */

    *rasun = atan(tan(lsun) * cecliptic);
    /* Restore to 4 quad, based on lsun */
    if(lsun <= M_PI_2)             /* no-op */ ;
    else if(lsun <= 1.5 * M_PI)    *rasun += M_PI;
    else                           *rasun += M_PI * 2;
      
    *decsun = asin(sin(lsun) * secliptic);
    return;
}

/* SKY obs. Heliocentric velocity correction, km/sec;
 *   The sky correction uses the position of the Sun, not the
 *   RA and Dec the telescope is pointing at.
 */
void
skybcv_vel(jdn, longitude, latitude, elevation, bcv, hcv)
     double jdn;               /* Julian day of observation */
     double longitude;         /* in radians west of 0 */
     double latitude;          /* in radians */
     double elevation;         /* observatory altitude, meters */
     double *bcv;              /* RETURNED = barycentric velocity */
     double *hcv;              /* RETURNED = heliocentric velocity */
{
    float bcvf, hcvf;
    double rasun, decsun;       /* Sun coordinate, radians */
    double epoch;
    void sun_coords(), bcvcalc_();

    sun_coords(jdn, &rasun, &decsun);

    /* epoch for solar coordinates is epoch of date */
    epoch =  (jdn - 2415020.0 - 0.313)/365.24219572 + 1900.0;

    bcvcalc_(&jdn, &longitude, &latitude, &elevation, 
	     &rasun, &decsun, &epoch, &bcvf, &hcvf);
    *bcv = (double)bcvf;
    *hcv = (double)hcvf;

    return;
}


/* OLD ROUTINE, WITH JULY 1990 CORRECTIONS, Good to about 30 meters/sec.
 * N.B. - removed the corrections for archival comparisons, August 1990.
 *    the routine is now offset by about +40 meters/sec, so is good to
 *    about 70 meters/sec.
 */
double 
hcvel(jdn, rastar, decstar, longitude, latitude, dvorb, dvrot)
     double jdn;               /* Julian day of observation */
     double rastar, decstar;   /* Star coordinate, radians, CURRENT EPOCH */
     double longitude;         /* in radians west of 0 */
     double latitude;          /* in radians */
     double *dvorb;   /* RETURN if non-NULL = orbital part */
     double *dvrot;   /* RETURN if non-NULL = rotational part */
{
    double rsidtime(), sun3long();
    double sra  = sin(rastar),  cra  = cos(rastar);
    double sdec = sin(decstar), cdec = cos(decstar);
    double lsun, mlsun, masun;
    double lperisun;
    double dvo, dvr;
    
    lsun = sun3long(jdn, &mlsun, &masun);
    lperisun = mlsun - masun + M_PI;   /* longitude of perihelion */

    /* NOTE - this constant should be 29.79, but since this routine has been
     * superseeded by skybcv_vel, the old value is left in as a check
     * against NOVA-derived files.
     */
    dvo = 29.75 * 
         (sin(lsun) * cdec * cra -
          (sdec * secliptic + cdec * cecliptic * sra) * cos(lsun))
       - 0.49736 *
         (sin(lperisun) * cdec * cra -
          (sdec * secliptic + cdec * sra * cecliptic) * cos(lperisun));

    dvr = -0.46509 * cos(latitude) * cdec * 
      sin(rsidtime(jdn, longitude, 0) - rastar);

    if(dvorb != (double *)NULL)  *dvorb = dvo;
    if(dvrot != (double *)NULL)  *dvrot = dvr;

    return(dvo + dvr);
}

/* OLD ROUTINE, WITH JULY 1990 CORRECTIONS, Good to about 30 meters/sec
 * N.B. - removed the corrections for archival comparisons, August 1990.
 *    the routine is now offset by about +40 meters/sec, so is good to
 *    about 70 meters/sec.
 */

/* SKY obs. Heliocentric velocity correction, km/sec;
 *   The sky correction uses the position of the Sun, not the
 *   RA and Dec the telescope is pointing at.
 */
double 
skyvel(jdn, longitude, latitude, dvorb, dvrot)
     double jdn;               /* Julian day of observation */
     double longitude;         /* in radians west of 0 */
     double latitude;          /* in radians */
     double *dvorb;   /* RETURN if non-NULL = orbital part */
     double *dvrot;   /* RETURN if non-NULL = rotational part */
{
    double rsidtime(), sun3long();
    double rasun, decsun;       /* Sun coordinate, radians */
    double sra,  cra;
    double sdec, cdec;
    double lsun, mlsun, masun;
    double lperisun;
    double dvo, dvr;
    
    lsun = sun3long(jdn, &mlsun, &masun);
    lperisun = mlsun - masun + M_PI;   /* longitude of perihelion */

    rasun = atan(tan(lsun) * cecliptic);
    /* Restore to 4 quad, based on lsun */
    if(lsun <= M_PI_2)       /* no-op */ ;
    else if(lsun <= 1.5 * M_PI)    rasun += M_PI;
    else                           rasun += M_PI * 2;
      
    decsun = asin(sin(lsun) * secliptic);
    
    sra  = sin(rasun);
    cra  = cos(rasun);
    sdec = sin(decsun);
    cdec = cos(decsun);

    /* NOTE - this constant should be 29.79, but since this routine has been
     * superseeded by skybcv_vel, the old value is left in as a check
     * against NOVA-derived files.
     */
    dvo = 29.75 * 
         (sin(lsun) * cdec * cra -
          (sdec * secliptic + cdec * cecliptic * sra) * cos(lsun))
       - 0.49736 *
         (sin(lperisun) * cdec * cra -
          (sdec * secliptic + cdec * sra * cecliptic) * cos(lperisun));

    dvr = -0.46509 * cos(latitude) * cdec * 
      sin(rsidtime(jdn, longitude, 0) - rasun);

    if(dvorb != (double *)NULL)  *dvorb = dvo;
    if(dvrot != (double *)NULL)  *dvrot = dvr;

    return(dvo + dvr);
}

