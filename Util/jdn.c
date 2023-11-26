static char SccsID[] = "@(#)jdn.c	2.2\t11/29/99";
/* Routines to return Julian day from various forms.
 *
 * CREATION 06/29/87 Bill Wyatt
 */

#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#include "trig.h"
#include "poly.h"

#define ABS(x) ((x) < 0 ? -(x) : (x))

void doxperror();

/* NOTE: UT 1 Jan 1970, 0000 hours (midnight) == 2440587.5 Julian date, 
 * and Ultrix time is kept in seconds (and, supposedly,  microseconds)
 *  of UT from UT 1 Jan 1970, 0000 hours.
 */
#define JD1970 2440587.5

/* INPUT is a time value structure. OUTPUT is the JDN (to a fraction
 * of a second. The fraction is limited by the true accuracy
 * of the system time (currently about 1/60th sec).
 */
 
double
jdn_tval(tp)
     struct timeval *tp;
{
    return(JD1970 + (((double)(tp->tv_sec))/8.6400e+4) + 
           (((double)(tp->tv_usec))/8.64e+10));
}

/* Input a Julian Day Number, return a ptr to a tm structure 
 * as in gmtime(3).
 */
struct tm *
tm_jdn(jdn)
     double jdn;
{
    struct tm *gmtime();
    double rint();
    long ut = rint((jdn - JD1970) * 86400.0);

    return(gmtime(&ut));
}

/* Input a JDN, return a timeval pointer, number of secs since 1970.0 
 * Return is negative if JDN is before 1970.0.
 */
tval_jdn(jdn, tval)
     double jdn;
     struct timeval *tval;
{
    if(jdn < JD1970) return(-1);
    tval->tv_sec  = (jdn - JD1970) * 86400.0;
    tval->tv_usec = rint(((jdn - JD1970) * 86400.0 - tval->tv_sec) * 1.e6);
    return(0);
}

/* Input is a Julian Day Number. Output is an ASCII string giving UT time
 * in the form as given by asctime(3).
 */
char *
ut_jdn(jdn)
     double jdn;
{
    struct tm *gmtime();
    char *asctime();

    long ut = (jdn - JD1970) * 86400.0 + 0.5;
    return(asctime(gmtime(&ut)));
}

/* OUTPUT the current JDN */
double 
jdn_current()
{
    int gettimeofday();
    struct timeval tp;
    struct timezone tz;
    double jdn_tval();

    if(gettimeofday(&tp, &tz) == -1) {
        doxperror("jdn_now");
        return(0.0);
    }
    return(jdn_tval(&tp));
}

/* convert from GMT Year, Month, Day, Hours, Minutes, Seconds to JDN.
 * 
 * ******* ONLY VALID FROM JAN 0, 1900 through FEB 28th, 2100 ***********
 *     (because no century non-leap-year corrections are made)
 * 
 * N.B.: midnight on Jan 0, 1900 is Julian day 2415019.5 
 */
static unsigned short daynum[] = {
    0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 
  };

double
jdn_date(year, month, day, hour, min, sec, usec)
     int year;  /* since 1900 (i.e. 0 to 99 for 20th cent. */
     int month; /* range is 0 to 11 */
     int day;   /* range is 1 to 31  - no limit checking */
     int hour;  /* range 0 to 23 */
     int min, sec;  /* ranges 0 to 59 */
     int usec;  /* microseconds (because Unix gives it to us! */
{
    int ndays;

    if(month < 0 || month > 11) month = 0;
    ndays = (year * 365) + (year / 4) + daynum[month] + day -
      		(((year % 4) == 0 && month < 2 && year > 0) ? 1 : 0);
    
    return(2415019.5 + (double)ndays + 
           (((double)((hour * 3600) + (min * 60) + sec) + 
             (((double)usec) / 1.e+6)) / 86400.0));
}

/* INPUT is in `struct tm' form as defined in <sys/time.h>. The number of
 * microseconds is passed separately since the tm structure doesn't
 * hold it.
 */    
double
jdn_tmdate(tmp, usecs)
     struct tm *tmp;
     int usecs;
{
    double jdn_date();

    return(jdn_date(tmp->tm_year, tmp->tm_mon, 
                    tmp->tm_mday, tmp->tm_hour,
                    tmp->tm_min, tmp->tm_sec, usecs));
}

/* return solar ecliptic longitude in RADIANS */

#ifdef oldcode
static double 
	TU1[] = { 0.000303, 36000.768920, 279.69668 },
	TU2[] = { -3.0e-6, -1.50e-4, 35999.04975, 358.47583 };

struct spoly 
	tu1 = { 3, TU1, 0.0, 1.0 },
	tu2 = { 4, TU2, 0.0, 1.0 };

double
sun3long(jdn, msun, masun)
     double jdn;
     double *msun;   /* RETURN mean long. in RADIANS if non-NULL */
     double *masun;  /* RETURN mean anomaly in RADIANS if non-NULL */
{
    double evalpoly(), fmod();
        /* tu == Julian centuries from 0 Jan 1900 12:00 noon */
    double tu = (jdn - 2415020.0) / 36525.0; 
    double mlsun, ma, lsun;

    mlsun = fmod(evalpoly(&tu1, tu), 360.0) - 0.005; /* mean long. in deg. */
    ma    = fmod(evalpoly(&tu2, tu), 360.0);    /* mean anomaly */
    lsun = (sin(DEGTORAD(2.0 * ma)) * 0.020) +
      (sin(DEGTORAD(ma)) * 1.916) + mlsun;

    if(msun != (double *)NULL)  *msun  = DEGTORAD(mlsun);
    if(masun != (double *)NULL) *masun = DEGTORAD(ma);
    return(DEGTORAD(lsun));
}
#endif
double
sun3long(jdn, msun, masun)
     double jdn;
     double *msun;   /* RETURN mean long. in RADIANS if non-NULL */
     double *masun;  /* RETURN mean anomaly in RADIANS if non-NULL */
{
    double fmod();
    /* days from 0 Jan 1900 12:00 noon */
    double Du = (jdn - 2415020.0) / 10000.0;
    double du = jdn - 2415020.0;
    double mlsun, ma, lsun;

    /* mean long. in deg. */
    mlsun =
      fmod(279.696678 + 0.9856473354 * du + 2.267e-5 * Du * Du, 360.0);

    /* mean anomaly */
    ma    = 
      fmod(358.475833 + 0.9856002670 * du - 
	   Du * (1.12e-5 * Du + 7.e-8 * Du * Du), 360.0);

    lsun = (sin(DEGTORAD(2.0 * ma)) * 0.020) +
      (sin(DEGTORAD(ma)) * 1.916) + mlsun;

    if(msun != (double *)NULL)  *msun  = DEGTORAD(mlsun);
    if(masun != (double *)NULL) *masun = DEGTORAD(ma);
    return(DEGTORAD(lsun));
}

double
sunlong(jdn)
     double jdn;
{
    double sun3long();

    return(sun3long(jdn, (double *)NULL, (double *)NULL));
}


/* static double secliptic = 0.3978166;  /* sine and cosine of obliquity  */
/* static double cecliptic = 0.9174650;  /* of ecliptic (23 deg. 27 min.) */

/* make correction for geocentric Julian day to heliocentric Julian day */
double 
hjdn(gjdn, rastar, decstar)
     double gjdn;     /* topocentric Julian Day */
     double rastar,
            decstar;  /* coordinates in RADIANS */
{
    void  obl1950_();
    double sunlong();
    double slong = sunlong(gjdn);  /* solar ecliptic longitude in RADIANS */
    double ecliptic;
    double ca = cos(rastar),
    	   sa = sin(rastar);
    double cd = cos(decstar), 
    	   sd = sin(decstar);
    double cls = cos(slong), 
    	   sls = sin(slong);
    double cm = cos(0.0172021242 * (gjdn - 2442416.01329));
    double delta;

    obl1950_(&gjdn, &ecliptic);
    delta = -0.0057755 *  (1.0 - 0.01673 * cm) *
      ((cls * ca * cd) + (sls * ((sin(ecliptic) * sd) + 
                                 (cos(ecliptic) * cd * sa))));
    return(gjdn + delta);
}

double 
hjdn_iraf(gjdn, rastar, decstar)
     double gjdn;     /* topocentric Julian Day */
     double rastar,
            decstar;  /* coordinates in RADIANS */
{
    double tanom, rsun, delta;
    double starlat, starlong;
    double slong;
    double ecliptic;
    double fmod();
    void  obl1950_(), ast_coord();

    double t = (gjdn - 2415020.0) / 36525.0;
    double manom = 358.47583 +
      t * (35999.04975 - t * (0.000150 + t * 0.000003));
    double lperi = 101.22083 +
      t * (1.7191733 + t * (0.000453 + t * 0.000003));
    double eccen = 0.01675104 - 
      t * (0.00004180 + t * 0.000000126);

    manom = DEGTORAD(fmod(manom, 360.0));
    lperi = DEGTORAD(fmod(lperi, 360.0));
    
    /* true anomaly, radians */
    tanom = manom + 
      eccen * (2. - 0.25 * eccen * eccen) * sin(manom) +
	1.25 * eccen * eccen * sin(2. * manom) +
	  13./12. * eccen * eccen * eccen * sin(3. * manom);

    /* solar ecliptic longitude in RADIANS */
    slong = lperi + tanom + M_PI;

    obl1950_(&gjdn, &ecliptic);

    ast_coord(0., 0., -(M_PI_2), M_PI_2 - ecliptic,
	      rastar, decstar, &starlong, &starlat);
    rsun = (1.0 - eccen * eccen) / (1. + eccen * cos(tanom));

    delta = -0.00577 * rsun * cos(starlat) * cos(starlong - slong);
    return(gjdn + delta);
}

/* AST_COORD -- Convert spherical coordinates to new system.
 *
 * This procedure converts the longitude-latitude coordinates (a1, b1)
 * of a point on a sphere into corresponding coordinates (a2, b2) in a
 * different coordinate system that is specified by the coordinates of its
 * origin (ao, bo).  The range of a2 will be from -pi to pi.
 */
void
ast_coord (ao, bo, ap, bp, a1, b1, a2, b2)
     double  ao, bo;         /* Origin of new coordinates (radians) */
     double  ap, bp;         /* Pole of new coordinates (radians)   */
     double  a1, b1;         /* Coordinates to be converted (radians)*/
     double  *a2, *b2;         /* Converted coordinates (radians)     */
{ 
    double  sao, cao, sbo, cbo, sbp, cbp;
    double  x, y, z, xp, yp, zp, temp;

    x = cos (a1) * cos (b1);
    y = sin (a1) * cos (b1);
    z = sin (b1);
    xp = cos (ap) * cos (bp);
    yp = sin (ap) * cos (bp);
    zp = sin (bp);
 
    /* Rotate the origin about z. */
    sao = sin (ao);
    cao = cos (ao);
    sbo = sin (bo);
    cbo = cos (bo);
    temp = -xp * sao + yp * cao;
    xp = xp * cao + yp * sao;
    yp = temp;
    temp = -x * sao + y * cao;
    x = x * cao + y * sao;
    y = temp;

    /* Rotate the origin about y. */
    temp = -xp * sbo + zp * cbo;
    xp = xp * cbo + zp * sbo;
    zp = temp;
    temp = -x * sbo + z * cbo;
    x = x * cbo + z * sbo;
    z = temp;
 
    /* Rotate pole around x. */
    sbp = zp;
    cbp = yp;
    temp = y * cbp + z * sbp;;
    y = y * sbp - z * cbp;
    z = temp;
 
    /* Final angular coordinates. */
    *a2 = atan2 (y, x);
    *b2 = asin (z);
    return;
}

