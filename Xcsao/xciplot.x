# File rvsao/Xcsao/xciplot.x
# February 20, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Stephen Levine and Jon Morse

# Copyright(c) 1990-2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Plot summary page for XCSAO

include	<gset.h>
include	<gio.h>
include	"rvsao.h"
include	"contin.h"
define  SZ_MONTH	 3
define  SZ_TIME		24
 
# Print additional information on XCSAO summary page

procedure xciplot (gp, specfile, image, mspec, strtwav, finwav)

pointer	gp		# IRAF graphics environment structure
char	specfile[ARB]	# spectrum image file name
pointer	image		# IRAF image descriptor structure
int	mspec		# Number of spectrum to read from multispec file
double	strtwav		# Starting wavelength of cross-correlation
double	finwav		# End wavelength of cross-correlation

int	j, i, mm, dd, yyyy, ntempx
real	lineno
double	ra, dec, ut, dj, hdj, epoch
double	tmpvel, fracpeak
double	dindef
char	text[SZ_LINE+1], text1[SZ_LINE+1], str[SZ_LINE+1]
char	colon
int	imaccf(), stridx()

string	month	"JanFebMarAprMayJunJulAugSepOctNovDec"

include "rvsao.com"               
include "results.com"

begin                                                          
	dindef = INDEFD
	colon = ':'

#  Set up graphics window
	call gseti (gp, G_WCS, 0)
	call gswind (gp, 0.0, 1.0, 0.0, 1.0)
	call gsview (gp, 0.0, 1.0, 0.0, 1.0)

#  Date and time of correlation
	call logtime (text1, SZ_TIME)
	call sprintf(text,SZ_LINE,"rvsao.xcsao %s %s")
	    call pargstr (VERSION)
	    call pargstr (text1)
#	call htsize (0.8)
	call gtext (gp, 0.01, 0.01, text, "f=r;q=m;s=0.6")
	call gflush (gp)

#  Type of peak used and the number of points
	if (pkfrac < 0.d0)
	    fracpeak = -pkfrac
	else
	    fracpeak = pkfrac
	if (fracpeak < 1.d0) {
	    call sprintf(text,SZ_LINE,"%4.2f-ht. peak fit, %d pts.")
	    call pargd (fracpeak)
	    call pargi (npmax)
	    }
	else {
	    call sprintf(text,SZ_LINE,"Peak fit with %d pts.")
	    call pargi (npmax)
	    }
	call gtext (gp, 0.35, 0.01, text, "f=r;q=m")
	call gflush (gp)

#  File name
	call sprintf (text, SZ_LINE, "File: %s")
	    call pargstr (specid)
#	call htsize (0.8)
	call gtext (gp, 0.025, 0.96, text, "f=r;q=m;s=0.8")
	call gflush (gp)
	
#  Object, RA, Dec
	ra = dindef
	dec = dindef
	epoch = 1950.0
	call imgdpar (image, "RA", ra)
	if (ra != dindef) {
	    call imgspar (image, "RA", str, SZ_LINE)
	    if (stridx (colon, str) == 0)
		ra = ra / 15.0
	    }
	call imgdpar (image, "DEC", dec)
	call imgdpar (image, "EPOCH", epoch)
	call imgdpar (image, "EQUINOX", epoch)
	if ((ra != dindef) && (dec != dindef)) {
	    call sprintf(text, SZ_LINE,
		"Object:%-10.10s RA: %011.2h DEC: %011.2h %.1f")
		call pargstr (specname)
		call pargd (ra)
		call pargd (dec)
		call pargd (epoch)
 	    }
	else {
	    call sprintf(text, SZ_LINE, "Object: %s  (No position)")
		call pargstr (specname)
	    }
#	call htsize (0.8)
 	call gtext (gp, 0.025, 0.91, text, "f=r;q=m;s=0.8")
	call gflush (gp)

#  Switch windows and put in numbers from program     
	call gseti (gp, G_WCS, 3)
	call gswind (gp, 0.0, 1.0, 0.0, 1.0)
 	call gsview (gp, 0.67, 1.0, 0.0, 1.0)
	lineno = 0.96

#  Date and time of observation
	dj = 0.0d0
	hdj = 0.0d0
	call juldate (image, ut, dj, hdj, debug)
	mm = INDEFI
	dd = INDEFI
	yyyy = INDEFI
	if (imaccf (image, "DATE-OBS") == YES)
	    call imgdate (image, "DATE-OBS", mm, dd, yyyy)
	else if (imaccf (image, "DATE") == YES)
	    call imgdate (image, "DATE", mm, dd, yyyy)
	yyyy = INDEFI
	if (imaccf (image, "DATE-OBS") == YES)
	    call imgdate (image, "DATE-OBS", mm, dd, yyyy)
	else if (imaccf (image, "DATE") == YES)
	    call imgdate (image, "DATE", mm, dd, yyyy)
	yyyy = INDEFI
	if (imaccf (image, "DATE-OBS") == YES)
	    call imgdate (image, "DATE-OBS", mm, dd, yyyy)
	else if (imaccf (image, "DATE") == YES)
	    call imgdate (image, "DATE", mm, dd, yyyy)
	if (mm != INDEFI && dd != INDEFI && yyyy != INDEFI) {
	    call sprintf(text,SZ_LINE,"%04d-%3.3s-%02d %011.2h")
		call pargi (yyyy)
		call pargstr (month[(mm - 1) * SZ_MONTH + 1])
		call pargi (dd)
		call pargd (ut)
	    }
	else
	    call sprintf(text,SZ_LINE,"No observation date")
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")
	call gflush (gp)

#  Template information
	call gsetr (gp, G_TXSIZE, 0.75)

# Heliocentric or geocentric Julian date 
	if (hdj != 0.d0) {
	    call sprintf (text, SZ_LINE, "HJD: %.5f")
		call pargd (hdj)
	    }
	else if (dj != 0.d0) {
	    call sprintf (text, SZ_LINE, "JulDate: %.5f")
		call pargd (dj)
	    }
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")
	call gflush (gp)

#  Heliocentric velocity correction
	switch (svcor) {
	    case HCV:
		call sprintf(text,SZ_LINE, "Object HCV: %7.3f")
		call pargd (spechcv)
	    case FBCV:
		if (specvb) {
		    call sprintf(text,SZ_LINE, "Object file BCV: %7.3f")
		    call pargd (spechcv)
		    }
		else {
		    call sprintf(text,SZ_LINE, "Object file HCV: %7.3f")
		    call pargd (spechcv)
		    }
	    case FHCV:
		if (specvb) {
		    call sprintf(text,SZ_LINE, "Object file BCV: %7.3f")
		    call pargd (spechcv)
		    }
		else {
		    call sprintf(text,SZ_LINE, "Object file HCV: %7.3f")
		    call pargd (spechcv)
		    }
	    case BCV: 
		call sprintf(text,SZ_LINE, "Object BCV: %7.3f")
		call pargd (spechcv)
	    default: 
		call sprintf(text,SZ_LINE, "Object BCV:    0.000")
	    }
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")

#  Wavelength limits
       	call sprintf(text, SZ_LINE, "lambda: %.1f %.1f")
	call pargd (strtwav)
	call pargd (finwav)
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")

#  Number of bins and tshift if nonzero
	call sprintf(text, SZ_LINE, "nbins: %d")
	    call pargi (npts)
	if (zpad)
	    call strcat ("x2",text,SZ_LINE)
	if (tshift != 0.0) {
	    call sprintf(text1, SZ_LINE, "  Tsh: %6.2f")
		call pargd (tshift)
	    call strcat (text1,text,SZ_LINE)
	    }
	if (tschop[itmax])
	    call strcat (" -em",text,SZ_LINE)
	if (tachop[itmax])
	    call strcat (" -ab",text,SZ_LINE)
	if (!tscont[itmax])
	    call strcat (" +tc",text,SZ_LINE)
	if (tconproc[itmax] == ZEROCONT)
	    call strcat (" divc",text,SZ_LINE)
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")

#  Transform filter limits
	call sprintf(text, SZ_LINE, "Filter: %d %d %d %d")
	call pargi(lo)	
	call pargi(toplo)
	call pargi(topnrn)
	call pargi(nrun)
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")

#  Apodization 
	call sprintf(text, SZ_LINE, "Frac. endmask: %.2f")        
	call pargr (han)
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")
	call gflush (gp)

#  Template information
	call gsetr (gp, G_TXSIZE, 0.65)

#  Template velocities and velocity shift from template
	if (ntemp < 6) {
	    lineno = lineno - 0.04
	    call sprintf (text, SZ_LINE, "Template   vel.    hcv   peak")
	    call gtext (gp, 0.0, lineno, text,"f=r;q=m")
	    do j = 1, ntemp {
		i = itr[j]
		tmpvel = tempvel[i] + tempshift[i] + tvshift[i]
		call sprintf(text, SZ_LINE, "%-7.7s %7.3f %7.3f %6.2f")
		    call pargstr (tempid[1,i])
		    call pargd (tmpvel)
		    call pargd (temphcv[i])
		    call pargd (cz[i])
		lineno = lineno - 0.04
		call gtext (gp, 0.0, lineno, text, "f=r;q=m")
		}
	    call gflush (gp)
	    }

#  Corrected heliocentric velocity and error
	lineno = lineno - 0.06
	if (correlate == COR_PIX || correlate == COR_WAV)
	    call sprintf (text, SZ_LINE, "Template Shift  error    R")
	else
	    call sprintf (text, SZ_LINE, "Template  CZ    error    R")
#	call htsize (0.7)
	call gtext (gp, 0.0, lineno, text, "f=r;q=m;s=0.7")
	if (ntemp > 10)
	    ntempx = 10
	else
	    ntempx = ntemp
	do j = 1, ntempx {
	    i = itr[j]
	    call sprintf (text, SZ_LINE, "%-7.7s %7.3f %7.3f %5.2f ")
		call pargstr (tempid[1,i])
		call pargd (zvel[i])
		call pargd (czerr[i])
		call pargd (czr[i])
	    if (toverlap[i])
		call strcat ("+w",text,SZ_LINE)
	    if (tschop[i])
		call strcat ("-e",text,SZ_LINE)
	    if (tachop[i])
		call strcat ("-a",text,SZ_LINE)
	    if (tempfilt[i] == 1)
		call strcat ("-f",text,SZ_LINE)
	    else if (tempfilt[i] == 2)
		call strcat ("+h",text,SZ_LINE)
	    else if (tempfilt[i] == 3)
		call strcat ("-f+h",text,SZ_LINE)
	    if (!tscont[i])
		call strcat ("+tc",text,SZ_LINE)
	    if (tconproc[i] == ZEROCONT)
		call strcat ("/c",text,SZ_LINE)
	    lineno = lineno - 0.04
	    call gtext (gp, 0.0, lineno, text, "f=r;q=m")
	    }                                                    
	call gflush(gp)
end                         
 
# June	1988	Stephen Levine (madraf::levine)
#               alter peak window for graph -- Morsey
# April	1990	Doug Mink (mink@cfa.harvard.edu)
#		alter display font sizes so everything fits
#		alter range for velocity plot
# Aug	1990	Add multiple peak fits and error calcs
# Oct	1990	Add file BCV or HCV
# June	1991	Simplify code
# July	1991	Scale data plot from 0
# Sept	1991	Make device a variable
# Sept	1991	Add template velocity shift
# Nov 20 1991	Shrink axis tick mark labels
# Dec  3 1991	Compute default width of velocity plot using peak width

# Feb 14 1992	Print mspec instead of eline for multispec spectra
# 		Print only result line if >5 templates
# Apr 20 1992	Pass mspec as argument
# Apr 21 1992	Change xcplot to xcorplot and spplot to specplot
#		Rename lsretup to openplot and lsrend to closeplot
# Aug 17 1992	Use htext to use Hershey fonts instead of ugly IRAF font
# Nov 24 1992	Go back to gtext instead of htext
# Dec  1 1992	If correlation error is near zero, set plot width to 300 km/s

# May 20 1993	Use different variable name for corrected template velocity
# Aug  9 1993	Plot spectrum within parameter limits if given

# Mar 23 1994	Print template id instead of object
# Apr 13 1994	Drop c0 as it is not used
# Jun 23 1994	Pass velocities in fquot insted of getim labelled common
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 15 1994	Move scaling and plot device handling to Util/plotutil.x
# Aug 16 1994	Fix zero-padding bug
# Aug 17 1994	Rename XCIPLOT for consistency
# Aug 19 1994	Add emission line chopping flag per template
# Dec 14 1994	Add correlation filter flag per template
# Dec 15 1994	Drop peak mode from summary line
# Dec 19 1994	Add high filter flag per template

# Jan 11 1995	Add additional filter flag value
# Mar 13 1995	Add notes for absorption line removal
# Mar 15 1995	Add note for template continuum left in
# May 10 1995	Print info for first 12 templates only
# May 10 1995	Sort templates by R-value for output
# May 15 1995	Note if continuum has been divided instead of subtracted
# Jun 19 1995	Use SPECID instead of making it up
# Jun 26 1995	Note if full overlap region is being used per template
# Mar 14 1997	Print HJD if available, else JD, else no date
# Apr  7 1997	Move R-value sorting to XCPLOT
# May  2 1997	Always test against dindef, not INDEFD
# May  6 1997	Add message if no observation date
# May  6 1997	Add message if no position
# Sep 25 1997	Flush graphics buffer more often
# Oct  1 1997	Print 4-digit year in time stamp
# Dec 17 1997   Change order of date from dd-mmm-yyyy to yyyy-mmm-dd
# Dec 17 1997   Use EQUINOX if it is present instead of EPOCH

# Apr  8 1999	Move Julian Date to under date on right side; shrink text
# May 21 1999	Read DATE if DATE-OBS is no found
# Aug 18 1999	Use int tconproc instead of bool tdivcon

# Sep 20 2000	Print Shift instead of CZ if wavelength or pixel correlation

# Feb 20 2009	Assume decimal RA keyword value to be in degrees
