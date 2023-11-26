# File rvsao/Emsao/emiplot.x
# November 19, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1994-2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Plot summary page for EMSAO

include	<gset.h>
include	<gio.h>
include <imhdr.h>
include	"rvsao.h"
include "emv.h"
define  SZ_MONTH	 3
define  SZ_TIME		24
 
procedure emiplot (gp, specfile, specim, mspec)

pointer	gp		# IRAF graphics environment structure
char	specfile[ARB]	# spectrum image file name
pointer	specim		# IRAF image descriptor structure
int	mspec		# Number of spectrum to read from multispec file

int	i, mm, dd, yyyy, itemp
real	lineno
double	ra, dec, ut, epoch
double	dindef
#double	tmpvel
char	text[SZ_LINE+1], text1[SZ_LINE+1], str[SZ_LINE+1]

char	qflag		# Quality flag character

char	wch		# Line weight indicator
char	colon
double	dj, hdj
int	imaccf(), stridx()

string	month	"JanFebMarAprMayJunJulAugSepOctNovDec"

include "rvsao.com"
include "results.com"
include "emv.com"

begin                                                          
	dindef = INDEFD
	colon = ':'

#  Set up graphics window for summary page
	call gseti (gp, G_WCS, 0)
	call gswind (gp, 0.0, 1.0, 0.0, 1.0)
	call gsview (gp, 0.0, 1.0, 0.0, 1.0)

#  Date and time of emission line search
	call logtime (text1, SZ_TIME)
	call sprintf(text,SZ_LINE,"%s %s %s")
	    call pargstr (taskname)
	    call pargstr (VERSION)
	    call pargstr (text1)
	call gtext (gp, 0.65, 0.01, text, "f=r;q=m;s=0.6")

#  File name
	call juldate (specim, ut, dj, hdj, debug)
	if (dj != 0.d0) {
	    call sprintf (text, SZ_LINE, "File: %s  JulDate: %.5f")
		call pargstr (specid)
		call pargd (dj)
	    }
	else {
	    call sprintf (text, SZ_LINE, "File: %s")
		call pargstr (specid)
	    }
	call gtext (gp, 0.025, 0.96, text, "f=r;q=m;s=0.8")
	
#  Object, RA, Dec
	ra = dindef
	dec = dindef
	epoch = 1950.0
	call imgdpar (specim, "RA", ra)
	if (ra != dindef) {
	    call imgspar (specim, "RA", str, SZ_LINE)
	    if (stridx (colon, str) == 0)
		ra = ra / 15.0
	    }
	call imgdpar (specim, "DEC", dec)
	call imgdpar (specim, "EPOCH", epoch)
	call imgdpar (specim, "EQUINOX", epoch)
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
 	call gtext (gp, 0.025, 0.91, text, "f=r;q=m;s=0.8")
                        
#  Switch windows and put in numbers from program     
	call gseti (gp, G_WCS, 3)
	call gswind (gp, 0.0, 1.0, 0.0, 1.0)
 	call gsview (gp, 0.67, 1.0, 0.0, 1.0)
	lineno = 0.96

#  Date of observation
	mm = INDEFI
	dd = INDEFI
	yyyy = INDEFI
	if (imaccf (specim, "DATE-OBS") == YES)
	    call imgdate (specim, "DATE-OBS", mm, dd, yyyy)
	else if (imaccf (specim, "DATE") == YES)
	    call imgdate (specim, "DATE", mm, dd, yyyy)

	if (mm != INDEFI && dd != INDEFI && yyyy != INDEFI) {
	    call sprintf(text,SZ_LINE,"%04d-%3.3s-%02d %011.2h")
		call pargi (yyyy)
		call pargstr (month[(mm - 1) * SZ_MONTH + 1])
		call pargi (dd)
		call pargd (ut)
	    }
	else
	    call sprintf(text,SZ_LINE,"(No observation date)")
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")

	call gsetr (gp, G_TXSIZE, 0.60)

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
		call sprintf(text,SZ_LINE, "Object file HCV: %7.3f")
		call pargd (spechcv)
	    case BCV: 
		call sprintf(text,SZ_LINE, "Object BCV: %7.3f")
		call pargd (spechcv)
	    default: 
		call sprintf(text,SZ_LINE, "Object BCV:    0.000")
	    }
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text, "f=r;q=m")

#  Initial guess at velocity, if this is what is plotted
	if (vplot == VGUESS) {
	    call sprintf(text, SZ_LINE,"*Guessed vel = %.2f")
		call pargd (gvel)
	    }

#  Combined velocity
	if (spvel != dindef) {
	    if (spvqual == 1)
		qflag = '?'
	    else if (spvqual == 2)
		qflag = '?'
	    else if (spvqual == 3)
		qflag = 'X'
	    else if (spvqual == 4)
		qflag = 'Q'
	    else
		qflag = '_'
	    call sprintf(text, SZ_LINE,
		" VELOCITY = %.2f +- %.2f km/sec %c")
		call pargd (spvel)
		call pargd (sperr)
		call pargc (qflag)
	    if (vplot == VCOMB || vplot == VGUESS)
		text[1] = '*'
	    lineno = lineno - 0.04
	    call gtext (gp, 0.0, lineno, text,"f=r;q=m")
	    }

#  Cross-correlation velocity
	if (spxvel != dindef) {
	    call sprintf(text, SZ_LINE,
	    " Corr vel = %.2f +- %.2f km/sec R= %.2f")
		call pargd (spxvel)
		call pargd (spxerr)
		call pargd (spxr)
	    if (vplot == VCORREL)
		text[1] = '*'
	    lineno = lineno - 0.04
	    call gtext (gp, 0.0, lineno, text,"f=r;q=m")
	    }

#  Emission line velocity
	call sprintf(text, SZ_LINE,
	" Emis vel = %.2f +- %.2f km/sec %d/%d lines")
	    call pargd (spevel)
	    call pargd (speerr)
	    call pargi (nfit)
	    call pargi (nfound)
	    if (vplot == VEMISS)
		text[1] = '*'
	lineno = lineno - 0.04
	call gtext (gp, 0.0, lineno, text,"f=r;q=m")

#  Template velocities and velocity shift from template
#	if (ntemp < 6) {
#	    lineno = lineno - 0.04
#	    call sprintf (text, SZ_LINE, "Template   vel.    hcv   peak")
#	    call gtext (gp, 0.0, lineno, text,"f=r;q=m")
#	    do itemp = 1, ntemp {
#		tmpvel = tempvel[itemp] + tempshift[itemp] + tvshift[itemp]
#		call sprintf(text, SZ_LINE, "%-7.7s %7.3f %7.3f %6.2f")
#		    call pargstr (tempid[1,itemp])
#		    call pargd (tmpvel)
#		    call pargd (temphcv[itemp])
#		    call pargd (cz[itemp])
#		lineno = lineno - 0.04
#		call gtext (gp, 0.0, lineno, text, "f=r;q=m")
#		}
#	    }

#  Template corrected heliocentric velocity and error
	if (ntemp > 0 && (ntemp + nfound < 17)) {
	    lineno = lineno - 0.06
	    call sprintf (text, SZ_LINE, "Template   CZ   error    R")
	    call gtext (gp, 0.0, lineno, text, "f=r;q=m;s=0.7")
	    do i = 1, ntemp {
		itemp = itr[i]
		call sprintf (text, SZ_LINE, "%-7.7s %9.3f %7.3f %6.2f")
		    call pargstr (tempid[1,itemp])
		    call pargd (zvel[itemp])
		    call pargd (czerr[itemp])
		    call pargd (czr[itemp])
		if (tschop[itemp])
		    call strcat (" -e",text,SZ_LINE)
		if (tempfilt[itemp] == 1)
		    call strcat (" -f",text,SZ_LINE)
		else if (tempfilt[itemp] == 2)
		    call strcat (" +h",text,SZ_LINE)
		else if (tempfilt[itemp] == 3)
		    call strcat (" -f +h",text,SZ_LINE)
		lineno = lineno - 0.04
		call gtext (gp, 0.0, lineno, text, "f=r;q=m")
		}                                                    
	    }                                                    

#  Emission line information if any were found
	if (nfound > 0) {
	    call sprintf(text,SZ_LINE,"Line  Rest   Obs.      CZ     error")
	    lineno = lineno - 0.06
	    call gtext (gp, 0.0, lineno, text,"f=r;q=m")
	    do i = 1, nfound {
		wtobs[i] = emparams[10,1,i]
		if (wtobs[i] > 0) {
		    if (override[i] < 0)
			wch = '-'
		    else if (override[i] > 0)
			wch = '+'
		    else
			wch = ' '
		    }
		else
		    wch = 'X'
		call sprintf(text,SZ_LINE,"%-4.4s %7.2f %7.2f %8.2f %5.2f %c")
		    call pargstr (nmobs[1,i])
		    call pargd (wlrest[i])
		    call pargd (wlobs[i])
		    call pargd ((c0*emparams[9,1,i]) + spechcv)
		    call pargd (c0*emparams[9,2,i])
		    call pargc (wch)
		if (linedrop[i] > 0) {
		    call sprintf(wch,1,"%d")
			call pargi (linedrop[i])
		    call strcat (wch,text,SZ_LINE)
		    }
		lineno = lineno - 0.04
		call gtext (gp, 0.0, lineno, text,"f=r;q=m")
		}
	    }
	else {
	    call sprintf(text,SZ_LINE,"No emission lines found\n")
	    lineno = lineno - 0.04
	    call gtext (gp, 0.0, lineno, text,"f=r;q=m")
	    }
	if (imaccf (specim, "VELSET") == YES) {
	    call imgspar (specim,"VELSET",text,SZ_LINE)
	    lineno = lineno - 0.04
	    call gtext (gp, 0.0, lineno, text,"f=r;q=m")
	    }

	call gflush(gp)
end                         
# Dec  7 1994	Read header info in GETSPEC, not here
# Dec 14 1994   Add emission line chopping flag per template
# Dec 14 1994   Add correlation filter flag per template

# Jan  6 1995	Add quality flag with combination velocity
# Jan 11 1995	Add filter flag 3 flags on correlations
# Feb 15 1995	Add heliocentric vel. correction to individual line velocities
# Jun 19 1995	Print SPECID instead of SPECFILE
# Jul 13 1995	Get Julian Date from JULDATE instead of JULDAY
# Jul 13 1995	Get UT from JULDATE
# Oct  2 1995	Keep reading of initial velocity GVEL in emfit
# Oct  3 1995	Sort templates by R-value; drop CLGET of DEVICE
# Oct  3 1995	Use PARGC for line weighting flag
# Oct 20 1995	Fix template flag printing

# Jan 24 1996	Write VELOCITY instead of Comb. Vel.
# Dec 13 1996	Check for indef, not 0 to see if velocities are set

# Jan 15 1997   Add HDJ to JULDATE arguments
# Feb 19 1997	Add flag to show why a line was rejected
# Apr  7 1997	Move template index sorting to EMPLOT
# May  2 1997	Always test against dindef, not INDEFD
# May  6 1997	Add message if no observation date
# May  6 1997	Add message if no position
# May  6 1997	Add epoch to position
# Oct  1 1997	Print 4-digit year in time stamp
# Dec 17 1997   Change order of date from dd-mmm-yyyy to yyyy-mmm-dd

# May 21 1999	Read DATE if DATE-OBS is no found

# Jan  7 2000	Don't print template information if lines + templates > 16

# Feb 20 2009	Assume decimal RA keyword value to be in degrees
# Nov 19 2009	Add * in front of VELOCITY if vplot is VGUESS
