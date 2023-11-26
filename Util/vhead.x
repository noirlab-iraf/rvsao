# File rvsao/Util/vhead.x
# May 14, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# vwhead (mspec,specim) saves combined velocity results in the image header

# Copyright(c) 1995-2008 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# VWHEAD saves combined velocity results in the object spectrum image header
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include "rvsao.h"
include	"emv.h"

procedure vwhead (mspec,specim)

int	mspec		# Object spectrum for which to write header in multispec file
pointer specim		# Object image header structure

char	qflag

char	keyword[SZ_HKWORD]	# header keyword
char	mline[SZ_LINE]		# one line character buffer
int	imaccf()

include	"rvsao.com"
include "results.com"
include "emv.com"
 
begin
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

# Save best Cz, error, and R to image header
	IM_UPDATE(specim) = YES
	if (mspec <= 0) {
	    call imaddd (specim,"VELOCITY",spvel, "cz (km/sec)")
	    call imaddd (specim, "Z", spvel / c0, "dwl / wl")
	    call imaddd (specim,"CZERR   ",sperr, "cz error (km/sec)")
	    call imaddd (specim, "ZERR", sperr / c0, "error in dwl / wl")
	    call sprintf (mline,SZ_LINE,"%c")
		call pargc (qflag)
	    call imastr (specim,"VELQUAL ",mline)
	    if (imaccf (specim, "BCV     ") == NO)
		call imaddd (specim, "BCV     ",spechcv)
	    else
		call imputd (specim,"BCV      ",spechcv)
	    if (specref > 0) {
		call sprintf (mline, SZ_HSTRING,
			      "Velocity from %.3fA %s emission line")
		    call pargd (waverest)
		    call pargstr (nmref[1,specref])
		call imastr (specim, "VELSET  ", mline)
		}
	    else if (specref < 0) {
		call sprintf (mline, SZ_HSTRING,
			      "Velocity from %.3fA %s absorption line")
		    call pargd (waverest)
		    call pargstr (nmabs[1,-specref])
		call imastr (specim, "VELSET  ", mline)
		}
	    else if (waverest != 0.d0) {
		call sprintf (mline, SZ_HSTRING,
			      "Velocity from line at %.3fA")
		    call pargd (waverest)
		call imastr (specim, "VELSET  ", mline)
		}
	    }
	else  if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		call pargi (mspec)
	    call sprintf (mline,SZ_HSTRING,"%.3f %.3f %c %.3f")
		call pargd (spvel)
		call pargd (sperr)
		call pargc (qflag)
		call pargd (spechcv)
	    call imastr (specim,keyword,mline)
	    }

end


procedure vrhead (mspec0, specim, sbcv, shcv)

int	mspec0		# Object spectrum to read from multispec file
pointer specim		# Object image header structure
double	sbcv		# Object spectrum barycentric velocity correction
double	shcv		# Object spectrum heliocentric velocity correction

char	keyword[SZ_HKWORD]	# header keyword
char	mline[SZ_HSTRING]	# one line character buffer
char	temp[SZ_LINE]
int	ip, ctod()
int	istat, sscan()
double	dtemp, dindef
char	qflag
int	imaccf()
int	mspec
bool	vfound, afound
 
include	"rvsao.com"
include "results.com"
include "emv.com"
include "oldres.com"
 
begin
	if (debug) call printf ("VRHEAD: Reading combined velocity info\n")
	mspec = mspec0
	dindef = INDEFD
	vfound = FALSE

# Initialize velocity parameters
	spvel = dindef
	sperr = dindef
	spvqual = 0
	spechcv = dindef

# Read solar system velocity correction(s) first
	sbcv = dindef
	shcv = dindef
	call imgdpar (specim,"BCV",sbcv)
	call imgdpar (specim,"HCV",shcv)

# Look for VELOCITY or APVELnn keywords in head and set flags if found
	vfound = FALSE
	afound = FALSE
	if (mspec == 0) {
	    call strcpy ("VELOCITY",keyword, SZ_HKWORD)
	    if (imaccf (specim, keyword) == YES)
		vfound = TRUE
	    else {
		mspec = 1
		call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		    call pargi (mspec)
		if (imaccf (specim, keyword) == YES)
		    afound = TRUE
		mspec = 0
		}
	    }
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		call pargi (mspec)
	    if (imaccf (specim, keyword) == YES)
		afound = TRUE
	    else {
		call strcpy ("VELOCITY",keyword, SZ_HKWORD)
		if (imaccf (specim, keyword) == YES) {
		    vfound = TRUE
		    mspec = 0
		    }
		}
	    }
	else {
	    call strcpy ("VELOCITY",keyword, SZ_HKWORD)
	    if (imaccf (specim, keyword) == YES)
		vfound = TRUE
	    }

# Read best Cz, error, and R from image header
	spr0 = 0.d0
	if (vfound) {
	    call imgdpar (specim,"VELOCITY",spvel)
	    call imgdpar (specim,"VELERR  ",sperr)
	    call imgdpar (specim,"CZERR   ",sperr)
	    call imgdpar (specim,"CZXCR   ",spr0)
	    qflag = '_'
	    call imgcpar (specim,"VELQUAL ",qflag)
	    }
	else if (afound) {
	    call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		call pargi (mspec)
	    call imgspar (specim,keyword,mline,SZ_HSTRING)
	    temp[1] = EOS
	    istat = sscan (mline)
		call gargd (spvel)
		call gargd (sperr)
		call gargc (qflag)
		call gargwrd (temp)
	    dtemp = dindef
	    ip = 1
	    if (ctod (temp,ip,dtemp) > 0) {
		sbcv = dtemp
		shcv = dtemp
		}
	    }
	if (vfound || afound) {
	    if (qflag == 'X')
		spvqual = 3
	    else if (qflag == '?')
		spvqual = 1
	    else if (qflag == 'Q')
		spvqual = 4
	    else
		spvqual = 0
	    }

	spvel0 = spvel
	sperr0 = sperr
	spechcv0 = sbcv
	spqual0 = qflag

	if (debug) {
	    call printf ("VRHEAD: velocity = %.4f, bcv = %.4f, hcv = %.4f\n")
		call pargd (spvel)
		call pargd (sbcv)
		call pargd (shcv)
	    }
	return
end


procedure vthead (mspec0, tempim, tvel, tbcv, thcv, debug)

int	mspec0		# Template spectrum to read from multispec file
pointer tempim		# Template spectrum image header structure
double	tvel		# Template spectrum velocity
double	tbcv		# Template spectrum barycentric velocity correction
double	thcv		# Template spectrum heliocentric velocity correction
bool	debug

char	keyword[SZ_HKWORD]	# header keyword
char	mline[SZ_HSTRING]	# one line character buffer
char	temp[SZ_LINE]
int	ip, ctod()
int	istat, sscan()
double	dtemp, dindef
double	terr
char	tflag
int	imaccf()
int	mspec
bool	vfound,afound
 
begin
	mspec = mspec0
	dindef = INDEFD
	vfound = FALSE

# Read solar system velocity correction(s) first
	tbcv = dindef
	thcv = dindef
	call imgdpar (tempim,"BCV",tbcv)
	call imgdpar (tempim,"HCV",thcv)

# Look for VELOCITY or APVELnn keywords in head and set flags if found
	vfound = FALSE
	afound = FALSE
	if (mspec == 0) {
	    call strcpy ("VELOCITY",keyword, SZ_HKWORD)
	    if (imaccf (tempim, keyword) == YES)
		vfound = TRUE
	    else {
		mspec = 1
		call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		    call pargi (mspec)
		if (imaccf (tempim, keyword) == YES)
		    afound = TRUE
		mspec = 0
		}
	    }
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		call pargi (mspec)
	    if (imaccf (tempim, keyword) == YES)
		afound = TRUE
	    else {
		call strcpy ("VELOCITY",keyword, SZ_HKWORD)
		if (imaccf (tempim, keyword) == YES) {
		    vfound = TRUE
		    mspec = 0
		    }
		}
	    }
	else {
	    call strcpy ("VELOCITY",keyword, SZ_HKWORD)
	    if (imaccf (tempim, keyword) == YES)
		vfound = TRUE
	    }

# Initialize velocity parameters
	tvel = dindef
	terr = dindef

# Read VELOCITY from header if aperture is 0 and keyword is present
	if (vfound)
	    call imgdpar (tempim,"VELOCITY",tvel)

# Read VELOCITY and BCV, if present, from image header if APVELn is present
	else if (afound) {
	    if (mspec == 0)
		mspec = 1
	    call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		call pargi (mspec)
	    call imgspar (tempim,keyword,mline,SZ_HSTRING)
	    temp[1] = EOS
	    istat = sscan (mline)
		call gargd (tvel)
		call gargd (terr)
		call gargc (tflag)
		call gargwrd (temp)
	    dtemp = dindef
	    ip = 1
	    if (ctod (temp,ip,dtemp) > 0) {
		tbcv = dtemp
		thcv = dtemp
		}
	    }

# If template velocity is not found, set it to zero
	if (tvel == dindef)
	    tvel = 0.d0

	if (debug) {
	    call printf ("VTHEAD: velocity = %.4f, bcv = %.4f, hcv = %.4f\n")
		call pargd (tvel)
		call pargd (tbcv)
		call pargd (thcv)
	    }
end

# QWHEAD writes qplot history into an IRAF image header

procedure qwhead (mspec, specim, prognam)

int	mspec		# Number of spectrum to read from multispec file
pointer	specim		# Object image header structure
char	prognam[ARB]	# Name of program doing qplotting

char	dstr[SZ_LINE]	# Date string
char	hstr[SZ_LINE]	# New history string
char	hstr0[SZ_LINE]	# Old history string
char	qflag
char	keyword[SZ_HKWORD]	# header keyword
char	mline[SZ_LINE]		# one line character buffer
int	imaccf()

include	"rvsao.com"
include "results.com"
include "emv.com"
 
begin
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

# Save best Cz, error, and R to image header
	IM_UPDATE(specim) = YES
	if (mspec <= 0) {
	    call imaddd (specim,"VELOCITY",spvel)
	    call imaddd (specim, "Z", spvel / c0)
	    call imaddd (specim,"CZERR   ",sperr)
	    call sprintf (mline,SZ_LINE,"%c")
		call pargc (qflag)
	    call imastr (specim,"VELQUAL ",mline)
	    if (imaccf (specim, "BCV     ") == NO)
		call imaddd (specim, "BCV     ",spechcv)
	    if (specref > 0) {
		call sprintf (mline, SZ_HSTRING,
			      "Velocity from %.3fA %s emission line")
		    call pargd (waverest)
		    call pargstr (nmref[1,specref])
		call imastr (specim, "VELSET  ", mline)
		}
	    else if (specref < 0) {
		call sprintf (mline, SZ_HSTRING,
			      "Velocity from %.3fA %s absorption line")
		    call pargd (waverest)
		    call pargstr (nmabs[1,-specref])
		call imastr (specim, "VELSET  ", mline)
		}
	    else if (waverest != 0.d0) {
		call sprintf (mline, SZ_HSTRING,
			      "Velocity from line at %.3fA")
		    call pargd (waverest)
		call imastr (specim, "VELSET  ", mline)
		}
	    }
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVEL%d")
		call pargi (mspec)
	    call sprintf (mline,SZ_HSTRING,"%.3f %.3f %c %.3f")
		call pargd (spvel)
		call pargd (sperr)
		call pargc (qflag)
		call pargd (spechcv)
	    call imastr (specim,keyword,mline)
	    }

#  Put date and program name into HISTORY
	call logtime (dstr, SZ_LINE)
	if (mspec > 0 && nmspec == 1) {
	    call sprintf (hstr,SZ_LINE,"rvsao.%s %s %s Quality[%s] set")
		call pargstr (prognam)
		call pargstr (VERSION)
		call pargstr (dstr)
		call pargstr (specnums)
	    }
	else {
	    call sprintf (hstr,SZ_LINE,"rvsao.%s %s %s VELOCITY = %.2f Q = %c")
		call pargstr (prognam)
		call pargstr (VERSION)
		call pargstr (dstr)
		call pargd (spvel)
		call pargc (qflag)
	    }
	hstr0[1] = EOS
	call imgspar (specim,"QPLOT",hstr0,SZ_HSTRING)
	if (hstr0[1] != EOS)
	    call imputh (specim, "HISTORY", hstr0)
	call imastr (specim, "QPLOT", hstr)
	return

end
# Jul 19 1995	New subroutines
# Sep 21 1995	Add QWHEAD to add qplot history line
# Oct  5 1995	Reformat history line to new standard
# Dec  5 1995	Initialize hcvel at start of VRHEAD

# Dec 18 1997	Set barycentric correction from BCV if APVELn not present

# Feb 13 1998	Do not read BCV if multispec file; do it in GETSPEC or GETTEMP
# Feb 13 1998	Add VTHEAD to read velocities for templates
# Feb 20 1998	Fix VTHEAD to deal with all conditions correctly
# Feb 23 1998	Fix VRHEAD to deal with all conditions correctly
# Apr  8 1998	Only set dtemp to INDEF when reading bcv and hcv from AVELnnn line

# Sep 16 2002	Set IM_UPDATE in vwhead() and qwhead() before calling imput()

# Aug 22 2005	Write Z to header

# Jan 27 2006	Write ZERR (Z error) to header
# Jul 13 2006	Do not try to create per-aperture keywords for apertures > 999

# Jan 31 2007	Write line used for velocity setting in vwhead() and qwhead()
# Jun 21 2007	In vwhead() always write BCV to header, even if there already

# May 14 2008	In vwhead(), do not write to header if ap > 999:1
