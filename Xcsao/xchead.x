# File rvsao/Xcsao/xchead.x
# July 13, 2006
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1994-2006 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# xcwhead (mspec, specim) saves cross-correlation results in a spectrum header
# xcrhead (mspec, specim) reads cross-correlation results from a spectrum header
 
# XCWHEAD saves cross-correlation results in the object spectrum image header
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include "rvsao.h"
include	"emv.h"

procedure xcwhead (mspec,specim)

int	mspec		# Object spectrum to read from multispec file
pointer specim		# Object image header structure

char	dstr[SZ_LINE]	# Date string
char	hstr[SZ_LINE]	# New history string
char	hstr0[SZ_LINE]	# Old history string

int	itemp, izpad
char	keyword[SZ_HKWORD]	# header keyword
char	mline[SZ_LINE]		# one line character buffer
double	width, tmpvel

include	"rvsao.com"
include "results.com"
include "emv.com"
 
begin
	if (zpad)
	    izpad = 2
	else
	    izpad = 1

# Save best Cz, error, and R to image header
	if (mspec <= 0) {
	    call imaddd (specim,"CZXC    ",spxvel,"Cross-correlation redshift velocity (km/sec)")
	    call imaddd (specim,"ZXC    ",spxvel / c0,"Cross-correlation z = dwl / wl")
	    call imaddd (specim,"CZXCERR ",spxerr,"Cross-correlation redshift velocity error (km/sec)")
	    call imaddd (specim,"CZXCR   ",spxr,"Cross-correlation R-value")
	    call imastr (specim,"BESTTEMP",tempid[1,itmax])
	    call sprintf (mline,SZ_HSTRING,"%d %d %d %d %d %d %d %.2f")
		call pargi (ntemp)
		call pargi (npts)
		call pargi (izpad)
		call pargi (lo)
		call pargi (toplo)
		call pargi (topnrn)
		call pargi (nrun)
		call pargr (han)
	    call imastr (specim,"CZPARAMS",mline)
	    }
	else {
	    call sprintf (keyword,SZ_HKWORD,"APVXC%d")
		call pargi (mspec)
	    mline[1] = EOS
	    call sprintf (mline,SZ_LINE,
		"%.3f %.3f %7g %d %s %d %d %d %d %d %d %d %.2f")
		call pargd (spxvel)
		call pargd (spxerr)
		call pargd (spxr)
		call pargi (itmax)
		call pargstr (tempid[1,itmax])
		call pargi (ntemp)
		call pargi (npts)
		call pargi (izpad)
		call pargi (lo)
		call pargi (toplo)
		call pargi (topnrn)
		call pargi (nrun)
		call pargr (han)
#	    call printf ("%s = %s\n")
#		call pargstr (keyword)
#		call pargstr (mline)
	    call imastr (specim,keyword,mline)
	    }
	    call imputh (specim,"COMMENT","Template information from rvsao.xcsao:")
	    call imputh (specim,"COMMENT","file name, cross-correlation velocity, error, R-value")
	    call imputh (specim,"COMMENT","correlation peak height (0-1) and width (km/sec)")
	    call imputh (specim,"COMMENT","velocity relative to template, template BCV")
	    call imputh (specim,"COMMENT","template em. chopped, template filtered")

# Template information
	do itemp = 1, ntemp {
	    if (mspec <= 0) {
		call sprintf (keyword,SZ_HKWORD,"XCTEMP%d")
		    call pargi (itemp)
		}
	    else {
		call sprintf (keyword,SZ_HKWORD,"XT%d_%d")
		    call pargi (mspec)
		    call pargi (itemp)
		}
	    if (taa[itemp] != 0) {
		dlogw =  1.d0 / (dlog (10.d0) * taa[itemp])
		width = c0 * ((10.d0 ** (twdth[itemp] * dlogw)) - 1.d0)
		}
	    else
		width = tvw[itemp]
	    tmpvel = tempvel[itemp] + tempshift[itemp] + tvshift[itemp]
	    mline[1] = EOS
	    call sprintf (mline,SZ_LINE,
		"%s %.3f %.3f %7g %.3f %.3f %.3f %.3f %1d %1d")
		call pargstr (tempid[1,itemp])
		call pargd (zvel[itemp])
		call pargd (czerr[itemp])
		call pargd (czr[itemp])
		call pargd (thght[itemp])
		call pargd (width)
		call pargd (tmpvel)
		call pargd (temphcv[itemp])
		call pargb (tschop[itemp])
		call pargi (tempfilt[itemp])
#	    call printf ("%s = %s\n")
#		call pargstr (keyword)
#		call pargstr (mline)
	    call imastr (specim,keyword,mline)
	    }

# Put date and program name into XCSAO, saving old XCSAO value as HISTORY
	call logtime (dstr, SZ_LINE)
	if (mspec > 0 && nmspec == 1) {
	    call sprintf (hstr,SZ_LINE,"rvsao.xcsao %s %s CZXC [%s]")
		call pargstr (VERSION)
		call pargstr (dstr)
		call pargstr (specnums)
	    }
	else {
	    call sprintf (hstr,SZ_LINE,"rvsao.xcsao %s %s CZXC = %.2f R = %.2f")
		call pargstr (VERSION)
		call pargstr (dstr)
		call pargd (spxvel)
		call pargd (spxr)
	    }
	hstr0[1] = EOS
	call imgspar (specim,"XCSAO",hstr0,SZ_HSTRING)
	if (hstr0[1] != EOS)
	    call imputh (specim, "HISTORY", hstr0)
	call imastr (specim, "XCSAO", hstr)
	IM_UPDATE(specim) = YES
end


# XCRHEAD reads cross-correlation results from a spectrum header

procedure xcrhead (mspec0,specim)

int	mspec0		# Object spectrum to read from multispec file
pointer specim		# Object image header structure

int	itemp
char	keyword[SZ_HKWORD]	# header keyword
char	mline[SZ_HSTRING]	# one line character buffer
char	bestid[SZ_FNAME]	# name of best template
int	npad, i
int	istat, sscan()
int	imaccf()
int	lbest, ltemp, lmax, strlen(), strncmp()
int	mspec, nchar
 
include	"rvsao.com"
include "results.com"
include "emv.com"
 
begin
	if (debug) call printf ("XCRHEAD: Reading cross-correlation info\n")
	mspec = mspec0

# Initialize velocity parameters
	spxvel = INDEFD
	spxerr = INDEFD
	spxr = 0.d0

# Return if cross-correlation velocity keyword not found in header
	ntemp = 0
	if (mspec == 0)
	    call strcpy ("CZXC    ",keyword, SZ_HKWORD)
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVXC%d")
	    call pargi (mspec)
	    }
	if (imaccf (specim, keyword) == NO)
	    return

	if (mspec <= 0) {
	    call imgdpar (specim,"CZXC    ",spxvel)
	    call imgdpar (specim,"CZXCERR ",spxerr)
	    call imgdpar (specim,"CZXCR   ",spxr)
	    call imgspar (specim,"BESTTEMP",bestid,SZ_FNAME)
	    lbest = strlen (bestid)
	    itmax = 0
	    mline[1] = '0'
	    mline[2] = EOS
	    if (imaccf (specim, "CZPARAMS") == YES) {
		call imgspar (specim,"CZPARAMS",mline,SZ_HSTRING)
		istat = sscan (mline)
		    call gargi (ntemp)
		if (ntemp > 0) {
		    call gargi (npts)
		    call gargi (npad)
		    if (npad == 2)
			zpad = TRUE
		    else
			zpad = FALSE
		    call gargi (lo)
		    call gargi (toplo)
		    call gargi (topnrn)
		    call gargi (nrun)
		    call gargr (han)
		    }
		}
	    else {
		ntemp = 0
		call imgipar (specim,"NTEMP   ",ntemp)
		}
	    }
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVXC%d")
		call pargi (mspec)
	    call imgspar (specim,keyword,mline,SZ_HSTRING)
	    lbest = 0
	    istat = sscan (mline)
		call gargd (spxvel)
		call gargd (spxerr)
		call gargd (spxr)
		call gargi (itmax)
		call gargwrd (tempid[1,itmax],SZ_FNAME)
		call gargi (ntemp)
		call gargi (npts)
		call gargi (npad)
		if (npad == 2)
		    zpad = TRUE
		else
		    zpad = FALSE
		call gargi (lo)
		call gargi (toplo)
		call gargi (topnrn)
		call gargi (nrun)
		call gargr (han)
	    }

# Template information
	if (ntemp > 0) {
	    do itemp = 1, ntemp {
		if (mspec <= 0) {
		    call sprintf (keyword,SZ_HKWORD,"XCTEMP%d")
			call pargi (itemp)
		    if (imaccf (specim, keyword) == NO) {
			call sprintf (keyword,SZ_LINE,"XT1_%d")
			    call pargi (itemp)
			}
		    }
		else if (mspec < 1000) {
		    call sprintf (keyword,SZ_LINE,"XT%d_%d")
			call pargi (mspec)
			call pargi (itemp)
		    if (mspec == 1 && imaccf (specim, keyword) == NO) {
			call sprintf (keyword,SZ_HKWORD,"XCTEMP%d")
			    call pargi (itemp)
			}
		    }
		if (imaccf (specim, keyword) == NO)
		    return

		call imgspar (specim,keyword,mline,SZ_HSTRING)
		istat = sscan (mline)
		    call gargwrd (tempid[1,itemp],SZ_FNAME)
		    call gargd (zvel[itemp])
		    call gargd (czerr[itemp])
		    call gargd (czr[itemp])
		    call gargd (thght[itemp])
		    call gargd (tvw[itemp])
		    call gargd (tempshift[itemp])
		    call gargd (temphcv[itemp])
		    call gargb (tschop[itemp])
		    call gargi (tempfilt[itemp])
		if (itmax == 0) {
		    ltemp = strlen (tempid[1,itemp])
		    lmax = ltemp
		    if (lbest > lmax) lmax = lbest
		    if (strncmp (tempid[1,itemp], bestid, lmax) == 0)
			itmax = itemp
		    }
		}

	    if (debug) {
		call printf ("XCRHEAD: %d templates:")
		    call pargi (ntemp)
		nchar = 20
		do i = 1, ntemp {
		    nchar = nchar + strlen (tempid[1,i]) + 2
		    if (nchar > 80) {
			nchar = strlen (tempid[1,i]) + 2
			call printf ("\n")
			}
		    call printf (" %s")
			call pargstr (tempid[1,i])
		    if (i < ntemp)
			call printf (",")
		    }
		call printf ("\n")
		}
	    call flush (STDOUT)
	    }

end
# Sep 12 1994	Set ntemp to 0 if correlation velocity is not found
# Nov 16 1994	Read and write quality flag
# Dec 14 1994	Read and write filter and chop flag per template
# Dec 19 1994	Read filter and template information only if keyword is found

# Jan  6 1995	Fix quality numbers
# May 30 1995	Read and write barycentric velocity correction
# Jun  6 1995	Debug multispec output
# Jun  7 1995	Read in MSPEC=0 fit information if no MSPEC=1 results
# Jul 19 1995	Move combined velocity to VWHEAD and VRHEAD
# Sep 19 1995	Set measured width when reading computed width
# Sep 19 1995	Use measured width if per template pixel size (TAA) not set
# Sep 22 1995	Find best template number if not given
# Oct  4 1995	Add XCSAO parameter in new standard format

# Feb  5 1998	Fix XCSAO parameter for multispec data

# Aug 19 1999	If reading template info for mspec = 1 or 0, try both keywords
# Aug 24 1999	If CZPARAMS is not found, try NTEMP

# Nov 21 2002	Put linefeeds in list of templates to avoid linewrap

# Aug 24 2004	Fix bug so tschop is read and written as boolean, not integer

# Aug 22 2005	Write Z to header

# May 19 2006	Use imputh() to write out COMMENT as well as HISTORY
# Jul 13 2006	Do not create per-aperture keywords for apertures > 999
