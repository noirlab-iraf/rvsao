# File rvsao/Emsao/emhead.x
# July 13, 2006
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1994-2006 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include <smw.h>
include "rvsao.h"
include "emv.h"

# EMWHEAD writes emission line velocity information into an IRAF image header

procedure emwhead (mspec, specim)

int	mspec		# Number of spectrum to read from multispec file
pointer	specim		# Object image header structure

char	dstr[SZ_LINE]	# Date string
char	hstr[SZ_LINE]	# New history string
char	hstr0[SZ_LINE]	# Old history string
char	keyword[SZ_HKWORD]
char	line[SZ_HSTRING]
int	iline
double	dwdp, wl1, width, twt, wt
double	wcs_p2w()
 
include	"rvsao.com"
include	"emv.com"
 
begin

# Save Cz and error to image header, if requested
	if (mspec == 0) {
	    call imaddd (specim,"CZEM    ",spevel, "Emission line redshift velocity (km/sec)")
	    call imaddd (specim,"ZEM    ",spevel/c0, "Emission line z = dwl/wl")
	    call imaddd (specim,"CZEMERR ",speerr, "Emission line redshift velocity error (km/sec)")
	    call imaddi (specim,"CZEMNL ",nfound, "Number of emission lines found")
	    call imaddi (specim,"CZEMNLF ",nfit, "Number of emission lines used")
	    }
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVEM%d")
		call pargi (mspec)
	    call sprintf (line,SZ_HSTRING,"%.3f %.3f %d %d")
		call pargd (spevel)
		call pargd (speerr)
		call pargi (nfound)
		call pargi (nfit)
	    call imastr (specim,keyword,line)
	    }

# Print emission line information if any were found
	if (spnlf > 0 && mspec < 1000) {
	    twt = 0.d0
	    do iline = 1, nfound {
		wt = emparams[10,1,iline]
		if (wt != 0.d0)
		    twt = twt + 1.d0 / (wt * wt)
		}
	    do iline = 1, nfound {
		if (mspec == 0) {
		    call sprintf (keyword,SZ_HKWORD,"EMLINE%d")
			call pargi (iline)
		    }
		else {
		    call sprintf (keyword,SZ_HKWORD,"EL%d_%d")
			call pargi (mspec)
			call pargi (iline)
		    }
		wt = emparams[10,1,iline]
		if (wt != 0.d0)
		    wt = (1.d0 / (wt * wt)) / twt
		call sprintf(line,67,"%s %.2f %.2f %6f %6f %.2f %.2f %.2f %.3f")
		    call pargstr (nmobs[1,iline])	# Name
		    call pargd (wlrest[iline])		# Rest wavelength
		    call pargd (wlobs[iline])		# Observed wavelength
		    wl1 = wcs_p2w (emparams[4,1,iline]-1.d0)
		    dwdp = wlobs[iline] - wl1
		    call pargd (emparams[5,1,iline])	# Observed height
		    width = emparams[6,1,iline] * dwdp
		    call pargd (width)			# Observed width
		    call pargd (c0*emparams[9,1,iline])	# CZ (km/sec)
		    call pargd (c0*emparams[9,2,iline])	# CZ error (km/sec)
		    call pargd (emparams[7,1,iline])	# Equivalent width
		if (override[iline] < 0)
		    call pargd (-wt)			# Weighting factor
		else
		    call pargd (wt)			# Weighting factor
		call imastr (specim,keyword,line)
		}
	    }

#  Put date and program name into HISTORY
	call logtime (dstr, SZ_LINE)
	if (mspec > 0 && nmspec == 1) {
	    call sprintf (hstr,SZ_LINE,"rvsao.emsao %s %s CZEM [%s]")
		call pargstr (VERSION)
		call pargstr (dstr)
		call pargstr (specnums)
	    }
	else {
	    call sprintf (hstr,SZ_LINE,"rvsao.emsao %s %s CZEM = %.2f LINES = %d")
		call pargstr (VERSION)
		call pargstr (dstr)
		call pargd (spevel)
		call pargi (nfit)
	    }
	hstr0[1] = EOS
	call imgspar (specim,"EMSAO",hstr0,SZ_HSTRING)
	if (hstr0[1] != EOS)
	    call imputh (specim, "HISTORY", hstr0)
	call imastr (specim, "EMSAO", hstr)
	IM_UPDATE(specim) = YES
	return

end


# EMRHEAD reads emission line velocity information from an IRAF image header

procedure emrhead (mspec0, specim)

int	mspec0		# Number of spectrum to read from multispec file
pointer	specim		# Object image header structure

int	mspec		# Number of spectrum read from multispec file
char	keyword[SZ_HKWORD]
char	line[SZ_HSTRING]
int	iline, istat, sscan()
double	cz, czerr, wl1, dwdp, verr, width, wt
double	wcs_w2p(), wcs_p2w()
int	imaccf()
 
include	"rvsao.com"
include	"emv.com"
 
begin
	mspec = mspec0

# Initialize parameters
	spevel = INDEFD
	speerr = INDEFD
	spnl = 0
	spnlf = 0
	nfound = 0
	nfit = 0

# Return if emission line velocity keyword not found in header
	if (mspec == 0)
	    call strcpy ("CZEM    ",keyword, SZ_HKWORD)
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVEM%d")
		call pargi (mspec)
	    }
	if (imaccf (specim, keyword) == NO)
	    return

# Read emission line Cz and error from image header
	if (mspec == 0) {
	    call imgdpar (specim,"CZEM    ",spevel)
	    call imgdpar (specim,"CZEMERR ",speerr)
	    call imgipar (specim,"CZEMNL ",spnl)
	    call imgipar (specim,"CZEMNLF ",spnlf)
	    }
	else if (mspec < 1000) {
	    call sprintf (keyword,SZ_HKWORD,"APVEM%d")
		call pargi (mspec)
	    call imgspar (specim,keyword,line,SZ_LINE)
	    istat = sscan (line)
		call gargd (spevel)
		call gargd (speerr)
		call gargi (spnl)
		call gargi (spnlf)
	    }
	nfound = spnl
	nfit = spnlf

# Read emission line information if any were found
	if (spnlf > 0 && mspec < 999) {
	    do iline = 1, nfound {
		if (mspec == 0) {
		    call sprintf (keyword,SZ_HKWORD,"EMLINE%d")
			call pargi (iline)
		    }
		else {
		    call sprintf (keyword,SZ_HKWORD,"EL%d_%d")
			call pargi (mspec)
			call pargi (iline)
		    }

		if (imaccf (specim, keyword) == NO)
		    return
		call imgspar (specim,keyword,line,SZ_LINE)
		istat =  sscan (line)
		    call gargwrd (nmobs[1,iline],SZ_ELINE)	# Name

		    call gargd (wlrest[iline])		# Rest wavelength

		    call gargd (wlobs[iline])		# Observed wavelength
		    emparams[4,1,iline] = wcs_w2p (wlobs[iline])
		    wl1 = wcs_p2w (emparams[4,1,iline]-1.d0)
		    dwdp = wlobs[iline] - wl1

		    call gargd (emparams[5,1,iline])	# Observed height

		    call gargd (width)			# Observed width
		    if (dwdp != 0.d0)
			emparams[6,1,iline] = width / dwdp
		    else
			emparams[6,1,iline] = 0.d0

		    call gargd (cz)			# CZ (km/sec)
		    emparams[9,1,iline] = cz / c0

		    call gargd (czerr)			# CZ error (km/sec)
		    verr = czerr / c0
		    emparams[9,2,iline] = verr
		    if (dwdp != 0.d0)
			emparams[4,2,iline] = verr * wlobs[iline] / (2.d0*dwdp)
		    else
			emparams[4,2,iline] = 0.d0

		    call gargd (emparams[7,1,iline])	# Equivalent width

		    call gargd (wt)			# Weighting factor

		if (wt > 0) {
		    override[iline] = 0
		    emparams[10,1,iline] = 1.d0 / sqrt (wt)
		    }
		else if (wt < 0) {
		    override[iline] = -1
		    emparams[10,1,iline] = 1.d0 / sqrt (-wt)
		    }
		else {
		    override[iline] = 0
		    emparams[10,1,iline] = 0.d0
		    }
		wtobs[iline] = wt
		emparams[8,2,iline] = 1.d0
		}
	    }
	if (debug) call printf ("EMRHEAD: Emission line info read\n")

	return

end

# Aug 22 1994	New subroutine
# Nov 16 1994	Read and write quality flag

# Jan  6 1995	Fix quality flag numbers
# May 26 1995	Read and write barycentric velocity correction
# Jul  7 1995	Fix bug reading barycentric velocity correction from file
# Jun  7 1995   Read in MSPEC=0 fit information if no MSPEC=1 results
# Jul 19 1995	Move combined velocity parameters to VWHEAD and VRHEAD
# Sep 27 1995	Always check for divide by zero
# Oct  5 1995	Add EMSAO parameter in new standard format

# Apr  4 1996	Initialize NFOUND and NFIT to zero
# Aug  7 1996	Use smw.h

# Aug 22 2005	Print emission line Z in header
# Jul 13 2006	Do not use per-aperture keywords for apertures > 999
