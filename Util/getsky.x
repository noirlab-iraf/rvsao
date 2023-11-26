# File rvsao/Util/getsky.x
# August 1, 2003
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1994-2003 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
#  GETSKY opens the image specified by specfile and returns pointers
#  to the sky data and the image descriptor.
 
#  Notes:	Shares data in "rvsao.com"
 
###########################################################################
 
include <imhdr.h>
include <smw.h>
include	"rvsao.h"
 
procedure getsky (specfile, sspec, sband, spectrum, specim, skysh)
 
char	specfile[ARB]	# Data file name
int	sspec		# Number of spectrum to read from multispec file
int	sband		# Band of spectrum to read from multispec file
pointer	spectrum	# Spectrum data (returned)
pointer	specim		# Image header structure (returned)
pointer	skysh		# Spectrum header structure (returned)
 
int	npix
int	nline
 
include	"rvsao.com"
 
begin
	nline = sspec
	if (debug) {
	    call printf ("GETSKY: %s ap %d band %d\n")
		call pargstr (specfile)
		call pargi (nline)
		call pargi (sband)
	    }
	call getimage (specfile,mext,nline,sband,spectrum,specim,skysh,npix,
		       skyname,READ_ONLY)
	if (specim == ERR) {
	    call printf ("GETSKY:  Error reading %s[%d,%d]\n")
		call pargstr (specfile)
		call pargi (nline)
		call pargi (sband)
	    return
	    }
 
# redefine common variables
	specpix = npix
	specdc = DC(skysh)

end

# Apr 22 1992	Modify to handle multispec spectra

# Jun  4 1993	Modify to handle MWCS spectra
# Jul  7 1993	Add spectrum header to getimage arguments

# Apr 12 1994	Keep sky name separate from spectrum name; fix debug code
# Apr 12 1994	Return MWCS header pointer
# Jun 15 1994	Set SPECDC from DC(SKYSH), not DCFLAG
# Jun 23 1994	Keep mWCS pointer in SHDR structure
# Aug  3 1994	Change common and header from fquot to rvsao

# Jul 13 1995	Add DEBUG to COMPBCV calls

# Aug  7 1996	Use smw instead of shdr

# Aug 27 1997	Add argument for multispec band

# Apr 22 1998	Drop use of getim.com; extracted needed parameters from  header locally
# Apr 22 1998	Do not set velocity corrections in this subroutine

# Aug  1 2003	Add image extension argument to getimage()

# Aug 22 2013	Drop local variable mext
