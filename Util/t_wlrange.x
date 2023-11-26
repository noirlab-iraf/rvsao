# File rvsao/Sumspec/t_wlrange.x
# June 15, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# WLRANGE is an IRAF task which discovers the overlap within a list of spectra
# For arguments, see parameter file wlrange.par.
# Information is shared in common blocks defined in "rvsao.com".
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include "rvsao.h"

procedure t_wlrange ()

char	spectra[SZ_PATHNAME]	# List of input spectra
int	npix
int	indefi
double	minwav, maxwav, pixwav
double	dindef
bool	verbose
 
bool	clgetb()
int	clgeti()
double	clgetd()

include	"rvsao.com"
include	"sum.com"
 
begin
	c0 = 299792.5
	dindef = INDEFD
	indefi = INDEFI
	npix = 0
	pixwav = dindef

# Get task parameters.

# Spectra to combine
	call clgstr ("spectra",spectra,SZ_PATHNAME)

# Print processing information
	debug  = clgetb ("debug")
	verbose = clgetb ("verbose")

# Optional wavelength limits for output spectrum
	minwav0 = clgetd ("st_lambda")
	maxwav0 = clgetd ("end_lambda")
	pixwav = clgetd ("pix_lambda")
	npts = clgeti ("npts")

# Get wavelength limits from input spectra
	call wlrange (spectra, minwav, maxwav, npix)
	if (minwav0 != dindef)
	    minwav = minwav0
	if (maxwav0 != dindef)
	    maxwav = maxwav0
	if (npts != indefi)
	    npix = npts
	if (pixwav == dindef || pixwav == 0.d0)
	    pixwav = (maxwav - minwav) / double (npix - 1)
	else if (npts == indefi)
	    npix = 1 + nint ((maxwav - minwav) / pixwav)

# Save wavelength overlap as parameters
	call clputd ("wl1", minwav)
	call clputd ("wl2", maxwav)
	call clputd ("dwl", pixwav)
	call clputi ("npix", npix)

# Optionally print the wavelength overlap
	if (verbose) {
	    call printf ("WLRANGE: %d-point spectra from %.3fA to %.3fA by %.3fA\n")
		call pargi (npix)
		call pargd (minwav)
		call pargd (maxwav)
		call pargd (pixwav)
	    call flush (STDOUT)
	    }

end

# Jun 15 2009	New program based on part of SUMSPEC
