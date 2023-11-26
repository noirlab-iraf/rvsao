# File rvsao/Util/wlcon.x
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# August 12, 2008

# Subroutines to convert between wavelength and pixel

# WCS_SET	Set spectrum header for succeeding conversions
# WCS_PIXSHIFT	Set spectrum pixel shift for succeeding conversions
# WCS_W2P	Convert wavelength to pixel
# WCS_L2P	Convert log wavelength to pixel
# WCS_P2W	Convert pixel to wavelength

include <smw.h>

procedure wcs_set (sh)

pointer	sh	# Spectrum header structure

double	shdr_lw()
int	sn

include "rvmwcs.com"

begin
	wsh = sh
	sn = SN(wsh)
	wclog = DC(wsh)
	px1 = double (NP1(wsh)) - 0.5d0
	px2 = double (NP2(sh)) + 0.5d0
	wl1 = shdr_lw (sh, px1)
	wl2 = shdr_lw (sh, px2)
	pxshift = 0.0
	if (wl1 > wl2) {
	    wl2 = shdr_lw (sh, px1)
	    wl1 = shdr_lw (sh, px2)
	    }
end

procedure wcs_pixshift (shift)

double	shift	# Spectrum shift in pixels

include "rvmwcs.com"

begin
	pxshift = shift
	return
end

double procedure wcs_getshift ()

include "rvmwcs.com"

begin
	return (pxshift)
end


# Convert spectrum wavelength to pixel

double procedure wcs_w2p (wl)

double	wl	# wavelength in angstroms

double	px	# pixel value (returned)
double  twl,shdr_wl()

include "rvmwcs.com"

begin
	
	if (wl > wl2)
	    twl = wl2
	else if (wl < wl1)
	    twl = wl1
	else
	    twl = wl
	px = shdr_wl (wsh, wl) + pxshift
	return (px)
end


# Convert spectrum log wavelength to pixel

double procedure wcs_l2p (wl)

double	wl	# log10 wavelength in log angstroms

double	px	# pixel value (returned)
double	twl
double  shdr_wl()

include "rvmwcs.com"

begin
	twl = 10.d0 ** wl
	if (twl > wl2)
	    twl = wl2
	else if (twl < wl1)
	    twl = wl1
	px = shdr_wl (wsh, twl) + pxshift
	return (px)
end


# Convert spectrum pixel to wavelength

double procedure wcs_p2w (px)

double	px	# pixel value

double	wl	# wavelength in angstroms (returned)
double	tpx
double  shdr_lw()

include "rvmwcs.com"

begin
	if (px > px2)
	    tpx = px2
	else if (px < px1)
	    tpx = px1
	else
	    tpx = px
	tpx = tpx - pxshift
	wl = shdr_lw (wsh,tpx)
	return (wl)
end


# Convert spectrum pixel to log wavelength

double procedure wcs_p2l (px)

double	px	# pixel value

double	wl	# wavelength in angstroms (returned)
double	tpx
double  shdr_lw()

include "rvmwcs.com"

begin
	if (px > px2)
	    tpx = px2
	else if (px < px1)
	    tpx = px1
	else
	    tpx = px
	tpx = tpx - pxshift
	wl = dlog10 (shdr_lw (wsh,tpx))
	return (wl)
end

# Apr 26 1994	Use tpx in wcsp2w, not px
# Sep  9 1994	Fix log conversions

# Oct  6 1995	Change SHDR_* calls with SPHDR_* calls

# Aug  7 1996	Use smw.h
# Aug  7 1996	Add pixel to log wavelength WCS_P2L

# Jun 11 1998	Set limiting pixels from spectrum header WCS

# Mar 22 2001	Set wavelength limits correctly if spectrum reversed

# May 23 2005	Add pixel shift, wcs_pixshift() to set it
# May 25 2005	Add wcs_getshift() to retrieve pixel shift

# Aug 12 2008	Set limits to start (0.5) of first and end (n+0.5) of last pixel
