# File rvsao/Makespec/mkgauss.x
# March 26, 1997
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

procedure mkgaus (spectrum, center, width, height, wlmin, wlpix, npix, debug)

real	spectrum[ARB]	# Spectrum, linear in wavelength
double	center		# Center wavelength for Gaussian in Angstroms
double	width		# Half-width of Gaussian in Angstroms
double	height		# Height of emission line or (depth of absorption line)
double	wlmin		# Wavelength at first pixel in Angstroms
double	wlpix		# Wavelength per pixel in Angstroms
int	npix		# Number of pixels in spectrum
bool	debug		# yes to print values

double	wlmax		# Wavelength of last pixel
double	wl, wl1, wl2, ds2, dw, dw2, dg
int	ipix, ipix1, ipix2

begin

	wlmax = wlmin + wlpix * double (npix - 1)

# Set pixel limits for this computation
	wl1 = center - (width * 10.d0)
	if (wl1 < wlmin)
	    ipix1 = 1
	else if (wl1 > wlmax)
	    ipix1 = 0
	else
	    ipix1 = 1 + Idnint ((wl1 - wlmin) / wlpix)
	wl2 = center + (width * 10.d0)
	if (wl2 < wlmin)
	    ipix2 = 0
	else if (wl2 > wlmax)
	    ipix2 = npix
	else
	    ipix2 = 1 + Idnint ((wl2 - wlmin) / wlpix)

# Compute Gaussian
	ds2 = 2.d0 * width * width
	if (debug) {
	    call printf ("MKGAUS: center= %.3f, height= %5.3f, width= %.3f\n")
		call pargd (center)
		call pargd (height)
		call pargd (width)
	    }
	do ipix = ipix1, ipix2 {
	    wl = wlmin + wlpix * double (ipix - 1)
	    dw = wl - center
	    dw2 = dw * dw
	    dg = height * dexp (-dw2 / ds2)
	    spectrum[ipix] = spectrum[ipix] + real (dg)
	    if (debug) {
		call printf ("pixel %5d, wl %.3f, dwl %.3f: %.5f -> %.5f\n")
		    call pargi (ipix)
		    call pargd (wl)
		    call pargd (dw)
		    call pargd (dg)
		    call pargr (spectrum[ipix])
		}
	    }
	return
end

# Mar 26 1997	New subroutine
