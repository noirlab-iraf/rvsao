# File rvsao/Makespec/wtgauss.x
# March 26, 1997
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

procedure wtgaus (spectrum, width, wlmin, wlpix, npix, debug)

real	spectrum[ARB]	# Spectrum, linear in wavelength
double	width		# Half-width of Gaussian in Angstroms
double	wlmin		# Wavelength at first pixel in Angstroms
double	wlpix		# Wavelength per pixel in Angstroms
int	npix		# Number of pixels in spectrum
bool	debug		# yes to print everything

double	wlmax		# Wavelength of last pixel
double	wl, ds2, dw, dw2, dg, center, wtj, twt, spix
int	ipix, np, jp, jp1, jp2, jpix
pointer	sp, newspec, wt

begin

	call smark (sp)
	call salloc (newspec, npix, TY_REAL)
	wlmax = wlmin + wlpix * double (npix - 1)
	np = 5 * (width / wlpix)
	call salloc (wt, (np*2)+1, TY_DOUBLE)

# Compute weighting factors for this Gaussian
	center = 0.d0
	ds2 = 2.d0 * width * width
	if (debug) {
	    call printf ("WTGAUS: width= %f, npoints= %d\n")
		call pargd (width)
		call pargi (np)
	    }
	wl = center - (double (np) * wlpix)
	twt = 0.d0
	do jp = -np, np {
	    dw = wl - center
	    dw2 = dw * dw
	    dg = 1.d0 * dexp (-dw2 / ds2)
	    Memd[wt + np - jp] = dg
	    twt = twt + dg
	    if (debug) {
		call printf ("WTGAUS: %d: wl= %7.3f, wt= %.5f total= %.5f\n")
		    call pargi (jp)
		    call pargd (wl)
		    call pargd (dg)
		    call pargd (twt)
		}
	    wl = wl + wlpix
	    }

# Normalize weights so they total 1.0
	do jp = -np, np {
	    Memd[wt + np - jp] = Memd[wt + np - jp] / twt
	    }
	    
# Intitialize spectrum to all zeroes
        call aclrr (Memr[newspec], npix)

# Run through spectrum
	do ipix = 1, npix {
	    twt = 0.d0
	    spix = 0.d0

#  Compute contributions from pixels in Gaussian spread
	    Jp1 = -np
	    if (ipix + jp1 < 1)
		jp1 = 1 - ipix
	    jp2 = np
	    if (ipix + jp2 > npix)
		jp2 = npix - ipix
	    do jp = jp1, jp2 {
		jpix = ipix + jp
		if (spectrum[jpix] != 0.0) {
		    wtj = Memd[wt + np - jp]
		    spix = spix + (double (spectrum[jpix]) * wtj)
		    twt = twt + wtj
		    }
		}

#  Set spectrum pixel to normalized sum of Gaussian elements
	    Memr[newspec+ipix-1] = real (spix)
	    if (debug && spix != 0.d0) {
		wl = wlmin + (wlpix * double (ipix - 1))
		call printf ("WTGAUS: %d: %9.4f: %.4f -> %.4f / %.3f = %.4f\n")
		    call pargi (ipix)
		    call pargd (wl)
		    call pargr (spectrum[ipix])
		    call pargd (spix)
		    call pargd (twt)
		    call pargr (Memr[newspec+ipix-1])
		}
	    }
	if (debug)
	    call flush (STDOUT)

# Replace original spectrum with filtered spectrum
	call amovr (Memr[newspec], spectrum, npix)

	call sfree (sp)
	return
end

# Mar 26 1997	New subroutine
