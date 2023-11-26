#*** File Util/gsmooth.x
#*** February 26, 1997
#*** By Doug Mink, Harvard-Smithsonian Center for Astrophysics

#--- Smooth one-dimensional data

procedure gsmooth (newdata, olddata, npix, mode, ndsamp, ghalf)

real	newdata[ARB]		# Smoothed data (returned)
real	olddata[ARB]		# Input data
int	npix			# Number of pixels to smooth
int	mode			# 0=none, 1=1-2-1, 2=Gaussian, 3=Hanning
int	ndsamp			# No. of spectrum points to smooth
double	ghalf			# Half width of Gaussian in pixels

int	ndsamp2			# Half of ndsamp
double	wt[100]			# Gaussian weighting function per pixel
double	wtpix			# Gaussian weighting function total

double	twopi, sig2, rad2, sumdata, dsamp
real	spj,data0,data1,data2
int	ipixn
int	isamp,i,ipix,jwt,jpix,jpix0,jpix1

begin

# 1-2-1 smoothing (SAO traditional method)
	if (mode == 1) {
	    do ipix = 1, npix {
		newdata[ipix] = olddata[ipix]
		}
	    do I =1, ndsamp {
		data0 = newdata[1]
		data1 = newdata[1]
		do ipix = 1, npix {
		    ipixn = ipix + 1
		    if (ipixn > npix) ipixn = npix
		    data2 = newdata[ipixn]
		    newdata[ipix] = (data0 + data1 + data1 + data2) / 4.0
		    data0 = data1
		    data1 = data2
		    }
		}
	    }

# Gaussian or Hanning smoothing
	else if (mode == 2 || mode == 3) {
	    ndsamp2 = ndsamp / 2
	    twopi = 6.2831853
	    sig2 = 2.d0 * ghalf * ghalf

#	Set up weighting functions for Gaussian and Hanning smoothing
	    ndsamp = 1 + (ndsamp2 * 2)
	    do isamp = -ndsamp2, ndsamp2 {
		i = isamp + ndsamp2 + 1
		dsamp = double (isamp)
		rad2 = dsamp * dsamp
		if (mode == 2)
		    wt[i] = dexp (-rad2 / sig2)
		else
		    wt[i] = 0.5d0*(1.d0-dcos(twopi*double(i)/double(ndsamp-1)))
		}

	    do ipix = 1, npix {
		jpix0 = ipix - ndsamp2
		jpix1 = ipix + ndsamp2
		sumdata = 0.d0
		wtpix = 0.d0
		jwt =0
		do jpix = jpix0, jpix1 {
		    jwt = jwt + 1
		    if ((jpix >= 1) && (jpix <= npix)) {
			sumdata = sumdata + (wt[jwt] * olddata[jpix])
			wtpix = wtpix + wt[jwt]
			}
		    }
		if (wtpix > 0)
		    spj = real (sumdata / wtpix)
		else
		    spj = 0.
		newdata[ipix] = spj
		}
	    }

# Otherwise return unsmoothed data
	else {
	    do ipix = 1, npix {
		newdata[ipix] = olddata[ipix]
		}
	    }

	return
end

# Feb 26 1997	New subroutine
