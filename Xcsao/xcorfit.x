# File rvsao/Xcsao/xcorfit.x
# May 9, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After John Tonry and Guillermo Torres

# Copyright(c) 1991-2008 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Perform a digital cross-correlation between two spectra

include "rvsao.h"

procedure xcorfit (spectrum,template,wltemp,ncor,xcor,xvel,itemp,
		   tcz,tczerr,tczr)

real	spectrum[ARB]	# Object spectrum (log wavelength)
real	template[ARB]	# Template spectrum (log wavelength)
real	wltemp[ARB]	# Wavelengths for object and template spectra
int	ncor		# Number of points in cross-correlation
real	xcor[ARB]	# Cross-correlation vector
real	xvel[ARB]	# Velocities for cross-correlation vector
int	itemp		# template sequence number
double	tcz		# Cz for object and this template (returned)
double	tczerr		# Cz error for object and this template (returned)
double	tczr		# R value for this cross-correlation (returned)
double	tcz1		# velocity at center minus peak half-width
#complex	cfn[ARB]

int	log2, i, j, pkindx, ncor2
int	numpts, index
real	xcormin, xcormax
double	arms, srms, fnorm, dindef
double	pc, pw, ph, fracpeak
char	title[SZ_LINE]
double	shift
double	epsilon
double	sigma1, sigma2
double	usigma1, usigma2
double	pwidth, pheight, pkhalf
int	tfilt, npts2
bool	tflag, sflag
char	wlab[SZ_LINE]	# Label for spectrum and template
char	clab[SZ_LINE]	# Lable for cross-correlation results

include "rvsao.com"
include "results.com"
include "xcorf.com"

begin
	dindef = INDEFD
	if (correlate == COR_PIX)
	    call strcpy ("Pixel", wlab, SZ_LINE)
	else
	    call strcpy ("Wavelength in Angstroms", wlab, SZ_LINE)
	if (correlate == COR_PIX)
	    call strcpy ("Shift in Pixels", clab, SZ_LINE)
	else if (correlate == COR_WAV)
	    call strcpy ("Shift in Wavelength (Angstroms)", clab, SZ_LINE)
	else
	    call strcpy ("Shift in Velocity (km/sec)", clab, SZ_LINE)

	sflag = FALSE
	tflag = TRUE
	tfilt = tempfilt[itemp]
	ncor2 = ncor / 2
	c0 = 299792.5d0
	log2 = alog (real (ncor)) / LN2 + 0.001


# Allocate transform vectors
	if (xind == NULL) {
	    npts2 = 2 * npts
	    call malloc (xind, npts2, TY_REAL)
	    call malloc (xifft, npts2, TY_REAL)
	    call malloc (pft, npts2, TY_REAL)
	    call malloc (tft, npts2, TY_COMPLEX)
	    call malloc (ftcfn, npts2, TY_COMPLEX)
	    call malloc (ft1, npts2, TY_COMPLEX)
	    call malloc (ft2, npts2, TY_COMPLEX)
	    call malloc (spexp, npts2, TY_REAL)
	    call malloc (xcont, npts, TY_REAL)
	    }

	if (debug) {
	    call printf ("XCORFIT: filter is %d %d %d %d\n")
		call pargi (lo)
		call pargi (toplo)
		call pargi (topnrn)
		call pargi (nrun)
	    }

# Set up X coordinates for cross-correlation function
	do i = 1, ncor {
	    Memr[xifft+i-1] = i
	    index = i - ncor2 - 1
	    Memr[xind+i-1] = index
	    if (correlate == COR_VEL || correlate == COR_YES)
		xvel[i] = c0 * (10.d0**(index*dlogw) - 1.d0) + tvel
	    else
		xvel[i] = index * delwav
	    }
	if (debug) {
	    if (correlate == COR_WAV)
		call printf ("XCORFIT: wavelengths from %.4f to %.4f, (+ %.4f)\n")
	    else if (correlate == COR_PIX)
		call printf ("XCORFIT: pixels from %.4f to %.4f, (+ %.4f)\n")
	    else
		call printf ("XCORFIT: velocities from %.4f to %.4f, (+ %.4f)\n")
		call pargr (xvel[1])
		call pargr (xvel[ncor])
		call pargd (tvel)
	    }

# Subtract continuum and chop unwanted lines from spectrum using ICFIT
	if (contfit)
	    call icsubcon (npts,spectrum,wltemp,specname,1,nsmooth,Memr[xcont])

# Apodize ends of spectrum
	call apodize (npts,spectrum,wltemp,han,"Spectrum apodized: ",specname)

# Fourier transform spectrum
	call amovr (spectrum,Memr[spexp],npts)
	if (zpad)
	    call aclrr (Memr[spexp+npts],npts]
	call rcvec (ncor,Memr[spexp],Memx[ft1])
	call fourm (Memx[ft1],log2,1)
	if (pltfft) {
	    call strcpy ("Spectrum FFT: ", title, SZ_LINE)
	    call strcat (specname, title, SZ_LINE)
	    do i = 1, ncor2 {
		call xcpower (Memx[ft1+i-1], Memr[pft+i-1]) 
		}
	    call plotspec (ncor2, Memr[pft], title, Memr[xifft],"",nsmooth)
	    }

# Subtract continuum and chop emission lines from template using ICFIT
	if (tcont) {
	    call icsubcon (npts,template,wltemp,tempname[1,itemp],2,nsmooth,
			   Memr[xcont])

# Apodize ends of template
	    call apodize (npts,template,wltemp,han,"Template apodized: ",tempname[1,itemp])
	    }

# Fourier transform template
	call amovr (template,Memr[spexp],npts)
	if (zpad)
	    call aclrr (Memr[spexp+npts],npts)
	call rcvec (ncor,Memr[spexp],Memx[ft2])
	call fourm (Memx[ft2],log2,1)
	if (pltfft) {
	    call strcpy ("Template FFT: ", title, SZ_LINE)
	    call strcat (tempname[1,itemp], title, SZ_LINE)
	    do i = 1, ncor2 {
		call xcpower (Memx[ft2+i-1], Memr[pft+i-1]) 
		}
	    call plotspec (ncor2, Memr[pft], title, Memr[xifft],"",nsmooth)
	    }

# Plot the unfiltered cross-correlation function if requested
	if (pltuc) {

	#  Compute the unfiltered normalization factor and sigmas
	    call rmsnorm (ncor,Memx[ft1],Memx[ft2],debug,usigma1, usigma2,fnorm)

	# Correlate the spectrum transform with the template transform
	    call correl (ncor,Memx[ft1],Memx[ft2],Memx[ftcfn])
	    call fourm (Memx[ftcfn],log2,-1)
	    call crvec (ncor,Memx[ftcfn],xcor)
	    call flip (ncor,xcor)

	# Normalize the unfiltered cross-correlation function if it's plotted
	    do i = 1, ncor {
		xcor[i] = xcor[i] * fnorm
		}
            call strcpy ("Unfiltered cross-correlation: ", title, SZ_LINE)
            call strcat (specname, title, SZ_LINE)
            call strcat (" X ", title, SZ_LINE)
            call strcat (tempname[1,itemp], title, SZ_LINE)
            call plotspec (ncor, xcor, title, xvel,clab,nsmooth)
	    call fourm (Memx[ftcfn],log2,1)
            }

#  Filter the transformed spectrum
	if (lo != nrun)
	    call flter2 (ncor,lo,toplo,topnrn,nrun,Memx[ft1],tfilt)

#  Plot the filtered transformed spectrum, if requested
	if (pltfft) {
	    call strcpy ("Filtered spectrum FFT: ", title, SZ_LINE)
	    call strcat (specname, title, SZ_LINE)
	    do i = 1, ncor2 {
		call xcpower (Memr[ft1+i-1], Memr[pft+i-1])
		}
	    call plotspec (ncor2, Memr[pft], title, Memr[xifft],"",nsmooth)
	    }

#  Filter the transformed template
	if (lo != nrun && tfilt != 1 && tfilt != 3)
	    call flter2 (ncor,lo,toplo,topnrn,nrun,Memx[ft2],tfilt)

#  Plot the filtered transformed template, if requested
	if (pltfft) {
	    call strcpy ("Filtered template FFT: ", title, SZ_LINE)
	    call strcat (tempname[1,itemp], title, SZ_LINE)
	    do i = 1, ncor2 {
		call xcpower (Memr[ft2+i-1], Memr[pft+i-1])
		}
	    call plotspec (ncor2, Memr[pft], title, Memr[xifft],"",nsmooth)
	    }

#  Plot the transform of the filtered transformed spectrum, if requested
	if (plttft) {
	    call strcpy ("Filtered spectrum: ", title, SZ_LINE)
	    call strcat (specname, title, SZ_LINE)
	    call amovx (Memx[ft1],Memx[tft],ncor)
	    call fourm (Memx[tft],log2,-1)
	    call crvec (ncor,Memx[tft],Memr[pft])
	    if (zpad)
		call plotspec (ncor2, Memr[pft], title, wltemp, wlab, nsmooth)
	    else
		call plotspec (ncor, Memr[pft], title, wltemp, wlab, nsmooth)
	    }

#  Plot the transform of the filtered transformed template, if requested
	if (plttft) {
	    call strcpy ("Filtered template: ", title, SZ_LINE)
	    call strcat (tempname[1,itemp], title, SZ_LINE)
	    call amovx (Memx[ft2],Memx[tft],ncor)
	    call fourm (Memx[tft],log2,-1)
	    call crvec (ncor,Memx[tft],Memr[pft])
	    if (zpad)
		call plotspec (ncor2, Memr[pft], title, wltemp, wlab, nsmooth)
	    else
		call plotspec (ncor, Memr[pft], title, wltemp, wlab, nsmooth)
	    }

#  Correlate the spectrum transform with the template transform
	call correl (ncor,Memx[ft1],Memx[ft2],Memx[ftcfn])
#	call amovx (Memx[ftcfn],Memx[cfn],ncor)

	call fourm (Memx[ftcfn],log2,-1)
	call crvec (ncor,Memx[ftcfn],xcor)
	call flip (ncor,xcor)

#  Compute the normalization factor
	call rmsnorm (ncor,Memx[ft1],Memx[ft2],debug,sigma1, sigma2, fnorm)

#  Normalize the cross-correlation function
	do i = 1, ncor {
	    xcor[i] = xcor[i] * fnorm
	    }
	if (debug) {
	    xcormin = xcor[1]
	    xcormax = xcor[1]
	    do i = 1, ncor-1 {
		if (xcor[i] < xcormin) xcormin = xcor[i]
		if (xcor[i] > xcormax) xcormax = xcor[i]
		}
	    call printf ("XCORFIT: cross correlation normalized: %.7f - %.7f\n")
		call pargr (xcormin)
		call pargr (xcormax)
	    call flush (STDOUT)
	    }

# Plot the filtered cross-correlation function if requested
# Note that a peak can be chosen by the cursor "p" command if plotted
	z[1] = 0
        if (pltcor) {
            call strcpy ("Filtered cross-correlation: ", title, SZ_LINE)
            call strcat (specname, title, SZ_LINE)
            call strcat (" X ", title, SZ_LINE)
            call strcat (tempname[1,itemp], title, SZ_LINE)
            call plotspec (ncor, xcor, title, xvel,clab,nsmooth)
            }
	if (z[1] > 0)
	    pkindx = z[1]
	else
	    pkindx = 0
	if (debug) {
	    call printf ("XCORFIT:  ready to fit peak, centered at %d\n")
		call pargi (pkindx)
	    call flush (STDOUT)
	    }

#  Calculate the redshift with a fit to the correlation peak
	if (pkfrac < 0.d0)
	    fracpeak = -pkfrac
	else
	    fracpeak = pkfrac
	if (pkmode0 == 2)
	    call pkfitq (ncor,Memr[xind],xcor,xvel,fracpeak,pkindx,numpts,pc,pw,ph,debug)
	else if (pkmode0 == 3)
	    call pkfitc (ncor,Memr[xind],xcor,xvel,fracpeak,pkindx,numpts,pc,pw,ph,debug)
	else {
	    pkmode0 = 1
	    call pkfitp (ncor,Memr[xind],xcor,xvel,fracpeak,pkindx,numpts,pc,pw,ph,debug)
	    }
	if (debug) {
	    if (correlate == COR_WAV) {
		call printf ("XCORFIT:  pc= %f, delwav=%f\n")
		    call pargd (pc)
		    call pargd (delwav)
		}
	    else if (correlate == COR_PIX) {
		call printf ("XCORFIT:  pc= %f, delpix=%f\n")
		    call pargd (pc)
		    call pargd (delwav)
		}
	    else {
		call printf ("XCORFIT:  pc= %f, dlogw=%f\n")
		    call pargd (pc)
		    call pargd (dlogw)
		}
	    call flush (STDOUT)
	    }

# Compute R and error in cz

	if (correlate == COR_WAV || correlate == COR_PIX)
	    tcz = (pc * delwav)
	else
	    tcz = c0 * (10.d0 ** (pc * dlogw) - 1.d0)
	shift = pc

	if (debug) {
	    call printf ("XCORFIT:  tcz: %9.3f  c: %8.3f  dlogw: %11.7g\n")
		call pargd (tcz)
		call pargd (pc)
		call pargd (dlogw)
	    call printf ("XCORFIT:  npts: %d w: %6.3f h: %5.3f\n")
		call pargi (numpts)
		call pargd (pw)
		call pargd (ph)
	    call flush (STDOUT)
	    }

#  Compute the anti-symmetric error
	call fourm (Memx[ftcfn],log2,1)
	call aspartf (ncor,lo,nrun,shift,Memx[ftcfn],arms,srms, tfilt)
	arms = arms * fnorm
	srms = srms * fnorm
	tcent[itemp] = pc
	thght[itemp] = ph
	twdth[itemp] = pw
	tarms[itemp] = arms
	tsrms[itemp] = srms
	tnpfit[itemp] = numpts
	tsig1[itemp] = sigma1
	tsig2[itemp] = sigma2

	if (arms > 0.) {
	    if (pw > 0.)
		tczr = ph / (2.d0 * arms)
	    else
		tczr = 0.d0
	    pkhalf = 0.5d0
	    if (z[1] > 0)
		pkindx = z[1]
	    else
		pkindx = 0
	    call pkwidth (ncor,xcor,xvel,pkhalf,pkindx,pheight,pwidth,i,j)
	    if (pwidth < 0.d0) pwidth = -pwidth
	    epsilon = 3.d0 * pwidth / 8.d0 / (1.d0 + tczr)
	    if (correlate == COR_WAV || correlate == COR_PIX)
		tczerr = epsilon * delwav
	    else
		tczerr = c0 * (10.d0 ** (epsilon * dlogw) - 1.d0)
	    }
	else {
	    pheight = 0.d0
	    tczr = 0.d0
	    tczerr = 0.d0
	    pwidth = 0.d0
	    }
	tpcent[itemp] = Memr[xind+pkindx-1]
	tphght[itemp] = pheight
	tpwdth[itemp] = 0.5d0 * (Memr[xind+pkindx+(pwidth*0.5d0)-1.d0] -
			Memr[xind+pkindx-(pwidth*0.5d0)-1.d0])

	if ((minvel != dindef && tcz+tvel < minvel) ||
	    (maxvel != dindef && tcz+tvel > maxvel)) {
	    tcz = xvel[pkindx]
	    tczr = 0.d0
	    tczerr = 0.d0
	    tcent[itemp] = Memr[xind+pkindx-1]
	    thght[itemp] = xcor[pkindx]
	    twdth[itemp] = 0.d0
	    }

	taa[itemp] = 1.d0 / dlogw / dlog (10.d0)
	if (debug) {
	    call printf ("XCORFIT:  R: %8.4f  epsilon: %8.6f  error: %8.4f\n")
		call pargd (tczr)
		call pargd (epsilon)
		call pargd (tczerr)
	    call printf ("XCORFIT:  arms: %8.6f  srms: %8.6f\n")
		call pargd (arms)
		call pargd (srms)
	    }

# Compute peak width in Cz for autoscaling of velocity plot
	if (correlate == COR_WAV || correlate == COR_PIX)
	    tcz1 = (pc - (pwidth * 0.5d0)) * delwav
	else
	    tcz1 = c0 * (10.d0 ** ((pc - (pwidth * 0.5d0)) * dlogw) - 1.d0)
	tvw[itemp] = tcz - tcz1
	if (debug) {
	    call printf ("XCORFIT:  vwidth: %8.3f km/sec\n")
		call pargd (tvw[itemp])
	    }
	
	call aclrr (Memr[spexp], npts2)
	call aclrr (Memr[xind], npts2)
	call aclrr (Memr[xifft], npts2)
	call aclrr (Memr[pft], npts2)
	call aclrx (Memx[tft], npts2)
	call aclrx (Memx[ftcfn], npts2)
	call aclrx (Memx[ft1], npts2)
	call aclrx (Memx[ft2], npts2)
	call aclrr (Memr[spexp], npts2)
	call aclrr (Memr[xcont], npts)

	return
end


# Compute power of point in Fourier transform

procedure xcpower (cx, rx)

complex	cx
real	rx

begin
	rx = sqrt ((real(cx) * real (cx)) + (aimag (cx) * aimag (cx)))
	if (rx > 0)
	    rx = alog10(rx)
	else
	    rx = 0.0
	return
end

# Sep 16 1991	Change ccvec call to amovx
# Sep 24 1991	Add results.com and double arms, srms
# Sep 25 1991	Double all peak fitting parameters
# Sep 25 1991	Pass number of points fit
# Oct 25 1991	Pass x axis values and labels to plotspec
# Oct 30 1991	Fix template name in plot headings
# Dec  2 1991	Compute peak width in velocity for plot autoscaling
# Dec 16 1991	Pass wavelength vector to ICSUBCON, EMCHOP, SUBCON, and APODIZE

# Mar 27 1992	Pass velocity vector to PKWIDTH (and change argument order)
# Mar 27 1993	Return velocity vector 
# May 27 1992	Negate pkindx if peak position set in plotspec
# May 28 1992	Do not negate pkindx if peak position set in plotspec
# Jun 16 1992	Drop old polynomial continuum fitting
# Aug 12 1992	Make all constants double precision
# Dec  1 1992	Move debug printout from after to before aspartf call
# Dec  1 1992	Pass debug to PKWIDTH and peak-fitting subroutines

# Jan 20 1994	Fix bug so different peak-fitting modes are really used
# Feb 11 1994	Save sigma1 and sigma2 for special type 6 report
# May  9 1994	Add number of times to smooth spectrum argument to PLOTSPEC call
# Jun 15 1994	Pass per-template filter flag to fltr
# Jun 29 1994	Use stack pointer instead of heap pointer
# Jun 29 1994	Return R=0 if template is all zeroes
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug  4 1994	Pad spectrum with zeroes to avoid wrap-around
# Aug  8 1994	Make zero-padding optional
# Aug 17 1994	Get ZPAD from labelled common
# Dec 15 1994	Add filter flag to RMSFIL, FLTER, and ASPARTF
# Dec 16 1994	Get number of points in RMS from RMSFIL
# Dec 19 1994	Filter transforms before correlating them

# Feb  9 1995	Return unfiltered sigmas, not filtered ones
# Feb 10 1995	Use square root of filter for equivalence with old program
# Mar 13 1995	Allow option of not removing continuum from template
# Mar 15 1995	Pass object and template names to ICSUBCON
# Mar 22 1995	Set FRACPEAK from RVSAO common instead of reading par file here
# Jul  3 1995	If peak width from fit is zero, set r-value to zero
# Jul  3 1995	Save measured, as well as fit, peak height
# Aug 18 1995	If fit velocity outside limits, use maximum peak within limits
# Aug 21 1995	Drop DEBUG from PKWIDTH call
# Aug 22 1995	Use absolute value of PKFRAC

# Jan 22 1997	Plot POWER of spectrum and template FFTs
# Jan 24 1997	Drop phase plot
# Jan 28 1997	Plot transformed transform to see what is being correlated
# Mar 12 1997	Drop delcarations for variables in common
# Apr 14 1997	Label unfiltered cross-correlation more explicitly
# May  6 1997	Add NSMOOTH argument to ICSUBCON
# Dec 22 1997	Label wavelength units in graphs

# May 15 1998	Allow INDEFD maxvel and minvel

# Sep 19 2000	Add options to cross-correlate in wavelength and pixel space

# Jul 31 2001	Add delpix to return cross-correlation shifts in original pixels

# Jun 20 2007	Put transform vector pointers in xcorf.com
# Jun 20 2007	Allocate all transform vectors instead of fixing dimension
# Jun 25 2007	Allocate transform vectors with 2*npts instead of npts pixels
# Jun 25 2007	Allocate xcont (formerly work) and spexp only once per task run

# Mar 10 2008	Pixels per pixel is passed in delwav if correlating in pixels
# May  9 2008	Removed object spectrum continuum only if contfit is true
