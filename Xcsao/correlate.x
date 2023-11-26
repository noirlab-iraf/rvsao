# File rvsao/Xcor/correlate.x
# January 30, 2007
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After John Tonry and Guillermo Torres

# Copyright(c) 1993-2007 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# Various subroutines to perform correlations and other things

procedure correl (n,ft1,ft2,cfn)

int	n
complex	ft1[ARB]
complex ft2[ARB]
complex cfn[ARB]

int	j

begin
	do j = 1, n {
	    cfn[j] = ft1[j] * conjg (ft2[j])
	    }

	return
end


procedure rmsfil (n, k1,k2,k3,k4, x, tflag, tfilt, debug, rms, sqnc)

# return the rms of a real function passed through a bandpass filter in
# fourier space. the bandpass is the square root of a cosine bell.

int	n		# Number of point in transform
int	k1,k2,k3,k4	# Filter limits
complex	x[ARB]		# Transform for which to compute RMS
bool	tflag		# If yes, template, else object
int	tfilt		# If =1, template has already been filtered
			# If =2, turn off high pass filter
			# If =3, turn off high pass filter on object
			#        and do not filter template
bool	debug		# True for debugging information
double	rms		# RMS (returned)
double	sqnc		# Square root of number of points in sigma

double pi,arg,factor,f
int	n1,n2,n3,n4	# Filter limits used to compute RMS
int	k, nc

begin

# Turn off filter on template altogether, if requested
	if ((tfilt == 1 || tfilt == 3) && tflag) {
	    n1 = 0
	    n2 = 0
	    n3 = n - 1
	    n4 = n - 1
	    }

# Turn off high end filtering, if requested
	else if (tfilt > 1) {
	    n1 = k1
	    n2 = k2
	    n3 = n
	    n4 = n
	    }

# Otherwise, keep filter
	else {
	    n1 = k1
	    n2 = k2
	    n3 = k3
	    n4 = k4
	    }
	nc = n4 - n1 + 1
	if (nc > 0)
	    sqnc = sqrt (double (nc))
	else
	    sqnc = 1.d0
	if (debug) {
	    call printf ("RMSFIL: filter is %d %d %d %d, nc = %f\n")
		call pargi (n1)
		call pargi (n2)
		call pargi (n3)
		call pargi (n4)
		call pargd (sqnc)
	    }

	pi = 3.1415926535897932d0

	rms = 0.d0
	do k = n1, n4 {
	    if (k == 0 || (k == n/2))
	        f = 1.d0
	    else
	        f = 2.d0

	    if (k < n2) {
	    	arg = pi * double (k-n1) / double (n2-n1)
	    	factor = .5d0 * (1.d0 - cos (arg))
		}
	    else if (k > n3) {
	    	arg = pi * double (k-n3) / double (n4-n3)
	    	factor = .5d0 * (1.d0 + cos (arg))
		}
	    else
	    	factor = 1.d0

	    rms = rms + f * factor *
		  (real(x[k+1])*real(x[k+1]) + aimag(x[k+1])*aimag(x[k+1]))
	    }

#  because it is the transform

	rms = sqrt (rms) / double (n)
	if (debug) {
	    call printf ("RMSFIL: rms is %f\n")
		call pargd (rms)
	    }

	return
end


procedure rmsnorm (n, x1, x2, debug, rms1, rms2, fnorm)

# return the correlation normalization based on two Fourier transforms
# which have already been filtered

int	n		# Number of points in each transform
complex	x1[ARB]		# Spectrum transform
complex	x2[ARB]		# Template transform
bool	debug		# True for debugging information
double	fnorm		# Normalization factor (returned)
double	rms1, rms2	# RMS for spectrum and template

double	f, dn, xr, xi, rms1n, rms2n
int	k, nc1

begin

	nc1 = n - 1
	dn = double (n)

	rms1 = 0.d0
	do k = 0, nc1 {
	    if (k == 0 || (k == n/2))
	        f = 1.d0
	    else
	        f = 2.d0
	    xr = real (x1[k+1])
	    xi = aimag (x1[k+1])
	    rms1 = rms1 + f * ((xr * xr) + (xi * xi))
	    }

#  because it is the transform
	rms1n = sqrt (rms1)
	rms1 = rms1n / dn

	rms2 = 0.d0
	do k = 0, nc1 {
	    if (k == 0 || (k == n/2))
	        f = 1.d0
	    else
	        f = 2.d0
	    xr = real (x2[k+1])
	    xi = aimag (x2[k+1])
	    rms2 = rms2 + f * ((xr * xr) + (xi * xi))
	    }

#  because it is the transform
	rms2n = sqrt (rms2)
	rms2 = rms2n / dn

#  Compute the normalization factor
	if (rms1 != 0. && rms2 != 0.)
	    fnorm = 2.d0 / (rms1 * rms2 * dn)
	else if (rms1 != 0.)
	    fnorm = 2.d0 / rms1n
	else if (rms2 != 0.)
	    fnorm = 2.d0 / rms2n
	else
	    fnorm = 2.d0 / dn

	if (debug) {
	    call printf ("RMSNORM: rms = %10g, %10g, n = %d, fnorm = %10g\n")
		call pargd (rms1)
		call pargd (rms2)
		call pargi (n)
		call pargd (fnorm)
	    }

	return
end


# Filter a Fourier tranform by multiplying by a cosine bell

procedure flter (n, k1,k2,k3,k4, x, tfilt)

int	n,k1,k2,k3,k4
complex	x[ARB]
int	tfilt		# If =1, do not filter template
			# If =2, turn off high pass filter
			# If =3, turn off high pass filter on object
			#        and do not filter template

int	j, numa ,number, n1, n2, n3, n4
double	arg, factor, pi

begin

	pi = 3.1415926535897932d0

# Turn off high end filtering, if requested
	if (tfilt > 1) {
	    n1 = k1
	    n2 = k2
	    n3 = n
	    n4 = n
	    }

# Otherwise, keep filter
	else {
	    n1 = k1
	    n2 = k2
	    n3 = k3
	    n4 = k4
	    }

	do j = 1, n {
	    if (j < n/2+1)
	    	number = j - 1
	    else
	    	number = j - n - 1

	    numa = abs (number)
	    if ((numa < n1) || (numa > n4))
		x[j] = 0.d0
	    else if (numa <  n2) {
	    	arg = pi * double (numa-n1) / double (n2-n1)
	    	factor = .5d0 * (1.d0 - cos (arg))
		x[j] = x[j] * factor
		}
	    else if (numa > n3) {
	    	arg = pi * double (numa-n3) / double (n4-n3)
	    	factor = .5d0 * (1.d0 + cos (arg))
		x[j] = x[j] * factor
		}
	    }

	return
end


procedure flter2 (n, k1,k2,k3,k4, x, tfilt)

# filter a fourier tranform by multiplying by square-root cosine bell

int	n,k1,k2,k3,k4
complex	x[ARB]
int	tfilt		# If =1, do not filter template
			# If =2, turn off high pass filter
			# If =3, turn off high pass filter on object
			#        and do not filter template

int	j, numa ,number
int	n1,n2,n3,n4	# Filter limits used
double	arg, factor, pi

begin

	pi = 3.1415926535897932d0

# Turn off high end filtering, if requested
	if (tfilt > 1) {
	    n1 = k1
	    n2 = k2
	    n3 = n
	    n4 = n
	    }

# Otherwise, keep filter
	else {
	    n1 = k1
	    n2 = k2
	    n3 = k3
	    n4 = k4
	    }

	do j = 1, n {
	    if (j < n/2+1)
	    	number = j - 1
	    else
	    	number = j - n - 1

	    numa = abs (number)
	    if ((numa < n1) || (numa > n4))
	    	x[j] = 0.d0
	    else if (numa < n2) {
	    	arg = pi * double (numa-n1) / double (n2-n1)
	    	factor = .5d0 * (1.d0 - cos (arg))
		x[j] = x[j] * sqrt (factor)
		}
	    else if (numa > n3) {
	    	arg = pi * double (numa-n3) / double (n4-n3)
	    	factor = .5d0 * (1.d0 + cos (arg))
		x[j] = x[j] * sqrt (factor)
		}
	    }

	return
end


# compute the phase of a complex function

procedure phband (n,shift,x,phase)

int	n
double	shift
complex	x[ARB]
real	phase[ARB]

double	pi,angle,npi
int	i,ncycle

begin
	pi = 3.1415926535897932d0

	do i = 1, n {
	    angle = 2.d0 * pi * double (i-1) * shift / double (n)
	    if ((aimag (x[i]) == 0.) || (real (x[i]) == 0.))
	    	phase[i] = 0.
	    else
	    	phase[i] = atan2 (aimag(x[i]), real(x[i]))

	    npi = (phase[i] - angle) / pi
	    if (npi > 0)
	        ncycle = (npi + 1) / 2
	    else
	        ncycle = (npi - 1) / 2

	    phase[i] = phase[i] - (2.d0 * pi * double (ncycle))
	    if (i < n/2) phase[i+n/2] = phase[i] - angle
	    }

	return
end

# Feb 25 1993	change .or. to || in line 134

# Jun 23 1994	In FLTER, use TFILT to tell whether template hass been filtered
# Dec  7 1994	In FLTER, also use TFILT to turn off high end of filter.
# Dec 15 1994	In RMSFIL, use TFILT to turn off all or high end of filter.
# Dec 16 1994	Return number of points used by RMSFIL
# Dec 19 1994	Filter transforms, not correlation

# Feb 10 1995	Really use square root in flter2
# May 10 1995	Limit length of diagnostic numbers in RMSNORM to 10
# May 11 1995	Print n in RMSNORM

# Jan 30 2007	Remove unnecessary parentheses
