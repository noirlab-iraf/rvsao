# File rvsao/Util/fourm.x
# September 16, 1991
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# From Press, et al., Numerical Recipes FOUR1
# Modified by Dan Fabricant, Center for Astrophysics

#--- Perform forward (isign=1) or reverse (isign=-1) Fourier transform on vector

procedure fourm (xdata,nin,isign)

real	xdata[ARB]	# Complex vector to transform (result returned)
int	nin		# Log base 2 of length of vector
int	isign		# 1=forward, -1=reverse

double	wr,wi,wpr,wpi,wtemp,theta
real	xr, xi, tempi, tempr
int	i,j,m,mmax,n,nn,istep

begin
	nn = 2 ** nin
	n = 2 * nn
	j = 1

#  Exchange the two complex numbers
	do i = 1 ,n, 2 {
	    if (j > i) {
		tempr = xdata[j]
		tempi = xdata[j+1]
		xdata[j] = xdata[i]
		xdata[j+1] = xdata[i+1]
		xdata[i] = tempr
		xdata[i+1] = tempi
		}
	    m = n / 2

	    while (m >= 2 && j > m) {
		j = j - m
		m = m / 2
		}
	    j = j+m
	    }

#  Danielson-Lanczos section executed nin times
	mmax = 2
	while (n > mmax) {
	    istep = 2 * mmax
	    theta = 6.28318530717959d0 / (isign * mmax)
	    wpr = -2.d0 * ((sin (0.5d0 * theta)) ** 2)
	    wpi = sin (theta)
	    wr = 1.d0
	    wi = 0.d0

	    do m = 1, mmax, 2 {
		do i = m, n, istep {
		    j = i + mmax

		#  Danielson-Lanczos formula
		    xr = real (wr)
		    xi = real (wi)
		    tempr = (xr * xdata[j]) - (xi * xdata[j+1])
		    tempi = (xr * xdata[j+1]) + (xi * xdata[j])
		    xdata[j] = xdata[i] - tempr
		    xdata[j+1] = xdata[i+1] - tempi

		    xdata[i] = xdata[i] + tempr
		    xdata[i+1] = xdata[i+1] + tempi
		    }
	    # Trigonometric recurrence
		wtemp = wr
		wr = (wr * wpr) - (wi * wpi) + wr
		wi = (wi * wpr) + (wtemp * wpi) + wi
		}
	    mmax = istep
	    }

# Added to make plug compatible with previously used subroutine (FFT2C)

	if (isign < 0) {
	    do i = 1, n {
		xdata[i] = xdata[i] / nn
		}
	    }

	return
end
