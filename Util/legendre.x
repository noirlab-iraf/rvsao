# File rvsao/Util/legendre.x
# April 13, 1994
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After John Tonry (11/19/82) and Guillermo Torres (Jan/1989)

# Copyright(c) 1994 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# All sorts of routines to fit and evaluate legendre polynomials, and to
# convert coefficients to those of plain polynomials

procedure fitlpo (npt,x,y,scale,ncoeff,coeff)

# x and y are data arrays with npt points to fit.
# scale receives xMin and xMax, and returns a and b.
# where z = a * (x - b)
# a polynomial of degree ncoeff-1 is fit and the coefficients are
# returned in coeff.

int	npt
real	x[npt]
real	y[npt]
real	scale[ARB]
int	ncoeff
double	coeff[ncoeff]

double	cov[8*8], vec[8]
double	a, b, z, lpi, det
int	i,j,n,maxfit, ij, ji, irst[2,8]

begin
	maxfit = 8

	a = 2.0d0 / (scale[2] - scale[1])
	b = 0.5d0 * (scale[2] + scale[1])

	do j = 1, ncoeff {
	    vec[J] = 0
	    }
	do i = 1, ncoeff*ncoeff {
	    cov[i] = 0
	    }

	do n = 1, npt {
	    z = a * (double (x[n]) - b)

	    do j = 1, ncoeff {
		coeff[j] = lpi (j,z)
		vec[j] = vec[j] + coeff[j] * y[n]
		do i = 1, j {
		    ji = (j - 1) * ncoeff + i
		    cov[ji] = cov[ji] + coeff[i] * coeff[j]
		   }
		}
	    }

	do j = 2, ncoeff {
	    do i = 1, j-1 {
		ij = (i - 1) * ncoeff + j
		ji = (j - 1) * ncoeff + i
		cov[ij] = cov[ji]
		}
	    }

	call invert (ncoeff,cov,irst,det)

	if (det == 0) {
	    call printf ("fitlpoly: singular covariance matrix\n")
	    return	
	    }

	do j = 1, ncoeff {
	    coeff[j] = 0
	    do i = 1, ncoeff {
		coeff[j] = coeff[j] + cov[(j-1)*ncoeff+i] * vec[i]
		}
	    }

	scale[1] = a
	scale[2] = b

	return
end


real procedure lpoly (x,scale,ncoeff,coeff)

real	x
double	coeff[ARB]
int	ncoeff
real	scale[ARB]

double	temp, z, lpi()
int	i

begin
	z = scale[1] * (x - scale[2])
	temp = 0
	do i = 1, ncoeff {
	    temp = temp + coeff[i] * lpi (i,z)
	    }
	lpoly = temp
	return
end


# legendre polynomial of order i-1

double procedure lpi (i,z)

int	i
double	z

begin
	switch (i) {
	    case 1:
		lpi = 1.d0
	    case 2:
		lpi = z
	    case 3:
		lpi =  -.5d0 + 1.5d0*z*z
	    case 4:
		lpi = z*(-1.5d0 + 2.5d0*z*z)
	    case 5:
		lpi = .375d0 + z*z*(-3.75d0 + 4.375d0*z*z)
	    case 6:
		lpi = z*(1.875d0 + z*z*(-8.75d0 + 7.875*z*z))
	    case 7:
		lpi = -.3125d0 + z*z*(6.5625d0 + z*z*(-19.6875d0 +
		      14.4375d0*z*z))
	    case 8:
		lpi = z*(2.1875d0 + z*z*(19.6875d0 + z*z*(-43.3125d0 +
		      26.8125d0*z*z)))
	    default:
	    }

	return
end


procedure polyc (nc,scale,coeff,pcoeff)

# routine to convert coefficients of a 1-d legendre polynomial to 
# coefficients of a plain polynomial
# The nc coefficients are stored with the low order coefficient first

int	nc
real	scale[2]
double	coeff[ARB]
double	pcoeff[ARB]

double	sc[2], temp
int	i,k,ibc()

# Legendre polynomials 0 - 7

double lcoeff[8,8]
data lcoeff/ 1d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,
	     0d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0,
	   -.5d0,0d0,1.5d0,0d0,0d0,0d0,0d0,0d0,
	     0d0,-1.5d0,0d0,2.5d0,0d0,0d0,0d0,0d0,
	  .375d0,0d0,-3.75d0,0d0,4.375d0,0d0,0d0,0d0,
	     0d0,1.875d0,0d0,-8.75d0,0d0,7.875,0d0,0d0,
	-.3125d0,0d0,6.5625d0,0d0,-19.6875d0,0d0,14.4375d0,0d0,
	     0d0,2.1875d0,0d0,19.6875d0,0d0,-43.3125d0,0d0,26.8125d0/

begin
	sc[1] = double (scale[1])
	sc[2] = double (scale[2])
	do k = 0, nc-1 {
	    temp = 0
	    do i = k, nc-1 {
	    	temp = temp + coeff[i+1] * lcoeff[k+1,i+1]
		}
	    pcoeff[k+1] = temp
	    }

	do k = 0, nc-1 {
	    temp = 0
	    do i = k, nc-1 {
	        if (sc(2) != 0)
	    	    temp = temp + pcoeff[i+1] * sc[1]**i * 
     			   (-sc[2])**(i-k) * double (ibc(k,i))
	        else
	    	    temp = temp + pcoeff[i+1] * sc[1]**i
		}
	    pcoeff[k+1] = temp
	    }

	return
end


#  return the (m n) binomial coefficient

int procedure ibc (m, n)

int	m,n

int	iret,i

begin
	iret = 1
	do i = m+1, n {
	    iret = iret * i
	    }
	do i = 2, n-m {
	    iret = iret / i
	    }
	return iret
end
# Sep 17 1991	New program

# Apr 13 1994	Add working array for invert
