# File rvsao/Xcor/xcmin.x
# October 31, 1991
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After John Tonry's MINI  (Revision 2.0, 11/17/82.)

# Copyright(c) 1991 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Black box minimization (and fitting) program.

procedure xcmin (npts,x,y,n,a,debug,covar,ierr)

int	npts		# Dimension of data vector
real	x[ARB]		# Independent variable vector
real	y[ARB]		# Data vector including data to be fit
int	n		# Number of parameters to be varied in searching
			#  for a minimum.
double	a[ARB]		# This vector serves three purposes:
			#  it passes the initial estimate for the parameters
			#  it passes arguments to the function to be minimized
			#  it returns the minimum parameter values
bool	debug		# Print diagnostic messages if TRUE
double	covar[ARB]	# RETURNED diagonal of covariance matrix 
int	ierr		# RETURNED error condition:
			#  0 for no error
			#  1 for maximum iteration exceeded
			#  2 for singular matrix

#     descriptions of some of the variables:
#     a - argument for the function
#     a0 - current guess for the minimum
#     ai - increments for a0 in computing derivatives
#     da - vector from a0 to new guess for the minimum
#     df - first derivatives of the function
#     d2f - second derivatives of the function
#     lambda - governs mix of gradient and analytic searches
#     iter - maximum number of iterations
#     qfrac - maximum fractional change for successful exit
#
#     The calling program should be as follows:
#********************************************************
#     real rdata[npts]
#     double a[4]
#     double covar[4]
#     ... (initialize a to the guess for the minimum)
#     call xcmin (npts,x,y,4,a,debug,covar,ier)
#     ...
#     procedure mfunk (npts,x,y,a,chi)
#     real x(npts),y(npts)
#     double a[4], chi
#     ... (define the function)
#************************************************************

double	df[16],a0[16],da[16],ai[16]
double	d2f[16,16],cov[16,16]
double	fnow, fthen, fminus, curve, qfrac, dfrac, places, vary, det
double	base, chi, err
int	lambda,itermax, i,j,k, ii, iter

begin

#  quit if n is greater than the space allocated
	if (n > 16) {
	    call printf ("EMMIN: too many parameters, maximum = 16\n")
	    return
	    }
	if (debug) {
	    call printf ("EMMIN: %d points fit with %d parameters\n")
		call pargi (npts)
		call pargi (n)
	    }

#  define a few parameters
	lambda = -3
	base = 10.d0
#	itermax = 20
	itermax = 40
	qfrac = 1.d-6
	dfrac = .02d0
	places = 1.d-7
	vary = 1.d-5

#  initialize a0
	call amovd (a, a0, n)
	call aclrd (cov,256)
	call aclrd (d2f,256)

	iter = 0
	call mfunk (npts,x,y,a,fnow)
#	if (debug) call mvprint(n,a0,fnow,iter,lambda)

#  initialize ai
	do i = 1, n {
	    ai[i] = dabs (vary * a0[i])
	    if (ai[i] == 0) ai[i] = 1d-6
	    }

#  begin iteration to find minimum
	do iter = 1, itermax {
	    fthen = fnow

#	Initialize parameters
	    call amovd (a0, a, n)

# 	Numerical derivatives

#	1st derivatives
	    do j = 1,n {
		a[j] = a0[j] + ai[j]
		call mfunk (npts,x,y,a,df[j])
		a[j] = a0[j]
		}

#	 off-diagonal 2nd derivatives.
	    do j = 2,n {
		do k = 1,j-1 {
	            a[k] = a0[k] + ai[k]
	            a[j] = a0[j] + ai[j]
		    call mfunk (npts,x,y,a,chi)
	            d2f[j,k] = (chi - df[j] - df[k] + fnow) / (ai[j]*ai[k])
	            d2f[k,j] = d2f[j,k]
	            a[j] = a0[j]
	            a[k] = a0[k]
		    }
		}

#	 on-diagonal 2nd derivatives, and fix the 1st ones.
	    do j = 1,n {
		a[j] = a0[j] - ai[j]
		call mfunk (npts,x,y,a,fminus)
		d2f[j,j] = (fminus + df[j] - 2*fnow) / (ai[j]*ai[j])
		df[j] = (df[j] - fminus) / (2*ai[j])
		a[j] = a0[j]
		}

#	compute better estimates for the increments.
	    do j = 1, n {
		curve = d2f[j,j]
		if (curve == 0) curve = 1.d-5
		ai[j] = sqrt ((df[j] * dfrac/curve)**2 + abs(fnow*places/curve))
		}

#     begin loop to find a direction along which function decreases

	    do ii = 1, 15 {

#	Get weight matrix
		do j = 1, n {
		    do k = 1, j-1 {
			cov[j,k] = d2f[j,k]
			cov[k,j] = d2f[j,k]
			}
		    cov[j,j] = dabs (d2f[j,j] * (1.d0 + base**lambda))
		    }

		call matinv (cov,n,det)
		if (det == 0) {
		    if (debug)
			call printf ("EMMIN:  singular matrix\n")
		    do j = 1, n {
			do k = 1, j-1 {
			    cov[j,k] = d2f[j,k]
			    cov[k,j] = d2f[j,k]
			    }
			cov[j,j] = dabs(d2f[j,j]*(1 + base**lambda))
			}
		    ierr = 2
		    return
		    }

#	    Multiply to get dA
		do j = 1, n {
		    da[j] = 0
		    do k = 1, n {
		        da[j] = da[j] - (cov[j,k] * df[k])
			}
		    }

#	    Get new function value
		do j = 1, n {
		    a[j] = a0[j] + da[j]
		    }
		call mfunk (npts,x,y,a,fnow)

#     test for whether the function has decreased
#     if so, adopt the new point and decrement lambda
#     else, increment lambda, and get a new weight matrix

		if (fnow < fthen) break
		lambda = lambda + 1
		}

#     normal exit, the function at a0 + da is less than at a0

	    call amovd (a,a0,n)
	    lambda = lambda - 1

#     print the current status and test to see whether the function
#     has varied fractionally less than qfrac.

#	    if (debug) call mvprint(n,a0,fnow,iter,lambda)
	    if (dabs ((fthen - fnow)/fnow) < qfrac) break
	    }

#     this is the final computation of the covariance matrix
#     quit if no minimum was found in the allowed number of iterations

	if (iter > itermax) {
	    if (debug)
		call printf ("EMMIN: maximum iteration exceeded\n")
	    ierr = 1
	    }
	else
	    ierr = 0

#     finally, compute the covariance matrix

	do j = 1, n {
	    do k = 1, n {
	        cov[j,k] = d2f[j,k] / 2.d0
		}
	    }
	call matinv (cov,n,det)
	do j = 1, n {
	    err = sqrt (abs (cov[j,j]))
	    if (cov[j,j] < 0) err = -err
	    cov[j,j] = err
	    }
	do j = 2, n {
	    do k = 1, j-1 {
	        cov[j,k] = cov[j,k] / (cov[j,j]*cov[k,k])
	        cov[k,j] = cov[j,k]
		}
	    }

	do j = 1, n {
	    covar[j] = cov[j,j]
	    }

	return
end


# simple subroutine to print the current parameters

procedure mvprint (n,a,fnow,niter,lambda)

int	n
double	a[ARB]
double	fnow
int	niter
int	lambda

int	i

begin
	call printf (" a(i) = ")
	do i = 1, n {
	    call printf ("%10.4f ")
		call pargd (a[i])
	    }

	call printf ("\n f = %15.7g    iter = %d   lambda = %d\n")
	    call pargd (fnow)
	    call pargi (niter)
	    call pargi (lambda)
	return
end
# Oct 31 1991	New program

# Apr 13 1994	Drop unused variable n2
