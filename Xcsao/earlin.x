# File rvsao/Xcor/earlin.x
# July 28, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1991-2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

define NMAX	1024

procedure earlin (x, y, n, v, meg)

real	x[ARB]		# Cross-correlation vector indices
real	y[ARB]		# Cross-correlation vector
int	n		# Number of points used in peak fit
double	v[ARB]		# Fit vector (returned)
int	meg		# Number of coefficients in fit

double	c[NMAX,NMAX],xx[NMAX],e[NMAX],det,xxx
int	m,mp1,i,j,k,jmi,nmx

begin
	nmx = NMAX - 1
	if (n > nmx) {
	    call eprintf ("*** Too many points in correlation, %d -> %d\n")
		call pargi (n)
		call pargi (nmx)
	    n = nmx
	    }
	m = abs (meg) 
	mp1 = m + 1 
	do i = 1, m {
	    do j = 1, mp1 {
		c[i,j] = 0.d0
		}
	    }

	do i = 1, n {
	    xxx = 1.d0
	    do j = 1, m {
		if (j != 1) xxx = xxx * x[i]
		xx[j] = xxx
		c[j,m+1] = c[j,m+1] + xx[j] * y[i]
		do k = 1,j {
		    c[j,k] = c[j,k] + xx[j] * xx[k]
		    }
		}
	    }

	do j = 2, m {
	    jmi = j - 1
	    do k = 1, jmi {
		c[k,j] = c[j,k]
		}
	    }

	call dsimul (m,c,v,0.d0,0,NMAX,det)

	if (meg < 0) return 

	do i = 1, m {
	    e[i] = sqrt (c[i,i])
	    do j = 1, i {
		c[i,j] = c[i,j] / (e[i] * e[j])
		}
	    }

	return 
end 


procedure dsimul (n,a,x,eps,indic,nrc, deter)
 
#--- This function returns the value of the determinant of a
#    matrix.  in addition, the inverse matrix may be calculated 
#    in place, and the solution vector of the corresponding linear 
#    system computed.  Gauss-jordan elimination with maximum
#    pivot strategy is employed, using double precision arithmetic.
#    if the matrix exceeds the maximum size (NMAX by NMAX), or if it is
#    singular, a true zero is returned. 
 
int	n		# Size of the matrix (n by n) 
double	a[nrc,nrc]	# Matrix 
double	x[n]		# Solution vector 
double	eps 		# Small number to be used as a test for singularity 
int	indic		# Control parameter:
			#  <0, the inverse is computed in place 
			#  =0, the matrix is assumed to be augmented
			#      and the solution and inverse are computed
			#  >0, only the solution is computed 
int	nrc		# Dimension of A in the calling program 
double	deter		# Determinant

double	y[NMAX],aijck,pivot
int	irow[NMAX],jcol[NMAX],jord[NMAX]
int	max,i,j,k,iscan,jscan,irowi,jcoli,intch,jtemp,km1,nm1,ip1
int	irowj,jcolj,irowk,jcolk
bool	bflag

begin
	max = n 
	if (indic >= 0) max = n + 1 
	if (n > NMAX) {
	    deter = 0.d0
	    return 
	    }
	deter = 1.0 d0
	do k = 1, n {
	    km1 = k - 1
	    pivot = 0.d0
	    do i = 1, n {
		do j = 1, n {
		    bflag = FALSE
		    if (k > 1) {
			do iscan = 1, km1 {
			    if (i == irow[iscan]) {
				bflag = TRUE
				break
				}
			    }
			do jscan = 1, km1 {
			    if (j == jcol[jscan]) {
				bflag = TRUE
				break
				}
			    }
			}
		    if (dabs(a[i,j]) > dabs(pivot) && !bflag) {
			pivot = a[i,j]
			irow[k] = i 
			jcol[k] = j 
			}
		    }
		}
	    if (abs (pivot) <= eps) {
		deter = 0.d0
		return
		}
	    irowk = irow[k]
	    jcolk = jcol[k]
	    deter = deter * pivot
	    do j = 1, max {
		a[irowk,j] = a[irowk,j] / pivot
		}
	    a[irowk,jcolk] = 1.d0 / pivot
	    do i = 1,n {
		aijck = a[i,jcolk]
		if (i != irowk) {
		    a[i,jcolk] = -aijck / pivot
		    do j = 1,max {
			if (j != jcolk)
			    a[i,j] = a[i,j] - (aijck * a[irowk,j])
			}
		    }
		}
	     }

	do i = 1,n {
	    irowi = irow[i]
	    jcoli = jcol[i]
	    jord[irowi] = jcoli
	    if (indic >= 0)
		x[jcoli] = a[irowi,max]
	    }

	intch = 0
	nm1 = n - 1
	do i = 1, nm1 {
	    ip1 = i + 1
	    do j = ip1, n {
		if (jord[j] < jord[i]) {
		    jtemp = jord[j]
		    jord[j] = jord[i]
		    jord[i] = jtemp
		    intch = intch + 1
		    }
		}
	    }

	if (intch/2*2 .le. intch) deter = -deter
	if (indic > 0) {
	    return
	    }
	do j = 1, n {
	    do i = 1, n {
		irowi = irow[i]
		jcoli = jcol[i]
		y[jcoli] = a[irowi,j]
		}
	    do i = 1, n {
		a[i,j] = y[i]
		}
	    }

	do i = 1, n {
	    do j = 1, n {
		irowj = irow[j]
		jcolj = jcol[j]
		y[irowj] = a[i,jcolj]
		}
	    do j = 1, n {
		a[i,j] = y[j] 
		}
	    }

	return 
end
# Sep 16 1991	Program written

# Feb 21 1995	Increase maximum dimension from 11 to 256

# Apr 21 1997	Use parameter NMAX for maximum dimension of 512

# Jul 28 2009	Increase maximum number of points, NMAX, from 512 to 1024
