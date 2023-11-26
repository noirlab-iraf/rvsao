# File rvsao/Util/polfit.x
# October 31, 1991
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Bevington, page 141

# Polynomial least squares fitting program, almost identical to the
# one in Bevington, "Data Reduction and Error Analysis for the
# Physical Sciences," page 141.  I changed the argument list and
# removed the weighting:  y = a(1) + a(2)*x + a(3)*x**2 + a(3)*x**3 + . . .

procedure polfit (x,y,npts,nterms,a,chisqr)

double	x[ARB]		# Array of independent variable points
double	y[ARB]		# Array of dependent variable points
int	npts		# Number of data points to fit
int	nterms		# Number of parameters to fit
double	a[ARB]		# Vector containing current fit values
double	chisqr

double xterm,yterm,delta,chisq,determ()
double sumx[19],sumy[10],array[10,10],xi,yi
int	i,j,k,l,n,nmax

begin
	nmax = (2 * nterms) - 1 
	call aclrd (sumx,nmax)
	call aclrd (sumy,nterms)

# accumulate weighted sums
	chisq = 0.d0
	do i = 1, npts {
	    xi = x[i] 
	    yi = y[i] 
	    xterm = 1.d0
	    do n = 1, nmax {
		sumx[n] = sumx[n] + xterm
		xterm = xterm * xi
		}
	    yterm = yi 
	    do n = 1, nterms {
		sumy[n] = sumy[n] + yterm
		yterm = yterm * xi
		}
	    chisq = chisq + yi*yi
	    }

# construct matrices and calculate coeffients
	do j = 1, nterms {
	    do k = 1, nterms {
		n = j + k - 1 
		array[j,k] = sumx[n]
		}
	    }
	delta =  determ (array,nterms)
	if (delta == 0.d0) {
	    chisqr = 0.d0
	    call aclrd (a, nterms)
	    return
	    }

	do l = 1, nterms {
	    do j = 1, nterms {
		do k = 1, nterms {
		    n = j + k - 1 
		    array[j,k] = sumx[n]
		    }
		array[j,l] = sumy[j]
		}
	    a[l] = determ (array,nterms) / delta
	    }

# calculate chi square
	do j = 1, nterms {
	    chisq = chisq - (2.d0 * a[j] * sumy[j])
	    do k = 1, nterms {
		n = j + k - 1
		chisq = chisq + (a[j] * a[k] * sumx[n])
		}
	    }
	chisqr = chisq / (npts - nterms)

	return
end


#--- calculate the determinant of a square matrix
#    this subprogram destroys the input matrix array
#    from bevington, page 294.

double procedure determ (array,norder)

double array[10,10]	# Input matrix array
int norder		# Order of determinant (degree of matrix)

double	save,det
int	i,j,k,k1

begin
	det = 1.d0
	do k = 1, norder {

	# Interchange columns if diagnoal element is zero
	    if (array[k,k] == 0) {
		do j = k, norder {
		    if (array[k,j] == 0) {
			det = 0.d0
			return det
			}
		    }
		do i = k, norder {
		    save = array[i,j] 
		    array[i,j] = array[i,k]
		    array[i,k] = save 
		    }
		det = -det
		}

	# Subtract row k from lower rows to get diagonal matrix 
	    det = det * array[k,k]
	    if (k < norder) {
		k1 = k+1
		do i = k1, norder {
		    do j = k1, norder {
			array[i,j] = array[i,j]-array[i,k]*array[k,j]/array[k,k]
			}
		    }
		}
	    }

	return det
end
