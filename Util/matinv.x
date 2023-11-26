# File rvsao/Util/matinv.x
# October 29, 1991
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Bevington, pages 302-303.

#  Invert a symmetric matrix up to 10x10 and calculate its determinant

procedure matinv (array,norder,det)

double	array[16,16]	# Input matrix which is replaced by its inverse
int	norder		# Degree of matrix (order of determinant)
double	det		# Determinant of input matrix

double	amax,save
int	ik[10],jk[10]
int	i,j,k,l

begin
	det = 1.d0
	do k = 1,norder {

# find largest element array[i,j] in rest of matrix
	    amax = 0. 
21	    do i = k, norder {
		do j = k, norder {
		    if (dabs(amax) <= dabs(array[i,j])) {
			amax = array[i,j] 
			ik[k] = i 
			jk[k] = j 
			}
		    }
		}

# interchange rows and columns to put amax in array[k,k]
	    if (amax == 0) {
		det = 0.
		return
		}
	    i = ik[k] 
	    if (i < k) go to 21
	    if (i > k) {
		do j = 1, norder {
		    save = array[k,j] 
		    array[k,j] = array[i,j]
		    array[i,j] = -save
		    }
		}
	    j = jk[k] 
	    if (j < k) go to 21
	    if (j > k) {
		do i = 1, norder {
		    save = array[i,k] 
		    array[i,k] = array[i,j]
		    array[i,j] = -save
		    }
		}

# accumulate elements of inverse matrix 
	    do i = 1,norder {
		if (i != k) array[i,k] = -array[i,k]/amax
		}
	    do i = 1, norder {
		do j = 1, norder {
		    if (i != k && j != k) 
			array[i,j] = array[i,j]+array[i,k]*array[k,j]
		    }
		}
	    do j = 1,norder {
		if (j != k) array[k,j] = array[k,j]/amax
		}
	    array[k,k] = 1./amax
	    det = det*amax
	    }

#  Restore ordering of matrix
	do l = 1, norder {
	    k = norder-l+1
	    j = ik[k] 
	    if (j > k) {
		do i = 1, norder {
		    save = array[i,k] 
		    array[i,k] = -array[i,j]
		    array[i,j] = save 
		    }
		}
	    i = jk[k] 
	    if (i > k) {
		do j = 1, norder {
		    save = array[k,j] 
		    array[k,j] = -array[i,j]
		    array[i,j] = save 
		    }
		}
	    }

	return
	end
