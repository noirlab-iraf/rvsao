# File rvsao/Util/invert.x
# December 2, 1991
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After John Tonry (9/2/80)

# Copyright(c) 1991 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

procedure invert (n,a,rst,det)

int	n		# Dimension of the square matrix
double	a[n,n]		# Matrix to be inverted.
			# Upon successful inversion, a contains the inverse.
int	rst[2,n]	# Scratch vector (the row status vector)
double	det		# Determinant of the matrix
			# Set to 0 for a singular matrix, a then is garbage.


double	save, pivot, onrow, cprev, cnow, decr
int	mrank,isign,i,j,k,l,nrow,ncol

begin
	mrank = 0
	isign = 1
	det = 0.d0
	do j = 1, n {
	    do i = 1, 2 {
		rst[i,j] = 0
		}
	    }

#	Loop over columns, reducing each
	do i = 1, n {

#	Find the pivot element
	    pivot = 0
	    nrow = 0
	    ncol = 0
	    do j = 1, n {
		if (rst[1,j] != 0) next
		do k = 1, n {
		    if (rst[1,k] != 0) next
		    if (pivot >= dabs (a[j,k])) next
		    pivot = dabs (a[j,k])
		    nrow = j
		    ncol = k
		    }
		}
	    pivot = a[nrow,ncol]
	    if (pivot == 0) {
		det = 0
		return
		}
	    rst[1,ncol] = nrow
	    rst[2,ncol] = i

#	Swap pivot element onto the diagonal
	    do k = 1, n {
		save = a[nrow,k]
		a[nrow,k] = a[ncol,k]
		a[ncol,k] = save
		}

#	Reduce pivot column
	    do j = 1, n {
		a[j,ncol] = -a[j,ncol] / pivot
		}
	    a[ncol,ncol] = 1 / pivot

#	Reduce other columns
	    do k = 1, n {
		if (k == ncol) next

#	Find maximum of column to check for singularity
		cprev = 0
		do j = 1, n {
		    cprev = dmax1(cprev,dabs(a[j,k]))
		    }

#	Reduce the column
		onrow = a[ncol,k]
		a[ncol,k] = 0
		do j = 1, n {
		    a[j,k] = a[j,k] + onrow * a[j,ncol]
		    }

#	Find the new maximum of the column
		cnow = 0
		do j = 1, n {
		    cnow = dmax1 (cnow,dabs(a[j,k]))
		    }

#	Quit if too many figures accuracy were lost (singular)
		if (cnow == 0)  {
		    det = 0
		    return
		    }
		decr = cprev / cnow
		if (decr > 1.e8) {
		    det = 0
		    return
		    }

		}

	    det = det + dlog10 (dabs (pivot))
	    if (pivot < 0) isign = -isign
	    mrank = mrank + 1
	    }

#     now untangle the mess
	do j = 1, n {
	    do k = 1, n {
		if (rst[2,k] != (n + 1 - j)) next
		ncol = rst[1,k]
		if(ncol == k) break
		do l = 1, n {
		    save = a[l,ncol]
		    a[l,ncol] = a[l,k]
		    a[l,k] = save
		    }
		break
		}
	    }
	if (abs(det) < 35) det = isign * (10. ** det)
	return
end
