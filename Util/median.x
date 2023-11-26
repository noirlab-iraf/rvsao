# File Util/median.x
# August 5,2 008
# By Doug Mink

procedure median (x, n, xmed, xm1, xm2)

double	x[ARB]		# Vector of numbers of which median is to be found
int	n		# Size of vector
double	xmed		# Median (returned)
double	xm1,xm2		# Quartile points (returned)

double	xx
int	l,ir,i,j,n2,n14,n34

define first_	10
define second_	20

begin

	if (n <= 0) {
	    xmed = 0.
	    return
	    }
	else if (n == 1) {
	    xmed = x[1]
	    return
	    }

	l = n / 2 + 1
	ir = n

first_
	if (l > 1) {
	    l = l - 1
	    xx = x[l]
	    }
	else {
	    xx = x[ir]
	    x[ir] = x[1]
	    ir = ir - 1
	    if (ir == 1) {
		x[1] = xx
		n2 = n / 2
		if (mod (n, 2) == 0) {
		    xmed = 0.5 * (x[n2] + x[n2+1])
		    }
		else {
		    xmed = x[n2+1]
		    }
		n14 = int (n * 0.25 + 0.5)
		n34 = int (n * 0.75 + 0.5)
		xm1 = x[n14]
		xm2 = x[n34]
		return
		}
	    }
	i = l
	j = l + l

second_
      if (j <= ir) {
	    if (j < ir) {
		if (x[j] < x[j+1])
		    j = j + 1
		}
	    if (xx < x[j]) {
		x[I] = x[j]
		i = j
		j = j + j
		}
	    else {
		j = ir + 1
		}
	    go to second_
	    }
	x[i] = xx
	go to first_

end

# Aug  5 2008	New subroutine separated from PXCSAO
