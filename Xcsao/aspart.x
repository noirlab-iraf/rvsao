# File rvsao/Xcor/aspart.x
# May 1, 1998
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After John Tonry

# Copyright(c) 1998 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.


procedure aspartf (n, k1,k4, shift, x, arms,srms, tfilt)

# Compute the anti- and symmetric parts of a complex buffer, shifted by shift
# Assume a transform of a real function, and filter cuttoffs k1 and k2.

int	n
int	k, k1, k4	# Filter
double	shift		# delta from zero of the correlation peak
			# in pixels of the log-wavelngth array. 
complex x[ARB]		# Cross-correlation vector		
double	arms		# anti-symmetric rms
double	srms		# symmetric rms
int	tfilt		# If =1, template has already been filtered
			# If =2, turn off high pass filter
			# If =3, turn off high pass filter on object
			#        and do not filter template

double	angle, f, aterm, sterm
complex phase
double	pi
int	n4

begin
	pi = 3.1415926535897932d0
	if (tfilt > 1) {
	    n4 = n
	    }
	else {
	    n4 = k4
	    }

	arms = 0.d0
	srms = 0.d0

	do k = k1, n4 {
	    angle = -2.d0 * pi * double (k) * shift / double (n)
	    phase = complex (cos(angle), sin(angle))
	    if (k == 0 || k == n/2)
	        f = 1.d0
	    else
	        f = 2.d0
	    aterm = aimag (phase*x[k+1])
	    sterm = real (phase*x[k+1])
	    arms = arms + f * aterm * aterm
	    srms = srms + f * sterm * sterm
	    }

# Divide by n since it is the transform

	arms = sqrt (arms) / n
	srms = sqrt (srms) / n

	return
end

# Jul 24 1990	New program

# Dec 15 1994	Make high frequency filter optional

# Sep 20 1995	Put in 16 digits of PI

# May  1 1998	Drop old aspart subroutine
