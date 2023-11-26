# File rvsao/Xcor/buffers.x
# September 16, 1991
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After John Tonry and Guillermo Torres

# Copyright(c) 1991 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# Various subroutines to copy and flip buffers

#--- Copy a real buffer to a complex one

procedure rcvec (n,buf,x)

int	n		# Length of vectors
real	buf[n]		# Real vector
complex	x[n]		# Complex vector (returned)

int	j

begin
	do j = n,1,-1 {
	    x[j] = complex (buf[j])
	    }

	return
end


#--- Copy a complex buffer to a real one

procedure crvec (n,x,buf)

int	n		# Length of vectors
complex	x[n]		# Complex vector
real	buf[n]	# Real*4 vector (returned)

int	j

begin
	do j = 1, n {
	    buf[j] = real (x[j])
	    }

	return
end


#--- Reverse the order of a real vector

procedure flip (n,buf)

int	n		# Length of vector
real	buf[n]

int	j
real	temp

begin
	do j = 1, n/2 {
	    temp = buf[j]
	    buf[j] = buf[j+n/2]
	    buf[j+n/2] = temp
	    }

	return
end


#--- Reverse the order of a complex vector

procedure cflip (n,buf)

int	n		# Length of vector
complex	buf[n]

complex	temp
int	j

begin

	do j = 1, n/2 {
	    temp = buf[j]
	    buf[j] = buf[j+n/2]
	    buf[j+n/2] = temp
	    }

	return
end
