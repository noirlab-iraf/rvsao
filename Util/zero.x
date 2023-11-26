# File hectospec/eigenres/zero.x
# May 4, 1999
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1999 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

procedure zerolist (npts, pix, wav, debug)
int	npts		# Number of points in spectrum
real	pix[ARB]	# Array of points to edit
real	wav[ARB]	# Array of wavelengths corresponding to points to edit
bool	debug		# TRUE print diagnostic messages

real	w1, w2		# Wavelengths between which to interpolate
int	fd, lpath, nlrep
char	linename[SZ_LINE]
char	dirpath[SZ_PATHNAME]
char	filename[SZ_PATHNAME]
char	filepath[SZ_PATHNAME]
int	fscan(), open(), stridx(), strlen()

begin

#  Directory containing line lists
	call clgstr ("linedir", dirpath, SZ_PATHNAME)
	lpath = strlen (dirpath)
	if (dirpath[lpath] != '/') {
	    dirpath[lpath+1] = '/'
	    dirpath[lpath+2] = EOS
	    }

#  Lines to replace with surrounding spectrum
	call clgstr ("badlines", filename, SZ_PATHNAME)
	if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
	    call strcpy (filename,filepath,SZ_PATHNAME)
	else {
	    call strcpy (dirpath,filepath,SZ_PATHNAME)
	    call strcat (filename,filepath,SZ_PATHNAME)
	    }
	fd = open (filepath, READ_ONLY, TEXT_FILE)
	if (fd >0) {
	    nlrep = 0
	    while (fscan(fd) != EOF) {
		call gargr (w1)
		call gargr (w2)
		call gargstr (linename,SZ_LINE)
		call zerowav (npts, pix, wav, w1, w2, linename, debug)
		nlrep = nlrep + 1
		}
	    call close (fd)
	    if (debug) {
		call printf ("ZERO_LIST: %d lines replaced from line file %s\n")
		    call pargi (nlrep)
		    call pargstr (filepath)
		}
	    }
	else {
	    call printf ("ZERO_LIST: Cannot find bad line file %s\n")
		call pargstr (filepath)
	    }
end

procedure zerowav (npts, pix, wav, w1, w2, linename, debug)
 
int	npts		# Number of points in spectrum
real	pix[ARB]	# Array of points to edit
real	wav[ARB]	# Array of wavelengths corresponding to points to edit
real	w1, w2		# Wavelengths between which to interpolate
char	linename[ARB]	# Name of line to be cut
bool	debug		# true for debugging information

int	p1, p2		# Pixels between which to interpolate
int	i, i1, i2
real	dw2

begin

# Find closest pixel to starting wavelength within limits
	if (w1 <= wav[1])
	    p1 = 1
	else if (w1 >= wav[npts])
	    p1 = npts
	else {
	    i1 = 1
	    i2 = npts - 1
	    do i = i1, i2 {
		if (w1 >= wav[i] && w1 < wav[i+1]) {
		    dw2 = (wav[i+1] - wav[i]) * 0.5
		    if (w1 < wav[i] + dw2)
			p1 = i
		    else
			p1 = i + 1
		    break
		    }
		}
	    }

# Find closest pixel to ending wavelength within limits
	if (w2 <= wav[1])
	    p2 = 1
	else if (w2 >= wav[npts])
	    p2 = npts
	else {
	    i1 = 1
	    i2 = npts - 1
	    do i = i1, i2 {
		if (w2 >= wav[i] && w2 < wav[i+1]) {
		    dw2 = (wav[i+1] - wav[i]) * 0.5
		    if (w2 < wav[i] + dw2)
			p2 = i
		    else
			p2 = i + 1
		    break
		    }
		}
	    }
	if (debug) {
	    call printf ("ZEROWAV: %s replaced from %.2dA to %.2dA = %d to %d\n")
		call pargstr (linename)
		call pargr (w1)
		call pargr (w2)
		call pargi (p1)
		call pargi (p2)
	    }

	call zeropix (npts, pix, p1, p2)
end

 
procedure zeropix (npts, pix, p1i, p2i)
 
int	npts		# Number of points in spectrum
real	pix[ARB]	# Array of points to edit
int	p1i, p2i	# Pixels between which to interpolate
 
int	i
real	p1, p2
 
begin
 
# If p2 < p1, switch the order.
	if (p2i < p1i) {
	    p1 = p2i
	    p2 = p1i
	    }
 
# If p2 = p1, increase it by one.
	else if (p2 == p1) {
	    p1 = p1i
	    p2 = p2i + 1
	    }

	else {
	    p1 = p1i
	    p2 = p2i
	    }
 
# Check limits
	if (p1 < 1)
	    p1 = 1
	if (p2 > npts)
	    p2 = npts
 
# Fill range of pixels with zeroes
	do i = p1, p2 {
	    pix[i] = 0.0
	    }
end

# May  4 1999	New file from rvsao/Util/fill.x
# Sep  8 1999	Fix bug interpolating wavelength limits
