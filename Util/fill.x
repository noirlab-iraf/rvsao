# File rvsao/Util/fill.x
# April 29, 2014
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1999-2014 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

procedure filllist (npts, pix, wav, order, pixfill, debug)
int	npts		# Number of points in spectrum
real	pix[ARB]	# Array of points to edit
real	wav[ARB]	# Array of wavelengths corresponding to points to edit
int	order		# Echelle order=aperture, 0=ignore
bool	pixfill		# TRUE if pixels filled, FALSE if wavelength
bool	debug		# TRUE print diagnostic messages

real	w1, w2		# Wavelengths between which to interpolate
int	p1, p2		# Pixels between which to interpolate
int	fd, lpath, nlrep, i, ispix
char	linename[SZ_LINE]
char	dirpath[SZ_PATHNAME]
char	filename[SZ_PATHNAME]
char	filepath[SZ_PATHNAME]
int	fscan(), open(), stridx(), strlen(), strmatch()

begin

#  Directory containing line lists
	call clgstr ("linedir", dirpath, SZ_PATHNAME)
	lpath = strlen (dirpath)
	if (dirpath[lpath] != '/') {
	    dirpath[lpath+1] = '/'
	    dirpath[lpath+2] = EOS
	    }

	pixfill = FALSE

#  Lines to replace with surrounding spectrum
	call clgstr ("badlines", filename, SZ_PATHNAME)
	if (strlen (filename) > 0) {
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
		    call gargwrd (linename,SZ_LINE)
#		    if (debug) {
#			call printf ("FILL_LIST: Line %d is %4.4s: %7.2f-%7.2f\n")
#			    call pargi (nlrep+1)
#			    call pargstr (linename)
#			    call pargr (w1)
#			    call pargr (w2)
#			}
		    ispix = strmatch (linename, "pix")
		    if (ispix > 0) {
			call gargi (i)
			if (order > 0) {
#			    if (debug) {
#				call printf ("FILL_LIST: order=%d, i=%d\n")
#				    call pargi (order)
#				    call pargi (i)
#				}
			    if (i < order)
				next
			    if (i > order)
				break
			    }
			p1 = int (w1 + 0.5)
			p2 = int (w2 + 0.5)
#			if (debug) {
#			    call printf ("FILL_LIST: About to fill pixels %d=%d: %d-%d\n")
#				call pargi (order)
#				call pargi (i)
#				call pargi (p1)
#				call pargi (p2)
#			    }
			call fillpix (npts, pix, p1, p2, debug)
			pixfill = TRUE
			}
		    else
			call fillwav (npts, pix, wav, w1, w2, linename, debug)
		    nlrep = nlrep + 1
		    }
		call close (fd)
		if (debug) {
		    call printf ("FILL_LIST: %d lines replaced from line file %s\n")
			call pargi (nlrep)
			call pargstr (filepath)
		    if (pixfill)
			call printf ("FILL_LIST: Lines filled in pixel space\n")
		    }
		}
	    else {
		call printf ("FILL_LIST: Cannot find bad line file %s\n")
		    call pargstr (filepath)
		}
	    }
	else
	    call printf ("FILL_LIST: No bad line file specified\n")
end

procedure fillwav (npts, pix, wav, w1, w2, linename, debug)
 
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
	    call printf ("FILL_WAV: %s line replaced from %.2dA to %.2dA\n")
		call pargstr (linename)
		call pargr (w1)
		call pargr (w2)
	    }

	call fillpix (npts, pix, p1, p2, debug)
end

 
procedure fillpix (npts, pix, ip1, ip2, debug)
 
int	npts		# Number of points in spectrum
real	pix[ARB]	# Array of points to edit
int	ip1, ip2	# Pixels between which to interpolate
bool	debug		# true for debugging information
 
int	i, temp
real	p1, p2, slope
 
begin
 
# If ip2 < ip1, switch the order.
	if (ip2 < ip1) {
	    temp = ip1
	    ip1 = ip2
	    ip2 = temp
	    }
 
	if (ip1 < 2 && ip2 < npts) {
	    p1 = pix[ip2 + 1]
	    slope = 0
	    }
	else if (ip2 > npts-1 && ip1 > 1) {
	    p1 = pix[ip1 - 1]
	    slope = 0
	    }
	else if ( ip1 > 1 && ip2 < npts ) {
	    p1 = pix[ip1 - 1]
	    p2 = pix [ip2 + 1]
	    slope = (p2 - p1) / real (ip2 -ip1 + 2)
	    }
 
	do  i = ip1, ip2 {
	    pix[i] = p1 + slope * real (i - ip1 + 1)
	    }
	if (debug) {
	    call printf ("FILL_PIX: pixels replaced from %d to %d\n")
		call pargi (ip1)
		call pargi (ip2)
	    }
end

 
procedure fillpt (npts, pix, ip1, val, debug)
 
int	npts		# Number of points in spectrum
real	pix[ARB]	# Array of points to edit
int	ip1		# Pixel of value to replace
real	val		# Replacement value for this pixel of spectrum
bool	debug		# true for debugging information
 
int	ip
 
begin
 
	ip = ip1
	if (ip > npts) {
	    ip = npts
	    }
	else if (ip < 1) {
	    ip = 1
	    }
 
	pix[ip] = val

	if (debug) {
	    call printf ("FILL_PT: pixel %d replaced with %f\n")
		call pargi (ip)
		call pargi (val)
	    }
end

# Jan 31 1997	New file renamed subroutine FIXIT to FILLPIX

# May 21 1999	Test for existence of bad pixel file
# Sep  8 1999	Fix bug interpolating wavelength limits

# Feb  6 2002	Fix bug filling single pixels and rewrite fillpix()

# Mar  6 2008	Add pixfill argument
# Mar  7 2008	If line name contains "pix", fill pixels, not wavelength
# Mar 10 2008	If order specified, use only lines with that order after "pix"
# May  6 2008	Add debug argument to fillpix()

# Apr 30 2014	Add fillpt() to set a single pixel value
