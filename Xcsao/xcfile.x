# File rvsao/Xcor/xcfile.x
# February 28, 2012
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1993-2012 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# Print a portion of the filtered cross-correlation function

include "rvsao.h"

procedure xcfile (specfile,specim,mspec,tempfile,tspec,itemp,ncor,xvel,
		  xcor,xcorfile)

char	specfile[SZ_PATHNAME]	# Object spectrum file name
pointer	specim		# Spectrum header structure
int	mspec		# Number of object aperture (ignored if < 1)
char	tempfile[SZ_PATHNAME]	# Template spectrum file name
int	tspec		# Number of template aperture (ignored if < 1)
int	itemp		# Index of template in template list
int	ncor		# Number of points in cross-correlation
real	xvel[ARB]	# Cross-correlation velocities from xcorfit
real	xcor[ARB]	# Cross-correlation returned from xcorfit
char	xcorfile[SZ_PATHNAME]	# Cross-correlation output file name (returned)

real	rindef
real	xcr		# Center velocity in km/sec
real	dxcr		# Half-width in velocity
char	title[SZ_PATHNAME]	# Title line in output file
char	inst[SZ_LINE]	# Instrument name from input file
int	xcfd		# File descriptor for correlation output
real	vmin,vmax	# Velocity limits for correlation output
int	i, ifile, open()
double	dj		# Julian Day of observation
int	strncmp()
char	aperture[8]

int	stridx(), strldx(), strmatch()

include "rvsao.com"
include "results.com"

begin
	rindef = INDEFR

#  If SAO echelle, prefix telescope code to file name
	call imgspar (specim, "INSTRUME",inst,SZ_LINE)
	if (strncmp (inst,"echelle",7) == 0) {
	    xcorfile[1] = 'T'
	    xcorfile[2] = EOS
	    }
	else if (strncmp (inst,"mmtech",6) == 0) {
	    xcorfile[1] = 'M'
	    xcorfile[2] = EOS
	    }
	else if (strncmp (inst,"oroech",6) == 0) {
	    xcorfile[1] = 'W'
	    xcorfile[2] = EOS
	    }
	else
	    xcorfile[1] = EOS

#  Make up a file name as <spectrum file>.<template file>
	if (stridx ('/',specfile[1]) > 0) {
	    ifile = strldx ('/',specfile) + 1
	    call strcat (specfile[ifile], xcorfile, SZ_PATHNAME)
	    }
	else {
	    call strcat (specfile, xcorfile, SZ_PATHNAME)
	    }

# Remove .ms.fits or .fits from end of object file name
	i = strmatch (xcorfile, ".ms.fits")
	if (i > 0) {
	    xcorfile[i-8] = EOS
	    }
	else {
	    i = strmatch (xcorfile, ".fits")
	    if (i > 0) {
		xcorfile[i-5] = EOS
		}
	    }

# Add object spectrum aperture if multispec file
	if (mspec > 0) {
	    call sprintf (aperture, 8, "_%d")
		call pargi (mspec)
	    call strcat (aperture, xcorfile, SZ_PATHNAME)
	    }

# Add template file name
	call strcat (".", xcorfile, SZ_PATHNAME)
	if (stridx ('/',tempfile[1]) > 0) {
	    ifile = strldx ('/',tempfile) + 1
	    call strcat (tempfile[ifile], xcorfile, SZ_PATHNAME)
	    }
	else {
	    call strcat (tempfile, xcorfile, SZ_PATHNAME)
	    }

# Remove .ms.fits or .fits from end of template file name
        i = strmatch (xcorfile, ".ms.fits")
        if (i > 0) {
            xcorfile[i-8] = EOS
            }
        else {
            i = strmatch (xcorfile, ".fits")
            if (i > 0) {
                xcorfile[i-5] = EOS
                }
            }

# Add template spectrum aperture if multispec file
	if (tspec > 0) {
	    call sprintf (aperture, 8, "_%d")
		call pargi (tspec)
	    call strcat (aperture, xcorfile, SZ_PATHNAME)
	    }

# Add .xcor file extension
	call strcat (".xcor", xcorfile)

# Open cross-correlation file
	iferr {xcfd = open (xcorfile, NEW_FILE, TEXT_FILE)} then {
	    call printf ("XCFILE:  Cannot write %s\n")
	    call pargstr (xcorfile)
	    return
	    }

#  Put together a title line and write it to the file
	call strcpy ("filtered cross-correlation: ", title, SZ_PATHNAME)
	call strcat (specname, title, SZ_PATHNAME)
	if (mspec > 0) {
	    call sprintf (aperture, 8, "[%d]")
		call pargi (mspec)
	    call strcat (aperture, title, SZ_PATHNAME)
	    }
	call strcat (" X ", title, SZ_LINE)
	call strcat (tempname[1,itemp], title, SZ_PATHNAME)
	if (tspec > 0) {
	    call sprintf (aperture, 8, "[%d]")
		call pargi (tspec)
	    call strcat (aperture, title, SZ_PATHNAME)
	    }
	dj = 0.d0
	call imgdpar (specim,"HJDN",dj)
	if (dj  == 0.d0) {
	    call imgdpar (specim,"GJDN",dj)
	    }
	call fprintf (xcfd,"%s %.5f\n")
	    call pargstr (title)
	    call pargd (dj)

# Set limits of correlation to print
	if (xcr0 == rindef)
	    xcr = xcrmax
	else
	    xcr = xcr0
	if (xcrdif == rindef) {
	    dxcr = 20. * tvw[itmax]
	    if (dxcr < 1.d0) dxcr = 300.d0
	    }
	else
	    dxcr = xcrdif
	vmin = xcr - dxcr
	vmax = xcr + dxcr

	do i = 1, ncor {
	    if (xvel[i] < vmin)
		next
	    else if (xvel[i] > vmax)
		break
	    call fprintf (xcfd," %f %f\n")
		call pargr (xvel[i])
		call pargr (xcor[i])
	    }
	call close (xcfd)
	return
end
# Jun 15 1993	Fix error message

# Feb 10 1994	Prefix observatory code to file name
# Feb 10 1994	Add Julian day to title
# Feb 11 1994	Return correlation file name
# Apr 13 1994	Drop unused variable c0
# Aug  3 1994	Change common and header from fquot to rvsao

# Jan 31 1995	Change lengths of file and directory names to sz_pathname
# Mar 13 1995	Make sure alloc'd memory is always freed
# Oct 13 1995	SALLOC instead of MALLOC temporary string

# Apr 30 1997	Use NCOR instead of NPTS to loop through cross-correlation
# May  2 1997	Always test against rindef, not INDEF
# Dec  9 1997	Add apertures for object and template spectra

# Feb  9 1998	Do not assume dimensions for XCOR and XVEL

# Aug 13 2007	No longer allocate instrument, just dimension it

# Feb 22 2012	If spectrum or template filename is paths, use file name
# Feb 28 2012	Add .xcor file extension to output file
