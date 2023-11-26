# File rvsao/Util/imgpar.x
# By Doug Mink, Center for Astrophysics
# June 14, 2007

include	<syserr.h>
include	<ctype.h>
define	IDB_ENDVALUE	30

# IMGPAR -- These subroutines do not modify the value argument if the
#           header parameter is not present
#
#	imgbpar (im, key, bval) parameter of type boolean
#	imgcpar (im, key, cval) parameter of type char
#	imgipar (im, key, ival) parameter of type int
#	imglpar (im, key, lval) parameter of type long
#	imgrpar (im, key, rval) parameter of type real
#	imgdpar (im, key, rval) parameter of type double
#	imgdate (im, key, mm, dd, yyyy) ISO or FITS format date parameter
#	imgdtim (im, key, ut) UT time from ISO or FITS format date parameter
#	imgxpar (im, key) delete single- or multi-line header parameter

# IMGBPAR -- Get an image header parameter of type boolean.  False is returned
# if the parameter cannot be found or if the value is not true.

procedure imgbpar (im, key, bval)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
bool	bval		# value returned

int	imaccf()
int	idb_findrecord()
pointer	rp

begin
	if (imaccf (im, key) == YES) {
	    if (idb_findrecord (im, key, rp) != 0)
		bval = (Memc[rp+IDB_ENDVALUE-1] == 'T')
	    }
	return
end


# IMGCPAR -- Get an image header parameter of type char.

procedure imgcpar (im, key, cval)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
char	cval		# parameter value returned

int	imaccf()
pointer	sp, vp
int	ip, cctoc()
errchk	syserrs, imgstr

begin
	if (imaccf (im, key) == YES) {
	    call smark (sp)
	    call salloc (vp, SZ_LINE, TY_CHAR)
	    ip = 1
	    call imgstr (im, key, Memc[vp], SZ_LINE)
	    if (cctoc (Memc[vp], ip, cval) == 0) {
		call printf ("IMGCPAR:  Bad header data %s = %s\n")
		    call pargstr (key)
		    call pargstr (Memc[vp])
		call flush (STDOUT)
		cval = EOS
		}
	    call sfree (sp)
	    }
	return
end


# IMGIPAR -- Get an image header parameter of type integer.

procedure imgipar (im, key, ival)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
int	ival		# parameter value returned

int	imaccf()
pointer	sp, vp
int	itemp, ip, ctoi()
errchk	syserrs, imgstr

begin
	if (imaccf (im, key) == YES) {
	    call smark (sp)
	    call salloc (vp, SZ_LINE, TY_CHAR)
	    ip = 1
	    call imgstr (im, key, Memc[vp], SZ_LINE)
	    if (ctoi (Memc[vp], ip, itemp) == 0) {
		call printf ("IMGIPAR:  Bad header data %s = %s\n")
		    call pargstr (key)
		    call pargstr (Memc[vp])
		call flush (STDOUT)
		}
	    else
		ival = itemp
	    call sfree (sp)
	    }
	return
end


# IMGLPAR -- Get an image header parameter of type long integer.

procedure imglpar (im, key, lval)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
long	lval		# parameter value returned

int	imaccf()
pointer	sp, vp
int	ip, ctol()
long	ltemp
errchk	syserrs, imgstr

begin
	if (imaccf (im, key) == YES) {
	    call smark (sp)
	    call salloc (vp, SZ_LINE, TY_CHAR)
	    ip = 1
	    call imgstr (im, key, Memc[vp], SZ_LINE)
	    if (ctol (Memc[vp], ip, ltemp) == 0) {
		call printf ("IMGLPAR:  Bad header data %s = %s\n")
		    call pargstr (key)
		    call pargstr (Memc[vp])
		call flush (STDOUT)
		}
	    else
		lval = ltemp
	    call sfree (sp)
	    }
	return
end


# IMGSPAR -- Get an image header string parameter

procedure imgspar (im, key, sval,lsval)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
char	sval[ARB]	# parameter value returned
int	lsval		# Maximum length for sval

int	imaccf()
errchk	imgstr

begin
	if (imaccf (im, key) == YES)
	    call imgstr (im, key, sval, lsval)
	return
end


# IMGRPAR -- Get an image header parameter of type real.

procedure imgrpar (im, key, rval)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
real	rval		# parameter value returned

int	imaccf()
pointer	sp, vp
int	ip, ctor()
real	rtemp
errchk	syserrs, imgstr

begin
	if (imaccf (im, key) == YES) {
	    call smark (sp)
	    call salloc (vp, SZ_LINE, TY_CHAR)
	    ip = 1
	    call imgstr (im, key, Memc[vp], SZ_LINE)
	    if (ctor (Memc[vp], ip, rtemp) == 0) {
		call printf ("IMGRPAR:  Bad header data %s = %s\n")
		    call pargstr (key)
		    call pargstr (Memc[vp])
		call flush (STDOUT)
		}
	    else
		rval = rval
	    call sfree (sp)
	    }
	return
end


# IMGDPAR -- Get an image header parameter of type double floating.  If the
# named parameter is a standard parameter return the value directly,
# else scan the user area for the named parameter and decode the value.

procedure imgdpar (im, key, dval)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
double	dval		# parameter value returned

int	imaccf()
pointer	sp, vp
double	dtemp
int	ip, ctod()
errchk	syserrs, imgstr

begin
	if (imaccf (im, key) == YES) {
	    call smark (sp)
	    call salloc (vp, SZ_LINE, TY_CHAR)
	    ip = 1
	    call imgstr (im, key, Memc[vp], SZ_LINE)
#	    call printf ("IMGDPAR: keyword = %s, value = %s\n")
#		call pargstr (key)
#		call pargstr (Memc[vp])
	    if (ctod (Memc[vp], ip, dtemp) == 0) {
		call printf ("IMGDPAR:  Bad header data %s = %s\n")
		    call pargstr (key)
		    call pargstr (Memc[vp])
		call flush (STDOUT)
		}
	    else
		dval = dtemp
	    call sfree (sp)
	    }
	else {
#	    call printf ("IMGDPAR: keyword %s not found\n")
#		call pargstr (key)
	    }
	return
end


# IMGDATE -- Extract month, day and year from a FITS image header string
#            parameter as day/month/year or as ISO yyyy-mm-dd

procedure imgdate (im, key, mm, dd, yyyy)

pointer	im		# image descriptor
char	key[ARB]	# parameter to be returned
int	mm, dd, yyyy

int	i, j, lstr, ix1, ix2, np
int	ctoi(), strlen(), stridx(), strldx()
int	imaccf()
char	slash, dash, tee
pointer	vp, sp, ip, tp

begin
	if (imaccf (im, key) == YES) {
	    call smark (sp)
	    call salloc (vp, SZ_LINE, TY_CHAR)
	    call salloc (tp, SZ_LINE, TY_CHAR)
	    call imgstr (im, key, Memc[vp], SZ_LINE)

#	String must be at least 8 characters long to be a date
	    lstr = strlen (Memc[vp])
	    if (lstr > 7) {

#	If there is a hyphen in the string, assume it is an ISO date
#	If the year is 31 or less, assume it to be to be dd-mm-yy[yy]
		dash = '-'
		tee = 'T'
		ix1 = stridx (dash, Memc[vp])
		if (ix1 > 0) {

#		Year (or day if less than 32)
		    dd = 0
		    ip = vp + ix1 - 1
		    np = ix1 - 1
		    call strcpy (Memc[vp], Memc[tp], np)
		    i = 1
		    if (ctoi(Memc[tp],i,j) > 0) {
			if (j > 31)
			    yyyy = j
			else
			    dd = j
			}
		    else
			yyyy = 0

#		Month
		    ix2 = strldx (dash, Memc[vp])
		    np = ix2 - ix1 - 1
		    call strcpy (Memc[ip+1], Memc[tp], np)
		    i = 1
		    if (ctoi(Memc[tp],i,j) > 0)
			mm = j
		    else
			mm = 0

#		Day (or Year if day is already set)
		    ix1 = ix2
		    ix2 = strldx (tee, Memc[vp])
		    if (ix2 > 0)
			np = ix2 - ix1 - 1
		    else
			np = lstr - ix1
		    ip = vp + ix1 - 1
		    call strcpy (Memc[ip+1],Memc[tp],np)
		    i = 1
		    
		    if (ctoi(Memc[tp],i,j) > 0) {
			if (dd > 0) {
			    yyyy = j
			    if (yyyy < 120)
				yyyy = 1900 + yyyy
			    }
			else
			    dd = j
			}
		    else
			dd = 0
		    }
	    
#	If there is a / in the string, assume is is dd/mm/yy[yy]
		else {
		    slash = '/'
		    ix1 = stridx (slash, Memc[vp])
		    if (ix1 > 0) {

#		Day
			ip = vp + ix1 - 1
			np = ix1 - 1
			call strcpy (Memc[vp], Memc[tp], np)
			i = 1
			if (ctoi(Memc[tp],i,j) > 0)
			    dd = j
			else
			    dd = 0

#		Month
			ix2 = strldx (slash, Memc[vp])
			if (ix2 < 1)
		            ix2 = strldx (dash, Memc[vp])
			np = ix2 - ix1 - 1
			call strcpy (Memc[ip+1], Memc[tp], np)
			i = 1
			if (ctoi(Memc[tp],i,j) > 0)
			    mm = j
			else
			    mm = 0
	
#		Year
			ip = vp + ix2 - 1
			call strcpy (Memc[ip+1],Memc[tp],4)
			i = 1
			if (ctoi(Memc[tp],i,j) > 0) {
			    yyyy = j
			    if (yyyy < 20)
				yyyy = 2000 + yyyy
			    else if (yyyy < 120)
				yyyy = 1900 + yyyy
			    }
			else
			    yyyy = 0
			}
		    }
		}
	    call sfree (sp)
	    }
	return
end                                     


# IMGDTIM -- Extract UT in hours from FITS image header string
#            as ISO yyyy-mm-ddThh:mm:ss.sss

procedure imgdtim (im, key, ut)

pointer	im		# image descriptor
char	key[ARB]	# parameter from which time will be extracted
double	ut		# UT in hours returned

int	lstr, ix1, ip
int	strlen(), stridx()
int	imaccf()
char	tee
pointer	vp, sp
double	dtemp
bool	debug
int	ctod()

begin
	debug = FALSE
#	debug = TRUE
	if (imaccf (im, key) == YES) {
	    call smark (sp)
	    call salloc (vp, SZ_LINE, TY_CHAR)

	    call imgstr (im, key, Memc[vp], SZ_LINE)

#	String must be more than 9 characters long for date plus T plus time
	    lstr = strlen (Memc[vp])
	    if (lstr > 9) {
		tee = 'T'
		ix1 = stridx (tee, Memc[vp])
		if (ix1 > 0) {
		    vp = vp + ix1
		    if (debug) {
			call printf ("IMGDTIM: %s UT is %s\n")
			    call pargstr (key)
			    call pargstr (Memc[vp])
			call flush (STDOUT)
			}
		    ip = 1
		    if (ctod (Memc[vp], ip, dtemp) == 0) {
			call printf ("IMGDTIM:  Bad header data %s = %s\n")
			    call pargstr (key)
			    call pargstr (Memc[vp])
			call flush (STDOUT)
			}
		    else {
			ut = dtemp
			}
		    }
		else if (debug) {
		    call printf ("IMGDTIM: keyword %s value has no time\n")
			call pargstr (key)
		    call flush (STDOUT)
		    }
		}
	    else if (debug) {
		call printf ("IMGDTIM: keyword %s value too short\n")
		    call pargstr (key)
		call flush (STDOUT)
		}
	    call sfree (sp)
	    }
	else if (debug) {
	    call printf ("IMGDTIM: keyword %s not found\n")
		call pargstr (key)
	    call flush (STDOUT)
	    }
	return
end

# IMGXPAR -- Delete single or multi-line image header parameter

procedure imgxpar (im, key)

pointer	im		# image descriptor
char	key[ARB]	# parameter or root to be deleted

int	imaccf()
char	keyword[16]
int	i
errchk	syserrs

begin
	if (imaccf (im, key) == YES) {
	    call imdelf (im, key)
#	    call printf ("IMGXPAR: keyword %s deleted\n")
#		call pargstr (key)
	    call flush (STDOUT)
	    }
	else {
	    call strcpy (key, keyword, 16)
	    call strcat ("_001", keyword, 16)
	    if (imaccf (im, keyword) == YES) {
		do i = 1, 1000 {
		    call sprintf (keyword, 16, "%s_%03d")
			call pargstr (key)
			call pargi (i)
		    if (imaccf (im, keyword) == YES) {
			call imdelf (im, keyword)
#			call printf ("IMGXPAR: keyword %s deleted\n")
#			    call pargstr (keyword)
			call flush (STDOUT)
			}
		    else
			break
		    }
		}
	    else {
#		call printf ("IMGXPAR: keyword %s not found\n")
#		    call pargstr (key)
		}
	    }
	return
end


# Jan 10 1995	Print error messages and return default values
# Apr 10 1995	Return null string from imgcpar if error
# Apr 10 1995	Return input value unchanged from other img$par if error
# Jul  7 1995	Use IMACCF instead of IDB_FINDRECORD for generality

# Sep 25 1997	Make IMGDATE accept ISO as well as FITS date format
# Oct  6 1997	Fix bug when there is no date in the header

# Apr 28 2000	Change years from 0-120, not just 0-99, to 1900+year

# Mar 19 2001	Add imgxpar() to delete single- or multi-line keywords

# May 22 2003	If year in old FITS date is < 20, it probably should be >2000

# Jan 30 2006	Add imgdtim() to extract time part of FITS ISO dates
# Feb 16 2006	Fix bug in imgdtim() declaring ctod()

# Jun 14 2007	Print keyword AND bad value if error
