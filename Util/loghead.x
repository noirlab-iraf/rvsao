# File rvsao/Utility/loghead.x
# October 1, 1997
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1997 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

include "rvsao.h"
 
# Create string containing IRAF, task, user, and time information
 
procedure loghead (taskname,str)

char	taskname	# Task name
char	str[ARB]	# Heading returned

char	tstr[SZ_LINE]

int	nchars
int	envfind()

begin
	call sprintf (str,SZ_LINE, "IRAF %s %s ")
	    call pargstr (taskname)
	    call pargstr (VERSION)
#	nchars = envfind ("version", tstr, SZ_LINE)
#	call strcat (tstr,str,SZ_LINE)
#	call strcat (" ",str,SZ_LINE)
	nchars = envfind ("userid", tstr, SZ_LINE)
	call strcat (tstr,str,SZ_LINE)
	call strcat ("@",str,SZ_LINE)
	call gethost (tstr, SZ_LINE)
	call strcat (tstr,str,SZ_LINE)
	call strcat (" ",str,SZ_LINE)
	call logtime (tstr, SZ_LINE)
	call strcat (tstr,str,SZ_LINE)

	return
end
# Aug 10 1994	New program

# Mar 14 1995	Drop IRAF version from header

# Oct  1 1997	Print 4-digit year in time stamp
