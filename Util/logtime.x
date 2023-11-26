# File rvsao/Util/logtime.x
# October 1, 1997
# Copyright(c) 1986 Association of Universities for Research in Astronomy Inc.
# Modified from sys/misc/cnvtime.x
# Modified by Doug Mink, Harvard-Smithsonian Center for Astrophysics

include	<time.h>

define	SZ_MONTH		3

# LOGTIME -- Return current time as a string, i.e., " 17-Mar-1982 16:30".
# The maximum length of the returned string is given by the parameter
# SZ_TIME in <time.h>.

procedure logtime (outstr, maxch)

char	outstr[maxch]
int	maxch

int	tm[LEN_TMSTRUCT]	# broken down time structure
string	month	"JanFebMarAprMayJunJulAugSepOctNovDec"
long	ltime			# seconds since 00:00:00 10-Jan-80
long	clktime()

begin
	ltime = 0
	ltime = clktime (ltime)
	call brktime (ltime+30, tm)
	
	call sprintf (outstr, maxch, " %02d-%3.3s-%04d %02d:%02d")
	    call pargi (TM_MDAY(tm))
	    call pargstr (month [(TM_MONTH(tm) - 1) * SZ_MONTH + 1])
	    call pargi (TM_YEAR(tm))
	    call pargi (TM_HOUR(tm))
	    call pargi (TM_MIN(tm))
end

# Oct  1 1997	Print 4-digit year
