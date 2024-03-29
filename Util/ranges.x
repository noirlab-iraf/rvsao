# File Util/ranges.x
# Modified September 28, 2004 by Doug Mink, Center for Astrophysics
# Copyright(c) 1986 Association of Universities for Research in Astronomy Inc.

include	<mach.h>
include	<ctype.h>

define	FIRST	0		# Default starting range
define	LAST	MAX_INT		# Default ending range
define	STEP	1		# Default step

# DECODE_RANGES -- Parse a string containing a list of integer numbers or
# ranges, delimited by either spaces or commas.  Return as output a list
# of ranges defining a list of numbers, and the count of list numbers.
# Range limits must be positive nonnegative integers.  ERR is returned as
# the function value if a conversion error occurs.  The list of ranges is
# delimited by a single NULL.

int procedure decode_ranges (range_string, ranges, max_ranges, nvalues)

char	range_string[ARB]	# Range string to be decoded
int	ranges[3, max_ranges]	# Range array
int	max_ranges		# Maximum number of ranges
int	nvalues			# The number of values in the ranges

int	i,ip, nrange, first, last, step, ctoi()

begin
	ip = 1
	nvalues = 0

#	call printf ("DECODE_RANGES: string is %s\n")
#	    call pargstr (range_string)

	do nrange = 1, max_ranges-1 {

	# Defaults to all positive integers
	    first = FIRST
	    last = LAST
	    step = STEP

	# Skip delimiters to start of range
	    while (IS_WHITE(range_string[ip]) || range_string[ip] == ',')
		ip = ip + 1

	# Get first limit.
	    # Must be a number, '-', ':', 'x', or EOS.  If not return ERR.
	    if (range_string[ip] == EOS) {			# end of list
		if (nrange == 1) {
		    # Null string defaults
		    ranges[1, 1] = first
		    if (first < 1)
			ranges[2, 1] = first
		    else
			ranges[2, 1] = last
		    ranges[3, 1] = step
		    ranges[3, 2] = NULL
	    	    nvalues = nvalues + (abs (ranges[2,1]-ranges[1,1])/step) + 1
		    return (OK)
		    }
		else {
		    ranges[3, nrange] = NULL
		    return (OK)
		    }
		}
	    else if (range_string[ip] == '-')
		;
	    else if (range_string[ip] == ':')
		;
	    else if (range_string[ip] == 'x')
		;
	    else if (IS_DIGIT(range_string[ip]))
		i = ctoi (range_string, ip, first)
	    else
		return (ERR)

	# Skip delimiters
	    while (IS_WHITE(range_string[ip]) || range_string[ip] == ',')
		ip = ip + 1

#	call printf ("last: %d %s\n")
#	call pargi (ip)
#	call pargstr (range_string[ip])

	# Get last limit
	    # Must be '-', or 'x' otherwise last = first.
	    if (range_string[ip] == 'x')
		;
	    else if (range_string[ip] == '-' || range_string[ip] == ':') {
		ip = ip + 1
	        while (IS_WHITE(range_string[ip]) || range_string[ip] == ',')
		    ip = ip + 1
		if (range_string[ip] == EOS)
		    ;
		else if (IS_DIGIT(range_string[ip]))
		    i = ctoi (range_string, ip, last)
		else if (range_string[ip] == 'x')
		    ;
		else
		    return (ERR)
		}
	    else
		last = first

	# Skip delimiters
	    while (IS_WHITE(range_string[ip]) || range_string[ip] == ',')
		ip = ip + 1

	# Get step.
	    # Must be 'x' or assume default step.
	    if (range_string[ip] == 'x') {
		ip = ip + 1
	        while (IS_WHITE(range_string[ip]) || range_string[ip] == ',')
		    ip = ip + 1
		if (range_string[ip] == EOS)
		    ;
		else if (IS_DIGIT(range_string[ip])) {
		    if (ctoi (range_string, ip, step) == 0)
		        ;
		} else if (range_string[ip] == '-')
		    ;
		else
		    return (ERR)
	    }

	# Output the range triple.
	    ranges[1, nrange] = first
	    ranges[2, nrange] = last
	    ranges[3, nrange] = step
	    nvalues = nvalues + abs (last-first) / step + 1
#	    call printf ("DECODE_RANGES: %d %d %d %d %d\n")
#		call pargi (nrange)
#		call pargi (first)
#		call pargi (last)
#		call pargi (step)
#		call pargi (nvalues)
	    }

	return (ERR)					# ran out of space
end


# GET_NEXT_NUMBER -- Given a list of ranges and the current file number,
# find and return the next file number.  Selection is done in such a way
# that list numbers are always returned in monotonically increasing order,
# regardless of the order in which the ranges are given.  Duplicate entries
# are ignored.  EOF is returned at the end of the list.

int procedure get_next_number (ranges, number)

int	ranges[ARB]		# Range array
int	number			# Both input and output parameter

int	ip, first, last, step, next_number, remainder

begin
	# If number+1 is anywhere in the list, that is the next number,
	# otherwise the next number is the smallest number in the list which
	# is greater than number+1.

	number = number + 1
	next_number = MAX_INT

	for (ip=1;  ranges[ip+2] != NULL;  ip=ip+3) {
	    first = min (ranges[ip], ranges[ip+1])
	    last = max (ranges[ip], ranges[ip+1])
	    step = ranges[ip+2]
	    if (number >= first && number <= last) {
		remainder = mod (number - first, step)
		if (remainder == 0)
		    return (number)
		if (number - remainder + step <= last)
		    next_number = number - remainder + step
	    } else if (first > number)
		next_number = min (next_number, first)
	}

	if (next_number == MAX_INT)
	    return (EOF)
	else {
	    number = next_number
	    return (number)
	}
end


# GET_PREVIOUS_NUMBER -- Given a list of ranges and the current file number,
# find and return the previous file number.  Selection is done in such a way
# that list numbers are always returned in monotonically decreasing order,
# regardless of the order in which the ranges are given.  Duplicate entries
# are ignored.  EOF is returned at the end of the list.

int procedure get_previous_number (ranges, number)

int	ranges[ARB]		# Range array
int	number			# Both input and output parameter

int	ip, first, last, step, next_number, remainder

begin
	# If number-1 is anywhere in the list, that is the previous number,
	# otherwise the previous number is the largest number in the list which
	# is less than number-1.

	number = number - 1
	next_number = 0

	for (ip=1;  ranges[ip+2] != NULL;  ip=ip+3) {
	    first = min (ranges[ip], ranges[ip+1])
	    last = max (ranges[ip], ranges[ip+1])
	    step = ranges[ip+2]
	    if (number >= first && number <= last) {
		remainder = mod (number - first, step)
		if (remainder == 0)
		    return (number)
		if (number - remainder >= first)
		    next_number = number - remainder
	    } else if (last < number) {
		remainder = mod (last - first, step)
		if (remainder == 0)
		    next_number = max (next_number, last)
		else if (last - remainder >= first)
		    next_number = max (next_number, last - remainder)
	    }
	}

	if (next_number == 0)
	    return (EOF)
	else {
	    number = next_number
	    return (number)
	}
end


# IS_IN_RANGE -- Test number to see if it is in range.

bool procedure is_in_range (ranges, number)

int	ranges[ARB]		# Range array
int	number			# Number to be tested against ranges

int	ip, first, last, step

begin
	for (ip=1;  ranges[ip+2] != NULL;  ip=ip+3) {
	    first = min (ranges[ip], ranges[ip+1])
	    last = max (ranges[ip], ranges[ip+1])
	    step = ranges[ip+2]
	    if (number >= first && number <= last)
		if (mod (number - first, step) == 0)
		    return (true)
	}

	return (false)
end

# May 22 1992	Allow first number in range to be zero

# Apr  9 1997	If first is 0 and no last, last is also 0

# Sep 28 2004	Allow : as well as - as separator
