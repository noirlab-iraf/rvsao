# File rvsao/Util/t_sindex.x
# February 6, 1995
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# SINDEX -- Use Fortran Index call to find one string in another string

procedure t_sindex()

char	istrng1[SZ_LINE]
char	istrng2[SZ_LINE]
char	fstrng1[SZ_LINE]
char	fstrng2[SZ_LINE]
int	loc

int	index

begin
	call clgstr ("string1",istrng1,SZ_LINE)
	call clgstr ("string2",istrng2,SZ_LINE)

	loc = index (fstrng1,fstrng2)
	CALL clputi ("location",loc)
end
