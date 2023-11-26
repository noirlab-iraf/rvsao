# File eqw.h
# May 14, 2009

# Band structure
define	LEN_BAND	9			# length of structure
define	BAND_ID		Memi[$1]		# ptr to band id string
define	BAND_FILTER	Memi[$1+1]		# ptr to filter string
define	BAND_W1		Memd[P2D($1+2)]		# lower wavelength limit
define	BAND_W2		Memd[P2D($1+4)]		# upper wavelength limit
define	BAND_FN		Memi[$1+6]		# no. of filter points
define	BAND_FW		Memi[$1+7]		# ptr to filter wavelengths
define	BAND_FR		Memi[$1+8]		# ptr to filter responses

# Multiple bands for indices and equivalent widths.
define	NBANDS		3			# maximum number of bands
define	BAND1		1
define	BAND2		2
define	BAND3		3
define	BAND		Memi[$1+($2-1)*NBANDS+($3-1)]

# May 14 2009	Change CW and DW to W1 and W2
