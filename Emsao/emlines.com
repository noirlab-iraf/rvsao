
#  Emission line record

common/emlin/ nfnd,nft

short nfnd		# Number of emission lines found
short nft		# Number of emission lines used in Cz fit

#  Emission line structure (one per line)

common/linrec/ lrest,lcent,lhght,lwidth,lcont,lslope,lgcent,lcente,
		lghght,lhghte,lgwidth,lwidthe,lcfit,leqw,leqwe,lchi2,
		ldegf,lwt

real lrest		# Rest wavelength in Angstroms
real lcent		# Center of line in pixels
real lhght		# Height of emission line
real lwidth		# Width of emission line
real lcont		# Continuum level
real lslope		# Slope of continuum
real lgcent		# Line center from Gaussian fit
real lcente		# Line center error from Gaussian fit
real lghght		# Line height from Gaussian fit
real lhghte		# Line height error from Gaussian fit
real lgwidth		# Line width from Gaussian fit
real lwidthe		# Line width error from Gaussian fit
real lcfit[3]		# Parabolic fit to continuum
real leqw		# Equivalent width in Angstroms
real leqwe		# Equivalent width error in Angstroms
real lchi2		# Chi ** 2
short ldegf		# Degrees of freedom
short lwt		# Weight in average
