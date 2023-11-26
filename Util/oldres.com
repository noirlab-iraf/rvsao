# Original header parameters that describe the object spectrum
 
common /old/ spvel0, sperr0, spr0, spechcv0, spqual0
 
double  spvel0		# Spectrum velocity
double  sperr0		# Spectrum velocity error
double  spr0		# Spectrum velocity R-value
double  spechcv0	# Spectrum heliocentric velocity
char	spqual0		# Spectrum quality flag
