# analysis summary record

int	qcstat		# Quality controll status
			#   0 - not reviewed
			#   1 - inconclusive velocity determination
			#   2 - insufficient wavelength coverage
			#   3 - incorrect redshift velocity
			#   4 - correct redshift velocity

real	cz0		# Overall radial velocity (Cz)
real	cz0err		# Overall Cz error
real	cz0conf		# Overall Cz confidence
real	czxc		# correlation radial velocity
real	czxcerr		# correlation Cz error
real	czxcr		# correlation R value
real	czem		# emission line radial velocity
real	czemerr		# emission line radial velocity error
real	czemscat	# emission line velocity scatter (chi-square)

common/ansum/ qcstat,cz0,cz0err,cz0conf, czxc,czxcerr,czxcr,
	      czem,czemerr,czemscat

