# Save velocity results in image header if savevel is true
# Save velocity results in parameter file if savevel is false

procedure printsum (specim, nap, dwdp)

pointer specim	# Object image header structure
int	nap	# Number of apertures for which shift has been computed
double	dwdp	# Wavelength per pixel

int	ip,jp,lfile	# Limits for multispec aperture decoding

char	str[SZ_LINE]
int	i,j,iline,ldir
int	fd
int	nsp, nflines, nap
double	avgvel, avgerr, medvel, avgpix, medpix, q1, q2
double	vl
 
include	"rvsao.com"
include	"emv.com"
include "contin.com"
include "results.com"
include "emsums.com"
 
begin
	c0 = 299792.5d0

# Save global results in the parameter or image file
	if (debug) {
	    if (savevel)
		call printf("PRINTSUMS: About to save results in header\n")
	    else
		call printf("PRINTSUMS: About to save results as parameters\n")
	    call flush (STDOUT)
	    }

# Mean velocity and error
	if (nap > 0)
	    avgvel = sumvel / double (nap)
	else
	    avgvel = 0.d0
	if (debug) {
	    call printf("PEMSAO: Mean velocity = %.3f saved as meanvel\n")
		call pargd (avgvel)
	    call flush (STDOUT)
	    }
	if (nap > 0)
	    avgerr = dsqrt (sumerr) / double (nap)
	else
	    avgerr = 0.d0
	if (debug) {
	    call printf("PEMSAO: Mean velocity error = %.3f saved as meanerr\n")
		call pargd (avgerr)
	    call flush (STDOUT)
	    }

# Median velocity and error
	call median (Memd[evel], nap, medvel, q1, q2)
	if (debug) {
	    call printf("PEMSAO: Median velocity = %.3f saved as medvel\n")
		call pargd (medvel)
	    call printf("PEMSAO: First quartile = %.3f saved as medq1\n")
		call pargd (q1)
	    call printf("PEMSAO: Second quartile = %.3f saved as medq2\n")
		    call pargd (q2)
	    call flush (STDOUT)
	    }

# If not saving in header, save in parameter file
	if (!savevel) {
	    call clputd ("meanvel", avgvel)
	    call clputd ("meanerr", avgerr)
	    call clputd ("medvel",medvel)
	    call clputd ("medq1",q1)
	    call clputd ("medq2",q2)
	    }

# Mean pixel shift
	if (nap > 0)
	    avgpix = sumpix / double (nap)
	else
	    avgpix = 0.d0
	if (debug) {
	    call printf("PEMSAO: Mean pixel shift = %.3f saved as meanpix\n")
		call pargd (avgpix)
	    call flush (STDOUT)
	    }

# Median pixel shift
	call median (Memd[epix], nap, medpix, q1, q2)
	if (debug) {
	    call printf("PEMSAO: Median pixel shift = %.3f saved as medpix\n")
		call pargd (medpix)
	    call flush (STDOUT)
	    }

	if (savevel) {
	    call imaddd (specim, "MEANPIX", avgpix, "Mean pixel shift")
	    call imaddd (specim, "MEDPIX", medpix, "Median pixel shift")
	    }
	else {
	    call clputd ("meanpix", avgpix)
	    call clputd ("medpix",medpix)
	    }

# Save some results for first line in the parameter file
	if (debug) {
	    call printf("PEMSAO: Number of lines fit = %d saved\n")
		call pargi (nflines)
	    call printf("PEMSAO: First line fit = %s saved as emline\n")
		call pargstr (nmobs[1,1])
	    call printf("PEMSAO: Rest wavelength = %.3f saved as wlrest\n")
		call pargd (wlrest[1])
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputi ("nlfit", nflines)
            call clpstr ("emline", nmobs[1,1])
            call clputd ("wlrest", wlrest[1])
	    }

# Line velocity
	linvel = (c0 * emparams[9,1,1]) + spechcv
	if (debug) {
	    call printf("PEMSAO: Velocity = %.4f saved as velocity\n")
		call pargd (linvel)
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputd ("velocity", linvel)
	    }

# Line velocity error
	linerr = c0 * emparams[9,2,1]
	if (debug) {
	    call printf("PEMSAO: Velocity error = %.4f saved as velerr\n")
		call pargd (linerr)
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputd ("velerr", linerr)
	    }

# Line height
	if (debug) {
	    call printf("PEMSAO: Line height = %.2f saved as lineheight\n")
		call pargd (emparams[5,1,1])
	    call flush (STDOUT)
	    }
	if (!savevel) {
	    call clputd ("lineheight", emparams[5,1,1])
	    }

# Line width
	vl = emparams[6,1,1] * abs (dwdp)
	if (debug) {
	    call printf("PEMSAO: Line width = %.4f saved as linewidth\n")
		call pargd (vl)
	    call printf("PEMSAO: Line equivalent width = %.4f saved as lineeqw\n")
		call pargd (emparams[7,1,1])
	    call flush (STDOUT)
	    }
	if (!savevel) {
            call clputd ("linewidth", vl)
	    call clputd ("lineeqw", emparams[7,1,1])
	    }

	return
end
