# File rvsao/Emvel/emarch.x
# November 16, 1994
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1994 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Produce an archive record with the results of emission line velocity fits

#  Input is in common in "fquot.com" and emv.com

include "rvsao.h"
include "emv.h"

procedure emarch (specim, specfile)

pointer	specim		# Spectrum file header
char	specfile[ARB]	# Name of archive file

char	header[48]
char	arcfile[SZ_FNAME]
char	osheader[48]
double	w2p(), pixel1, pixel2, center, width, wcs_p2w()
double	twt, velline, wpp, velerr, zline, bcz
int	il, nw, nh
int	open(), strlen()
int	arc_fd
int	nbrec,lhead
int	i

include "rvsao.com"
include "emlines.com"
include "ansum.com"
include "emv.com"

begin

#  Set up analysis summary output file name
	call strcpy (specfile,arcfile,SZ_FNAME)
	call strcat (".ansum",arcfile,SZ_FNAME)
	if (debug) {
	    call printf ("EMARCH:  writing %s\n")
		call pargstr (arcfile)
	    }
	iferr {arc_fd = open (arcfile, NEW_FILE, BINARY_FILE)} then {
	    call printf ("EMARCH:  Cannot write %s\n")
		call pargstr (arcfile)
	    return
	    }

#  Write analysis summary record header
	do i = 1, 48, 8 {
	    call strcpy ("        ", header[i], 8)
	    }
	call strcpy ("ANALYSIS_SUMMARY 40 ",header,20)
	if (debug) {
	    call printf ("EMARCH: header is %s\n")
		call pargstr (header)
	    }
	nh = 48
	call strpak (header,osheader,nh)
	nw = 24
	call write (arc_fd, osheader, nw)

#  Write analysis summary record
	qcstat = spvqual
	nw = 2
	call write (arc_fd, qcstat, nw)
	czxc = spxvel
	czxcerr = spxerr
	czxcr = spxr
	cz0 = spvel
	cz0err = sperr
	czem = spevel
	czemerr = speerr
	czemscat = 0.
	nw = 18
	call write (arc_fd, cz0, nw)
	call close (arc_fd)

#  Set up emission line output file name
	call strcpy (specfile,arcfile,SZ_FNAME)
	call strcat (".emlines",arcfile,SZ_FNAME)
	if (debug) {
	    call printf ("EMARCH:  writing %s\n")
		call pargstr (arcfile)
	    }
	iferr {arc_fd = open (arcfile, NEW_FILE, BINARY_FILE)} then {
	    call printf ("EMARCH:  Cannot write %s\n")
		call pargstr (arcfile)
	    return
	    }

#  Write emission line record header
	nbrec = 4 + (76 * nfound)
	call sprintf (header, 48, "EMISSION_LINES %d ")
	    call pargi (nbrec)
	if (debug) {
	    call printf ("EMARCH: header is %s\n")
		call pargstr (header)
	    }
	lhead = strlen (header)
	do i = lhead+1, 48 {
	    header[i] = ' '
	    }
	nh = 48
	call strpak (header,osheader,nh)
	nw = 24
	call write (arc_fd, osheader, nw)

#  Write number of emission lines found and fit
	nw = 2
	nfnd = nfound
	nft = nfit
	call write (arc_fd, nfnd, nw)

#  Set up wavelength to pixel conversion
	call wpinit (specim)

#  Compute total weight
	twt = 0.d0
        do il = 1, nfound {
	    twt = twt + emparams[10,1,il]
	    }

#  Save emission line parameters for each line found
	do il = 1, nfound {

	    # Rest wavelength
	    lrest = real (wlrest[il])

	    # Line center in pixels
	    lcent = real (w2p (wlobs[il]))

	    # Line height
	    lhght = real (emparams[5,1,il])

	    # Line center and width in wavelength
	    center = wcs_p2w (emparams[4,1,il])
	    width = wcs_p2w (emparams[4,1,il] + emparams[6,1,il]) -
	            wcs_p2w (emparams[4,1,il] - emparams[6,1,il])

	    # Line width in reticon pixels
	    pixel1 = w2p (center - width)
	    pixel2 = w2p (center + width)
	    if (pixel1 < pixel2)
		lwidth = real (pixel2 - pixel1)
	    else
		lwidth = real (pixel1 - pixel2)

	    # Continuum level and slope
	    lcont = real (emparams[1,1,il])
	    lslope = real (emparams[2,1,il])

	    # Emission line Gaussian coefficients
	    lgcent = real (w2p (center))
	    lcente = emparams[4,2,il]
	    lghght = emparams[5,1,il]
	    lhghte = emparams[5,2,il]
	    lgwidth = lwidth
	    lwidthe = emparams[6,2,il]

	    # Continuum polynomial coefficients
	    lcfit[1] = emparams[1,1,il]
	    lcfit[2] = emparams[2,1,il]
	    lcfit[3] = emparams[3,1,il]

	    # Equivalent width
	    leqw = emparams[7,1,il]
	    leqwe = emparams[7,2,il]

	    # Fit chi^2 and degrees of freedom
	    lchi2 = emparams[8,1,il]
	    ldegf = emparams[8,2,il]

	    # Lime weight in velocity fit
	    if (emparams[10,1,il] > 0)
		lwt = 10000.d0 * emparams[10,1,il] / twt
	    else
		lwt = 0

	    if (debug) {
		c0 = 299792.5d0
		bcz = 1.d0 + (spechcv / c0)
		zline = (wlobs[il] / wlrest[il]) / bcz
		velline = c0 * (zline - 1.d0)
		If (pixel1 != pixel2)
		    wpp = (2.d0 * width) / (pixel2 - pixel1)
		else {
		    pixel1 = w2p (center - 5.d0)
		    pixel2 = w2p (center + 5.d0)
		    wpp = 10.d0 / (pixel2 - pixel1)
		    }
		velerr = c0 * 2.d0 * emparams[4,2,il] * wpp / wlrest[il]
		call printf ("%8.2f %8.2f %7.4f %9.2f %7.2f\n")
		    call pargr (lrest)
		    call pargr (lcent)
		    call pargd (zline)
		    call pargd (velline)
		    call pargd (velerr)
		}
	    nw = 38
	    call write (arc_fd, lrest, nw)
	    }
	call close (arc_fd)
end

# Aug 12 1992	Convert wavelength to pixels for archive
# Nov 19 1992	Add continuum slope

# Aug 12 1993	Fix for mwcsg

# Apr 21 1994	Call WRITE as procedure rather than function
# Apr 25 1994	Change for loop to do loop
# Jun 23 1994	Pass velocities in fquot, not getim labelled common
# Aug  3 1994	Change common and header from fquot to rvsao
# Nov 16 1994	Set QCSTAT from SPVQUAL in rvsao common
