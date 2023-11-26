# File rvsao/Subpix/t_subpix.x
# June 9, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 2008 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
  
# SUBPIX is an IRAF task for obtaining sub-pixel relative radial velocities
# using an algorithm in "Fundamental photon noise limit to radial velocity
# measurements," by F. Bouchy, F. Pepe, and D. Queloz, A&A 374, 733-739 (2001)

# Information is shared in common blocks defined in "rvsao.com".

include <imhdr.h> 
include <imio.h>
include <fset.h>
include <smw.h>
include "rvsao.h"
include "xcv.h"
include "emv.h"

procedure t_subpix ()

int	i
char	specfile[SZ_PATHNAME]	# Object spectrum file name
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	tempfiles[SZ_PATHNAME]	# List of template spectra
char	specdir[SZ_PATHNAME]	# Directory for object spectra
#char	filt_type[SZ_LINE]	# Filter for transform
				# (ramp | welch | hanning | cos-bell)
char	svel_corr[SZ_LINE]	# Type of velocity correction for spectrum
				# (none | file | heliocentric | barycentric)
char	tvel_corr[SZ_LINE]	# Type of velocity correction for template
				# (none | file | heliocentric | barycentric)
bool	savevel0		# Save velocity, error, and R in data file header
int	rmode		# Report format (1=normal,2=one-line)
int	logfiles	# List of log files
char	logfile[SZ_PATHNAME] # Log file name
char	wtitle[20]	# Title for wavelength plots of spectrum

int	mspec		# Object aperture to read from multispec file
int	mband		# Object band to read from multispec file

pointer	speclist	# List of spectrum files
char	str[SZ_LINE]
int	fd
char	vel_plot[SZ_LINE]	# type of velocity for redshifting plot
				# correlation|emission|combination|search

int	nmspec0		# Number of object multispec spectra
int	mspec_range[3,MAX_RANGES]
int	ip,jp,lfile	# Limits for multispec aperture decoding
char	lbracket[3]	# "[({"
char	rbracket[3]	# "])}"
bool	echelle		# TRUE if echelle parameter is set
double	sumvel, sumerr, sumr, avgvel, avgerr, avgr, dnap
 
bool	clgetb()
int	clpopnu(), clgeti(), clgfil(), open()
int	strdic(), stridx(), stridxs()
real	clgetr()
int	decode_ranges(),get_next_number()
char	vel_init[SZ_LINE]	# type of velocity for initial value
int	imtgetim(), imaccess(), strlen(), ldir, clscan()
pointer	imtopenp()

define	newspec_ 10
define	newap_	 20
define	endxc_	 90

include	"rvsao.com"
include "results.com"
include	"emv.com"
include	"xcor.com"
include	"xcorf.com"
include	"xplt.com"
 
begin
	c0 = 299792.5
	qplot = FALSE
	nfound = 0
	maxpix = 0
	maxpts4 = 0
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	call sprintf (wtitle,20,"Wavelength")
	ntmp = 0
	xcor = NULL
	xvel = NULL
	shspec = NULL
	shtemp = NULL
	wltemp = NULL
	xind = NULL
	xifft = NULL
	ft1 = NULL
	ft2 = NULL
	ftcfn = NULL
	tft = NULL
	pft = NULL
	spexp = NULL
	xcont = NULL
	waverest = 0.d0
	specref = 0

# Get task parameters.

# Spectra for which to find shifts
	speclist = imtopenp ("spectra")

# Multispec spectrum numbers (use only first if multiple files)
	call clgstr ("specnum",specnums,SZ_LINE)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR){
	    call sprintf (str, SZ_LINE, "T_XCSAO: Illegal multispec list <%s>")
		call pargstr (specnums)
	    call error (1, str)
	    }
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }
	mband = clgeti ("specband")
	echelle = clgetb ("echelle")

# Templates against which to correlate spectra
	call clgstr ("templates",tempfiles,SZ_PATHNAME)

# Optional intermediate data plot switches
	pltspec = clgetb ("obj_plot")
	plttemp = clgetb ("temp_plot")

# Print processing information
	debug  = clgetb ("debug")

# Number of times to smooth (1-2-1) final data plot
	nsmooth = clgeti ("nsmooth")

# Velocity center and width of summary page shift graph
	xcr0 = clgetr ("cvel")
	xcrdif = clgetr ("dvel")

# Type of heliocentric velocity correction to be used
	call clgstr ("svel_corr",svel_corr,SZ_LINE)
	svcor = strdic (svel_corr,svel_corr,SZ_LINE, HC_VTYPES)
	call clgstr ("tvel_corr",tvel_corr,SZ_LINE)
	tvcor = strdic (tvel_corr,tvel_corr,SZ_LINE, HC_VTYPES)

# Type of velocity for initial redshift
        call clgstr ("vel_init",vel_init,SZ_LINE)
        vinit = strdic (vel_init,vel_init,SZ_LINE,XC_VTYPES)

# Image header result flag
	savevel0 = FALSE
	savevel0 = clgetb ("save_vel")

# Report mode for log file
	rmode = 1
	rmode = clgeti ("report_mode")

# Initialize emission and absorption lines for labelling
	call eminit (FALSE)
 
# Open log files and write a header.
	logfiles = clpopnu ("logfiles")
	call fseti (STDOUT, F_FLUSHNL, YES)
	i = 0
	call strcpy ("rvsao.xcsao",taskname,SZ_LINE)
	while (clgfil (logfiles, logfile, SZ_PATHNAME) != EOF) {
	    fd = open (logfile, APPEND, TEXT_FILE)
	    if (fd == ERR) break
	    if (rmode == 1) {
		call loghead (taskname,str)
		call fprintf (fd, "%s\n")
		    call pargstr (str)
		}
	    i = i + 1
	    logfd[i] = fd
	    }
	nlogfd = i
	call clpcls (logfiles)

# Type of velocity for plotting emission and absorption lines
	if (clscan("vel_plot") != EOF) {
	    call clgstr ("vel_plot",vel_plot,SZ_LINE)
	    vplot = strdic (vel_plot,vel_plot,SZ_LINE,PL_VTYPES)
	    }
	else
	    vplot = VCORREL

# Print tab table headers
	if (rmode < 0)
	    call xcrshead (rmode)

# Get next object spectrum file name from the list
newspec_
	if (imtgetim (speclist, specfile, SZ_PATHNAME) == EOF)
	   go to endxc_

# Check for specified apertures in multispec spectrum file
	ip = stridxs (lbracket,specfile)
	if (ip > 0) {
	    lfile = strlen (specfile)
	    specfile[ip] = EOS
	    jp = 0
	    ip = ip + 1
	    while (stridx (specfile[ip],rbracket) == 0 && ip <= lfile) {
		jp = jp + 1
		specnums[jp] = specfile[ip]
		specfile[ip] = EOS
		ip = ip + 1
		}
	    if (jp > 0)
		specnums[jp+1] = EOS
	    else
		call strcpy ("0",specnums,SZ_LINE)
	    if (decode_ranges (specnums,mspec_range,MAX_RANGES,nmspec) == ERR){
		call sprintf (str, SZ_LINE, "T_XCSAO: Illegal multispec list <%s>")
		    call pargstr (specnums)
		call error (1, str)
		}
	    }
	else
	    nmspec = nmspec0
	if (debug) {
	    call printf ("XCSAO: next file is %s [%s] = %d aps\n")
		call pargstr (specfile)
		call pargstr (specnums)
		call pargi (nmspec)
	    }

# Check for readability of object spectrum
	call strcpy (specdir,specpath,SZ_PATHNAME)
	call strcat (specfile,specpath,SZ_PATHNAME)
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("XCSAO: cannot read spectrum file %s \n")
		call pargstr (specpath)
	    go to newspec_
	    }

# Get next multispec number from list
	mspec = -1
	sumvel = 0.d0
	sumerr = 0.d0
	sumr = 0.d0
	dnap = 0.d0
newap_
	savevel = savevel0
	if (nmspec <= 0)
	    go to newspec_
	if (get_next_number (mspec_range, mspec) == EOF)
	    go to newspec_

	call pixshift (specfile, specdir, mspec, mband)

	if (echelle) {
	    sumvel = sumvel + zvel[itmax]
	    sumerr = sumerr + czerr[itmax]
	    sumr = sumr + czr[itmax]
	    dnap = dnap + 1.d0
	    }

# Move on to next aperture or next image
	nmspec = nmspec - 1
	if (nmspec > 0)
	    go to newap_
	go to newspec_
 
# Close the log files
endxc_	do i = 1, nlogfd {
	    call close (logfd[i])
	    }                                              

	if (debug) {
	    call printf ("Best Template(%d) %s: %.3f (%.3f) R = %.4f\n")
		call pargi (itmax)
		call pargstr (tempid[1,itmax])
		call pargd (zvel[itmax])
		call pargd (czerr[itmax])
		call pargd (czr[itmax])
	    }

# Save some results in parameter file
	call clpstr ("besttemp", tempid[1,itmax])
	if (echelle && dnap > 0.d0) {
	    avgvel = sumvel / dnap
	    call clputd ("velocity", avgvel)
	    avgerr = sumerr / dnap
	    call clputd ("velerr", avgerr)
	    }
	else {
	    call clputd ("velocity", zvel[itmax])
	    call clputd ("velerr", czerr[itmax])
	    }

#  Close spectrum list
	call imtclose (speclist)

# Free processing vectors allocated in xcfit()
	call mfree (xcor, TY_REAL)
	call mfree (xvel, TY_REAL)
	call mfree (shspec, TY_REAL)
	call mfree (shtemp, TY_REAL)
	call mfree (wltemp, TY_REAL)

# Free processing vectors allocated in xcorfit()
	call mfree (xind, TY_REAL)
	call mfree (xifft, TY_REAL)
	call mfree (pft, TY_REAL)
	call mfree (tft, TY_COMPLEX)
	call mfree (ftcfn, TY_COMPLEX)
	call mfree (ft1, TY_COMPLEX)
	call mfree (ft2, TY_COMPLEX)
	call mfree (spexp, TY_REAL)
	call mfree (xcont, TY_REAL)

# Free plotting vectors allocated in xcplot
	if (maxpix > 0) {
	    call mfree (scont, TY_REAL)
	    call mfree (smspec, TY_REAL)
	    call mfree (cspec, TY_REAL)
	    call mfree (smcspec, TY_REAL)
	    }

# Free plotting vectors allocated in xcorplot
	if (maxpts4 > 0) {
	    call mfree (xlev, TY_REAL)
	    call mfree (fraclev, TY_REAL)
	    }
end
 
# Jun  9 2008	New task based on pxcsao to compute sub-pixel shifts
