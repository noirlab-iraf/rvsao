# File lib/rvsao.h
# March 15, 2013

# RVSAO 2.8.2 header file

define	VERSION		"2.8.2"
define	C		2.9979e5
define	LN2		0.69314718
define	CLN10		6.90292e5
define	TWOPI		6.283185

define  MIN_PIXEL_VALUE -1000000        # Minimum legal pixel value

define	MAX_RANGES	50
define	MAXTEMPS	2000
define	MAXCACHE	0
define	MAXLOG		5
#define MAXPTS		16384
#define MAXPTS2		8192
define	NPAR		4
define	SZ_HSTRING	68
define	SZ_HKWORD	9
 
define	CZMIN		-500.
define	CZMAX		50000.
define	SIGMIN		25.
define	SIGMAX		750.
define	GAMMIN		0.25
define	GAMMAX		1.50
define	CHIMIN		0.2
define	CHIMAX		10.0
define	ERRMIN		0.0001
define	ERRMAX		10.

#  Velocity correction flags
define	HC_VTYPES	"|none|file|heliocentric|barycentric|hfile|"
define  NONE		1	# no correction
define  FBCV		2	# Barycentric correction from file BCV
define	HCV		3	# Heliocentric correction
define	BCV		4	# Solar system barycentric correction
define  FHCV		5	# Heliocentric correction from file HCV

#  Velocity types
define  VCORREL		1		# Velocity from cross-correlation
define	VEMISS		2		# Velocity from emission line fit
define	VCOMB		3		# Velocity from emission and correlation
define	VSEARCH		4		# Velocity guessed from spectrum
define	VGUESS		5		# Velocity guess from parameter file
define	ZGUESS		6		# Velocity/c guess from parameter file
define	VCORTEMP	7		# Velocity from correlation of template

#  Velocity plotting flags
define	PL_VTYPES	"|correlation|emission|combination|search|guess|"

#  Initial velocity flags
define	IEM_VTYPES	"|correlation|emission|combination|search|guess|zguess|cortemp|"

# Correlation flags
define  COR_TYPES       "|no|velocity|pixel|wavelength|yes|"
define  COR_NO          1       # Do not cross-correlate
define  COR_VEL         2       # Cross-correlate velocity shift
define  COR_PIX         3       # Cross-correlate pixel shift
define  COR_WAV         4       # Cross-correlate wavelength shift
define  COR_YES         5       # Cross-correlate velocity shift

# Jul  3 1995	Last undocumented revision

# Feb  5 1997	Add velocity plotting and initial velocity flags

# Mar 20 1997	Version 2.0
# Apr 25 1997	Add zguess velocity flag
# May 19 1997	Add vcortemp velocity flat for initial velocities

# Feb 13 1998	Increase maximum number of correlation points to 16,384
# Mar  4 1998	Fix BCVCORR
# Apr  7 1998	Fix SUMSPEC
# May 15 1998	Allow minvel and maxvel to be INDEFD
# Jun 12 1998	2.0b17: Use pixel limits of WCS
# Jul 31 1998	2.0: Add toggling of heading in emsao mode 3 summary graph
# Dec 30 1998	2.1b1: Fix bug in XCSAO template pathname list handling
# Dec 30 1998	2.1b1: Increase maximum templates to 256

# Mar 18 1999	2.1b3: Multiply by 1000 after renormalizing
# Mar 30 1999	2.1b5: Allow use of OBSERVAT keyword for observatory location
# Apr  6 1999	2.1b6: Read templates only once
# May 11 1999	2.1b7: Fix SUMSPEC spectrum stacking
# Jun 10 1999	2.1b9: Add more options to SUMSPEC
# Jun 30 1999	2.1b11: Fix bug using template vel and HCV in XCSAO
# Jul 16 1999	2.1b12: Fix bug so QPLOT updates headers
# Jul 23 1999	2.1b13: Fix null report entries and add EMSAO tab table reports
# Jul 29 1999	2.1b14: Add XCSAO tab table reports; SUMSPEC constant normalize
# Aug 18 1999	2.1b15: SUMSPEC divide continuum no longer subtracts 1
# Sep  7 1999	2.1b16: Set CD matrix in SUMSPEC; fix dependency bugs
# Sep 15 1999	2.1b17: Fix tab table output from XCSAO
# Sep 24 1999	2.1.18: Fix EMSAO and XCSAO to update if quality flag is set
# Dec  6 1999	2.1.19: Free image WCS structure correctly

# Jan 25 2000	2.1.20: Add option to not remove object spectrum continuum
# Mar  6 2000	2.1.21: Fix bug in LINESPEC when using redshift as Z
# Apr 22 2000	2.1.22: Add 1900 to years from 0-120, not just 0-99
# Jul  5 2000	2.1.24: Automate renormalization
# Jul 21 2000	2.1.25: Add BCV adjustment without target velocity to SUMSPEC
# Aug  2 2000	2.1.26: Fix echelle spectrum number tracking
# Aug 10 2000	2.1.27: Fix bug in EMSAO which sorted lines poorly
# Sep 13 2000	Version 2.2: Add wavelength correlation and partial templates

# Feb  9 2001	2.2.1: Fix XCSAO to work with reversed spectra
# Mar 28 2001	2.2.2: Fix EMSAO to work with reversed spectra

# Feb  5 2002	2.2.5: Drop ends of spectra with values below MIN_PIXEL_VALUE
# Mar 29 2002	2.2.6: Add per spectrum wavelength limits to SUMSPEC
# May 29 2002	2.2.7: Add per spectrum velocity shifts to SUMSPEC
# Aug  7 2002	2.2.8: Add report mode 16 to XCSAO for per template wavelenghs
# Sep 17 2002	2.2.9: Allocate large vectors only once
# Sep 30 2002	2.3.0: Add task to compute equivalent widths including redshift

# Jun  2 2003	2.3.2: Add smoothing to SUMSPEC; fix BCVCORR
# Aug  1 2003	2.3.3: Add ability to read first extension of multi-ext. FITS

# Jan 16 2004	2.3.4: Print nonspec disp message only in debug mode in getimage()
# May 13 2004	2.3.5: Return fit from SUMSPEC; add forced emission labels to XCSAO
# Jul 20 2004	2.3.6: Read only first token from filename list files
# Aug  3 2004	2.3.7: Fix bug in SUMSPEC: use contin, not contout for cont rem.
# Sep 28 2004	2.3.8: Add parameters to EQWIDTH so it can be used on skies
# Oct 19 2004	2.3.9: Fix equivalent width computation in EMSAO so units match

# Jan 28 2005	2.4.0: Finally make equivalent width task EQWIDTH usable
# Mar 15 2005	2.4.1: Fix bug in RA setting in BCVCORR and compbcv()
# Mar 23 2005	2.4.2: Add constant graph scaling to SUMSPEC
# May 19 2005	2.4.3: Add pixel shift report mode, 3 decimal places to EMSAO
# Jul 14 2005:	2.4.4: Fix bugs for Linux port
# Jul 28 2005:	2.4.5: Make EQWIDTH output zeroes if band is off spectrum
# Aug 22 2005:	2.4.6: Print Z, ZEM, ZXC in header
# Sep 19 2005:	2.4.7: Fix bug in EQWIDTH flux computation
# Nov  3 2005:	2.4.7: Force renormalization in XCSAO if -1 < spectrum < 1
# Dec 20 2005:	2.4.8: Free all malloc'ed arrays in XCSAO

# Jan 12 2006:	2.4.8: Do not try to write file unless enabled by parameter
# May 31 2006:  2.5.0: Read name and position from APID, if present
# May 31 2006:  2.5.0: Add report mode 17 to xcsao for Hectochelle
# Jul 13 2006:	2.5.1: Do not use per-aperture keywords for apertures > 999
# Aug 17 2006:	2.5.2: Do not use zero-valued pixels at ends of spectrum
# Sep 27 2006:	2.5.3: Fix bug renormalizing very small numbers in XCSAO
# Oct 24 2006:	2.5.4: Increase max. number of templates to 512 from 256

# Jan 17 2007:	2.5.5: Plot all emission lines if any, not just those found
# Feb 14 2007	2.5.6: Preserve APNUMi in SUMSPEC from WCS strings
# Mar 30 2007	2.5.7: SUMSPEC exposure time addition and EMSAO reversed spectra
# Jun 13 2007	2.5.8: Read RA and DEC as strings not numbers from APIDi in getimage.x
# Jun 22 2007	2.5.9: Allocate transform vectors, so maximum number of bins is unfixed
# Aug  2 2007	2.6.0: Increase MAXTEMPS from 512 to 1024; drop MAXPTS, MAXPTS2
# Aug 13 2007	2.6.0: Increase MAXTEMPS to 2000
# Oct 29 2007	2.6.2: Save mean shift, err, R if echelle

# Feb 11 2008	2.6.3: Increase MAXTEMPS to 3000
# Mar  5 2008	2.6.4: Bad lines are in pixels if names start with "pix"
# May 13 2008	2.6.5: Add xplot and eplot to plot spectra without fitting
# May 28 2008	2.6.6: Add error computation and more report modes to EQWIDTH
# Nov 18 2008	2.6.9: Label emission lines with wavelengths in abs not labelled

# Jan  7 2009	2.7.0: Fix bug in SUMSPEC which put back BCV when it should not
# Jul 16 2009	2.7.7: Fix error in EQWIDTH and renormalize divided continuum
# Nov 20 2009	2.7.8: Increase file sizes from 256 to 511 characters

# Mar  1 2010	2.7.9: Write pixel shifts to image header in PEMSAO
# Jul  8 2010	2.7.9: Add FITS extension parameter mext to PXCSAO and PEMSAO

# Jun 27 2011	2.8.0: Use HJD or BJD from header even if no GJD

# May 25 2012	2.8.1: Add maxlength to PXCSAO and xcrslts.x for long file names
