# File rvsao/Xcor/xcarch.x
# August 25, 2004
# By Doug Mink, Center for Astrophysics

# Copyright(c) 1994-2004 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Produce an archive record with the results of cross-correlations
#  of templates and the object spectrum.

#  Most input is in common in "rvsao.com"

include "rvsao.h"

procedure xcarch (specfile)

char	specfile[ARB]		# name of archive file

char	header[48]
char	arcfile[SZ_FNAME]
char	osheader[48]
double	x,tx,y,erfy
char	tname[16]
int	it, nw, nh
#int	id
int	open(),strlen()
int	arc_fd
int	nbrec,lhead,npeaks
real	tmpwl1,tmpwl2
int	i
bool	arcwrite

include "rvsao.com"
include "results.com"
include "corr.com"
include "ansum.com"

#int     outchop                 # number of emmision lines chopped
#int     choplist[2,10]          # center pixel, +/- range in pixels
#common/emch/ outchop,choplist

begin

#  Set up analysis summary output file name
	call strcpy (specfile,arcfile,SZ_FNAME)
	call strcat (".ansum",arcfile,SZ_FNAME)
	if (debug) {
	    call printf ("XCARCH:  writing %s\n")
		call pargstr (arcfile)
	    }
	arcwrite = TRUE
	iferr {arc_fd = open (arcfile, NEW_FILE, BINARY_FILE)} then {
	    call printf ("XCARCH:  Cannot write %s\n")
		call pargstr (arcfile)
	    arcwrite = FALSE
	    }
	if (arcwrite) {

#  Write analysis summary record header
	for (i=1; i<=48; i=i+8) {
	    call strcpy ("        ", header[i], 8)
	    }
	call strcpy ("ANALYSIS_SUMMARY 40 ",header,20)
	if (debug) {
	    call printf ("XCARCH: header is %s\n")
		call pargstr (header)
	    }
	nh = 48
	call strpak (header,osheader,nh)
	nw = 24
	call write (arc_fd, osheader, nw)

#  Compute confidence of correlation velocity
        x = czr[itmax] / dsqrt (2.d0)
        tx = 1.d0 + 0.3275911d0 * x
        y = 1.d0 / tx
        erfy = 0.254829592d0*y - 0.284496736d0*y*y + 1.421413741d0*y*y*y -
               1.453152027d0*y*y*y*y + 1.061405429*y*y*y*y*y
        npeaks = 20
        cz0conf = (1.d0 - erfy) ** npeaks

#  Write analysis summary record
	nw = 2
	qcstat = spvqual
	call write (arc_fd, qcstat, nw)
	spxvel = zvel[itmax]
	spxerr = czerr[itmax]
	spxr = czr[itmax]
	czxc = spxvel
	czxcerr = spxerr
	czxcr = spxr
	cz0 = spvel
	cz0err = sperr
	if (spevel != INDEFD)
	    czem = spevel
	else
	    czem = 0.
	if (speerr != INDEFD)
	    czemerr = speerr
	else
	    czemerr = 0.
	czemscat = 0.
	nw = 18
	call write (arc_fd, cz0, nw)
	call close (arc_fd)
	}

#  Set up correlation output file name
	call strcpy (specfile,arcfile,SZ_FNAME)
	call strcat (".corr",arcfile,SZ_FNAME)
	if (debug) {
	    call printf ("XCARCH:  writing %s\n")
		call pargstr (arcfile)
	    }
	arcwrite = TRUE
	iferr {arc_fd = open (arcfile, NEW_FILE, BINARY_FILE)} then {
	    call printf ("XCARCH:  Cannot write %s\n")
		call pargstr (arcfile)
	    arcwrite = FALSE
	    }
	if (arcwrite) {

#  Write correlation record header
	nchop = 0
	nbrec = 16 + (nchop * 4) + (ntemp * 44)
	call sprintf (header, 48, "CORRELATION %d ")
	    call pargi (nbrec)
	if (debug) {
	    call printf ("XCARCH: header is %s\n")
		call pargstr (header)
	    }
	lhead = strlen (header)
	for (i=lhead+1; i<=48; i=i+1) {
	    header[i] = ' '
	    }
	nh = 48
	call strpak (header,osheader,nh)
	nw = 24
	call write (arc_fd, osheader, nw)

#  Set global information
	objrms = 0.
	ntmpl = ntemp
	nchop = 0
	tmpwl1 = twl1[itmax]
	tmpwl2 = twl2[itmax]
	objrms = 0.
	call printf ("%fA - %fA %d templates %d lines\n")
	    call pargr (tmpwl1)
	    call pargr (tmpwl2)
	    call pargs (ntmpl)
	    call pargs (nchop)

	nw = 6
	call write (arc_fd, tmpwl1, nw)
	nw = 2
	call write (arc_fd, ntmpl, nw)

#  Save information for each emission line chopped
#	do il = 1, outchop {
#	    if (specdc)
#		cchop = 10 ** (logw0 + (choplist[2,il] - 1.) * dlogw)
#	    else {
#		dw = (tmpwl2 - tmpwl1) / (specpix - 1.d0)
#		cchop = tmpwl1 + (choplist[2,il] - 1.) * dw
#		}
#	    nw = 1
#	    call write (arc_fd, cchop, nw)
#	    rchop = choplist[2,il]
#	    call write (arc_fd, rchop, nw)
#	    }

#  Save cross-correlation results for each template used
	do it = 1, ntemp {
	    call strpak (tempname[1,it], tname, 16)
	    nw = 8
	    call write (arc_fd, tname, nw)
	    tcenter = tcent[it]
	    theight = thght[it]
	    twidth = twdth[it]
	    trmsa = tarms[it]
	    trmss = tsrms[it]
	    tshft = tempshift[it] - tempvel[it] - temphcv[it]
	    tpw = taa[it]
	    nw = 14
	    call write (arc_fd, tcenter, nw)
	    }
	call close (arc_fd)
	}
	return
end
# Nov 14 1991	Move vcombine into t_xcsao.x

# Aug 11 1992	Fix analysis summary velocity writing error

# Apr 21 1994	Call WRITE as procedure, not function
# Jun 23 1994	Pass velocities in fquot, not getim labelled common
# Aug  3 1994	Change common and header from fquot to rvsao
# Nov 16 1994	Set QCSTAT from SPVQUAL in rvsao common

# Apr  6 1999	Change ntmp to ntmpl to avoid conflict
# Nov 29 1999	Deal with INDEF values for emission line vel and err
# Nov 29 1999	Print error message for each file unwritten

# Aug 25 2005	Change tname to char
