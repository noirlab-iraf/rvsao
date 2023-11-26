# File rvsao/Emvel/eminit.x
# March 11, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1992-2008 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  EMINIT reads emission and absorption line information from files:
#  emlines contains emission lines and continuum regions
#  ablines contains absorption lines
#  emcombine contains combinations of emission lines to fit simultaneously
#  emsearch contains emission lines used to guess velocity

include "emv.h"
include "rvsao.h"

procedure eminit (all)
 
bool	all		# initialize combinations and guesses if true

int	fd		# File descriptor
int	lpath, lfil, i
int	fscan(), open(), stridx(), strlen()
char	filename[SZ_PATHNAME]
char	dirpath[SZ_PATHNAME]
char	filepath[SZ_PATHNAME]
int	icombo

include	"emv.com"
include "rvsao.com"

begin

#  Directory containing line lists
	call clgstr ("linedir", dirpath, SZ_PATHNAME)
	lpath = strlen (dirpath)
	if (lpath < 1)
	    dirpath[0] = EOS
	else if (dirpath[lpath] != '/') {
	    dirpath[lpath+1] = '/'
	    dirpath[lpath+2] = EOS
	    }

#  Absorption lines
	nabs = 0
	call clgstr ("ablines", filename, SZ_PATHNAME)
	lfil = strlen (filename)
	if (lfil > 0) {
	    if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
		call strcpy (filename,filepath,SZ_PATHNAME)
	    else {
		call strcpy (dirpath,filepath,SZ_PATHNAME)
		call strcat (filename,filepath,SZ_PATHNAME)
		}
	    fd = open (filepath, READ_ONLY, TEXT_FILE)
	    if (fd >0) {
		while (fscan(fd) != EOF && nabs < MAXABS) {
		    nabs = nabs + 1
		    call gargd (wlabs[nabs])
		    call gargwrd (nmabs[1,nabs],SZ_ELINE)
		    }
		call close (fd)
	        }
	    }

#  Emission lines
	nref = 0
	call clgstr ("emlines", filename, SZ_PATHNAME)
	lfil = strlen (filename)
	if (lfil > 0) {
	    plotem = FALSE
	    if (filename[1] == '+') {
		plotem = TRUE
		do i = 1, lfil-1 {
		    filename[i] = filename[i+1]
		    }
		filename[lfil] = EOS
		}
	    if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
		call strcpy (filename,filepath,SZ_PATHNAME)
	    else {
		call strcpy (dirpath,filepath,SZ_PATHNAME)
		call strcat (filename,filepath,SZ_PATHNAME)
		}
	    fd = open (filepath, READ_ONLY, TEXT_FILE)
	    if (fd >0) {
		while (fscan(fd) != EOF && nref < MAXREF) {
		    nref = nref + 1
		    call gargd (wlref[nref])
		    bcont[nref] = 0.d0
		    call gargd (bcont[nref])
		    if (bcont[nref] > 0.d0) {
			call gargd (rcont[nref])
			call gargd (wgfit[nref])
			call gargwrd (nmref[1,nref],SZ_ELINE)
			}
		    else {
			call sprintf (nmref[1,nref], SZ_LINE, "%9.4f")
			    call pargd (wlref[nref])
			rcont[nref] = 0.d0
			wgfit[nref] = 1.d0
			}
#		    if (debug) {
#			call printf ("EMINIT: %d: %s %7.2f +-%6.2f\n")
#			    call pargi (nref)
#			    call pargstr (nmref[1,nref])
#			    call pargd (wlref[nref])
#			    call pargd (wgfit[nref])
#			}
		    }
		call close (fd)
		if (debug && plotem) {
		    call printf ("EMINIT: Always labelling %d emission lines from %s\n")
			call pargi (nref)
			call pargstr (filename)
		    }
		}
	    }
	if (!all) return

#  Emission line combinations
	ncombo = 0
	call clgstr ("emcombine",filename,SZ_PATHNAME)
	lfil = strlen (filename)
	if (lfil > 0) {
	    if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
		call strcpy (filename,filepath,SZ_PATHNAME)
	    else {
		call strcpy (dirpath,filepath,SZ_PATHNAME)
		call strcat (filename,filepath,SZ_PATHNAME)
		}
	    fd = open (filepath, READ_ONLY,TEXT_FILE)
	    if (fd >0) {
		while (fscan (fd) != EOF && ncombo < MAXCOMB) {
		    ncombo = ncombo + 1
		    call gargi (numcom[ncombo])
		    call gargd (edwl[ncombo])
		    do icombo = 1, numcom[ncombo] {
			call gargd (elines[icombo, ncombo])
			call gargd (eht[icombo, ncombo])
			if (eht[icombo, ncombo] < 0) {
			    edrop[icombo, ncombo] = TRUE
		    	    eht[icombo, ncombo] = -eht[icombo, ncombo]
			    }
			else
			    edrop[icombo, ncombo] = FALSE
			}
		    }
		call close (fd)
		}
	    }

#  Emission lines for velocity from one line
	nsearch = 0
	call clgstr ("emsearch",filename,SZ_PATHNAME)
	lfil = strlen (filename)
	if (lfil > 0) {
	    if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
		call strcpy (filename,filepath,SZ_PATHNAME)
	    else {
		call strcpy (dirpath,filepath,SZ_PATHNAME)
		call strcat (filename,filepath,SZ_PATHNAME)
		}
	    fd = open (filepath, READ_ONLY,TEXT_FILE)
	    if (fd >0) {
		while (fscan (fd) != EOF && nsearch < MAXSEARCH) {
		    nsearch = nsearch + 1
		    call gargd (restwave[nsearch])
		    call gargd (bluelim[nsearch])
		    call gargd (redlim[nsearch])
		    call gargwrd (nmsearch[1,nsearch],SZ_ELINE)
		    }
		call close (fd)
		}
	    }

end
# Dec  1 1992	Add fit width to emission lines

# Jun 15 1994	Use separate directory parameter
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 22 1994	Consistently use sz_eline for length of line names

# Jan 31 1995	Change lengths of file and directory names to sz_pathname

# Mar  7 1997	Set switches to allow lines in combos to be dropped
# Aug  1 1997	For any null file name, set that number of lines to zero

# Jan 26 2001	Initialize line directory correctly

# May 12 2004	Set plotem = true if emission line filename starts with a +

# Mar 11 2008	If no continuum is given, use center wavelength as name
