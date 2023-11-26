# File rvsao/Emvel/eminitx.x
# October 3, 1995
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1995 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  EMINITX reads emission and absorption line information from files:
#  emlines contains emission lines and continuum regions
#  ablines contains absorption lines
#  emcombine contains combinations of emission lines to fit simultaneously
#  emsearch contains emission lines used to guess velocity

include "emv.h"
include "rvsao.h"

procedure eminitx (emset,all)
 
pointer	emset		# Emission line pset from which to read filenames
bool	all		# initialize combinations and guesses if true

int	fd		# File descriptor
int	lpath
int	fscan(), open(), stridx(), strlen()
char	filename[SZ_PATHNAME]
char	dirpath[SZ_PATHNAME]
char	filepath[SZ_PATHNAME]

include	"emv.com"
include "rvsao.com"

begin

#  Directory containing line lists
	call clgpseta (emset, "linedir", dirpath, SZ_PATHNAME)
	lpath = strlen (dirpath)
	if (dirpath[lpath] != '/') {
	    dirpath[lpath+1] = '/'
	    dirpath[lpath+2] = EOS
	    }

#  Absorption lines
	call clgpseta (emset, "ablines", filename, SZ_PATHNAME)
	if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
	    call strcpy (filename,filepath,SZ_PATHNAME)
	else {
	    call strcpy (dirpath,filepath,SZ_PATHNAME)
	    call strcat (filename,filepath,SZ_PATHNAME)
	    }
	fd = open (filepath, READ_ONLY, TEXT_FILE)
	nabs = 0
	if (fd >0) {
	    while (fscan(fd) != EOF && nabs < MAXABS) {
		nabs = nabs + 1
		call gargd (wlabs[nabs])
		call gargwrd (nmabs[1,nabs],SZ_ELINE)
		}
	    call close (fd)
	    }

#  Emission lines
	call clgpseta (emset, "emlines", filename, SZ_PATHNAME)
	if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
	    call strcpy (filename,filepath,SZ_PATHNAME)
	else {
	    call strcpy (dirpath,filepath,SZ_PATHNAME)
	    call strcat (filename,filepath,SZ_PATHNAME)
	    }
	fd = open (filepath, READ_ONLY, TEXT_FILE)
	nref = 0
	if (fd >0) {
	    while (fscan(fd) != EOF && nref < MAXREF) {
		nref = nref + 1
		call gargd (wlref[nref])
		call gargd (bcont[nref])
		call gargd (rcont[nref])
		call gargd (wgfit[nref])
		call gargwrd (nmref[1,nref],SZ_ELINE)
#		if (debug) {
#		    call printf ("EMINIT: %d: %s %7.2f +-%6.2f\n")
#			call pargi (nref)
#			call pargstr (nmref[1,nref])
#			call pargd (wlref[nref])
#			call pargd (wgfit[nref])
#		    }
		}
	    call close (fd)
	    }
	if (!all) return

#  Emission line combinations
	call clgpseta (emset, "emcombine",filename,SZ_PATHNAME)
	if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
	    call strcpy (filename,filepath,SZ_PATHNAME)
	else {
	    call strcpy (dirpath,filepath,SZ_PATHNAME)
	    call strcat (filename,filepath,SZ_PATHNAME)
	    }
	fd = open (filepath, READ_ONLY,TEXT_FILE)
	ncombo = 0
	if (fd >0) {
	    while (fscan (fd) != EOF && ncombo < MAXCOMB) {
		ncombo = ncombo + 1
		call gargi (numcom[ncombo])
		call gargd (edwl[ncombo])
		call gargd (elines[1,ncombo])
		call gargd (eht[1,ncombo])
		call gargd (elines[2,ncombo])
		call gargd (eht[2,ncombo])
		call gargd (elines[3,ncombo])
		call gargd (eht[3,ncombo])
		}
	    call close (fd)
	    }

#  Emission lines for velocity from one line
	call clgpseta (emset, "emsearch",filename,SZ_PATHNAME)
	if (stridx ("/",filename) > 0 || stridx ("$",filename) > 0)
	    call strcpy (filename,filepath,SZ_PATHNAME)
	else {
	    call strcpy (dirpath,filepath,SZ_PATHNAME)
	    call strcat (filename,filepath,SZ_PATHNAME)
	    }
	fd = open (filepath, READ_ONLY,TEXT_FILE)
	nsearch = 0
	if (fd >0) {
	    nsearch = 0
	    while (fscan (fd) != EOF && nsearch < MAXSEARCH) {
		nsearch = nsearch + 1
		call gargd (restwave[nsearch])
		call gargd (bluelim[nsearch])
		call gargd (redlim[nsearch])
		call gargwrd (nmsearch[1,nsearch],SZ_ELINE)
		}
	    call close (fd)
	    }

end
# Dec  1 1992	Add fit width to emission lines

# Jun 15 1994	Use separate directory parameter
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 22 1994	Consistently use sz_eline for length of line names

# Jan 31 1995	Change lengths of file and directory names to sz_pathname
