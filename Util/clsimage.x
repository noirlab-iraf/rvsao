#*** File rvsao/Util/clsimage.x
#*** July 9, 2010
#*** By Doug Mink

#  CLOSE_IMAGE -- Close IRAF image opened by GETIMAGE
 
include	<imhdr.h>
include	<smw.h>
 
procedure close_image (im, sh)
 
pointer	im		# Image header structure [returned]
pointer	sh		# Spectrum header structure [returned]

begin

	if (sh == ERR)
	    sh = NULL

#  Close spectrum WCS structure
	if (sh != NULL) {
	    if (MW(sh) != NULL) {
		call smw_close (MW(sh))
		MW(sh) = NULL
		}

#  Close spectrum header
	    call shdr_close (sh)
	    sh = NULL
	    }

#  Unmap the image
	if (im == ERR)
	    im = NULL
	if (im != NULL) {
	    call imunmap (im)
	    im = NULL
	    }

	return
end
# Mar 13 1995	New subroutine
# Mar 29 1995	Close header and image if not NULL or ERR
# Oct  5 1995	Change SHDR_CLOSE call to SPHD_CLOSE

# Aug  7 1996	Use smw.h; close MWCS structure separately

# Mar 20 1998	Drop separate closing of MWCS structure

# Nov 12 1999	Close spectral MWCS explicitly; shdr_close() does not do it

# Jul 09 2010	Don't close spectrum WCS structure if spectrum header is NULL
