#*** File arc_w2p.x
#*** October 23, 1992
#*** By Doug Mink, Center for Astrophysics

# Convert wavelength to pixel based on polynomial dispersion function

double procedure w2p (wavelength)

double wavelength

double pixel
double wl
int i

#  Common for appropriate order wavelength fit
common /wpconv/ pdum,pdim,pmid,pscale,pcoeff,dispoff
long	pdum		# Padding for RISC machines
long	pdim		# Dimension of wavelength to pixel polynomial
double	pmid		# Midpoint of independent variable in angstroms
double	pscale		# Scaling factor for wavelength
double	pcoeff[8]	# Coefficients for wavelength to pixel polynomial
int	dispoff		# Pixel offset for dispersion (1=new, 0=old)

begin

	pixel = 0.d0
	wl = (wavelength - pmid)
	if (pscale != 0.) wl = wl / pscale
	do i = pdim, 1, -1 {
	    pixel = pcoeff[i] + (wl * pixel)
	    }
	if (dispoff == 0)
	    pixel = pixel + 1

	return (pixel)
end


# Set up polynomial dispersion function for w2p from image header

procedure wpinit (im)

pointer im	# image pointer

int i

char	kstring[10]

#  Common for appropriate order wavelength fit
common /wpconv/ pdum,pdim,pmid,pscale,pcoeff,dispoff
long	pdum		# Padding for RISC machines
long	pdim		# Dimension of wavelength to pixel polynomial
double	pmid		# Midpoint of independent variable in angstroms
double	pscale		# Scaling factor for wavelength
double	pcoeff[8]	# Coefficients for wavelength to pixel polynomial
int	dispoff		# Pixel offset for dispersion (1=new, 0=old)

begin

	call imgipar (im, "WPNCOEFF", pdim)
	call imgdpar (im, "WPSCALE", pscale)
	call imgdpar (im, "WPMID", pmid)
	dispoff = 1
	call imgipar (im, "WPOFFSET", dispoff)
	do i = 1, pdim {
	    call sprintf (kstring,8,"WPCOEFF%d")
		call pargi (i)
	    call imgdpar (im, kstring, pcoeff[i])
	    }
end

# Oct 23 1992	Add dispersion function zero pixel offset
