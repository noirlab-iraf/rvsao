# Common for pixel <-> wavelength conversion

common/rvmwcs/ wl1, wl2, px1, px2, pxshift, wsh, wclog

double	wl1,wl2		# Blue and red wavelength limits
double	px1,px2		# Blue and red pixel limits
double	pxshift		# number of pixels to shift when converting
pointer	wsh		# Spectrum header structure
int	wclog		# Log wavelength flag
