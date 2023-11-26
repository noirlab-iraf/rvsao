common/emg/ pixl,pixr,npar,ngau
int	pixl
int	pixr
int	npar
int	ngau

common/emd/ ndata, nderiv,deriv,diffbuf
int	ndata		# number of data points max 
int	nderiv		# number of params 
pointer	deriv		# derivative storage, ndata * nderiv  
	                           # array storage order is deriv[npix][npar] 
pointer	diffbuf		# Observed-Fit, per data point 
