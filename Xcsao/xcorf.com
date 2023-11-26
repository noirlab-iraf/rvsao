# File Xcsao/xcorf.com
# June 25, 2007
# By Doug Mink

# Vectors for cross-correlation results in XCFIT

common /xcorf/ xcont, spexp, xind, xifft, ft1, ft2, ftcfn, pft, tft

pointer xcont	# Returned fit of spectrum continuum
pointer spexp	# Spectrum to be transformed (may be half zeroes)
pointer xind    # Indexes for cross-correlation vector (real)
pointer xifft   # Indexes for Fourier transforms (real)
pointer ft1     # Fourier tranform of spectrum (complex)
pointer ft2     # Fourier tranform of template (complex)
pointer ftcfn	# Cross-correlation (complex)
pointer tft     # Transform of transform (complex)
pointer pft     # Power spectrum (real)

# Jun 20 2007	New file
# Jun 25 2007	Add xcont (formerly work) and spexp from xcorfit task
