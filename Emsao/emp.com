# File Emsao/emsmv.com

common/emsmv/ smplot, smspec, smcont, sfspec, sfcont
pointer smplot          # Smoothed object spectrum (wavelength-binned)
pointer smspec          # Spectrum smoothed for searching (continuum removed)
pointer smcont          # Continuum smoothed for searching and plotting
pointer sfspec          # Spectrum smoothed for fitting
pointer sfcont          # Continuum smoothed for fitting
