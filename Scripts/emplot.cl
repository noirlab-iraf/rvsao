# File rvsao/emplot.cl
# April 11, 2012
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# EMPLOT -- Use EMSAO to plot emission and absorption lines already found
#           in a spectrum file or list of spectra

procedure emplot (spectra)

string	spectra=""	{prompt="RFN or file"}
int     specband=1	{prompt="Spectrum band to plot (in IRAF multispec)"}
string  specnum=""	{prompt="Spectrum aperture range to plot"}
int     specext=0	{prompt="Spectrum FITS extension to plot"}
int     nsmooth=0	{prompt="Number of time sto smooth plotted spectrum"}
string  specdir=""	{prompt="Directory for spectrum to plot"}
string  ablines="ablines.dat"	{prompt="File of absorption lines to label"}
string  emlines="emlines.dat"	{prompt="File of emission lines to label"}
string  linedir="rvsao$lib/"	{prompt="Directory for lines to label"}
string	velplot="combination"	{prompt="Velocity to plot",
				 enum="combination|emission|correlation"}
bool	curmode=no	{prompt="Wait for cursor commands after plotting (yes or no)"}
string	device="stdgraph"	{prompt="Display device"}

begin
string	eml
bool	cm

eml = "+" // emlines
cm = curmode

pemsao (spectra,specband=specband,specnum=specnum,specext=specext,specdir=specdir,nsmooth=nsmooth,vel_init="combination",linefit=no,curmode=cm,ablines=ablines,emlines=eml,emcomb="",linedir=linedir,logfiles="",save_vel=no,vel_plot=velplot,displot=yes,dispmode=3,device=device)
 
gflush

end

# May 13 2008	New script based on qplot
# May 29 2008	Change name from eplot to emplot for consistency with xcplot
# Jul  7 2008	Add specband, specnum, specdir, ablines, emlines, linedir
# Jul  7 2008	Prepend + to emission line file to force line labelling
# Jul 28 2008	Set all emsao parameters from input parameters

# Aug 13 2009	Add device parameter

# Mar 11 2010	Add curmode parameter so labelled graph can be processed
# Mar 11 2010	Add gflush call to clear graphics
# Aug 09 2010	Add specext argument and call pemsao

# Apr 11 2012	Add nsmooth parameter
