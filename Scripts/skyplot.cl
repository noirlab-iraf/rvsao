# File archive/Scripts/skyplot.cl
# April 13, 2012
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

# SKYPLOT -- Use EMSAO to plot emission lines in a night sky spectrum
#           from a spectrum file or list of spectra

procedure skyplot (spectra)

string	spectra=""	{prompt="Spectrum RFN or file"}
int	specband=3	{prompt="Spectrum band to plot (3 in IRAF multispec)"}
string	specnum=""	{prompt="Spectrum aperture range to plot"}
string	specdir=""	{prompt="Directory for spectrum to plot"}
int	nsmooth=0	{prompt="Number of times to smooth plotted spectrum"}
string	skylines="mmtsky.dat"	{prompt="File of emission lines to label"}
string	linedir="rvsao$lib/"	{prompt="Directory for lines to label"}
bool	curmode=no	{prompt="Wait for cursor commands after plotting (yes or no)"}
string	device="stdgraph"	{prompt="Display device"}

begin

string	skl
bool	cm

skl = "+" // skylines
cm = curmode

pemsao (spectra,specband=specband,specnum=specnum,specext=0,specdir=specdir,nsmooth=nsmooth,vel_init="combination",linefit=no,curmode=cm,ablines="",emlines=skl,emcomb="",linedir=linedir,logfiles="",save_vel=no,displot=yes,dispmode=3,device=device)
 
gflush

end

# Jul  7 2008	New script based on xcplot

# Apr 12 2012	Base on emplot instead of xcplot
