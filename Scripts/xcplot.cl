# File archive/Scripts/xcplot.cl
# July 7, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# XCPLOT -- Use XCSAO to plot emission and absorption lines already found
#           in a spectrum file or list of spectra

procedure xcplot (spectra)

string	spectra=""	{prompt="Spectrum RFN or file"}
string	specnum=""	{prompt="Spectrum aperture range to plot"}
int	specband=0	{prompt="Spectrum band to plot"}
string	specdir=""	{prompt="Directory for spectrum to plot"}
string	velplot="combination"	{prompt="Velocity to plot",
				 enum="combination|emission|correlation"}

begin

xcsao (spectra,specnum=specnum,specband=specband,specdir=specdir,correlate=no,logfiles="",curmode=yes,dispmode=2,save_vel=no,vel_plot=velplot,displot=yes)
 
end

# May 13 2008	New script based on qplot
# May 29 2008	Rename from xplot to xcplot to avoid name conflict
# Jul  7 2008	Add specnum, specband, and specdir parameters
