# File rvsao/Scripts/qplotc.cl
# May 31, 2007
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# QPLOTC -- Use XCSAO or EMSAO to plot emission and absorption lines already found
#            in a spectrum file or list of spectra

procedure qplotc (spectra)

string	spectra=""	{prompt="RFN or file"}
string	qtask="xcsao"	{prompt="program to run (xcsao or emsao)"}
string	velplot="combination"	{prompt="Velocity to plot",
				 enum="combination|emission|correlation"}
int	dispmode=4	{prompt="Display mode (2 or 4-continuum, bad lines)"}

begin

if (qtask == "emsao") {
    emsao (spectra,vel_init="combination",linefit=no,dispmode=dispmode,curmode=yes,save_vel=yes,vel_plot=velplot)
    }
else {
    xcsao (spectra,correlate=no,curmode=yes,dispmode=dispmode,save_vel=yes,vel_plot=velplot)
    }
 
end

# Sep 28 1995	New script

# Jul 16 1999	Fix so headers are updated properly
# Aug 19 1999	Add velplot parameter

# May 31 2007	New script which plots continuum-subtracted, line-chopped spectrum
